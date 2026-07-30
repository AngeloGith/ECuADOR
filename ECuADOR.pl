#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use File::Spec;
use File::Path qw(make_path);
use Cwd;

# -----------------------------------------------------------------
# Script header: Print start date/time.
# -----------------------------------------------------------------
my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900;
$mon  += 1;
printf "\nECuADOR run started at %02d/%02d/%04d %02d:%02d:%02d\n",
       $mday, $mon, $year, $hour, $min, $sec;
print "------------------------------------------------------------\n\n";

# Banner
print q(

▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀

     ███████╗░█████╗░██╗░░░██╗░█████╗░██████╗░░█████╗░██████╗░
     ██╔════╝██╔══██╗██║░░░██║██╔══██╗██╔══██╗██╔══██╗██╔══██╗
     █████╗░░██║░░╚═╝██║░░░██║███████║██║░░██║██║░░██║██████╔╝
     ██╔══╝░░██║░░██╗██║░░░██║██╔══██║██║░░██║██║░░██║██╔══██╗
     ███████╗╚█████╔╝╚██████╔╝██║░░██║██████╔╝╚█████╔╝██║░░██║
     ╚══════╝░╚════╝░░╚═════╝░╚═╝░░╚═╝╚═════╝░░╚════╝░╚═╝░░╚═╝

                     [ VERSION 2.0 BETA ]
              Sequence Curator for Plant Genomics
               © 2025 Biodiversity Genomics Team 
▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
);
print "\n\n";

# Global data structures
my %global_warnings;         # Warnings by category (Quality, Processing, Orientation, IR, etc.)
my @non_processed;           # Sequences that could not be processed
my %rotation_map;            # Mapping of sequence ID => rotation value
my @non_oriented_annotations;# Annotation coordinates for non-oriented sequences
my @oriented_annotations;    # Annotation coordinates for oriented sequences

# For oriented output (reference-based orientation)
my $ref_id;                  # Reference sequence ID used for orientation
my %ref_regions;             # Reference regions for LSC, IRa, SSC, IRb
my %oriented_regions;        # For each region type, an array of hashes { id => ..., seq => ... }
my @oriented_summary_lines;  # Summary messages specific for orientation

# Arrays to store summary info and notes for each sequence.
my @summary_lines;
my @notes_lines;

# Hashes to store global (concatenated) sequences per cpDNA.
my %global_sequences;             # For non-oriented sequences
my %global_sequences_oriented;    # For oriented sequences

# Sequence counter (global variable)
my $seqCounter = 0;

# -----------------------------------------------------------------
# Command-line parameters (with defaults)
# -----------------------------------------------------------------
my $input_dir;               # Directory with cpDNA sequence files.
my $file_format;             # Determined automatically if not provided.
my $winsize     = 1000;      # Default sliding window size.
my $out_param   = '';        # Regions to export (e.g., ALL or comma-separated list: LSC,SSC,IRa,IRb).
my $ir_thresh   = 500;       # Default IR fusion threshold (in bp).
my $rotate      = 0;         # Global rotation value (in bp); default is 0 (none).
my $map_file    = "";        # Mapping file (ID and rotation value per line).
my $orient_flag = "FALSE";   # Orientation flag; if TRUE, perform reference-based orientation.
my $reference_id;            # ID of the sequence to use as reference
my $reference_file;          # Path to a separate reference file
my $help;

GetOptions(
    'input_dir=s'     => \$input_dir,
    'format=s'        => \$file_format,
    'window=i'        => \$winsize,
    'out=s'           => \$out_param,
    'ir_thresh=i'     => \$ir_thresh,
    'rotate=i'        => \$rotate,
    'map=s'           => \$map_file,
    'orient=s'        => \$orient_flag,
    'reference=s'     => \$reference_id,        # New: sequence ID as reference
    'reference_file=s'=> \$reference_file,      # New: separate reference file
    'help'            => \$help,
) or die "Error in command-line arguments\n";

if ($help or not defined $input_dir) {
    print_usage();
    exit;
}

# -----------------------------------------------------------------
# Create output directories.
# -----------------------------------------------------------------
my $output_dir = "ECuADOR_output";
unless (-d $output_dir) {
    make_path($output_dir) or die "Failed to create output directory $output_dir: $!";
}
my $oriented_output_dir = $output_dir;
if (uc($orient_flag) eq "TRUE") {
    $oriented_output_dir = "ECuADOR_output_oriented";
    unless (-d $oriented_output_dir) {
        make_path($oriented_output_dir) or die "Failed to create oriented output directory $oriented_output_dir: $!";
    }
    print "INFO: Orientation option activated. Oriented outputs will be saved in $oriented_output_dir\n\n";
    
    # Print reference information
    if (defined $reference_file) {
        print "INFO: Using sequence from file '$reference_file' as reference for orientation.\n\n";
    } elsif (defined $reference_id) {
        print "INFO: Will use sequence with ID '$reference_id' as reference for orientation.\n\n";
    } else {
        print "WARNING: No reference specified (use --reference or --reference_file). The first processed sequence will be used as reference.\n\n";
    }
}

if ($rotate > 0) {
    print "WARNING: Global rotation specified as $rotate bp. It will be applied to all sequences not covered by the mapping file.\n\n";
}
if ($map_file) {
    print "INFO: Mapping file '$map_file' provided; individual rotation values will be used when available.\n\n";
}

# -----------------------------------------------------------------
# Function: print_usage
# -----------------------------------------------------------------
sub print_usage {
    print <<"END_USAGE";
Usage: perl $0 --input_dir <directory> [--window <size>] [--out <regions>] [--ir_thresh <bp>] [--rotate <bp>] [--map <file>] [--orient <TRUE|FALSE>] [--reference <ID>] [--reference_file <file>]
Examples:
  Basic (no orientation):
    perl $0 --input_dir /path/to/sequences --window 2000 --out ALL
  Global rotation:
    perl $0 --input_dir /path/to/sequences --window 2000 --out ALL --rotate 85668
  Rotation using a mapping file:
    perl $0 --input_dir /path/to/sequences --window 2000 --out ALL --map rotation.txt
  Combined with orientation:
    perl $0 --input_dir /path/to/sequences --window 2000 --out ALL --rotate 50000 --map rotation.txt --orient TRUE
  With specific reference sequence ID:
    perl $0 --input_dir /path/to/sequences --window 2000 --out ALL --orient TRUE --reference NC_015803
  With external reference file:
    perl $0 --input_dir /path/to/sequences --window 2000 --out ALL --orient TRUE --reference_file /path/to/reference.fasta
Note: Assumes all chloroplast sequences are circular.
END_USAGE
}

# -----------------------------------------------------------------
# Function: detect_file_format
# -----------------------------------------------------------------
sub detect_file_format {
    my ($filename) = @_;
    if ($filename =~ /\.(gb|gbk|genbank)$/i) {
        return 'genbank';
    } elsif ($filename =~ /\.(fa|fasta|fna)$/i) {
        return 'fasta';
    }
    return $file_format || 'fasta';
}

# -----------------------------------------------------------------
# Function: reverse_complement
# -----------------------------------------------------------------
sub reverse_complement {
    my ($seq) = @_;
    my $rev = reverse($seq);
    $rev =~ tr/ACGTacgt/TGCAtgca/;
    return $rev;
}

# -----------------------------------------------------------------
# Function: merge_intervals
# -----------------------------------------------------------------
sub merge_intervals {
    my ($nums_ref, $threshold) = @_;
    my @nums = sort { $a <=> $b } @$nums_ref;
    my @intervals;
    return () unless @nums;
    my $start = shift @nums;
    my $prev  = $start;
    foreach my $n (@nums) {
        if ($n - $prev <= $threshold) {
            $prev = $n;
        } else {
            push @intervals, [$start, $prev];
            $start = $n;
            $prev = $n;
        }
    }
    push @intervals, [$start, $prev];
    return @intervals;
}

# -----------------------------------------------------------------
# Function: merge_interval_pairs
# -----------------------------------------------------------------
sub merge_interval_pairs {
    my ($intervals_ref, $gap_thresh) = @_;
    return () unless @$intervals_ref;
    my @sorted = sort { $a->[0] <=> $b->[0] } @$intervals_ref;
    my @merged;
    my $current = [ $sorted[0]->[0], $sorted[0]->[1] ];
    foreach my $interval (@sorted) {
        if ($interval->[0] - $current->[1] <= $gap_thresh) {
            $current->[1] = $interval->[1] if $interval->[1] > $current->[1];
        } else {
            push @merged, $current;
            $current = [ $interval->[0], $interval->[1] ];
        }
    }
    push @merged, $current;
    return @merged;
}

# -----------------------------------------------------------------
# Function: choose_largest_interval
# -----------------------------------------------------------------
sub choose_largest_interval {
    my @intervals = @_;
    return unless @intervals;
    my $largest = $intervals[0];
    foreach my $int (@intervals) {
        if (($int->[1] - $int->[0]) > ($largest->[1] - $largest->[0])) {
            $largest = $int;
        }
    }
    return $largest;
}

# -----------------------------------------------------------------
# Function: check_sequence_quality
# -----------------------------------------------------------------
sub check_sequence_quality {
    my ($seq, $id) = @_;
    my $length = length($seq);
    my $n_count = () = $seq =~ /N/gi;
    my $gap_count = () = $seq =~ /-/g;
    my $ambiguous_count = () = $seq =~ /[RYWSKMBDHV]/gi;
    my $n_thresh = 0.10;
    my $gap_thresh = 0.05;
    my $ambiguous_thresh = 0.10;
    if ($n_count / $length > $n_thresh) {
       push @{ $global_warnings{'Quality'} },
         sprintf("Sequence %s: High proportion of N's: %.2f%%", $id, 100 * $n_count / $length);
    }
    if ($gap_count / $length > $gap_thresh) {
       push @{ $global_warnings{'Quality'} },
         sprintf("Sequence %s: High proportion of gaps (-): %.2f%%", $id, 100 * $gap_count / $length);
    }
    if ($ambiguous_count / $length > $ambiguous_thresh) {
       push @{ $global_warnings{'Quality'} },
         sprintf("Sequence %s: High proportion of ambiguous bases: %.2f%%", $id, 100 * $ambiguous_count / $length);
    }
}

# -----------------------------------------------------------------
# Function: extract_region
# -----------------------------------------------------------------
sub extract_region {
    my ($seqobj, $start, $end, $len) = @_;
    return "" if ($start >= $end);
    if ($end <= $len) {
       return $seqobj->subseq($start, $end);
    } else {
       my $part1 = $seqobj->subseq($start, $len);
       my $part2 = $seqobj->subseq(1, $end - $len);
       return $part1 . $part2;
    }
}

# -----------------------------------------------------------------
# Function: rotate_sequence
# -----------------------------------------------------------------
sub rotate_sequence {
    my ($seqobj, $rotate, $len) = @_;
    my $part1 = $seqobj->subseq($rotate + 1, $len);
    my $part2 = $seqobj->subseq(1, $rotate);
    my $newseq = $part1 . $part2;
    $seqobj->seq($newseq);
    return $seqobj;
}

# -----------------------------------------------------------------
# Function: compute_IRa_start
# -----------------------------------------------------------------
sub compute_IRa_start {
    my ($seq, $winsize) = @_;
    my $len = length($seq);
    my @ir_a_starts;
    my %fragments;
    for (my $i = 1; $i <= $len; $i++) {
         my $frag;
         if ($i + $winsize - 1 > $len) {
             my $part1 = substr($seq, $i - 1);
             my $part2 = substr($seq, 0, ($i + $winsize - 1) - $len);
             $frag = $part1 . $part2;
         } else {
             $frag = substr($seq, $i - 1, $winsize);
         }
         my $frag_rc = reverse_complement($frag);
         push @{ $fragments{$frag} }, $i;
         if (exists $fragments{$frag_rc}) {
             foreach my $prev (@{ $fragments{$frag_rc} }) {
                 push @ir_a_starts, $prev;
             }
         }
    }
    my @merged = merge_intervals(\@ir_a_starts, 5);
    return @merged ? $merged[0]->[0] : 0;
}

# -----------------------------------------------------------------
# Function: suggest_rotation
# -----------------------------------------------------------------
sub suggest_rotation {
    my ($seqobj, $winsize) = @_;
    my $orig_seq = $seqobj->seq;
    my $len = length($orig_seq);
    my $best_candidate = 0;
    my $best_IRa_start = 0;
    for (my $cand = 1000; $cand <= 10000; $cand += 100) {
         next if ($cand >= $len);
         my $rotated = substr($orig_seq, $cand) . substr($orig_seq, 0, $cand);
         my $IRa_start = compute_IRa_start($rotated, $winsize);
         if ($IRa_start > 1 and $IRa_start > $best_IRa_start) {
              $best_IRa_start = $IRa_start;
              $best_candidate = $cand;
         }
    }
    return $best_candidate;
}

# -----------------------------------------------------------------
# Function: orient_region
# Compares the first 20 bp of a region with the reference region and returns
# the region oriented to match the reference.
# -----------------------------------------------------------------
sub orient_region {
    my ($region, $ref) = @_;
    my $sub_len = 20;
    $sub_len = length($ref) if length($ref) < 20;
    my $ref_sub = substr($ref, 0, $sub_len);
    my $region_sub = substr($region, 0, $sub_len);
    my $region_rc = reverse_complement($region);
    my $region_rc_sub = substr($region_rc, 0, $sub_len);
    my $dist_normal = 0;
    my $dist_rc = 0;
    for (my $i = 0; $i < $sub_len; $i++) {
         $dist_normal++ if substr($region_sub, $i, 1) ne substr($ref_sub, $i, 1);
         $dist_rc++ if substr($region_rc_sub, $i, 1) ne substr($ref_sub, $i, 1);
    }
    return ($dist_rc < $dist_normal) ? $region_rc : $region;
}

# Variables para los arrays/hashes de exportación de regiones
my @regions_to_export;
my %exported_regions;

# -----------------------------------------------------------------
# Configure export of regions.
# -----------------------------------------------------------------
if ($out_param eq '' or uc($out_param) eq 'NONE') {
    @regions_to_export = ();
} elsif (uc($out_param) eq 'ALL') {
    @regions_to_export = qw(LSC IRa SSC IRb);
} else {
    @regions_to_export = split /,/, $out_param;
    @regions_to_export = map { uc($_) } map { s/^\s+|\s+$//gr } @regions_to_export;
}

# -----------------------------------------------------------------
# Function: process_sequence
# Processes a single sequence, with reference detection
# -----------------------------------------------------------------
sub process_sequence {
    my ($seqobj, $is_reference) = @_;
    
    $seqCounter++;  # Incrementa el contador global de secuencias
    my $id = $seqobj->display_id;
    my $len = $seqobj->length;
    
    # Pre-process: Clean and uppercase.
    my $raw = $seqobj->seq;
    $raw =~ s/[^ACGTacgt]//g;
    $seqobj->seq(uc($raw));
    $len = length($seqobj->seq);
    if ($len < $winsize) {
        my $msg = "Sequence $id is too short (length $len < $winsize). Skipping.";
        push @{ $global_warnings{'Processing'} }, $msg;
        push @non_processed, $msg;
        return 0; # Failed to process
    }
    
    # Adjust parameters for long cpDNA.
    my $local_winsize = $winsize;
    my $local_ir_thresh = $ir_thresh;
    if ($len > 180000) {
        $local_winsize = 2000;
        $local_ir_thresh = 1000;
        push @{ $global_warnings{'LongCP'} },
          "Sequence $id is long ($len bp). Adjusted window size to $local_winsize and IR threshold to $local_ir_thresh.";
    }
    
    # Determine rotation value.
    my $current_rotate = 0;
    if ($map_file and exists $rotation_map{$id}) {
        $current_rotate = $rotation_map{$id};
    } elsif ($rotate > 0 and $rotate < $len) {
        $current_rotate = $rotate;
    }
    if ($current_rotate) {
        print "INFO: Rotating sequence $id by $current_rotate bp.\n";
        $seqobj = rotate_sequence($seqobj, $current_rotate, $len);
    }
    
    # Quality Check.
    check_sequence_quality($seqobj->seq, $id);
    
    # IR Detection via Sliding Window.
    my (@ir_a_starts, @ir_a_ends, @ir_b_starts, @ir_b_ends);
    my (@ir_a_intervals, @ir_b_intervals);
    my %fragments;
    for (my $i = 1; $i <= $len; $i++) {
        my ($frag, $win_end);
        if ($i + $local_winsize - 1 > $len) {
            my $part1 = $seqobj->subseq($i, $len);
            my $part2 = $seqobj->subseq(1, ($i + $local_winsize - 1) - $len);
            $frag = $part1 . $part2;
            $win_end = $i + $local_winsize - 1;
        } else {
            $frag = $seqobj->subseq($i, $i + $local_winsize - 1);
            $win_end = $i + $local_winsize - 1;
        }
        my $frag_rc = reverse_complement($frag);
        push @{ $fragments{$frag} }, $i;
        if (exists $fragments{$frag_rc}) {
            foreach my $prev (@{ $fragments{$frag_rc} }) {
                my $prev_end = $prev + $local_winsize - 1;
                push @ir_a_starts, $prev;
                push @ir_a_ends,   $prev_end;
                push @ir_b_starts, $i;
                push @ir_b_ends,   $win_end;
                push @ir_a_intervals, [$prev, $prev_end];
                push @ir_b_intervals, [$i, $win_end];
            }
        }
    }
    unless (@ir_a_starts && @ir_a_ends && @ir_b_starts && @ir_b_ends) {
        my $msg = "Could not identify IR regions in $id.";
        push @{ $global_warnings{'IR'} }, $msg;
        push @non_processed, $msg;
        return 0; # Failed to process
    }
    
    # Strict merging.
    my $strict_thresh = 5;
    my @merged_ir_a_starts = merge_intervals(\@ir_a_starts, $strict_thresh);
    my @merged_ir_a_ends   = merge_intervals(\@ir_a_ends,   $strict_thresh);
    my @merged_ir_b_starts = merge_intervals(\@ir_b_starts, $strict_thresh);
    my @merged_ir_b_ends   = merge_intervals(\@ir_b_ends,   $strict_thresh);
    my $len_a = scalar(@merged_ir_a_ends);
    my $len_b = scalar(@merged_ir_b_ends);
    if ($len_a == 0 or $len_b == 0
        or !defined($merged_ir_a_ends[$len_a-1]) or !defined($merged_ir_a_ends[$len_a-1]->[1])
        or !defined($merged_ir_b_ends[$len_b-1]) or !defined($merged_ir_b_ends[$len_b-1]->[1])) {
        my $msg = "Could not merge IR regions in $id.";
        push @{ $global_warnings{'IR'} }, $msg;
        push @non_processed, $msg;
        return 0; # Failed to process
    }
    my $IRa_start = $merged_ir_a_starts[0]->[0];
    my $last_ir_a_end_ref = $merged_ir_a_ends[$len_a-1];
    my $IRa_end   = $last_ir_a_end_ref->[1];
    my $IRb_start = $merged_ir_b_starts[0]->[0];
    my $last_ir_b_end_ref = $merged_ir_b_ends[$len_b-1];
    my $IRb_end   = $last_ir_b_end_ref->[1];
    if ($IRa_end == 0 or $IRb_end == 0) {
        my $msg = "Merged IR intervals not defined for $id.";
        push @{ $global_warnings{'IR'} }, $msg;
        push @non_processed, $msg;
        return 0; # Failed to process
    }
    
    # Irregular reconstruction.
    my @merged_ir_a_irregular = merge_interval_pairs(\@ir_a_intervals, $local_ir_thresh);
    my @merged_ir_b_irregular = merge_interval_pairs(\@ir_b_intervals, $local_ir_thresh);
    my $reconstruction_warning = "";
    my $reconstructed_flag = 0;
    if (scalar(@merged_ir_a_irregular) == 1 && scalar(@merged_ir_b_irregular) == 1) {
        my $selected_ir_a = $merged_ir_a_irregular[0];
        my $selected_ir_b = $merged_ir_b_irregular[0];
        my $len_ir_a = $selected_ir_a->[1] - $selected_ir_a->[0] + 1;
        my $len_ir_b = $selected_ir_b->[1] - $selected_ir_b->[0] + 1;
        if (abs($len_ir_a - $len_ir_b) > ($len_ir_a * 0.1)) {
            $reconstruction_warning = "Sequence $id: IR regions differ significantly in length (IRa: $len_ir_a bp, IRb: $len_ir_b bp)";
            push @{ $global_warnings{'IR'} }, $reconstruction_warning;
        } else {
            $IRa_start = $selected_ir_a->[0];
            $IRa_end   = $selected_ir_a->[1];
            $IRb_start = $selected_ir_b->[0];
            $IRb_end   = $selected_ir_b->[1];
            $reconstructed_flag = 1;
        }
    } else {
        $reconstruction_warning = "Sequence $id: Multiple IR fragments detected (IRa fragments: " .
            scalar(@merged_ir_a_irregular) . ", IRb fragments: " . scalar(@merged_ir_b_irregular) . ")";
        push @{ $global_warnings{'IR'} }, $reconstruction_warning;
    }
    
    # Orientation and consistency checks.
    if ($IRa_start > 100) {
        push @{ $global_warnings{'Orientation'} },
          "Sequence $id: Does not appear to be normally oriented (LSC at beginning). IRa starts at $IRa_start.";
    }
    if ($IRa_end >= $IRb_start) {
        my $msg = "In $id, IRa_end ($IRa_end) is greater or equal to IRb_start ($IRb_start).";
        push @{ $global_warnings{'Orientation'} }, $msg;
        push @non_processed, $msg;
        if ($IRa_start == 1) {
            my $suggested = suggest_rotation($seqobj, $local_winsize);
            if ($suggested > 0) {
                push @{ $global_warnings{'Orientation'} },
                  "Sequence $id: Suggested rotation = $suggested (automatically computed).";
                $rotation_map{$id} = $suggested;
            }
        } else {
            my $candidate = $IRb_start - 1;
            push @{ $global_warnings{'Orientation'} },
              "Sequence $id: Suggested rotation = $candidate.";
            $rotation_map{$id} = $candidate;
        }
        return 0; # Failed to process
    }
    
    # Extract regions using circular extraction.
    my $LSC = extract_region($seqobj, 1, $IRa_start - 1, $len);
    my $IRa = extract_region($seqobj, $IRa_start, $IRa_end, $len);
    my $SSC = extract_region($seqobj, $IRa_end + 1, $IRb_start - 1, $len);
    my $IRb = extract_region($seqobj, $IRb_start, $IRb_end, $len);
    
    # Store annotation coordinates for non-oriented sequences.
    my %anno = (
        id        => $id,
        LSC_start => 1,
        LSC_end   => $IRa_start - 1,
        IRa_start => $IRa_start,
        IRa_end   => $IRa_end,
        SSC_start => $IRa_end + 1,
        SSC_end   => $IRb_start - 1,
        IRb_start => $IRb_start,
        IRb_end   => $IRb_end
    );
    push @non_oriented_annotations, \%anno;
    
    # Build global (concatenated) sequence for this cpDNA (non-oriented).
    my $global_seq = $LSC . $IRa . $SSC . $IRb;
    $global_sequences{$id} = $global_seq;
    
    # If orientation is active, perform reference-based orientation and store oriented regions.
    my $orientation_note = "";
    if (uc($orient_flag) eq "TRUE") {
        if (!defined $ref_id) {
            # This is the first sequence processed and should be the reference
            $ref_id = $id;
            $ref_regions{LSC} = $LSC;
            $ref_regions{IRa} = $IRa;
            $ref_regions{SSC} = $SSC;
            $ref_regions{IRb} = $IRb;
            $orientation_note = "Reference used for orientation: $ref_id";
            push @oriented_summary_lines, $orientation_note;
        } else {
            # Orient this sequence against the reference
            $LSC = orient_region($LSC, $ref_regions{LSC});
            $IRa = orient_region($IRa, $ref_regions{IRa});
            $SSC = orient_region($SSC, $ref_regions{SSC});
            $IRb = orient_region($IRb, $ref_regions{IRb});
            $orientation_note = "Oriented using reference: $ref_id";
            push @oriented_summary_lines, "Oriented $id using reference: $ref_id";
        }
        
        # Store oriented regions.
        push @{ $oriented_regions{LSC} }, { id => $id, seq => $LSC };
        push @{ $oriented_regions{IRa} }, { id => $id, seq => $IRa };
        push @{ $oriented_regions{SSC} }, { id => $id, seq => $SSC };
        push @{ $oriented_regions{IRb} }, { id => $id, seq => $IRb };
        
        # Build global oriented sequence.
        my $global_seq_oriented = $LSC . $IRa . $SSC . $IRb;
        $global_sequences_oriented{$id} = $global_seq_oriented;
        
        # Store annotation coordinates for oriented sequences (using same format as non-oriented)
        my $lsc_len = length($LSC);
        my $ira_len = length($IRa);
        my $ssc_len = length($SSC);
        my $irb_len = length($IRb);
        
        my %oriented_anno = (
            id        => $id,
            LSC_start => 1,
            LSC_end   => $lsc_len,
            IRa_start => $lsc_len + 1,
            IRa_end   => $lsc_len + $ira_len,
            SSC_start => $lsc_len + $ira_len + 1,
            SSC_end   => $lsc_len + $ira_len + $ssc_len,
            IRb_start => $lsc_len + $ira_len + $ssc_len + 1,
            IRb_end   => $lsc_len + $ira_len + $ssc_len + $irb_len
        );
        push @oriented_annotations, \%oriented_anno;
    }
    
    # Build main summary (without notes) and store note separately.
    my $main_summary = sprintf("Done Seq#%d\t%s\tTotal length: %d\tLSC: 1-%d\tIRa: %d-%d\tSSC: %d-%d\tIRb: %d-%d",
        $seqCounter, $id, $len, ($IRa_start - 1), $IRa_start, $IRa_end, ($IRa_end + 1), ($IRb_start - 1), $IRb_start, $IRb_end);
    push @summary_lines, $main_summary;
    push @notes_lines, ($orientation_note ? "$id: $orientation_note" : "$id: (no note)");
    
    # For non-oriented FASTA export, modify header to include region name for uniqueness.
    if (@regions_to_export) {  # Usa la variable global @regions_to_export
        my %region_map_local = (
            'LSC' => $LSC,
            'IRA' => $IRa,
            'IRa' => $IRa,
            'IRB' => $IRb,
            'IRb' => $IRb,
            'SSC' => $SSC
        );
        foreach my $region (@regions_to_export) {
            my $seq = $region_map_local{$region} || next;
            $exported_regions{$region} .= ">$id\_$region\n$seq\n";  # Usa la variable global %exported_regions
        }
    }
    print "$main_summary\n";
    
    return 1; # Successfully processed
}

# -----------------------------------------------------------------
# Load mapping file if provided.
# -----------------------------------------------------------------
if ($map_file) {
    open(my $mf, '<', $map_file) or die "Cannot open mapping file $map_file: $!\n";
    while (<$mf>) {
        chomp;
        next if /^\s*$/;
        my ($id_key, $rot_val) = split(/\s+/, $_);
        if (defined $id_key and defined $rot_val and $rot_val =~ /^\d+$/) {
            $rotation_map{$id_key} = $rot_val;
        }
    }
    close($mf);
}

# -----------------------------------------------------------------
# First, process reference file if provided
# -----------------------------------------------------------------
if (uc($orient_flag) eq "TRUE" && defined $reference_file) {
    if (-e $reference_file) {
        my $ref_format = detect_file_format($reference_file);
        my $ref_seqio = Bio::SeqIO->new(-file => $reference_file, -format => $ref_format);
        
        my $ref_seq = $ref_seqio->next_seq();
        if ($ref_seq) {
            print "INFO: Processing reference sequence from file: $reference_file\n";
            process_sequence($ref_seq, 1); # Process as reference
        } else {
            die "ERROR: Could not extract a valid sequence from reference file: $reference_file\n";
        }
    } else {
        die "ERROR: Reference file not found: $reference_file\n";
    }
}

# -----------------------------------------------------------------
# Track processed IDs and sequences for reference-based processing
# -----------------------------------------------------------------
my %processed_ids;            # Keep track of processed sequence IDs

# -----------------------------------------------------------------
# Main Processing.
# -----------------------------------------------------------------
$seqCounter = 0;  # Resetear contador de secuencias
opendir(my $dh, $input_dir) or die "Cannot open '$input_dir': $!\n";
my @files = grep { /\.(fa|fasta|fna|gb|gbk|genbank)$/i } readdir($dh);
closedir($dh);
die "No sequence files found in $input_dir\n" unless @files;

# First pass - find and process reference sequence if specified by ID
if (uc($orient_flag) eq "TRUE" && defined $reference_id && !defined $ref_id) {
    print "INFO: Looking for reference sequence with ID: $reference_id\n";
    my $found_reference = 0;
    
    FILE_REF: foreach my $file (@files) {
        my $file_path = File::Spec->catfile($input_dir, $file);
        my $format = detect_file_format($file);
        my $seqio = Bio::SeqIO->new(-file => $file_path, -format => $format);
        
        while (my $seqobj = $seqio->next_seq) {
            my $id = $seqobj->display_id;
            $processed_ids{$id} = 1;
            
            if ($id eq $reference_id) {
                # Found the reference sequence
                print "INFO: Found reference sequence $id in file $file\n";
                $found_reference = 1;
                process_sequence($seqobj, 1); # Process as reference
                last FILE_REF;
            }
        }
    }
    
    if (!$found_reference) {
        die "ERROR: Reference sequence ID '$reference_id' not found in any input file.\n";
    }
}

# Second pass - process all other sequences
FILE: foreach my $file (@files) {
    my $file_path = File::Spec->catfile($input_dir, $file);
    my $format = detect_file_format($file);
    my $seqio = Bio::SeqIO->new(-file => $file_path, -format => $format);
    
    while (my $seqobj = $seqio->next_seq) {
        my $id = $seqobj->display_id;
        
        # Skip if already processed (e.g., reference sequence)
        next if exists $processed_ids{$id};
        $processed_ids{$id} = 1;
        
        # If orientation is on but no reference yet, this becomes the reference
        if (uc($orient_flag) eq "TRUE" && !defined $ref_id) {
            print "INFO: Using sequence $id as reference (first sequence processed)\n";
            process_sequence($seqobj, 1); # Process as reference
        } else {
            # Regular sequence processing
            process_sequence($seqobj, 0); # Not a reference
        }
    }
}

print "\n\n";

# Export non-oriented region FASTA files.
if (@regions_to_export) {
    foreach my $region (@regions_to_export) {
        next unless defined $exported_regions{$region};
        my $outfile = File::Spec->catfile($output_dir, "Export_$region.fasta");
        open(my $fh, '>', $outfile) or die "Cannot open $outfile: $!\n";
        print $fh $exported_regions{$region};
        close($fh);
        print "Region $region exported to $outfile\n";
    }
}

print "\n";

# If orientation is active, export oriented region FASTA files.
if (uc($orient_flag) eq "TRUE") {
    foreach my $region (qw(LSC IRa SSC IRb)) {
         my $outfile = File::Spec->catfile($oriented_output_dir, "Export_${region}_oriented.fasta");
         open(my $fh, '>', $outfile) or die "Cannot open $outfile: $!\n";
         if (exists $oriented_regions{$region}) {
             foreach my $entry (@{ $oriented_regions{$region} }) {
                 print $fh ">$entry->{id}\_$region\n$entry->{seq}\n";
             }
         }
         close($fh);
         print "Oriented region $region exported to $outfile\n";
    }
}

print "\n";

# Export GFF3 annotations for non-oriented sequences.
my $non_oriented_gff3 = File::Spec->catfile($output_dir, "non_oriented_annotations.gff3");
open(my $ngff, '>', $non_oriented_gff3) or die "Cannot open $non_oriented_gff3: $!\n";
print $ngff "##gff-version 3\n\n";
foreach my $anno (@non_oriented_annotations) {
    my $id = $anno->{id};
    print $ngff join("\t", $id, "ECuADOR", "LSC", 1, $anno->{LSC_end}, ".", "+", ".", "ID=LSC;Name=LSC"), "\n";
    print $ngff join("\t", $id, "ECuADOR", "IRa", $anno->{IRa_start}, $anno->{IRa_end}, ".", "+", ".", "ID=IRa;Name=IRa"), "\n";
    print $ngff join("\t", $id, "ECuADOR", "SSC", $anno->{SSC_start}, $anno->{SSC_end}, ".", "+", ".", "ID=SSC;Name=SSC"), "\n";
    print $ngff join("\t", $id, "ECuADOR", "IRb", $anno->{IRb_start}, $anno->{IRb_end}, ".", "+", ".", "ID=IRb;Name=IRb"), "\n";
}
# Add FASTA section to non-oriented GFF3.
print $ngff "##FASTA\n";
foreach my $id (sort keys %global_sequences) {
    print $ngff ">$id\n$global_sequences{$id}\n";
}
close($ngff);
print "GFF3 annotations for non-oriented cpDNA exported to $non_oriented_gff3\n";

# If orientation is active, export a combined oriented GFF3 file (with FASTA section).
if (uc($orient_flag) eq "TRUE") {
    my $oriented_gff3 = File::Spec->catfile($oriented_output_dir, "oriented_annotations.gff3");
    open(my $ogff, '>', $oriented_gff3) or die "Cannot open $oriented_gff3: $!\n";
    print $ogff "##gff-version 3\n\n";
    
    # Use the same format as non-oriented GFF3
    foreach my $anno (@oriented_annotations) {
        my $id = $anno->{id};
        print $ogff join("\t", $id, "ECuADOR", "LSC", $anno->{LSC_start}, $anno->{LSC_end}, ".", "+", ".", "ID=LSC;Name=LSC"), "\n";
        print $ogff join("\t", $id, "ECuADOR", "IRa", $anno->{IRa_start}, $anno->{IRa_end}, ".", "+", ".", "ID=IRa;Name=IRa"), "\n";
        print $ogff join("\t", $id, "ECuADOR", "SSC", $anno->{SSC_start}, $anno->{SSC_end}, ".", "+", ".", "ID=SSC;Name=SSC"), "\n";
        print $ogff join("\t", $id, "ECuADOR", "IRb", $anno->{IRb_start}, $anno->{IRb_end}, ".", "+", ".", "ID=IRb;Name=IRb"), "\n";
    }
    
    # Add FASTA section for oriented sequences.
    print $ogff "##FASTA\n";
    foreach my $id (sort keys %global_sequences_oriented) {
        print $ogff ">$id\n$global_sequences_oriented{$id}\n";
    }
    close($ogff);
    print "GFF3 annotations for oriented cpDNA exported to $oriented_gff3\n\n";
}

# Export summary file for non-oriented output.
my $summary_file = File::Spec->catfile($output_dir, "Summary_Regions.txt");
open(my $sfh, '>', $summary_file) or die "Cannot open $summary_file: $!\n";
print $sfh join("\n", @summary_lines), "\n\n";
print $sfh "Warnings:\n";
foreach my $category (sort keys %global_warnings) {
    print $sfh "\n-- $category Warnings --\n";
    foreach my $warn_msg (@{ $global_warnings{$category} }) {
         print $sfh "$warn_msg\n";
    }
}
print $sfh "\n\nNon-processed Sequences:\n";
foreach my $np (@non_processed) {
    print $sfh "$np\n";
}
if (%rotation_map) {
    print $sfh "\n\nCandidate Rotations:\n";
    foreach my $seq_id (sort keys %rotation_map) {
         print $sfh "Sequence $seq_id: Suggested rotation = $rotation_map{$seq_id}\n";
    }
}
close($sfh);
print "Summary exported to $summary_file\n";

# If orientation is active, export a separate oriented summary file.
if (uc($orient_flag) eq "TRUE") {
    my $oriented_summary_file = File::Spec->catfile($oriented_output_dir, "Summary_Regions_oriented.txt");
    open(my $osfh, '>', $oriented_summary_file) or die "Cannot open $oriented_summary_file: $!\n";
    print $osfh join("\n", @summary_lines), "\n\n";
    print $osfh "Orientation Summary:\n";
    print $osfh join("\n", @oriented_summary_lines), "\n\n";
    print $osfh "Notes for each sequence:\n";
    print $osfh join("\n", @notes_lines), "\n\n";
    print $osfh "Warnings:\n";
    foreach my $category (sort keys %global_warnings) {
        print $osfh "\n-- $category Warnings --\n";
        foreach my $warn_msg (@{ $global_warnings{$category} }) {
             print $osfh "$warn_msg\n";
        }
    }
    print $osfh "\n\nNon-processed Sequences:\n";
    foreach my $np (@non_processed) {
        print $osfh "$np\n";
    }
    if (%rotation_map) {
        print $osfh "\n\nCandidate Rotations:\n";
        foreach my $seq_id (sort keys %rotation_map) {
             print $osfh "Sequence $seq_id: Suggested rotation = $rotation_map{$seq_id}\n";
        }
    }
    close($osfh);
    print "Oriented summary exported to $oriented_summary_file\n\n";
}

# Final message.
($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900; 
$mon  += 1;
printf "\nECuADOR run finished at %02d/%02d/%04d %02d:%02d:%02d\n",
       $mday, $mon, $year, $hour, $min, $sec;
print "------------------------------------------------------------\n";