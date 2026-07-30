# ECuADOR — *v2.0 Beta*
**E**xtraction and **Cu**ration of **A**rchitectural **D**omains for **O**rganelle **R**esearch

> A Perl tool for the automatic identification, extraction, and orientation of the four structural regions of circular chloroplast genomes: LSC, IRa, SSC, and IRb.

Originally developed as a Master's thesis project. The original version is also available at  
**[BiodivGenomic/ECuADOR](https://github.com/BiodivGenomic/ECuADOR)**.

---

## Table of Contents
- [Overview](#overview)
- [What's New in v2.0](#whats-new-in-v20)
- [Installation](#installation)
- [Quick Test](#quick-test)
- [Usage](#usage)
- [Options](#options)
- [Examples](#examples)
- [Output Files](#output-files)
- [Citation](#citation)

---

## Overview

Chloroplast genomes share a conserved quadripartite architecture: two inverted repeats (IRa and IRb) flanking a large single-copy region (LSC) and a small single-copy region (SSC). ECuADOR automates the detection and extraction of these four regions from one or multiple FASTA/GenBank sequence files, and can optionally orient all sequences consistently relative to a reference.

---

## What's New in v2.0

### New Options

#### `--reference`
Specify which sequence should serve as the reference for orientation:

```
--reference ID                # Select a specific sequence ID from your input directory
--reference_file path/to/ref  # Use an external reference file
```

#### `--rotate`
Apply a global rotation (in bp) to reposition the origin before the four-region search. This is especially useful when ECuADOR cannot detect the four expected regions due to an unconventional starting position. In those cases ECuADOR will automatically suggest a tentative rotation value to help reinitialize the search.

```bash
perl ECuADOR.pl --input_dir /path/to/sequences --window 800 --out ALL --rotate 50000 --orient TRUE
perl ECuADOR.pl --input_dir /path/to/sequences --window 800 --out ALL --rotate 50000 --orient TRUE --reference NC_015803
```

#### `--map`
Provide per-sequence rotation values via a two-column tab-separated file (ID and rotation in bp). Useful when a first run identifies different problematic offsets for individual sequences:

```
NC_015830    50000
NC_015826    52000
```

```bash
perl ECuADOR.pl --input_dir /path/to/sequences --window 800 --out ALL --map rotation.txt --orient TRUE
perl ECuADOR.pl --input_dir /path/to/sequences --window 800 --out ALL --map rotation.txt --orient TRUE --reference NC_015803
```

---

### Speed Improvements (~100×)

The orientation function has been completely rewritten. Instead of performing a full alignment of each region, it now examines only the first ~20 bp to determine strand orientation — like checking the first page of a book to decide whether to flip it. This achieves equivalent results in the vast majority of cases at a fraction of the cost.

---

### Bug Fixes & Code Optimization

- Modularized key processes into reusable functions for readability and maintainability.
- Fixed variable scope issues that caused problems with newer Perl installations.
- Shorter, cleaner, and more efficient codebase overall.

---

## Installation

Install required Perl modules:

```bash
sudo cpan Getopt::Long Bio::SeqIO File::Spec File::Path Cwd
```

Make the script executable:

```bash
chmod +x ECuADOR.pl
```

---

## Quick Test

```bash
perl ECuADOR.pl --input_dir TEST --window 800 --out ALL --orient TRUE
```

---

## Usage

```
perl ECuADOR.pl --input_dir <directory> [options]
```

---

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--input_dir` | path | *(required)* | Directory with sequence files (.fasta / .fa / .fna / .gb / .gbk) |
| `--window` | int | 1000 | Sliding window size (bp) for IR detection |
| `--out` | string | *(none)* | Regions to export: `ALL` or comma-separated `LSC,IRa,SSC,IRb` |
| `--ir_thresh` | int | 500 | Gap threshold (bp) for merging IR fragments |
| `--rotate` | int | 0 | Global rotation value (bp) applied to all sequences |
| `--map` | file | *(none)* | Tab-delimited file with per-sequence rotation values |
| `--orient` | TRUE/FALSE | FALSE | Enable reference-based strand orientation |
| `--reference` | string | *(none)* | Sequence ID from the input directory to use as reference |
| `--reference_file` | path | *(none)* | External file containing the reference sequence |
| `--help` | flag | — | Print usage information |

---

## Examples

**Basic extraction (no orientation):**
```bash
perl ECuADOR.pl --input_dir /path/to/sequences --window 2000 --out ALL
```

**Global rotation:**
```bash
perl ECuADOR.pl --input_dir /path/to/sequences --window 2000 --out ALL --rotate 85668
```

**Per-sequence rotation from a mapping file:**
```bash
perl ECuADOR.pl --input_dir /path/to/sequences --window 2000 --out ALL --map rotation.txt
```

**Rotation + orientation with a reference ID:**
```bash
perl ECuADOR.pl --input_dir /path/to/sequences --window 2000 --out ALL \
  --rotate 50000 --map rotation.txt --orient TRUE --reference NC_015803
```

**Orientation with an external reference file:**
```bash
perl ECuADOR.pl --input_dir /path/to/sequences --window 2000 --out ALL \
  --orient TRUE --reference_file /path/to/reference.fasta
```

---

## Output Files

ECuADOR generates two output directories depending on the options used:

### `ECuADOR_output/`
| File | Description |
|---|---|
| `Export_LSC.fasta` | LSC sequences (one per input sequence) |
| `Export_IRa.fasta` | IRa sequences |
| `Export_SSC.fasta` | SSC sequences |
| `Export_IRb.fasta` | IRb sequences |
| `non_oriented_annotations.gff3` | GFF3 annotation with embedded FASTA section |
| `Summary_Regions.txt` | Summary table, warnings, and candidate rotations |

### `ECuADOR_output_oriented/` *(when `--orient TRUE`)*
| File | Description |
|---|---|
| `Export_*_oriented.fasta` | Strand-oriented region FASTA files |
| `oriented_annotations.gff3` | GFF3 for oriented sequences with embedded FASTA section |
| `Summary_Regions_oriented.txt` | Full summary with per-sequence orientation notes |

---

## Citation

If you use ECuADOR in your research, please cite the original work:

> Original version: [https://github.com/BiodivGenomic/ECuADOR](https://github.com/BiodivGenomic/ECuADOR)

---

*© 2025 Biodiversity Genomics Team — Angelo*
