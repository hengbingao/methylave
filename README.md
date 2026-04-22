# methylave

Calculate CG methylation levels from **allc** files over genomic regions defined in BED files.

## Features

| | |
|---|---|
| **Two modes** | `summary` (one row per BED file) · `region` (one row per region) |
| **Two methods** | `method1` coverage-weighted · `method2` equal-weight per CpG site |
| **Any BED input** | single file or multiple files at once |
| **Parallel** | multiple BED files processed concurrently |
| **No tabix needed** | uses bedtools only |

---

## Requirements

| Tool | Notes |
|------|-------|
| Python ≥ 3.6 | standard library only |
| [bedtools](https://bedtools.readthedocs.io) | must be on PATH or passed via `--bedtools` |

---

## Installation

```bash
git clone https://github.com/hengbingao/methylave.git
cd methylave

# Option A: install as a command (recommended)
pip install -e .
methylave --help

# Option B: run directly without installing
python methylave_bin summary --help
```

After `pip install -e .` you can run `methylave` from anywhere.

---

## Input formats

### allc file (tab-separated, 1-based position)

```
chr1    3000827    +    CGT    3    10    1
chr1    3001006    +    CHH    0     5    0
chr1    3001011    +    CGA    8    10    1
```

| col | field | description |
|-----|-------|-------------|
| 1 | chrom | chromosome |
| 2 | pos | 1-based position |
| 3 | strand | +/- |
| 4 | mc_class | methylation context (CG*, CH*, …) |
| 5 | mc_count | methylated reads |
| 6 | total_count | total reads covering this site |
| 7 | methylated | 0/1 flag (optional) |

Supports plain text or `.gz`.

### BED file (tab-separated, 0-based half-open)

```
chr1    3000000    3001000
chr1    3005000    3010000
```

At minimum 3 columns (chrom, start, end). Extra columns are ignored.

---

## Modes & methods

### Method 1 — coverage-weighted (default)

$$\text{mean} = \frac{\sum mc\_count}{\sum total\_count}$$

Recommended for most use cases. High-coverage sites contribute more weight.

### Method 2 — equal-weight per CpG site

$$\text{mean} = \frac{1}{N} \sum_{i=1}^{N} \frac{mc\_count_i}{total\_count_i}$$

Every CpG site contributes equally regardless of coverage.

---

## Usage

### `summary` — one row per BED file

```bash
methylave summary \
    --allc  sample.allc.gz \
    --bed   E1.bed E2.bed E3.bed \
    --method method1 \
    --threads 4 \
    --output summary_results.tsv
```

**Output columns**

| column | description |
|--------|-------------|
| `bed_file` | input BED filename |
| `num_regions` | number of regions in the BED file |
| `total_region_length_bp` | sum of all region lengths (bp) |
| `num_CpG_sites` | CpG sites with coverage ≥ `--min_coverage` |
| `total_mc` | sum of methylated reads |
| `total_coverage` | sum of total reads |
| `mean_methylation` | aggregated methylation level |
| `method` | method used |

---

### `region` — one row per region

```bash
methylave region \
    --allc      sample.allc.gz \
    --bed       E1.bed E2.bed \
    --method    method2 \
    --output_dir region_results/ \
    --threads   4
```

Each BED file produces its own output TSV:
```
region_results/
  E1_region_methylation.tsv
  E2_region_methylation.tsv
```

```bash
methylave region \
    --allc      sample.allc.gz \
    --bed       E1.bed \
    --method    method2 \
    --suffix E1_sample_CG_methylation_method2.tsv \
    --output_dir region_results/ \
    --threads   4
```

**Output columns**

| column | description |
|--------|-------------|
| `chrom` | chromosome |
| `start` | region start (0-based) |
| `end` | region end |
| `num_CpG_sites` | covered CpG sites in this region |
| `total_mc` | methylated reads |
| `total_coverage` | total reads |
| `mean_methylation` | methylation level (`NA` if no coverage) |
| `method` | method used |

---

## All options

```
usage: methylave {summary,region} [options]

Modes:
  summary             One row per BED file
  region              One row per region

Shared options:
  --allc FILE         allc or allc.gz file (required)
  --bed FILE [...]    One or more BED files (required)
  --method            method1 (default) or method2
  --context STR       Methylation context prefix [default: CG]
  --bedtools PATH     Path to bedtools binary
  --threads N         Parallel threads [default: 4]
  --min_coverage INT  Minimum coverage per CpG site [default: 1]

summary options:
  --output FILE       Output TSV [default: summary_methylation.tsv]

region options:
  --output_dir DIR    Output directory [default: .]
  --suffix STR        Output filename suffix
                      [default: _region_methylation.tsv]
```

---

## Examples

```bash
# All 11 BEDs, summary, method1
methylave summary \
    --allc data/sample.allc.gz \
    --bed  beds/E{1..11}.bed \
    --threads 11 \
    --output results/all_summary.tsv

# Per-region, method2, filter low coverage
methylave region \
    --allc         data/sample.allc.gz \
    --bed          beds/E1.bed beds/E2.bed \
    --method       method2 \
    --min_coverage 5 \
    --output_dir   results/regions/

# Custom bedtools path
methylave summary \
    --allc     data/sample.allc.gz \
    --bed      beds/E1.bed \
    --bedtools /home/hgao/apps/bedtools2-2.25.0/bin/bedtools \
    --output   results/E1_summary.tsv
```

---

## Project structure

```
methylave/
├── methylave_bin           # executable entry point
├── methylave/
│   ├── __init__.py
│   ├── utils.py            # logging, tool resolution
│   ├── allc.py             # allc reading & sorting
│   ├── intersect.py        # bedtools intersect + aggregation
│   ├── summary.py          # summary mode
│   └── region.py           # region mode
├── setup.py                # pip install support
├── README.md
└── LICENSE
```
