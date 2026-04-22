#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

from methylave.utils   import log, resolve_bedtools
from methylave.summary import run_summary
from methylave.region  import run_region


# ── custom formatter: wider help column, preserves newlines ──────────────────

class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    def __init__(self, prog, **kwargs):
        kwargs.setdefault("max_help_position", 32)
        kwargs.setdefault("width", 90)
        super().__init__(prog, **kwargs)


# ── top-level description ─────────────────────────────────────────────────────

_TOP_DESC = """\
methylave — methylation level calculator from allc files over BED regions

Calculate average methylation levels for genomic regions defined in BED files,
using allc-format methylation data as input. Supports any methylation context
(CG, CH, CA, CHH, CHG, ...) and two aggregation methods.
"""

_TOP_EPILOG = """\
modes:
  summary   Aggregate all regions in each BED into one methylation value.
            One output row per BED file.

  region    Calculate methylation for each region individually.
            One output row per region; each BED file gets its own output file.

methods:
  method1   coverage-weighted  :  mean = sum(mc) / sum(coverage)
            Recommended for most cases. High-coverage sites carry more weight.

  method2   equal-weight       :  mean = (1/N) * sum(mc_i / cov_i)
            Each site contributes equally regardless of coverage.

context (--context):
  Uses prefix matching against the mc_class column of the allc file.
    CG   → matches CGA, CGT, CGC, CGG  (default)
    CH   → matches CHH, CHG
    CHH  → matches CHH only
    CHG  → matches CHG only
    CA   → matches CAA, CAT, CAC, CAG
    C    → matches all cytosines

examples:
  # Summary, multiple BEDs, default method
  methylave summary --allc sample.allc.gz --bed E1.bed E2.bed --output out.tsv

  # Summary, method2, filter low-coverage sites
  methylave summary --allc sample.allc.gz --bed E1.bed \\
      --method method2 --min_coverage 5 --output out.tsv

  # Region mode, CH context, custom output directory and suffix
  methylave region --allc sample.allc.gz --bed E1.bed E2.bed \\
      --context CH --output_dir results/ --suffix _mCH.tsv

  # All 11 BEDs in parallel, custom bedtools path
  methylave summary --allc sample.allc.gz --bed E{1..11}.bed \\
      --threads 11 --bedtools /path/to/bedtools --output summary.tsv
"""

# ── subcommand descriptions ───────────────────────────────────────────────────

_SUMMARY_DESC = """\
summary mode — one output row per BED file

Aggregates methylation across ALL regions within each BED file into a single
value. Useful for comparing overall methylation levels between different sets
of genomic regions.

output columns:
  bed_file               input BED filename
  num_regions            number of regions in the BED file
  total_region_length_bp total length of all regions combined (bp)
  num_CpG_sites          sites with coverage >= --min_coverage
  total_mc               sum of methylated reads
  total_coverage         sum of total reads
  mean_methylation       aggregated methylation level (or NA)
  method                 method1 or method2
"""

_SUMMARY_EPILOG = """\
examples:
  methylave summary --allc sample.allc.gz --bed E1.bed --output E1_summary.tsv

  methylave summary --allc sample.allc.gz --bed E1.bed E2.bed E3.bed \\
      --method method2 --threads 3 --output all_summary.tsv

  methylave summary --allc sample.allc.gz --bed E1.bed \\
      --context CH --min_coverage 3 --output E1_CH.tsv
"""

_REGION_DESC = """\
region mode — one output row per region

Calculates methylation for each region individually. Each input BED file
produces its own output TSV. Regions with no covered sites are still reported
with zero counts and NA methylation.

output columns:
  chrom            chromosome
  start            region start (0-based)
  end              region end
  num_CpG_sites    covered sites in this region (>= --min_coverage)
  total_mc         methylated reads in this region
  total_coverage   total reads in this region
  mean_methylation methylation level (NA if no covered site)
  method           method1 or method2

output filenames:
  Each BED file produces: <bed_basename><suffix>
  e.g.  E1.bed  +  --suffix _mCG.tsv  →  E1_mCG.tsv
"""

_REGION_EPILOG = """\
examples:
  methylave region --allc sample.allc.gz --bed E1.bed --output_dir results/

  methylave region --allc sample.allc.gz --bed E1.bed E2.bed \\
      --method method2 --suffix _mCG.tsv --output_dir results/

  methylave region --allc sample.allc.gz --bed E1.bed \\
      --context CHH --min_coverage 5 --output_dir results/
"""


# ── parser builder ────────────────────────────────────────────────────────────

def build_parser():

    # ── top-level parser ──
    p = argparse.ArgumentParser(
        prog="methylave",
        description=_TOP_DESC,
        epilog=_TOP_EPILOG,
        formatter_class=HelpFormatter,
    )
    p.add_argument("--version", action="version", version="methylave 1.0.0")

    # ── shared arguments (inherited by both subcommands) ──
    shared = argparse.ArgumentParser(add_help=False)

    shared.add_argument(
        "--allc",
        required=True,
        metavar="FILE",
        help=(
            "Path to allc or allc.gz file.\n"
            "Tab-separated columns (1-based position):\n"
            "  chrom  pos  strand  mc_class  mc_count  total_count  [methylated]"
        ),
    )
    shared.add_argument(
        "--bed",
        required=True,
        nargs="+",
        metavar="FILE",
        help=(
            "One or more BED files (space-separated).\n"
            "Required columns: chrom, start (0-based), end.\n"
            "Extra columns are ignored."
        ),
    )
    shared.add_argument(
        "--method",
        choices=["method1", "method2"],
        default="method1",
        metavar="{method1,method2}",
        help=(
            "Aggregation method [default: method1]\n"
            "  method1  coverage-weighted: sum(mc) / sum(coverage)\n"
            "  method2  equal-weight:      mean of per-site fractions"
        ),
    )
    shared.add_argument(
        "--context",
        default="CG",
        metavar="STR",
        help=(
            "Methylation context prefix [default: CG]\n"
            "Uses prefix matching on the mc_class column.\n"
            "  CG   → CGA/CGT/CGC/CGG\n"
            "  CH   → CHH/CHG\n"
            "  CHH  → CHH only\n"
            "  CA   → CA*\n"
            "  C    → all cytosines"
        ),
    )
    shared.add_argument(
        "--min_coverage",
        type=int,
        default=1,
        metavar="INT",
        help=(
            "Minimum read coverage per site [default: 1]\n"
            "Sites with coverage below this value are excluded."
        ),
    )
    shared.add_argument(
        "--threads",
        type=int,
        default=4,
        metavar="N",
        help=(
            "Number of parallel threads [default: 4]\n"
            "Multiple BED files are processed concurrently.\n"
            "Set to the number of BED files for maximum speed."
        ),
    )
    shared.add_argument(
        "--bedtools",
        default=None,
        metavar="PATH",
        help=(
            "Full path to bedtools binary [default: auto-detect]\n"
            "Auto-detected from PATH if not specified.\n"
            "e.g. /home/user/apps/bedtools2/bin/bedtools"
        ),
    )

    # ── subcommands ──
    sub = p.add_subparsers(dest="mode", required=True, title="modes")

    # summary
    sp = sub.add_parser(
        "summary",
        parents=[shared],
        help="One row per BED file — aggregate methylation across all regions",
        description=_SUMMARY_DESC,
        epilog=_SUMMARY_EPILOG,
        formatter_class=HelpFormatter,
    )
    sp.add_argument(
        "--output",
        default="summary_methylation.tsv",
        metavar="FILE",
        help="Output TSV file path [default: summary_methylation.tsv]",
    )

    # region
    rp = sub.add_parser(
        "region",
        parents=[shared],
        help="One row per region — per-region methylation breakdown",
        description=_REGION_DESC,
        epilog=_REGION_EPILOG,
        formatter_class=HelpFormatter,
    )
    rp.add_argument(
        "--output_dir",
        default=".",
        metavar="DIR",
        help=(
            "Directory for output files [default: current directory]\n"
            "Created automatically if it does not exist."
        ),
    )
    rp.add_argument(
        "--suffix",
        default="_region_methylation.tsv",
        metavar="STR",
        help=(
            "Suffix for output filenames [default: _region_methylation.tsv]\n"
            "Output name = <bed_basename><suffix>\n"
            "e.g. E1.bed + '_mCG.tsv' → E1_mCG.tsv"
        ),
    )

    return p


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    parser = build_parser()
    args   = parser.parse_args()

    bedtools_bin = resolve_bedtools(args.bedtools)
    log("methylave  |  mode={}  method={}  context={}  threads={}".format(
        args.mode, args.method, args.context, args.threads))
    log("bedtools: {}".format(bedtools_bin))

    bed_paths = []
    for b in args.bed:
        p = Path(b)
        if not p.exists():
            log("ERROR: BED file not found: {}".format(b))
            sys.exit(1)
        bed_paths.append(p)
    log("BED files ({}): {}".format(
        len(bed_paths), ", ".join(str(p) for p in bed_paths)))

    if args.mode == "summary":
        run_summary(args, bedtools_bin, bed_paths)
    else:
        run_region(args, bedtools_bin, bed_paths)


if __name__ == "__main__":
    main()
