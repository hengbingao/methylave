"""
allc I/O helpers.

allc column layout (1-based pos, tab-separated):
  1  chrom
  2  pos          (1-based)
  3  strand
  4  mc_class     e.g. CGT, CHH …
  5  mc_count
  6  total_count
  7  methylated   (0/1 flag, optional)
"""

import gzip
import os
import subprocess
import sys
import time

from methylave.utils import log


def make_allc_bed(allc_path, context, out_path):
    """
    Stream-convert allc → 5-column BED, keeping only sites whose
    mc_class starts with *context* (e.g. 'CG').

    Output columns (0-based BED, tab-separated):
        chrom  start  end  mc_count  total_count
    """
    log("Reading allc: {}".format(allc_path))
    opener  = gzip.open if allc_path.endswith(".gz") else open
    count   = 0
    skipped = 0
    report  = 1_000_000

    with opener(allc_path, "rt") as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                skipped += 1
                continue
            if not parts[3].startswith(context):
                continue
            pos = int(parts[1])
            fout.write("{}\t{}\t{}\t{}\t{}\n".format(
                parts[0], pos - 1, pos, parts[4], parts[5]))
            count += 1
            if count % report == 0:
                log("  ... {:,} {} sites written".format(count, context))

    log("  Done: {:,} {} sites extracted  ({} lines skipped)".format(
        count, context, skipped))
    return count


def sort_bed(bedtools_bin, in_path, out_path, label=""):
    """Sort a BED file with bedtools sort (required for -sorted intersect)."""
    t0 = time.time()
    label = label or os.path.basename(in_path)
    log("  Sorting {} ...".format(label))
    with open(out_path, "w") as fout:
        result = subprocess.run(
            [bedtools_bin, "sort", "-i", in_path],
            stdout=fout,
            stderr=subprocess.PIPE,
        )
    if result.returncode != 0:
        log("  ERROR sorting {}: {}".format(label, result.stderr.decode().strip()))
        sys.exit(1)
    log("  Sorted  {} ({:.1f}s)".format(label, time.time() - t0))
