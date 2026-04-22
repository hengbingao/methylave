"""
Core bedtools intersect wrapper.

Returns raw hits as a list of (mc_count, total_count) per overlapping
allc site.  The caller decides how to aggregate.
"""

import subprocess
import sys

from methylave.utils import log


def bedtools_intersect(bedtools_bin, allc_sorted, bed_sorted, extra_flags=None):
    """
    Run bedtools intersect -a allc_sorted -b bed_sorted -sorted -wo.

    Yields tuples:
        (chrom, start, end,          <- allc site (BED coords, 0-based)
         mc_count, total_count,      <- from allc
         reg_chrom, reg_start, reg_end)   <- from region BED
    """
    cmd = [
        bedtools_bin, "intersect",
        "-a", allc_sorted,
        "-b", bed_sorted,
        "-sorted",
        "-wo",          # write both A and B entries + overlap length
    ]
    if extra_flags:
        cmd.extend(extra_flags)

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        log("ERROR bedtools intersect: {}".format(proc.stderr.decode().strip()))
        sys.exit(1)

    for line in proc.stdout.decode("utf-8").splitlines():
        parts = line.split("\t")
        # allc BED:  chrom start end mc total  (5 cols)
        # region BED: chrom start end [extra…] (≥3 cols)
        # last col: overlap length (bedtools -wo)
        if len(parts) < 9:
            continue
        chrom      = parts[0]
        start      = int(parts[1])
        end        = int(parts[2])
        mc_count   = int(parts[3])
        total_count= int(parts[4])
        reg_chrom  = parts[5]
        reg_start  = int(parts[6])
        reg_end    = int(parts[7])
        yield chrom, start, end, mc_count, total_count, reg_chrom, reg_start, reg_end


def aggregate_method1(hits, min_coverage=1):
    """
    Method 1: total_mc / total_coverage  (coverage-weighted).

    Returns (total_mc, total_cov, mean_meth, num_sites).
    """
    total_mc  = 0
    total_cov = 0
    num_sites = 0
    for _, _s, _e, mc, cov, *_ in hits:
        if cov < min_coverage:
            continue
        total_mc  += mc
        total_cov += cov
        num_sites += 1
    mean = (total_mc / total_cov) if total_cov > 0 else float("nan")
    return total_mc, total_cov, mean, num_sites


def aggregate_method2(hits, min_coverage=1):
    """
    Method 2: mean of per-CpG fractions (equal-weight per site).

    Returns (total_mc, total_cov, mean_meth, num_sites).
    total_mc / total_cov are kept for reference output.
    """
    frac_sum  = 0.0
    total_mc  = 0
    total_cov = 0
    num_sites = 0
    for _, _s, _e, mc, cov, *_ in hits:
        if cov < min_coverage:
            continue
        frac_sum  += mc / cov
        total_mc  += mc
        total_cov += cov
        num_sites += 1
    mean = (frac_sum / num_sites) if num_sites > 0 else float("nan")
    return total_mc, total_cov, mean, num_sites
