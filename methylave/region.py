"""
region mode
-----------
One output row per region per BED file.
Each BED file produces its own output TSV in --output_dir.

Output columns per file:
    chrom  start  end  num_CpG_sites  total_mc  total_coverage
    mean_methylation  method
"""

import os
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from methylave.utils     import log
from methylave.allc      import make_allc_bed, sort_bed
from methylave.intersect import bedtools_intersect


# ── per-region aggregation helpers ───────────────────────────────────────────

def _region_key(reg_chrom, reg_start, reg_end):
    return (reg_chrom, reg_start, reg_end)


def _aggregate_regions_method1(hits, min_coverage):
    """Coverage-weighted: sum mc / sum cov per region."""
    mc_map  = defaultdict(int)
    cov_map = defaultdict(int)
    cnt_map = defaultdict(int)

    for _, _s, _e, mc, cov, rc, rs, re in hits:
        if cov < min_coverage:
            continue
        key = _region_key(rc, rs, re)
        mc_map[key]  += mc
        cov_map[key] += cov
        cnt_map[key] += 1

    rows = []
    for key in sorted(mc_map):
        chrom, start, end = key
        mc    = mc_map[key]
        cov   = cov_map[key]
        sites = cnt_map[key]
        mean  = mc / cov if cov > 0 else float("nan")
        rows.append((chrom, start, end, sites, mc, cov, mean))
    return rows


def _aggregate_regions_method2(hits, min_coverage):
    """Equal-weight per site: mean of fractions per region."""
    frac_map = defaultdict(float)
    mc_map   = defaultdict(int)
    cov_map  = defaultdict(int)
    cnt_map  = defaultdict(int)

    for _, _s, _e, mc, cov, rc, rs, re in hits:
        if cov < min_coverage:
            continue
        key = _region_key(rc, rs, re)
        frac_map[key] += mc / cov
        mc_map[key]   += mc
        cov_map[key]  += cov
        cnt_map[key]  += 1

    rows = []
    for key in sorted(frac_map):
        chrom, start, end = key
        sites = cnt_map[key]
        mc    = mc_map[key]
        cov   = cov_map[key]
        mean  = frac_map[key] / sites if sites > 0 else float("nan")
        rows.append((chrom, start, end, sites, mc, cov, mean))
    return rows


# ── zero-coverage regions (in BED but no CpG hits) ───────────────────────────

def _read_all_regions(bed_path):
    """Return sorted list of (chrom, start, end) from a BED file."""
    regions = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            regions.append((cols[0], int(cols[1]), int(cols[2])))
    return sorted(regions)


# ── per-BED worker ────────────────────────────────────────────────────────────

def _process_one_bed(bedtools_bin, allc_sorted, bed_path, sorted_bed,
                     method, min_coverage, output_path):
    t0    = time.time()
    label = Path(bed_path).name
    log("  [{}] intersecting (region mode) ...".format(label))

    all_hits = list(bedtools_intersect(bedtools_bin, allc_sorted, sorted_bed))

    if method == "method1":
        hit_rows = _aggregate_regions_method1(all_hits, min_coverage)
    else:
        hit_rows = _aggregate_regions_method2(all_hits, min_coverage)

    # build a dict keyed by region for fast lookup
    hit_dict = {(r[0], r[1], r[2]): r for r in hit_rows}

    # merge with full region list so zero-coverage regions also appear
    all_regions = _read_all_regions(bed_path)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    header = [
        "chrom", "start", "end",
        "num_CpG_sites", "total_mc", "total_coverage",
        "mean_methylation", "method",
    ]
    with open(output_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for reg in all_regions:
            if reg in hit_dict:
                chrom, start, end, sites, mc, cov, mean = hit_dict[reg]
                mean_str = "{:.6f}".format(mean) if sites > 0 else "NA"
            else:
                chrom, start, end = reg
                sites, mc, cov = 0, 0, 0
                mean_str = "NA"
            out.write("\t".join([
                chrom, str(start), str(end),
                str(sites), str(mc), str(cov),
                mean_str, method,
            ]) + "\n")

    elapsed = time.time() - t0
    log("  [{}] done ({:.1f}s): {} regions written → {}".format(
        label, elapsed, len(all_regions), output_path))
    return label, len(all_regions)


# ── public entry point ────────────────────────────────────────────────────────

def run_region(args, bedtools_bin, bed_paths):
    t_start = time.time()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        # Step 1
        log("=== Step 1/3: Extract {} sites from allc ===".format(args.context))
        raw_allc = os.path.join(tmp, "allc_raw.bed")
        srt_allc = os.path.join(tmp, "allc_sorted.bed")
        n_sites  = make_allc_bed(args.allc, args.context, raw_allc)
        if n_sites == 0:
            log("WARNING: no {} sites in allc.".format(args.context))

        # Step 2
        log("=== Step 2/3: Sort BED files ===")
        sort_bed(bedtools_bin, raw_allc, srt_allc,
                 label="allc_{}".format(args.context))
        os.remove(raw_allc)

        sorted_beds = {}
        for bp in bed_paths:
            sb = os.path.join(tmp, bp.stem + "_sorted.bed")
            sort_bed(bedtools_bin, str(bp), sb, label=bp.name)
            sorted_beds[bp] = sb

        # Step 3
        log("=== Step 3/3: Region intersect (threads={}) ===".format(args.threads))
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = {
                pool.submit(
                    _process_one_bed,
                    bedtools_bin, srt_allc, str(bp), sorted_beds[bp],
                    args.method, args.min_coverage,
                    str(out_dir / (bp.stem + args.suffix))
                ): bp
                for bp in bed_paths
            }
            done = 0
            for fut in as_completed(futures):
                fut.result()
                done += 1
                log("Progress: {}/{} BED files done".format(done, len(bed_paths)))

    log("=== REGION DONE in {:.1f}s  →  {} ===".format(
        time.time() - t_start, out_dir))
