"""
summary mode
------------
One output row per BED file.

Output columns:
    bed_file  num_regions  total_region_length_bp
    num_CpG_sites  total_mc  total_coverage
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
from methylave.intersect import bedtools_intersect, aggregate_method1, aggregate_method2


# ── per-BED worker ────────────────────────────────────────────────────────────

def _process_one_bed(bedtools_bin, allc_sorted, bed_path, sorted_bed,
                     method, min_coverage):
    """Worker function: intersect + aggregate for a single BED file."""
    t0 = time.time()
    label = Path(bed_path).name
    log("  [{}] intersecting ...".format(label))

    agg_fn = aggregate_method1 if method == "method1" else aggregate_method2

    # collect all hits, then aggregate
    hits = list(bedtools_intersect(bedtools_bin, allc_sorted, sorted_bed))

    total_mc, total_cov, mean_meth, num_sites = agg_fn(
        ((c, s, e, mc, cov) for c, s, e, mc, cov, *_ in hits),
        min_coverage=min_coverage,
    )

    # region count and total length from BED
    num_regions = 0
    total_length = 0
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            num_regions  += 1
            total_length += int(cols[2]) - int(cols[1])

    elapsed = time.time() - t0
    log("  [{}] done ({:.1f}s): {} regions, {} CpG sites, mean={:.4f}".format(
        label, elapsed,
        num_regions, num_sites,
        mean_meth if num_sites > 0 else float("nan")))

    return {
        "bed_file":            label,
        "num_regions":         num_regions,
        "total_region_length_bp": total_length,
        "num_CpG_sites":       num_sites,
        "total_mc":            total_mc,
        "total_coverage":      total_cov,
        "mean_methylation":    "{:.6f}".format(mean_meth) if num_sites > 0 else "NA",
        "method":              method,
    }


# ── public entry point ────────────────────────────────────────────────────────

def run_summary(args, bedtools_bin, bed_paths):
    t_start = time.time()

    with tempfile.TemporaryDirectory() as tmp:
        # Step 1: allc → context BED (once)
        log("=== Step 1/3: Extract {} sites from allc ===".format(args.context))
        raw_allc  = os.path.join(tmp, "allc_raw.bed")
        srt_allc  = os.path.join(tmp, "allc_sorted.bed")
        n_sites   = make_allc_bed(args.allc, args.context, raw_allc)
        if n_sites == 0:
            log("WARNING: no {} sites found in allc — output will be empty.".format(
                args.context))

        # Step 2: sort allc BED + all region BEDs
        log("=== Step 2/3: Sort BED files ===")
        sort_bed(bedtools_bin, raw_allc, srt_allc, label="allc_{}".format(args.context))
        os.remove(raw_allc)

        sorted_beds = {}
        for bp in bed_paths:
            sb = os.path.join(tmp, bp.stem + "_sorted.bed")
            sort_bed(bedtools_bin, str(bp), sb, label=bp.name)
            sorted_beds[bp] = sb

        # Step 3: parallel intersect + aggregate
        log("=== Step 3/3: Summary intersect (threads={}) ===".format(args.threads))
        results_map = {}
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = {
                pool.submit(
                    _process_one_bed,
                    bedtools_bin, srt_allc, str(bp), sorted_beds[bp],
                    args.method, args.min_coverage
                ): bp
                for bp in bed_paths
            }
            done = 0
            for fut in as_completed(futures):
                bp = futures[fut]
                results_map[bp] = fut.result()
                done += 1
                log("Progress: {}/{} BED files done".format(done, len(bed_paths)))

    # Write output (preserve input order)
    header = [
        "bed_file", "num_regions", "total_region_length_bp",
        "num_CpG_sites", "total_mc", "total_coverage",
        "mean_methylation", "method",
    ]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for bp in bed_paths:
            r = results_map[bp]
            out.write("\t".join(str(r[h]) for h in header) + "\n")

    log("=== SUMMARY DONE in {:.1f}s  →  {} ===".format(
        time.time() - t_start, args.output))

    # preview
    print()
    with open(args.output) as f:
        for line in f:
            print(line, end="")
