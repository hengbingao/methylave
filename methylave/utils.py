"""Shared utilities: logging, tool resolution."""

import os
import shutil
import sys
import time


def log(msg):
    """Timestamped log line, flushed immediately."""
    ts = time.strftime("%H:%M:%S")
    print("[{}] {}".format(ts, msg), flush=True)


def resolve_bedtools(user_path):
    """Return a usable bedtools binary path or exit with a clear error."""
    if user_path:
        if os.path.isfile(user_path):
            return user_path
        log("ERROR: specified bedtools path does not exist: {}".format(user_path))
        sys.exit(1)

    found = shutil.which("bedtools")
    if found:
        return found

    # common HPC / conda locations
    hints = [
        os.path.expanduser("~/apps/bedtools2-2.25.0/bin/bedtools"),
        os.path.expanduser("~/miniconda3/bin/bedtools"),
        os.path.expanduser("~/anaconda3/bin/bedtools"),
        "/home/hgao/apps/bedtools2-2.25.0/bin/bedtools",
    ]
    for h in hints:
        if os.path.isfile(h):
            return h

    log("ERROR: bedtools not found. Install it or pass --bedtools /path/to/bedtools")
    sys.exit(1)
