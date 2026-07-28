#!/usr/bin/env python3
"""
Cluster SimRNA frames by exact `.ss_detected` dot-bracket string.

Given a directory containing per-frame `*.ss_detected` files (produced by
`SimRNA_trafl2pdbs`), groups frames by the detected 2D pattern and emits a
TSV on stdout with: rank, ss, count, pct. Sorted by population descending.

Usage: ss_cluster.py <dir>
"""
import sys
from collections import Counter
from pathlib import Path

VALID_SS_CHARS = set(".()[]{}<>-&")


def read_ss(path: Path):
    with path.open() as f:
        for line in f:
            s = line.strip()
            if s and set(s) <= VALID_SS_CHARS:
                return s
    return None


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <dir>", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    files = sorted(root.glob("*.ss_detected"))
    if not files:
        print("no .ss_detected files found", file=sys.stderr)
        return 1

    counter: Counter = Counter()
    for path in files:
        ss = read_ss(path)
        if ss is not None:
            counter[ss] += 1

    total = sum(counter.values())
    print("rank\tss\tcount\tpct")
    for rank, (ss, count) in enumerate(counter.most_common(), start=1):
        pct = 100.0 * count / total if total else 0.0
        print(f"{rank}\t{ss}\t{count}\t{pct:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
