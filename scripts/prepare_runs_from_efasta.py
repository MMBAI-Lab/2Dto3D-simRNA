#!/usr/bin/env python3
"""
Split an extended-FASTA (`.efasta`) file into one SimRNA input directory per
record, ready for `scripts/run_simrna_set.sh`.

The extended-FASTA format used across this repo is a 3-line record:

    > free-text description
    ACGUACGU...
    ((((....))))...

Records sharing a sequence but carrying different dot-brackets are normal —
that is how alternative 2D restraints on the same aptamer are compared (see
the APT-PF1 reference set).

For each record this writes `<out_dir>/<run_name>/`:

    example        the sequence (SimRNA `-s`)
    example.str    the dot-bracket restraint (SimRNA `-S`)
    config.in      copy of the set-level config
    _source.txt    the original header, for traceability

Every record is validated before anything is written: sequence and structure
must be the same length, the sequence must be pure ACGU (SimRNA is an RNA
engine — transcribe DNA to RNA upstream), and the brackets must balance. A
single bad record aborts the whole run rather than half-staging a set.

Usage:
    python3 scripts/prepare_runs_from_efasta.py inputs/comercialApt/comercialApt_structures.fa \\
        --config inputs/comercialApt/config.in \\
        --names AN6_consensus AN6_IPknot M_apt_16_consensus MC_J3T2_consensus
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

VALID_BASES = set("ACGU")
RUN_SEQ_FILE = "example"      # SimRNA is invoked as `-s example -S example.str`
RUN_STR_FILE = "example.str"


class Record:
    def __init__(self, header: str, seq: str, ss: str):
        self.header = header
        self.seq = seq
        self.ss = ss

    @property
    def ident(self) -> str:
        """The token before the first dash/em-dash in the header."""
        return re.split(r"\s+[—–-]\s+", self.header, maxsplit=1)[0].strip()


def parse_efasta(path: Path) -> list[Record]:
    records: list[Record] = []
    header: str | None = None
    body: list[str] = []

    def flush(lineno: int) -> None:
        if header is None:
            return
        if len(body) != 2:
            die(f"{path}:{lineno}: record {header!r} has {len(body)} data line(s), expected 2 "
                f"(sequence then dot-bracket)")
        records.append(Record(header, body[0], body[1]))

    with path.open() as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush(lineno)
                header = line[1:].strip()
                body = []
            elif header is None:
                die(f"{path}:{lineno}: data line before any '>' header")
            else:
                body.append(line)
    flush(lineno)

    if not records:
        die(f"{path}: no records found")
    return records


def validate(records: list[Record]) -> None:
    problems: list[str] = []
    for i, r in enumerate(records, start=1):
        tag = f"record {i} ({r.ident!r})"
        if len(r.seq) != len(r.ss):
            problems.append(f"{tag}: sequence is {len(r.seq)} nt but dot-bracket is {len(r.ss)}")
        bad = sorted(set(r.seq) - VALID_BASES)
        if bad:
            problems.append(f"{tag}: non-ACGU character(s) {bad} in sequence — "
                            f"SimRNA needs RNA; transcribe T→U upstream")
        depth = 0
        for c in r.ss:
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth < 0:
                    problems.append(f"{tag}: unbalanced dot-bracket (closes before it opens)")
                    break
            elif c != ".":
                problems.append(f"{tag}: unsupported dot-bracket character {c!r} — "
                                f"SimRNA's -S takes plain ()/. only, no pseudoknot brackets")
                break
        else:
            if depth != 0:
                problems.append(f"{tag}: unbalanced dot-bracket ({depth} unclosed)")
    if problems:
        die("input validation failed:\n  - " + "\n  - ".join(problems))


def slugify(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return slug or "run"


def resolve_names(records: list[Record], explicit: list[str] | None) -> list[str]:
    if explicit:
        if len(explicit) != len(records):
            die(f"--names got {len(explicit)} name(s) for {len(records)} record(s)")
        names = [slugify(n) for n in explicit]
    else:
        names = [slugify(r.ident) for r in records]
        # Disambiguate repeated identifiers (same aptamer, alternative 2D).
        seen: dict[str, int] = {}
        for i, n in enumerate(names):
            if names.count(n) > 1:
                seen[n] = seen.get(n, 0) + 1
                names[i] = f"{n}_{seen[n]:02d}"
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        die(f"run names are not unique: {sorted(dupes)}")
    return names


def die(msg: str) -> None:
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("efasta", type=Path, help="Extended-FASTA file (header / sequence / dot-bracket)")
    ap.add_argument("--config", type=Path, required=True, help="Set-level SimRNA config.in to copy per run")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Where to write the per-run directories (default: the .efasta's own directory)")
    ap.add_argument("--names", nargs="*", default=None,
                    help="Explicit run names, one per record, in file order. "
                         "Defaults to the header identifier, de-duplicated.")
    ap.add_argument("--force", action="store_true",
                    help="Overwrite existing per-run input files")
    args = ap.parse_args()

    if not args.config.exists():
        die(f"config not found: {args.config}")

    records = parse_efasta(args.efasta)
    validate(records)
    names = resolve_names(records, args.names)
    out_dir = (args.out_dir or args.efasta.parent).resolve()

    for name, rec in zip(names, records):
        run_dir = out_dir / name
        seq_path = run_dir / RUN_SEQ_FILE
        if seq_path.exists() and not args.force:
            die(f"{seq_path} already exists — pass --force to overwrite")
        run_dir.mkdir(parents=True, exist_ok=True)
        seq_path.write_text(rec.seq + "\n")
        (run_dir / RUN_STR_FILE).write_text(rec.ss + "\n")
        (run_dir / "_source.txt").write_text(
            f"{rec.header}\n\nsource: {args.efasta}\nlength: {len(rec.seq)} nt\n"
            f"base pairs: {rec.ss.count('(')}\n"
        )
        shutil.copyfile(args.config, run_dir / "config.in")
        print(f"{name:<24} {len(rec.seq):>3} nt  {rec.ss.count('('):>2} bp  {rec.header}")

    print(f"\nstaged {len(records)} run(s) under {out_dir}")


if __name__ == "__main__":
    main()
