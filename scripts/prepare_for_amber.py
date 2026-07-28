#!/usr/bin/env python3
"""
Prepare a SimRNA all-atom PDB for AMBER MD (tleap-ready).

Transformations applied:
  1. Residues renamed from single-letter A/C/G/U to AMBER 3-letter RNA names
     RA/RC/RG/RU, with `5`/`3` suffix on the terminal residues.
  2. The 5'-terminal phosphate atoms (P/OP1/OP2/OP3) are stripped from
     residue 1 so tleap caps the 5' end with H5T on O5' per
     leaprc.RNA.OL3.
  3. Atom serial numbers renumbered 1..N.
  4. TER + END records appended (tleap expects them).

Usage: prepare_for_amber.py <in.pdb> <out.pdb>
"""
import argparse
import sys
from pathlib import Path

RNA_MAP = {"A": "RA", "C": "RC", "G": "RG", "U": "RU"}
STRIP_FIRST = {"P", "OP1", "OP2", "OP3"}


def write_atom(serial: int, atom: str, altloc: str, resn: str, chain: str,
               resnum: int, icode: str, rest: str) -> str:
    """Build an ATOM record obeying the PDB column layout.

    PDB fixed columns (1-based): 1-6 record, 7-11 serial, 13-16 atom name,
    17 altLoc, 18-20 resName, 22 chainID, 23-26 resSeq, 27 iCode, then the
    coordinates and occ/b-factor block (columns 31-66+) carried over verbatim.
    """
    atom_field = atom if len(atom) == 4 else f" {atom:<3}"
    return (
        f"ATOM  "                # 1-6
        f"{serial:>5d} "          # 7-11 + space
        f"{atom_field:<4}"        # 13-16
        f"{altloc or ' '}"        # 17
        f"{resn:<3}"              # 18-20
        f" "                      # 21
        f"{chain or ' '}"         # 22
        f"{resnum:>4d}"           # 23-26
        f"{icode or ' '}"         # 27
        f"   "                    # 28-30
        f"{rest}"                 # 31-66+ (coords, occ, b-factor, element, charge)
    )


def prepare(in_pdb: Path, out_pdb: Path) -> None:
    lines = in_pdb.read_text().splitlines()

    residues: list[int] = []
    seen = set()
    for L in lines:
        if L.startswith(("ATOM", "HETATM")):
            n = int(L[22:26])
            if n not in seen:
                seen.add(n)
                residues.append(n)
    if not residues:
        raise RuntimeError(f"no ATOM records found in {in_pdb}")
    first, last = residues[0], residues[-1]

    out: list[str] = []
    serial = 0
    for L in lines:
        if not L.startswith(("ATOM", "HETATM")):
            continue
        atom = L[12:16].strip()
        resn = L[17:20].strip()
        resnum = int(L[22:26])
        altloc = L[16:17]
        chain = L[21:22]
        icode = L[26:27]
        rest = L[30:]

        if resnum == first and atom in STRIP_FIRST:
            continue

        base = RNA_MAP.get(resn)
        if base is None:
            raise RuntimeError(f"unexpected residue '{resn}' at pos {resnum} in {in_pdb}")
        if resnum == first:
            new_resn = base + "5"
        elif resnum == last:
            new_resn = base + "3"
        else:
            new_resn = base

        serial += 1
        out.append(write_atom(serial, atom, altloc, new_resn, chain, resnum, icode, rest))

    out.append("TER")
    out.append("END")
    out_pdb.write_text("\n".join(out) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("in_pdb", type=Path)
    ap.add_argument("out_pdb", type=Path)
    args = ap.parse_args()
    prepare(args.in_pdb, args.out_pdb)
    print(f"wrote {args.out_pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
