#!/usr/bin/env python3
"""
Assemble the AMBER-ready starting structures for a result set: pick the
populated substates, convert their centroids to tleap-compatible PDBs, and
write `tleap.in` plus the bilingual branded README.

Output goes to results/<set>/global_analysis/PDBtoMD/, mirroring the APT-PF1
layout that was built by hand.

Selection: every cluster whose population exceeds --min-pop (default 15 %) at
the run's reference cutoff. The cutoff is per-run because RMSD is an absolute
distance -- a set mixing chain lengths needs different cutoffs per run, and
those must match the ones its REPORT reports in detail, or the structures fed
to MD are not the substates the report describes.

Files are named `<run>_c<cutoff>_clust<NN>.pdb`. APT-PF1's hand-built directory
used short labels (`VFold2D_c12_clust01.pdb`) from a mapping table that only
existed for that set; the full run name needs no table and cannot collide.

Selection can also be given explicitly with --select, which overrides the
automatic criterion and may mix both kinds of substate:

    RUN:rmsd:<cutoff>:<NN>   centroid of RMSD cluster NN at that cutoff
    RUN:ss:<rank>            representative of the rank-th 2D pattern of the
                             low-T pool (lowest-energy frame carrying it)

The two are not equivalent. An RMSD cluster is a group of similar 3D structures
with a designated representative; a 2D-pattern group only shares base-pairing
and can span diverse folds, so its representative stands for the pairing, not
for a conformational cluster. Mixing them is a deliberate choice for sets where
RMSD clustering does not resolve every aptamer -- say so in the report.

Usage:
    python3 scripts/build_pdbtomd.py ApF20053R
    python3 scripts/build_pdbtomd.py comercialApt \\
        --select MC_J3T2_consensus:ss:1 --select AN6_consensus:rmsd:12:01
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_ss_representative import extract_frame, pick  # noqa: E402
from prepare_for_amber import prepare  # noqa: E402
from report_brand import LANGS, brand_footer, brand_header, out_path  # noqa: E402

BASE = "README"
DEFAULT_CUTOFF = 12.0
DEFAULT_MIN_POP = 15.0

TLEAP = """\
# Minimal AMBER tleap input to solvate + neutralize one of the {n} prepared PDBs.
# Swap the loadpdb line for each starting structure, then run:
#   tleap -f tleap.in
#
# Requires AmberTools >= 22 (for leaprc.RNA.OL3 / leaprc.water.tip3p).

source leaprc.RNA.OL3
source leaprc.water.tip3p

# === starting structure (edit for each of the {n}) ===
{loadpdb}

# quick sanity check
check rna

# Neutralize with Na+ counterions, solvate in a 10 A TIP3P octahedron.
addions rna Na+ 0
solvateOct rna TIP3PBOX 10.0

# Save topology + inpcrd ready for sander/pmemd.
saveamberparm rna system.prmtop system.inpcrd

quit
"""

T = {
    "en": {
        "title": "PDBtoMD — AMBER-ready starting structures",
        "intro": "{n} centroids selected from the RMSD clustering of the `{set}` set, converted to "
                 "AMBER-compatible PDBs (`leaprc.RNA.OL3`).",
        "hdr": "| file | substate | source | nt | population |",
        "crit": "**Selection criterion**: substates with population > {pop:g} % of the low-T pool "
                "({frames} frames per run) at each run's reference cutoff.",
        "crit_note": "Cutoffs used: {cuts}. RMSD is an absolute distance, so a set spanning several "
                     "chain lengths needs a cutoff per run — the same ones its `REPORT_EN.md` reports "
                     "in detail.",
        "did": "## What `prepare_for_amber.py` did to the SimRNA AA output",
        "d1": "- Renamed residues A/C/G/U → RA/RC/RG/RU; first residue → `RA5`/`RC5`/`RG5`/`RU5`, "
              "last residue → `RA3`/`RC3`/`RG3`/`RU3`.",
        "d2": "- Stripped the 5'-terminal phosphate (`P`, `OP1`, `OP2`, `OP3`) from residue 1 so "
              "tleap caps with `H5T`.",
        "d3": "- Renumbered atom serials 1..N and appended `TER` + `END`.",
        "run": "## Running tleap",
        "run_text": "The bundled [tleap.in](tleap.in) loads one PDB, neutralizes with Na⁺, solvates in a "
                    "10 Å TIP3P truncated octahedron, and saves `system.prmtop` / `system.inpcrd`. Edit "
                    "the `loadpdb` line to switch starting structure, or loop them in a shell wrapper.",
        "align": "Centroids of the same aptamer share a coordinate frame (aligned by C1' upstream), but "
                 "alignment has no effect on MD — you can disregard it for production runs.",
        "mixed": "This selection mixes two kinds of substate, which are **not** equivalent. An *RMSD* "
                 "entry is the representative of a 3D cluster: a group of structurally similar frames. "
                 "A *2D pattern* entry is the lowest-energy frame carrying one exact dot-bracket across "
                 "the low-T pool — it stands for the base-pairing, and the frames sharing that pairing "
                 "may still span diverse folds. 2D-pattern substates were used where RMSD clustering "
                 "did not resolve the aptamer: on the shortest chains a single cluster absorbs almost "
                 "the whole pool, so the conformational signal lives in the base-pairing instead.",
        "target": "**Next step**: explicit-solvent molecular dynamics (AMBER). Each structure below is "
                  "an independent starting point.",
    },
    "es": {
        "title": "PDBtoMD — estructuras de partida listas para AMBER",
        "intro": "{n} centroides seleccionados del clustering RMSD del set `{set}`, convertidos a PDBs "
                 "compatibles con AMBER (`leaprc.RNA.OL3`).",
        "hdr": "| archivo | subestado | origen | nt | población |",
        "crit": "**Criterio de selección**: subestados con población > {pop:g} % del pool low-T "
                "({frames} frames por run) al cutoff de referencia de cada corrida.",
        "crit_note": "Cutoffs usados: {cuts}. El RMSD es una distancia absoluta, así que un set que "
                     "abarca varias longitudes de cadena necesita un cutoff por corrida — los mismos "
                     "que su `REPORT_ES.md` reporta en detalle.",
        "did": "## Qué le hizo `prepare_for_amber.py` a la salida AA de SimRNA",
        "d1": "- Renombró los residuos A/C/G/U → RA/RC/RG/RU; primer residuo → `RA5`/`RC5`/`RG5`/`RU5`, "
              "último residuo → `RA3`/`RC3`/`RG3`/`RU3`.",
        "d2": "- Removió el fosfato 5' terminal (`P`, `OP1`, `OP2`, `OP3`) del residuo 1 para que tleap "
              "cierre con `H5T`.",
        "d3": "- Renumeró los atom serials 1..N y agregó `TER` + `END` al final.",
        "run": "## Correr tleap",
        "run_text": "El [tleap.in](tleap.in) incluido carga un PDB, neutraliza con Na⁺, solvata en un "
                    "octaedro truncado TIP3P de 10 Å y guarda `system.prmtop` / `system.inpcrd`. Editar "
                    "la línea `loadpdb` para cambiar de estructura, o iterarlas con un wrapper de shell.",
        "align": "Los centroides de un mismo aptámero comparten marco de coordenadas (alineados por C1' "
                 "aguas arriba), pero el alineamiento no afecta a la MD — es descartable en producción.",
        "mixed": "Esta selección mezcla dos tipos de subestado, que **no** son equivalentes. Una entrada "
                 "*RMSD* es el representante de un cluster 3D: un grupo de frames estructuralmente "
                 "similares. Una entrada *patrón 2D* es el frame de menor energía que porta un "
                 "dot-bracket exacto en el pool low-T — representa el apareamiento de bases, y los "
                 "frames que lo comparten pueden abarcar plegamientos diversos. Se usaron subestados de "
                 "patrón 2D donde el clustering RMSD no resolvía el aptámero: en las cadenas más cortas "
                 "un único cluster absorbe casi todo el pool, así que la señal conformacional está en el "
                 "apareamiento de bases.",
        "target": "**Paso siguiente**: dinámica molecular en solvente explícito (AMBER). Cada estructura "
                  "de abajo es un punto de partida independiente.",
    },
}


def read_centroids(tsv: Path) -> list[dict]:
    rows = []
    if not tsv.exists():
        return rows
    for line in tsv.read_text().splitlines()[1:]:
        p = line.split("\t")
        if len(p) >= 4:
            rows.append({"rank": int(p[0]), "cluster": p[1],
                         "frames": int(p[2]), "pct": float(p[3])})
    return rows


def first_line(path: Path) -> str:
    if not path.exists():
        return ""
    for ln in path.read_text().splitlines():
        if ln.strip():
            return ln.strip()
    return ""


def fmt(c: float) -> str:
    return f"{c:g}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("set_name")
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--cutoff", action="append", default=[], metavar="RUN=CUTOFF",
                    help="Reference cutoff in A for one run; repeatable.")
    ap.add_argument("--default-cutoff", type=float, default=DEFAULT_CUTOFF)
    ap.add_argument("--min-pop", type=float, default=DEFAULT_MIN_POP)
    ap.add_argument("--select", action="append", default=[],
                    metavar="RUN:rmsd:CUTOFF:NN | RUN:ss:RANK",
                    help="Explicit substate to include; repeatable. Overrides --min-pop.")
    args = ap.parse_args()

    set_root = args.repo / "results" / args.set_name
    ga = set_root / "global_analysis"
    if not ga.is_dir():
        print(f"no global_analysis for {args.set_name}: {ga}", file=sys.stderr)
        return 1

    overrides: dict[str, float] = {}
    for spec in args.cutoff:
        if "=" not in spec:
            print(f"--cutoff expects RUN=CUTOFF, got {spec!r}", file=sys.stderr)
            return 1
        k, v = spec.split("=", 1)
        overrides[k] = float(v)

    out_dir = ga / "PDBtoMD"
    out_dir.mkdir(parents=True, exist_ok=True)

    selected: list[dict] = []
    cutoffs_used: set[float] = set()
    simrna_dir = args.repo / "SimRNA_64bitIntel_Linux"

    if args.select:
        for spec in args.select:
            bits = spec.split(":")
            run_name, kind = bits[0], bits[1] if len(bits) > 1 else ""
            run_res = set_root / run_name
            if not run_res.is_dir():
                print(f"--select: no such run {run_name}", file=sys.stderr)
                return 1
            seq = first_line(run_res / "example")

            if kind == "rmsd" and len(bits) == 4:
                cut, cid = float(bits[2]), bits[3].zfill(2)
                src = ga / run_name / f"cutoff_{int(round(cut))}A" / f"clust{cid}.pdb"
                if not src.exists():
                    print(f"--select: missing {src}", file=sys.stderr)
                    return 1
                pct = next((r["pct"] for r in read_centroids(
                    ga / run_name / f"cutoff_{int(round(cut))}A" / "_centroids.tsv")
                    if r["cluster"] == cid), 0.0)
                name = f"{run_name}_c{int(round(cut))}_clust{cid}.pdb"
                prepare(src, out_dir / name)
                cutoffs_used.add(cut)
                selected.append({"file": name, "run": run_name,
                                 "origin": f"{run_name}/cutoff_{int(round(cut))}A/clust{cid}",
                                 "kind": "rmsd", "detail": f"RMSD {fmt(cut)} Å",
                                 "pct": pct, "nt": len(seq), "ss": ""})

            elif kind == "ss" and len(bits) == 3:
                rank = int(bits[2])
                pattern, pct, frame, energy = pick(run_res, rank)
                raw = out_dir / f".raw_{run_name}_ss{rank:02d}.pdb"
                extract_frame(run_res, frame, raw, simrna_dir)
                name = f"{run_name}_ss{rank:02d}.pdb"
                prepare(raw, out_dir / name)
                raw.unlink(missing_ok=True)
                selected.append({"file": name, "run": run_name,
                                 "origin": f"{run_name} 2D pattern #{rank} (frame {frame}, E={energy:.1f})",
                                 "kind": "ss", "detail": f"2D pattern #{rank}",
                                 "pct": pct, "nt": len(seq), "ss": pattern})
            else:
                print(f"--select: cannot parse {spec!r}", file=sys.stderr)
                return 1
            print(f"{selected[-1]['file']:<44} {selected[-1]['pct']:5.1f} %  "
                  f"{selected[-1]['detail']}")
        return finish(args, out_dir, selected, cutoffs_used)

    for run in sorted(p for p in ga.iterdir() if p.is_dir() and p.name != "PDBtoMD"):
        cut = overrides.get(run.name, args.default_cutoff)
        tsv = run / f"cutoff_{int(round(cut))}A" / "_centroids.tsv"
        rows = read_centroids(tsv)
        if not rows:
            print(f"warning: no centroids for {run.name} at {fmt(cut)} A ({tsv})", file=sys.stderr)
            continue
        cutoffs_used.add(cut)
        seq = first_line(set_root / run.name / "example")
        for r in rows:
            if r["pct"] <= args.min_pop:
                continue
            src = run / f"cutoff_{int(round(cut))}A" / f"clust{r['cluster']}.pdb"
            if not src.exists():
                print(f"warning: missing centroid PDB {src}", file=sys.stderr)
                continue
            name = f"{run.name}_c{int(round(cut))}_clust{r['cluster']}.pdb"
            prepare(src, out_dir / name)
            selected.append({"file": name, "run": run.name,
                             "origin": f"{run.name}/cutoff_{int(round(cut))}A/clust{r['cluster']}",
                             "kind": "rmsd", "detail": f"RMSD {fmt(cut)} Å",
                             "pct": r["pct"], "nt": len(seq), "ss": ""})
            print(f"{name:<48} {r['pct']:5.1f} %")

    return finish(args, out_dir, selected, cutoffs_used)


def finish(args, out_dir: Path, selected: list[dict], cutoffs_used: set[float]) -> int:
    if not selected:
        print("nothing selected", file=sys.stderr)
        return 1

    n = len(selected)
    load = "\n".join(
        (f"rna = loadpdb {s['file']}" if i == 0 else f"# rna = loadpdb {s['file']}")
        for i, s in enumerate(selected))
    (out_dir / "tleap.in").write_text(TLEAP.format(n=n, loadpdb=load))

    mixed = len({s["kind"] for s in selected}) > 1
    cuts_txt = ", ".join(f"{fmt(c)} Å" for c in sorted(cutoffs_used)) or "—"

    for lang in LANGS:
        t = T[lang]
        parts = brand_header(out_dir, lang, BASE)
        parts.append(f"# {t['title']}")
        parts.append("")
        parts.append(t["intro"].format(n=n, set=args.set_name))
        parts.append("")
        parts.append(t["target"])
        parts.append("")
        parts.append(t["hdr"])
        parts.append("|:---|:---|:---|---:|---:|")
        for s in selected:
            parts.append(f"| `{s['file']}` | {s['detail']} | `{s['origin']}` | "
                         f"{s['nt']} | {s['pct']:.1f} % |")
        parts.append("")
        if mixed:
            parts.append(t["mixed"])
            parts.append("")
        elif getattr(args, "select", None):
            pass
        else:
            parts.append(t["crit"].format(pop=args.min_pop, frames="30 003"))
            parts.append("")
            parts.append(t["crit_note"].format(cuts=cuts_txt))
            parts.append("")
        ss_rows = [s for s in selected if s["ss"]]
        if ss_rows:
            parts.append("| " + ("run" if lang == "en" else "run")
                         + " | 2D (dot-bracket) |")
            parts.append("|:---|:---|")
            for s in ss_rows:
                parts.append(f"| `{s['file']}` | `{s['ss']}` |")
            parts.append("")
        parts.append(t["did"])
        parts.append("")
        parts += [t["d1"], t["d2"], t["d3"]]
        parts.append("")
        parts.append(t["run"])
        parts.append("")
        parts.append("```")
        parts.append("tleap -f tleap.in")
        parts.append("```")
        parts.append("")
        parts.append(t["run_text"])
        parts.append("")
        parts.append(t["align"])
        parts.append("")
        parts += brand_footer(lang)
        p = out_path(out_dir, lang, BASE)
        p.write_text("\n".join(parts))
        print(f"wrote {p}")

    print(f"wrote {out_dir}/tleap.in  ({n} structures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
