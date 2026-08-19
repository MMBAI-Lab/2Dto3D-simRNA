#!/usr/bin/env python3
"""
Build the cross-run global report for a result set, as a bilingual branded pair
`REPORT_EN.md` + `REPORT_ES.md`.

This is the set-aware counterpart of `build_global_readme.py`, which is
hardwired to APT-PF1 (three runs over a single 40-nt sequence, one shared
alignment anchor). A generic set does not have that shape: `comercialApt` holds
four runs over *three* different sequences, so only runs sharing a sequence can
be structurally superposed or compared cluster-by-cluster.

Everything is derived from files already on disk — nothing is hardcoded per set:

  inputs/<set>/<run>/example        sequence
  inputs/<set>/<run>/example.str    2D restraint
  inputs/<set>/<run>/_source.txt    provenance header (written by
                                    prepare_runs_from_efasta.py); its first
                                    line's leading token groups runs by aptamer
  results/<set>/<run>/rmsd_clusters/          per-cutoff cluster populations
  results/<set>/<run>/ss_clusters/ss_clusters.tsv   2D-pattern populations
  results/<set>/global_analysis/<run>/cutoff_<N>A/  centroid PDB/PNG (optional)

The **reference cutoff** is the one whose clustering is reported in detail and
carried into the next stage. It is per-run because RMSD is an absolute distance:
a cutoff that resolves substates on a 60-nt chain puts every frame of a 13-nt
hairpin in one cluster. Pass it explicitly:

    python3 scripts/build_set_report.py comercialApt \\
        --ref-cutoff AN6_consensus=12 --ref-cutoff AN6_IPknot=12 \\
        --ref-cutoff M_apt_16_consensus=4 --ref-cutoff MC_J3T2_consensus=2

Runs without an explicit value fall back to --default-ref-cutoff.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_cluster_readme import count_frames, rmsd_table, ss_table  # noqa: E402
from report_brand import LANGS, brand_footer, brand_header, out_path  # noqa: E402

BASE = "REPORT"
TOP_N = 10


# ---------- localized strings ------------------------------------------------

T = {
    "en": {
        "title": "{set} — conformational substates (global analysis)",
        "intro": "Cross-run summary of the {nruns} SimRNA REMC simulations in the `{set}` set, "
                 "covering {napt} aptamer(s). Pipeline: `scripts/run_simrna_set.sh` to sample, "
                 "`scripts/run_postprocess_set.sh` to pool and cluster, "
                 "`scripts/extract_global_centroids.py` for all-atom centroids.",
        "b_pool": "- **Frame pool**: REMC levels 2–4 (T≈1 on the relative scale), {frames} frames per run.",
        "b_ladder": "- **Temperature ladder**: 16 replicas, {lo}–{hi} on the relative scale.",
        "b_clust": "- **RMSD clustering**: SimRNA's `clustering` binary, cutoffs {cutoffs}.",
        "b_back": "- **Coarse-grain → all-atom back-mapping**: `SimRNA_trafl2pdbs <ref> <cluster.trafl> 1 AA` "
                  "— frame 1 of each cluster is the representative (centroid).",
        "design": "## Experimental design",
        "design_note": "Runs sharing an aptamer differ only in the 2D restraint — that is the comparison "
                       "the set is built to make. Runs on different aptamers are independent and are "
                       "**not** structurally comparable to each other.",
        "d_hdr": "| aptamer | run | nt | bp | 2D restraint (dot-bracket) |",
        "refcut": "## Reference cutoff per run",
        "refcut_text": "RMSD is an absolute distance, so one cutoff does not serve every chain length: "
                       "a threshold that resolves substates on a 60-nt chain collapses a 13-nt hairpin "
                       "into a single cluster. The cutoff chosen per run is the one that resolves "
                       "conformational heterogeneity without either collapsing to a single mega-class or "
                       "fragmenting into thermal noise.",
        "r_hdr": "| run | nt | reference cutoff | top-1 population | verdict |",
        "sweep": "## Cutoff sweep — dominant-cluster population",
        "sweep_text": "How the top-1 cluster population responds to the cutoff. A run is *saturated* when "
                      "the value sits at ~100 % (the cutoff resolves nothing) and *fragmented* when the "
                      "top cluster is a few per cent (the cutoff is below the thermal fluctuation scale).",
        "ss_cmp": "## Reference 2D (restraint) vs. detected 2D of the top-1 cluster",
        "ss_cmp_text": "The centroid's detected 2D rarely matches the restraint exactly. The gap tells us "
                       "how hard SimRNA is pushing off the target 2D.",
        "s_hdr": "| run | 2D restraint | 2D of top-1 @ reference cutoff | match? |",
        "tables": "## Top-{n} clusters at the reference cutoff",
        "t_hdr": "| rank | cluster | frames | % | centroid 2D |",
        "sspat": "## Most populated 2D patterns per run",
        "sspat_text": "Population of each **exact dot-bracket pattern** across the low-T pool, independent "
                      "of 3D clustering. The restraint's own pattern is marked **=**.",
        "p_hdr": "| rank | % | 2D (dot-bracket) | |",
        "pair": "## Restraint comparison — {apt}",
        "pair_text": "These runs share a sequence, so their clusterings are directly comparable and their "
                     "centroids can be superposed. This is the head-to-head the set exists to answer.",
        "gallery": "## Gallery — top-1 centroids",
        "gallery_cap": "_All-atom centroid back-mapped from frame 1 of the dominant cluster at each run's "
                       "reference cutoff, rendered in ChimeraX (tube/slab muffler, rainbow by residue "
                       "5' blue → 3' red). Centroids of the same aptamer share an alignment anchor and "
                       "camera, so they can be compared directly; different aptamers cannot._",
        "no_gallery": "_(no centroid renders yet — run `scripts/extract_global_centroids.py --set {set}` "
                      "then `scripts/render_centroids_chimerax.py --set {set}`)_",
        "note": "_Cluster `.trafl`, `.pdb` and `.png` files live under `global_analysis/`. Only the REPORTs "
                "and the per-run `README_{{EN,ES}}.md` are versioned; trajectories and cluster data are gitignored._",
        "yes": "yes", "no": "no",
        "v_ok": "resolves", "v_sat": "saturated", "v_frag": "fragmented",
        "md": "## Structures selected for molecular dynamics",
        "md_text": "The substates below are the starting points for the next stage: **explicit-solvent "
                   "molecular dynamics** with AMBER. AMBER-ready PDBs and a starter `tleap.in` are in "
                   "[`global_analysis/PDBtoMD/`](global_analysis/PDBtoMD) — see its "
                   "[README](global_analysis/PDBtoMD/README_EN.md) for the conversion details.",
        "md_hdr": "| structure | aptamer | substate | population |",
        "md_mixed": "The selection deliberately mixes two kinds of substate. An **RMSD** entry represents "
                    "a 3D cluster of structurally similar frames. A **2D pattern** entry is the "
                    "lowest-energy frame carrying one exact dot-bracket — it represents the base-pairing, "
                    "not a conformational cluster. 2D-pattern substates were taken where RMSD clustering "
                    "does not resolve the aptamer: on the shortest chains one cluster absorbs nearly the "
                    "whole pool, so what varies is the base-pairing rather than the fold.",
        "per_run": "## Per-run detail",
        "per_run_line": "- [`{run}` →]({run}/README_EN.md) — full RMSD and 2D-pattern clustering, all cutoffs.",
    },
    "es": {
        "title": "{set} — subestados conformacionales (análisis global)",
        "intro": "Resumen cruzado de las {nruns} simulaciones SimRNA REMC del set `{set}`, "
                 "sobre {napt} aptámero(s). Pipeline: `scripts/run_simrna_set.sh` para muestrear, "
                 "`scripts/run_postprocess_set.sh` para poolear y clusterizar, "
                 "`scripts/extract_global_centroids.py` para los centroides all-atom.",
        "b_pool": "- **Pool de frames**: niveles REMC 2–4 (T≈1 en la escala relativa), {frames} frames por run.",
        "b_ladder": "- **Escalera de temperaturas**: 16 réplicas, {lo}–{hi} en la escala relativa.",
        "b_clust": "- **Clustering RMSD**: binario `clustering` de SimRNA, cutoffs {cutoffs}.",
        "b_back": "- **Back-mapping a all-atom**: `SimRNA_trafl2pdbs <ref> <cluster.trafl> 1 AA` "
                  "— frame 1 de cada cluster = estructura representativa (centroide).",
        "design": "## Diseño experimental",
        "design_note": "Las corridas que comparten aptámero difieren únicamente en la 2D restraint — esa es "
                       "la comparación para la que está armado el set. Las corridas sobre aptámeros "
                       "distintos son independientes y **no** son comparables estructuralmente entre sí.",
        "d_hdr": "| aptámero | run | nt | pb | 2D restraint (dot-bracket) |",
        "refcut": "## Cutoff de referencia por run",
        "refcut_text": "El RMSD es una distancia absoluta, así que un mismo cutoff no sirve para toda "
                       "longitud de cadena: un umbral que resuelve subestados en 60 nt colapsa una "
                       "horquilla de 13 nt en un solo cluster. El cutoff elegido por run es el que resuelve "
                       "la heterogeneidad conformacional sin colapsar en una megaclase única ni "
                       "fragmentarse en ruido térmico.",
        "r_hdr": "| run | nt | cutoff de referencia | población top-1 | veredicto |",
        "sweep": "## Barrido de cutoffs — población del cluster dominante",
        "sweep_text": "Cómo responde la población del cluster top-1 al cutoff. Un run está *saturado* cuando "
                      "el valor se sienta en ~100 % (el cutoff no resuelve nada) y *fragmentado* cuando el "
                      "cluster dominante es de unos pocos por ciento (el cutoff está por debajo de la "
                      "escala de fluctuación térmica).",
        "ss_cmp": "## 2D de referencia (restraint) vs 2D detectada del cluster top-1",
        "ss_cmp_text": "La 2D detectada del centroide rara vez coincide exactamente con el restraint. La "
                       "diferencia informa sobre cuánto está \"empujando\" SimRNA fuera de la 2D objetivo.",
        "s_hdr": "| run | 2D restraint | 2D top-1 @ cutoff de referencia | ¿coinciden? |",
        "tables": "## Top-{n} clusters al cutoff de referencia",
        "t_hdr": "| rank | cluster | frames | % | 2D centroide |",
        "sspat": "## Patrones 2D más poblados por run",
        "sspat_text": "Población de cada **patrón dot-bracket exacto** en el pool low-T, independiente del "
                      "clustering 3D. El patrón del propio restraint va marcado con **=**.",
        "p_hdr": "| rank | % | 2D (dot-bracket) | |",
        "pair": "## Comparación de restraints — {apt}",
        "pair_text": "Estas corridas comparten secuencia, así que sus clusterings son directamente "
                     "comparables y sus centroides se pueden superponer. Es el mano a mano que el set "
                     "existe para responder.",
        "gallery": "## Galería — centroides top-1",
        "gallery_cap": "_Centroide all-atom back-mapeado desde el frame 1 del cluster dominante al cutoff de "
                       "referencia de cada run, renderizado en ChimeraX (tube/slab muffler, rainbow por "
                       "residuo 5' azul → 3' rojo). Los centroides de un mismo aptámero comparten ancla de "
                       "alineamiento y cámara, así que son comparables entre sí; los de aptámeros distintos no._",
        "no_gallery": "_(todavía no hay renders de centroides — correr `scripts/extract_global_centroids.py --set {set}` "
                      "y luego `scripts/render_centroids_chimerax.py --set {set}`)_",
        "note": "_Los `.trafl`, `.pdb` y `.png` de clusters viven en `global_analysis/`. Sólo se versionan los "
                "REPORTs y los `README_{{EN,ES}}.md` por run; trayectorias y datos de clustering están gitignored._",
        "yes": "sí", "no": "no",
        "v_ok": "resuelve", "v_sat": "saturado", "v_frag": "fragmentado",
        "md": "## Estructuras seleccionadas para dinámica molecular",
        "md_text": "Los subestados de abajo son los puntos de partida de la etapa siguiente: **dinámica "
                   "molecular en solvente explícito** con AMBER. Los PDBs listos para AMBER y un "
                   "`tleap.in` de arranque están en "
                   "[`global_analysis/PDBtoMD/`](global_analysis/PDBtoMD) — ver su "
                   "[README](global_analysis/PDBtoMD/README_ES.md) para el detalle de la conversión.",
        "md_hdr": "| estructura | aptámero | subestado | población |",
        "md_mixed": "La selección mezcla deliberadamente dos tipos de subestado. Una entrada **RMSD** "
                    "representa un cluster 3D de frames estructuralmente similares. Una entrada **patrón "
                    "2D** es el frame de menor energía que porta un dot-bracket exacto — representa el "
                    "apareamiento de bases, no un cluster conformacional. Se tomaron subestados de patrón "
                    "2D donde el clustering RMSD no resuelve el aptámero: en las cadenas más cortas un "
                    "cluster absorbe casi todo el pool, así que lo que varía es el apareamiento y no el "
                    "plegamiento.",
        "per_run": "## Detalle por run",
        "per_run_line": "- [`{run}` →]({run}/README_ES.md) — clustering RMSD y por patrón 2D completo, todos los cutoffs.",
    },
}

SATURATED_PCT = 97.0    # top-1 at or above this: the cutoff resolves nothing
FRAGMENTED_PCT = 15.0   # top-1 below this: no dominant substate at this cutoff


# ---------- data gathering ---------------------------------------------------

class Run:
    def __init__(self, set_root: Path, in_root: Path, name: str):
        self.name = name
        self.run_dir = set_root / name
        self.in_dir = in_root / name
        self.seq = first_line(self.in_dir / "example") or first_line(self.run_dir / "example")
        self.ss = first_line(self.in_dir / "example.str") or first_line(self.run_dir / "example.str")
        self.header = ""
        src = self.in_dir / "_source.txt"
        if src.exists():
            self.header = first_line(src)
        self.aptamer = re.split(r"\s+[—–-]\s+", self.header, maxsplit=1)[0].strip() or name
        self.descriptor = ""
        if " " in self.header and self.aptamer != self.header:
            self.descriptor = self.header[len(self.aptamer):].lstrip(" —–-").strip()

        pool = self.run_dir / "low_temp.trafl"
        self.frames = count_frames(pool) if pool.exists() else 0
        rmsd_dir = self.run_dir / "rmsd_clusters"
        self.tables = rmsd_table(rmsd_dir, top_n=TOP_N) if rmsd_dir.is_dir() else {}
        ss_tsv = self.run_dir / "ss_clusters" / "ss_clusters.tsv"
        self.ss_rows = ss_table(ss_tsv, top_n=TOP_N) if ss_tsv.exists() else []
        self.ref_cutoff: float | None = None
        self.centroids: dict[float, list[dict]] = {}

    @property
    def nt(self) -> int:
        return len(self.seq)

    @property
    def bp(self) -> int:
        return self.ss.count("(")

    def top1_pct(self, cutoff: float) -> float | None:
        rows = self.tables.get(cutoff)
        return rows[0][2] if rows else None

    def verdict(self, lang: str, cutoff: float) -> str:
        pct = self.top1_pct(cutoff)
        t = T[lang]
        if pct is None:
            return "—"
        if pct >= SATURATED_PCT:
            return t["v_sat"]
        if pct < FRAGMENTED_PCT:
            return t["v_frag"]
        return t["v_ok"]

    def ref_rows(self) -> list[tuple[int, int, float]]:
        return self.tables.get(self.ref_cutoff, []) if self.ref_cutoff else []

    def centroid_ss(self, cutoff: float, cluster_id: int) -> str:
        for row in self.centroids.get(cutoff, []):
            if int(row["cluster"]) == cluster_id:
                return row["ss"]
        return ""

    def centroid_png(self, global_root: Path, cutoff: float, cluster_id: int) -> Path | None:
        p = global_root / self.name / f"cutoff_{int(round(cutoff))}A" / f"clust{cluster_id:02d}.png"
        return p if p.exists() else None


def first_line(path: Path) -> str:
    if not path.exists():
        return ""
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if s:
            return s
    return ""


def load_centroids(run: Run, global_root: Path) -> None:
    base = global_root / run.name
    if not base.is_dir():
        return
    for cut_dir in sorted(base.glob("cutoff_*A")):
        m = re.match(r"cutoff_(\d+)A", cut_dir.name)
        tsv = cut_dir / "_centroids.tsv"
        if not m or not tsv.exists():
            continue
        rows = []
        with tsv.open() as f:
            f.readline()
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                rows.append({"rank": int(parts[0]), "cluster": parts[1],
                             "frames": int(parts[2]), "pct": float(parts[3]), "ss": parts[4]})
        run.centroids[float(m.group(1))] = rows


def read_ladder(run: Run) -> tuple[str, str]:
    log = run.run_dir / "example.log"
    if not log.exists():
        return "", ""
    temps = re.findall(r"replica Temperature = ([\d.]+)", log.read_text(errors="ignore"))
    return (temps[0], temps[-1]) if temps else ("", "")


def fmt_cut(c: float | None) -> str:
    return "—" if c is None else f"{c:g}"


# ---------- section builders -------------------------------------------------

def design_section(runs: list[Run], lang: str) -> list[str]:
    t = T[lang]
    lines = [t["design"], "", t["design_note"], "", t["d_hdr"], "|:---|:---|---:|---:|:---|"]
    for r in runs:
        apt = f"**{r.aptamer}**"
        desc = f"`{r.name}`" + (f"<br>{r.descriptor}" if r.descriptor else "")
        lines.append(f"| {apt} | {desc} | {r.nt} | {r.bp} | `{r.ss}` |")
    lines.append("")
    # One line per aptamer, not per run: runs of the same aptamer share a sequence.
    seen: set[str] = set()
    for r in runs:
        if r.seq in seen:
            continue
        seen.add(r.seq)
        lines.append(f"`{r.aptamer}` ({r.nt} nt): `{r.seq}`  ")
    lines.append("")
    return lines


def refcut_section(runs: list[Run], lang: str) -> list[str]:
    t = T[lang]
    lines = [t["refcut"], "", t["refcut_text"], "", t["r_hdr"], "|:---|---:|---:|---:|:---|"]
    for r in runs:
        pct = r.top1_pct(r.ref_cutoff) if r.ref_cutoff else None
        lines.append(
            f"| `{r.name}` | {r.nt} | **{fmt_cut(r.ref_cutoff)} Å** | "
            f"{f'{pct:.1f} %' if pct is not None else '—'} | {r.verdict(lang, r.ref_cutoff)} |"
        )
    lines.append("")
    return lines


def sweep_section(runs: list[Run], lang: str) -> list[str]:
    t = T[lang]
    cutoffs = sorted({c for r in runs for c in r.tables})
    lines = [t["sweep"], "", t["sweep_text"], ""]
    lines.append("| run | " + " | ".join(f"{fmt_cut(c)} Å" for c in cutoffs) + " |")
    lines.append("|:---|" + "---:|" * len(cutoffs))
    for r in runs:
        cells = []
        for c in cutoffs:
            pct = r.top1_pct(c)
            if pct is None:
                cells.append("—")
            elif c == r.ref_cutoff:
                cells.append(f"**{pct:.1f} %**")
            else:
                cells.append(f"{pct:.1f} %")
        lines.append(f"| `{r.name}` | " + " | ".join(cells) + " |")
    lines.append("")
    return lines


def ss_compare_section(runs: list[Run], lang: str) -> list[str]:
    t = T[lang]
    lines = [t["ss_cmp"], "", t["ss_cmp_text"], "", t["s_hdr"], "|:---|:---|:---|:---:|"]
    for r in runs:
        rows = r.ref_rows()
        top2d = r.centroid_ss(r.ref_cutoff, rows[0][0]) if rows and r.ref_cutoff else ""
        match = "✓" if top2d and top2d == r.ss else "✗"
        lines.append(f"| `{r.name}` | `{r.ss}` | `{top2d or '—'}` | {match if top2d else '—'} |")
    lines.append("")
    return lines


def cluster_tables(runs: list[Run], lang: str) -> list[str]:
    t = T[lang]
    lines = [t["tables"].format(n=TOP_N), ""]
    for r in runs:
        rows = r.ref_rows()
        lines.append(f"### `{r.name}` — {fmt_cut(r.ref_cutoff)} Å")
        lines.append("")
        if not rows:
            lines += ["_—_", ""]
            continue
        lines.append(t["t_hdr"])
        lines.append("|-----:|--------:|-------:|---:|:---|")
        for rank, (cid, cnt, pct) in enumerate(rows, start=1):
            ss = r.centroid_ss(r.ref_cutoff, cid)
            lines.append(f"| {rank} | {cid:02d} | {cnt} | {pct:.2f} | `{ss}` |" if ss
                         else f"| {rank} | {cid:02d} | {cnt} | {pct:.2f} | — |")
        lines.append("")
    return lines


def ss_pattern_section(runs: list[Run], lang: str) -> list[str]:
    t = T[lang]
    lines = [t["sspat"], "", t["sspat_text"], ""]
    for r in runs:
        lines.append(f"### `{r.name}`")
        lines.append("")
        if not r.ss_rows:
            lines += ["_—_", ""]
            continue
        lines.append(t["p_hdr"])
        lines.append("|-----:|---:|:-----------------|:--|")
        for rank, ss, cnt, pct in r.ss_rows:
            mark = "**=**" if ss == r.ss else ""
            lines.append(f"| {rank} | {pct:.2f} | `{ss}` | {mark} |")
        lines.append("")
    return lines


def pair_sections(runs: list[Run], lang: str) -> list[str]:
    """One head-to-head section per aptamer that has more than one run."""
    t = T[lang]
    by_apt: dict[str, list[Run]] = {}
    for r in runs:
        by_apt.setdefault(r.aptamer, []).append(r)

    lines: list[str] = []
    for apt, group in by_apt.items():
        if len(group) < 2:
            continue
        if len({r.seq for r in group}) != 1:
            continue  # not actually the same sequence — nothing comparable
        lines += [t["pair"].format(apt=apt), "", t["pair_text"], ""]
        cutoffs = sorted({c for r in group for c in r.tables})
        lines.append("| " + " | ".join(
            ["run"] +
            [f"{fmt_cut(c)} Å" for c in cutoffs] +
            [("restraint bp" if lang == "en" else "pb restraint"),
             ("restraint recovered" if lang == "en" else "restraint recuperado")]
        ) + " |")
        lines.append("|:---|" + "---:|" * len(cutoffs) + "---:|:---|")
        for r in group:
            cells = [f"`{r.name}`"]
            for c in cutoffs:
                pct = r.top1_pct(c)
                cells.append(f"{pct:.1f} %" if pct is not None else "—")
            cells.append(str(r.bp))
            hit = next((f"{pct:.2f} %" for _, ss, _, pct in r.ss_rows if ss == r.ss), None)
            cells.append(hit if hit else (T[lang]["no"] + " (top-10)"))
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return lines


def md_selection_section(runs: list[Run], set_root: Path, lang: str) -> list[str]:
    """Read back what build_pdbtomd.py selected, so the two never disagree."""
    t = T[lang]
    readme = set_root / "global_analysis" / "PDBtoMD" / f"README_{lang.upper()}.md"
    if not readme.exists():
        return []
    rows = []
    for line in readme.read_text().splitlines():
        if not line.startswith("| `") or ".pdb`" not in line:
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 5:
            continue
        rows.append(cells)
    if not rows:
        return []

    by_name = {r.name: r for r in runs}
    lines = [t["md"], "", t["md_text"], "", t["md_hdr"], "|:---|:---|:---|---:|"]
    for c in rows:
        fname = c[0].strip("`")
        run = next((n for n in by_name if fname.startswith(n)), "")
        apt = by_name[run].aptamer if run else ""
        lines.append(f"| {c[0]} | **{apt}** | {c[1]} | {c[4]} |")
    lines.append("")
    if len({c[1] for c in rows}) > 1:
        lines.append(t["md_mixed"])
        lines.append("")
    return lines


def gallery_section(runs: list[Run], global_root: Path, set_name: str, lang: str,
                    report_dir: Path) -> list[str]:
    t = T[lang]
    cells: list[str] = []
    for r in runs:
        rows = r.ref_rows()
        if not rows or not r.ref_cutoff:
            continue
        png = r.centroid_png(global_root, r.ref_cutoff, rows[0][0])
        if png:
            rel = png.relative_to(report_dir).as_posix()
            cells.append(f"| `{r.name}` | {rows[0][2]:.1f} % | ![{r.name}]({rel}) |")
    lines = [t["gallery"], ""]
    if not cells:
        lines += [t["no_gallery"].format(set=set_name), ""]
        return lines
    lines.append("| run | top-1 | centroid |" if lang == "en" else "| run | top-1 | centroide |")
    lines.append("|:---|---:|:---:|")
    lines += cells
    lines += ["", t["gallery_cap"], ""]
    return lines


# ---------- entry point -----------------------------------------------------

def build(set_name: str, runs: list[Run], set_root: Path, global_root: Path, lang: str) -> Path:
    t = T[lang]
    napt = len({r.aptamer for r in runs})
    frames = runs[0].frames if runs else 0
    lo, hi = read_ladder(runs[0]) if runs else ("", "")
    cutoffs = sorted({c for r in runs for c in r.tables})

    parts = brand_header(set_root, lang, BASE)
    parts.append(f"# {t['title'].format(set=set_name)}")
    parts.append("")
    parts.append(t["intro"].format(set=set_name, nruns=len(runs), napt=napt))
    parts.append("")
    parts.append(t["b_pool"].format(frames=f"{frames:,}".replace(",", " ")))
    if lo and hi:
        parts.append(t["b_ladder"].format(lo=lo, hi=hi))
    parts.append(t["b_clust"].format(cutoffs=", ".join(f"{fmt_cut(c)}" for c in cutoffs) + " Å"))
    parts.append(t["b_back"])
    parts.append("")
    parts += design_section(runs, lang)
    parts += refcut_section(runs, lang)
    parts += sweep_section(runs, lang)
    parts += ss_compare_section(runs, lang)
    parts += cluster_tables(runs, lang)
    parts += ss_pattern_section(runs, lang)
    parts += pair_sections(runs, lang)
    parts += gallery_section(runs, global_root, set_name, lang, set_root)
    parts += md_selection_section(runs, set_root, lang)
    parts += [t["per_run"], ""]
    for r in runs:
        parts.append(t["per_run_line"].format(run=r.name))
    parts.append("")
    parts.append(t["note"].format())
    parts.append("")
    parts += brand_footer(lang)

    path = out_path(set_root, lang, BASE)
    path.write_text("\n".join(parts))
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("set_name")
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--ref-cutoff", action="append", default=[], metavar="RUN=CUTOFF",
                    help="Reference cutoff in Å for one run; repeatable.")
    ap.add_argument("--default-ref-cutoff", type=float, default=12.0)
    args = ap.parse_args()

    set_root = args.repo / "results" / args.set_name
    in_root = args.repo / "inputs" / args.set_name
    global_root = set_root / "global_analysis"
    if not set_root.is_dir():
        print(f"no such result set: {set_root}", file=sys.stderr)
        return 1

    overrides: dict[str, float] = {}
    for spec in args.ref_cutoff:
        if "=" not in spec:
            print(f"--ref-cutoff expects RUN=CUTOFF, got {spec!r}", file=sys.stderr)
            return 1
        k, v = spec.split("=", 1)
        overrides[k] = float(v)

    names = sorted(p.name for p in set_root.iterdir()
                   if p.is_dir() and (p / "rmsd_clusters").is_dir())
    if not names:
        print(f"no runs with rmsd_clusters/ under {set_root}", file=sys.stderr)
        return 1

    runs = [Run(set_root, in_root, n) for n in names]
    runs.sort(key=lambda r: (r.nt, r.name))
    unknown = set(overrides) - {r.name for r in runs}
    if unknown:
        print(f"--ref-cutoff names not in the set: {sorted(unknown)}", file=sys.stderr)
        return 1

    for r in runs:
        load_centroids(r, global_root)
        want = overrides.get(r.name, args.default_ref_cutoff)
        if want in r.tables:
            r.ref_cutoff = want
        elif r.tables:
            # Fall back to the nearest clustered cutoff rather than emitting blanks.
            r.ref_cutoff = min(r.tables, key=lambda c: abs(c - want))
            print(f"warning: {r.name} has no {fmt_cut(want)} Å clustering; "
                  f"using {fmt_cut(r.ref_cutoff)} Å", file=sys.stderr)
        print(f"{r.name:<24} {r.nt:>3} nt  ref={fmt_cut(r.ref_cutoff or 0):>4} Å  "
              f"cutoffs={[fmt_cut(c) for c in sorted(r.tables)]}")

    for lang in LANGS:
        print(f"wrote {build(args.set_name, runs, set_root, global_root, lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
