#!/usr/bin/env python3
"""
Build the global APT-PF1 report as a bilingual pair, `REPORT_EN.md` and
`REPORT_ES.md`, both carrying the DansLab brand header (see
`scripts/report_brand.py`). The English version drops the per-row `[[pdb]]` /
`[[png]]` links requested for the cleaner docx conversion.

Reads:
  results/APT-PF1/<run>/{<name>,<name>.str}
  results/APT-PF1/global_analysis/<run>/cutoff_<N>A/_centroids.tsv

Writes:
  results/APT-PF1/REPORT_EN.md    (en, no links)
  results/APT-PF1/REPORT_ES.md    (es, with per-row links)
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_brand import LANGS, brand_footer, brand_header, out_path  # noqa: E402

BASE = "REPORT"

#: Per-language switch for the per-row `[[pdb]]` / `[[png]]` links. Spanish is
#: the working copy read on screen; English is the one converted to DOCX, where
#: the inline links clutter the tables.
WITH_LINKS = {"en": False, "es": True}

RUNS = [
    ("DNA_as_RNA_VFold2D", "VFold2D"),
    ("DNA_as_RNA_NUPACK4", "NUPACK4"),
    ("NA_as_DNA_NUPACK4", "NA_as_DNA"),
]

# Source attributions for the three 2D restraints (identical for both langs).
SOURCE_ES = {
    "DNA_as_RNA_VFold2D": "VFold2D (tratando la secuencia como RNA)",
    "DNA_as_RNA_NUPACK4": "NUPACK v4 (tratando la secuencia como RNA)",
    "NA_as_DNA_NUPACK4": "NUPACK v4 (tratando la secuencia como DNA)",
}
SOURCE_EN = {
    "DNA_as_RNA_VFold2D": "VFold2D (sequence treated as RNA)",
    "DNA_as_RNA_NUPACK4": "NUPACK v4 (sequence treated as RNA)",
    "NA_as_DNA_NUPACK4": "NUPACK v4 (sequence treated as DNA)",
}

MD_SELECTION = [
    # (run_dir, cluster file suffix without .pdb, short label)
    ("DNA_as_RNA_VFold2D", "cutoff_12A/clust01", "VFold2D c12 clust01"),
    ("DNA_as_RNA_VFold2D", "cutoff_12A/clust02", "VFold2D c12 clust02"),
    ("DNA_as_RNA_NUPACK4", "cutoff_12A/clust01", "NUPACK4 c12 clust01"),
    ("NA_as_DNA_NUPACK4",  "cutoff_12A/clust01", "NA_as_DNA c12 clust01"),
]


# ---------- helpers ----------------------------------------------------------

def read_seq(apt_root: Path, run_dir: str, name: str) -> str:
    f = apt_root / run_dir / name
    if f.exists():
        for ln in f.read_text().splitlines():
            s = ln.strip()
            if s:
                return s
    return ""


def read_str(apt_root: Path, run_dir: str, name: str) -> str:
    f = apt_root / run_dir / f"{name}.str"
    if f.exists():
        for ln in f.read_text().splitlines():
            s = ln.strip()
            if s:
                return s
    return ""


def read_centroids(apt_root: Path, run_dir: str, cutoff: int) -> list[dict]:
    tsv = apt_root / "global_analysis" / run_dir / f"cutoff_{cutoff}A" / "_centroids.tsv"
    rows: list[dict] = []
    if not tsv.exists():
        return rows
    with tsv.open() as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rows.append({
                "rank": int(parts[0]), "cluster": parts[1],
                "frames": int(parts[2]), "pct": float(parts[3]), "ss": parts[4],
            })
    return rows


# ---------- section builders ------------------------------------------------

def header(lang: str, seq: str) -> list[str]:
    if lang == "es":
        return [
            "# APT-PF1 — subestados conformacionales (análisis global)",
            "",
            "Resumen cruzado de las tres corridas SimRNA REMC del aptámero APT-PF1 sobre la misma secuencia con tres 2D restraints distintas. Pipeline: `scripts/run_cluster_populations.sh` por run, luego `scripts/extract_global_centroids.py` para producir centroides AA + PNG en [`global_analysis/`](global_analysis).",
            "",
            "- **Pool de frames**: niveles REMC 2–4 (T≈1 en la escala relativa), 30 003 frames por run.",
            "- **Clustering RMSD**: cutoffs 8, 12, 16 Å con el binario `clustering` de SimRNA.",
            "- **Back-mapping a all-atom**: `SimRNA_trafl2pdbs <ref> <cluster.trafl> 1 AA` — frame 1 de cada cluster = estructura representativa (centroide).",
            "- **Alineamiento**: los 86 PDBs fueron superpuestos por `align … @C1'` a `DNA_as_RNA_VFold2D/cutoff_12A/clust01` (cluster dominante al 60.4 %).",
            "- **Snapshot PNG**: estilo *tube/slab muffler* con cartoon de sección redonda, rainbow por residuo, vía ChimeraX 1.11.1 (modo headless `--nogui --offscreen`, OSMesa). Cámara fija compartida por los 86 renders.",
            "",
        ]
    return [
        "# APT-PF1 — conformational substates (global analysis)",
        "",
        "Cross-run summary of the three SimRNA REMC simulations of aptamer APT-PF1 sharing the same sequence but using three different 2D restraints. Pipeline: per-run `scripts/run_cluster_populations.sh`, followed by `scripts/extract_global_centroids.py` which produces all-atom centroids + PNGs in `global_analysis/`.",
        "",
        "- **Frame pool**: REMC levels 2–4 (T≈1 on the relative scale), 30,003 frames per run.",
        "- **RMSD clustering**: cutoffs 8, 12, 16 Å using SimRNA's `clustering` binary.",
        "- **Coarse-grain → all-atom back-mapping**: `SimRNA_trafl2pdbs <ref> <cluster.trafl> 1 AA` — frame 1 of each cluster is the representative (centroid).",
        "- **Alignment**: all 86 PDBs superposed by `align … @C1'` onto `DNA_as_RNA_VFold2D/cutoff_12A/clust01` (the 60.4 % dominant cluster).",
        "- **PNG snapshot**: *tube/slab muffler* nucleotide depiction with round cartoon cross-section, rainbow by residue, rendered via ChimeraX 1.11.1 (headless `--nogui --offscreen`, OSMesa). A single camera is shared across all 86 renders.",
        "",
    ]


def pipeline_image(lang: str) -> list[str]:
    alt = "Pipeline simRNA: muestreo REMC + reconstrucción all-atom" if lang == "es" \
          else "simRNA pipeline: REMC sampling + all-atom reconstruction"
    return [f"![{alt}](../../docs/simRNA-pipeline.png)", ""]


def design_section(apt_root: Path, lang: str, seq: str) -> list[str]:
    if lang == "es":
        lines = [
            "## Diseño experimental",
            "",
            f"Secuencia compartida (40 nt): `{seq}`",
            "",
            "| aptámero | 2D restraint (dot-bracket) | fuente del 2D |",
            "|:---|:---|:---|",
        ]
        src = SOURCE_ES
    else:
        lines = [
            "## Experimental design",
            "",
            f"Shared sequence (40 nt): `{seq}`",
            "",
            "| aptamer | 2D restraint (dot-bracket) | 2D source |",
            "|:---|:---|:---|",
        ]
        src = SOURCE_EN
    for run_dir, short in RUNS:
        s2d = read_str(apt_root, run_dir, "example")
        lines.append(f"| `{short}` | `{s2d}` | {src[run_dir]} |")
    lines.append("")
    return lines


def topcluster_summary(data: dict, lang: str) -> list[str]:
    if lang == "es":
        lines = [
            "## Resumen — población del cluster dominante por cutoff",
            "",
            "| run | top-1 @ 8Å | top-1 @ 12Å | top-1 @ 16Å |",
            "|:---|---:|---:|---:|",
        ]
    else:
        lines = [
            "## Summary — dominant-cluster population by cutoff",
            "",
            "| run | top-1 @ 8Å | top-1 @ 12Å | top-1 @ 16Å |",
            "|:---|---:|---:|---:|",
        ]
    for run_dir, short in RUNS:
        pcts = []
        for c in (8, 12, 16):
            rows = data[run_dir].get(c, [])
            pcts.append(f"{rows[0]['pct']:.1f} %" if rows else "—")
        lines.append(f"| `{short}` | {pcts[0]} | {pcts[1]} | {pcts[2]} |")
    lines.append("")
    return lines


def gallery(data: dict, lang: str) -> list[str]:
    if lang == "es":
        title = "## Galería — centroides top-1 por cutoff"
        header_row = "| cutoff | VFold2D | NUPACK4 | NA_as_DNA |"
        caption = "_Imágenes: centroide AA alineado al top-1 @ 12Å de `VFold2D` (ancla común, por `align @C1'`), en estilo **tube/slab muffler** (backbone tube + láminas muffler), rainbow por residuo 5' (azul) → 3' (rojo). Renderizado vía ChimeraX 1.11.1 headless (OSMesa)._"
    else:
        title = "## Gallery — top-1 centroids by cutoff"
        header_row = "| cutoff | VFold2D | NUPACK4 | NA_as_DNA |"
        caption = "_Images: AA centroid aligned onto the top-1 @ 12Å of `VFold2D` (shared anchor via `align @C1'`), rendered in **tube/slab muffler** style (tube backbone + muffler slabs), rainbow by residue from 5' (blue) to 3' (red). ChimeraX 1.11.1 headless (OSMesa)._"

    lines = [title, "", header_row, "|:---:|:---:|:---:|:---:|"]
    for cutoff in (8, 12, 16):
        cells = [f"**{cutoff} Å**"]
        for run_dir, short in RUNS:
            rows = data[run_dir].get(cutoff, [])
            if not rows:
                cells.append("—")
                continue
            top = rows[0]
            rel = f"global_analysis/{run_dir}/cutoff_{cutoff}A/clust{top['cluster']}.png"
            cells.append(f"![clust{top['cluster']} — {top['pct']:.1f}%]({rel})")
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    lines.append(caption)
    lines.append("")
    return lines


def ref_vs_top1(data: dict, apt_root: Path, lang: str) -> list[str]:
    if lang == "es":
        lines = [
            "## Referencia (2D restraint) vs 2D detectada del top-1 @ 12 Å",
            "",
            "La 2D detectada del centroide rara vez coincide exactamente con el restraint. La diferencia informa sobre cuánto está \"empujando\" SimRNA fuera de la 2D objetivo.",
            "",
            "| run | 2D restraint | 2D top-1 @ 12 Å | ¿coinciden? |",
            "|:---|:---|:---|:---:|",
        ]
    else:
        lines = [
            "## Reference (2D restraint) vs. detected 2D of the top-1 @ 12 Å",
            "",
            "The centroid's detected 2D rarely matches the restraint exactly. The gap tells us how hard SimRNA is pushing off the target 2D.",
            "",
            "| run | 2D restraint | 2D of top-1 @ 12 Å | match? |",
            "|:---|:---|:---|:---:|",
        ]
    for run_dir, short in RUNS:
        ref2d = read_str(apt_root, run_dir, "example")
        rows = data[run_dir].get(12, [])
        top2d = rows[0]["ss"] if rows else ""
        match = "✓" if ref2d and top2d and ref2d == top2d else "✗"
        lines.append(f"| `{short}` | `{ref2d}` | `{top2d or '—'}` | {match} |")
    lines.append("")
    return lines


def cluster_comparison_tables(data: dict, lang: str, with_links: bool) -> list[str]:
    if lang == "es":
        title = "## Comparación de los top-10 clusters por cutoff"
        col_headers = "| rank | VFold2D (% · 2D centroide) | NUPACK4 (% · 2D centroide) | NA_as_DNA (% · 2D centroide) |"
    else:
        title = "## Top-10 cluster comparison per cutoff"
        col_headers = "| rank | VFold2D (% · centroid 2D) | NUPACK4 (% · centroid 2D) | NA_as_DNA (% · centroid 2D) |"

    lines = [title, ""]
    for cutoff in (8, 12, 16):
        lines.append(f"### Cutoff {cutoff} Å")
        lines.append("")
        lines.append(col_headers)
        lines.append("|:---:|:---|:---|:---|")
        max_n = max(len(data[r].get(cutoff, [])) for r, _ in RUNS)
        for i in range(max_n):
            cells = [str(i + 1)]
            for run_dir, _ in RUNS:
                rows = data[run_dir].get(cutoff, [])
                if i >= len(rows):
                    cells.append("—")
                    continue
                r = rows[i]
                ss = r["ss"] or "—"
                cell = f"**{r['pct']:.2f} %** · `{ss}`"
                if with_links:
                    pdb_rel = f"global_analysis/{run_dir}/cutoff_{cutoff}A/clust{r['cluster']}.pdb"
                    png_rel = f"global_analysis/{run_dir}/cutoff_{cutoff}A/clust{r['cluster']}.png"
                    cell += f" [[pdb]]({pdb_rel}) [[png]]({png_rel})"
                cells.append(cell)
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return lines


def md_selection_section(data: dict, lang: str) -> list[str]:
    """The 'selected for MD' section that appears in every report."""
    # Sanity check populations (criterion > 15 %).
    pops = {}
    for run_dir, subpath, label in MD_SELECTION:
        _, cluster_part = subpath.split("/")
        cutoff = 12
        clust_id = int(cluster_part.replace("clust", ""))
        rows = data[run_dir].get(cutoff, [])
        for r in rows:
            if int(r["cluster"]) == clust_id:
                pops[label] = r["pct"]
                break

    if lang == "es":
        lines = [
            "## Selección de subestados para dinámica molecular (AMBER)",
            "",
            "Para el siguiente paso — dinámica atomística con AMBER — nos **quedamos con los resultados del clustering RMSD al cutoff 12 Å**. Ese nivel es el que mejor resuelve la heterogeneidad conformacional sin caer ni en ruido térmico (típico del 8 Å sobre 40 nt flexibles) ni en un única megaclase (16 Å).",
            "",
            "**Criterio**: mantener todos los subestados con **población > 15 %** del pool low-T (30 003 frames por run). Eso deja cuatro estructuras:",
            "",
            "| run | cluster | población (%) |",
            "|:---|:---|---:|",
        ]
        for _, _, label in MD_SELECTION:
            lines.append(f"| `{label.split(' c12 ')[0]}` | `{label.split(' c12 ')[1]}` | {pops.get(label, 0):.1f} |")
        lines.append("")
        lines.append("Los cuatro PDBs all-atom, listos para el pipeline `tleap`→`sander/pmemd`, viven en [`global_analysis/PDBtoMD/`](global_analysis/PDBtoMD). La conversión (residuos renombrados a `RA/RC/RG/RU` con cap 5'/3', fosfato 5' terminal removido, atom serials renumerados, `TER`/`END` finales) la hace `scripts/prepare_for_amber.py`; también se incluye un [`tleap.in`](global_analysis/PDBtoMD/tleap.in) de arranque.")
        lines.append("")
        lines.append("![4 centroides seleccionados para MD](global_analysis/summary_md_selection.png)")
        lines.append("")
        lines.append("_Los cuatro centroides alineados al mismo frame, renderizados en ChimeraX (tube/slab muffler, rainbow por residuo)._")
        lines.append("")
    else:
        lines = [
            "## Substates selected for atomistic dynamics (AMBER)",
            "",
            "For the next step — atomistic MD with AMBER — we **keep the RMSD clustering results at the 12 Å cutoff**. That cutoff resolves conformational heterogeneity best without collapsing into thermal noise (what 8 Å does on 40-nt flexible RNA) nor into a single mega-class (what 16 Å does).",
            "",
            "**Criterion**: keep every substate with **population > 15 %** of the low-T pool (30,003 frames per run). That yields four structures:",
            "",
            "| run | cluster | population (%) |",
            "|:---|:---|---:|",
        ]
        for _, _, label in MD_SELECTION:
            lines.append(f"| `{label.split(' c12 ')[0]}` | `{label.split(' c12 ')[1]}` | {pops.get(label, 0):.1f} |")
        lines.append("")
        lines.append("The four all-atom PDBs, ready for the `tleap`→`sander/pmemd` pipeline, live in `global_analysis/PDBtoMD/`. Conversion (residues renamed to `RA/RC/RG/RU` with 5'/3' caps, 5'-terminal phosphate stripped, atom serials renumbered, trailing `TER`/`END`) is performed by `scripts/prepare_for_amber.py`; a starter `tleap.in` is bundled with the PDBs.")
        lines.append("")
        lines.append("![4 centroids selected for MD](global_analysis/summary_md_selection.png)")
        lines.append("")
        lines.append("_The four centroids aligned onto the same frame, rendered in ChimeraX (tube/slab muffler, rainbow by residue)._")
        lines.append("")
    return lines


def per_run_links(lang: str) -> list[str]:
    suffix = "ES" if lang == "es" else "EN"
    if lang == "es":
        lines = ["## Detalle por run", ""]
        for run_dir, short in RUNS:
            lines.append(f"- [`{short}` →]({run_dir}/README_{suffix}.md) — RMSD + clustering por patrón 2D completo (incluye cutoff 8/12/16 Å).")
    else:
        lines = ["## Per-run detail", ""]
        for run_dir, short in RUNS:
            lines.append(f"- [`{short}` →]({run_dir}/README_{suffix}.md) — RMSD + 2D-pattern clustering (all cutoffs).")
    lines.append("")
    return lines


def data_note(lang: str) -> list[str]:
    if lang == "es":
        return [
            "_Los `.pdb` y `.png` de los centroides viven en [`global_analysis/`](global_analysis). Sólo los REPORTs se versionan; los `rmsd_clusters/`, `ss_clusters/`, `low_temp.trafl`, y el contenido de `global_analysis/` están gitignored._",
            "",
        ]
    return [
        "_Centroid `.pdb` and `.png` files live under `global_analysis/`. Only the REPORTs are versioned; `rmsd_clusters/`, `ss_clusters/`, `low_temp.trafl` and `global_analysis/` contents are gitignored._",
        "",
    ]


# ---------- entry point -----------------------------------------------------

def build(apt_root: Path, lang: str) -> None:
    seq = read_seq(apt_root, RUNS[0][0], "example")
    data = {run_dir: {c: read_centroids(apt_root, run_dir, c) for c in (8, 12, 16)} for run_dir, _ in RUNS}

    parts: list[str] = brand_header(apt_root, lang, BASE)
    parts += header(lang, seq)
    parts += pipeline_image(lang)
    parts += design_section(apt_root, lang, seq)
    parts += topcluster_summary(data, lang)
    parts += gallery(data, lang)
    parts += ref_vs_top1(data, apt_root, lang)
    parts += cluster_comparison_tables(data, lang, WITH_LINKS[lang])
    parts += md_selection_section(data, lang)
    parts += per_run_links(lang)
    parts += data_note(lang)
    parts += brand_footer(lang)

    path = out_path(apt_root, lang, BASE)
    path.write_text("\n".join(parts))
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    args = ap.parse_args()
    apt_root = args.repo / "results" / "APT-PF1"

    for lang in LANGS:
        build(apt_root, lang)


if __name__ == "__main__":
    main()
