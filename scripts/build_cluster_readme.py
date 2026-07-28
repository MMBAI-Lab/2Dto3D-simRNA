#!/usr/bin/env python3
"""
Build the per-run conformational-substate report inside a SimRNA run directory,
summarizing cluster populations for the low-T pool at every RMSD cutoff and by
exact 2D (`.ss_detected`) pattern.

Reports ship as a bilingual pair carrying the DansLab brand header (see
`scripts/report_brand.py`).

Inputs (under <run_dir>):
  low_temp.trafl                     — pooled low-T frames (for total count)
  <name>.str                         — reference dot-bracket (shown for context)
  rmsd_clusters/low_temp_thrs*A_clust*.trafl
                                     — per-cutoff cluster trafl files
  ss_clusters/ss_clusters.tsv        — TSV from ss_cluster.py

Writes:
  <run_dir>/README_EN.md
  <run_dir>/README_ES.md
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report_brand import LANGS, brand_footer, brand_header, out_path  # noqa: E402

FILENAME_RE = re.compile(r"low_temp_thrs(?P<cutoff>[\d.]+)A_clust(?P<clust>\d+)\.trafl")

BASE = "README"

# ---------- localized strings ------------------------------------------------

T = {
    "en": {
        "title": "Conformational substates — `{run}`",
        "intro": "Clustering analysis over the low-temperature frame pool (REMC levels 2–4, "
                 "T≈1 on the relative scale) of the SimRNA run in this directory.",
        "frames": "- **Frames in the low-T pool**: {n}",
        "seq": "- **Sequence**: `{seq}`",
        "ref2d": "- **Reference 2D** (SimRNA restraint): `{ss}`",
        "pipeline": "- **Pipeline**: `scripts/run_cluster_populations.sh` → `extract_low_temp_frames.py` "
                    "(--min-temp 2 --max-temp 4) → `clustering` ({cutoffs}) + `SimRNA_trafl2pdbs` + `ss_cluster.py`.",
        "rmsd_title": "## RMSD clustering — top 10 per cutoff",
        "rmsd_head": "| rank | cluster | frames | % |",
        "rmsd_note": "_Top 10 cover {n} frames of the pool ({k} clusters reported)._",
        "rmsd_empty": "## RMSD clustering\n\n_(no results — `rmsd_clusters/` is empty)_\n",
        "ss_title": "## Clustering by 2D pattern (`.ss_detected`) — top 10",
        "ss_intro": "Each substate is an **exact dot-bracket pattern**. The population is the fraction "
                    "of low-T pool frames exhibiting that 2D.",
        "ss_head": "| rank | 2D (dot-bracket) | frames | % |",
        "ss_empty": "## Clustering by 2D pattern\n\n_(no results — `ss_clusters/ss_clusters.tsv` missing)_\n",
        "raw": "_Raw files: `low_temp.trafl`, `rmsd_clusters/`, `ss_clusters/`. All excluded by "
               "`.gitignore` except these reports._",
    },
    "es": {
        "title": "Subestados conformacionales — `{run}`",
        "intro": "Análisis de clustering sobre el pool de frames de baja temperatura (niveles REMC 2–4, "
                 "T≈1 en la escala relativa) de la corrida SimRNA de esta carpeta.",
        "frames": "- **Frames en el pool low-T**: {n}",
        "seq": "- **Secuencia**: `{seq}`",
        "ref2d": "- **2D de referencia** (restricción SimRNA): `{ss}`",
        "pipeline": "- **Pipeline**: `scripts/run_cluster_populations.sh` → `extract_low_temp_frames.py` "
                    "(--min-temp 2 --max-temp 4) → `clustering` ({cutoffs}) + `SimRNA_trafl2pdbs` + `ss_cluster.py`.",
        "rmsd_title": "## RMSD clustering — top 10 por cutoff",
        "rmsd_head": "| rank | cluster | frames | % |",
        "rmsd_note": "_Top 10 cubren {n} frames del pool ({k} clusters reportados)._",
        "rmsd_empty": "## RMSD clustering\n\n_(no hay resultados — `rmsd_clusters/` vacío)_\n",
        "ss_title": "## Clustering por patrón 2D (`.ss_detected`) — top 10",
        "ss_intro": "Cada subestado corresponde a un **patrón dot-bracket exacto**. La población es la "
                    "fracción de frames del pool low-T que exhiben esa 2D.",
        "ss_head": "| rank | 2D (dot-bracket) | frames | % |",
        "ss_empty": "## Clustering por patrón 2D\n\n_(no hay resultados — `ss_clusters/ss_clusters.tsv` ausente)_\n",
        "raw": "_Archivos crudos: `low_temp.trafl`, `rmsd_clusters/`, `ss_clusters/`. Todos excluidos por "
               "`.gitignore` salvo estos reportes._",
    },
}


# ---------- data -------------------------------------------------------------

def count_frames(trafl: Path) -> int:
    n = 0
    with trafl.open() as f:
        for line in f:
            # headers start with a digit; coord lines start with whitespace
            if line and line[0].isdigit():
                n += 1
    return n


def rmsd_table(rmsd_dir: Path, top_n: int = 10) -> dict[float, list[tuple[int, int, float]]]:
    """Return {cutoff -> sorted list of (cluster_id, count, pct)} (desc by count)."""
    per_cutoff: dict[float, list[tuple[int, int]]] = {}
    for path in sorted(rmsd_dir.glob("low_temp_thrs*A_clust*.trafl")):
        m = FILENAME_RE.match(path.name)
        if not m:
            continue
        cutoff = float(m.group("cutoff"))
        clust_id = int(m.group("clust"))
        per_cutoff.setdefault(cutoff, []).append((clust_id, count_frames(path)))

    out: dict[float, list[tuple[int, int, float]]] = {}
    for cutoff, clusters in per_cutoff.items():
        total = sum(c for _, c in clusters) or 1
        clusters_sorted = sorted(clusters, key=lambda x: (-x[1], x[0]))
        out[cutoff] = [(cid, cnt, 100.0 * cnt / total) for cid, cnt in clusters_sorted[:top_n]]
    return out


def ss_table(ss_tsv: Path, top_n: int = 10) -> list[tuple[int, str, int, float]]:
    rows = []
    with ss_tsv.open() as f:
        f.readline()  # skip TSV header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            rank, ss, count, pct = parts
            rows.append((int(rank), ss, int(count), float(pct)))
            if len(rows) >= top_n:
                break
    return rows


# ---------- sections ---------------------------------------------------------

def fmt_rmsd_section(tables: dict[float, list[tuple[int, int, float]]], lang: str) -> str:
    t = T[lang]
    lines = [t["rmsd_title"], ""]
    for cutoff in sorted(tables.keys()):
        rows = tables[cutoff]
        total_in_top = sum(r[1] for r in rows)
        lines.append(f"### Cutoff {cutoff:.1f} Å")
        lines.append("")
        lines.append(t["rmsd_head"])
        lines.append("|-----:|--------:|-------:|---:|")
        for rank, (cid, cnt, pct) in enumerate(rows, start=1):
            lines.append(f"| {rank} | {cid:02d} | {cnt} | {pct:.2f} |")
        lines.append("")
        lines.append(t["rmsd_note"].format(n=total_in_top, k=len(rows)))
        lines.append("")
    return "\n".join(lines)


def fmt_ss_section(rows: list[tuple[int, str, int, float]], lang: str) -> str:
    t = T[lang]
    lines = [t["ss_title"], ""]
    lines.append(t["ss_intro"])
    lines.append("")
    lines.append(t["ss_head"])
    lines.append("|-----:|:-----------------|-------:|---:|")
    for rank, ss, cnt, pct in rows:
        lines.append(f"| {rank} | `{ss}` | {cnt} | {pct:.2f} |")
    lines.append("")
    return "\n".join(lines)


def build(run_dir: Path, lang: str, total_pool: int, ref_seq: str, ref_2d: str,
          rmsd_tables: dict, ss_rows: list, cutoffs_str: str) -> None:
    t = T[lang]

    parts: list[str] = brand_header(run_dir, lang, BASE)
    parts.append(f"# {t['title'].format(run=run_dir.name)}")
    parts.append("")
    parts.append(t["intro"])
    parts.append("")
    parts.append(t["frames"].format(n=total_pool))
    if ref_seq:
        parts.append(t["seq"].format(seq=ref_seq))
    if ref_2d:
        parts.append(t["ref2d"].format(ss=ref_2d))
    parts.append(t["pipeline"].format(cutoffs=cutoffs_str))
    parts.append("")
    parts.append(fmt_rmsd_section(rmsd_tables, lang) if rmsd_tables else t["rmsd_empty"])
    parts.append(fmt_ss_section(ss_rows, lang) if ss_rows else t["ss_empty"])
    parts.append(t["raw"])
    parts.append("")
    parts += brand_footer(lang)

    path = out_path(run_dir, lang, BASE)
    path.write_text("\n".join(parts))
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--name", required=True, help="Base name of the run (e.g. 'example')")
    args = ap.parse_args()

    run_dir: Path = args.run_dir.resolve()
    name: str = args.name

    pool = run_dir / "low_temp.trafl"
    seq_file = run_dir / name
    ref_str = run_dir / f"{name}.str"
    rmsd_dir = run_dir / "rmsd_clusters"
    ss_tsv = run_dir / "ss_clusters" / "ss_clusters.tsv"

    total_pool = count_frames(pool) if pool.exists() else 0

    ref_seq = ref_2d = ""
    if seq_file.exists():
        for ln in seq_file.read_text().splitlines():
            s = ln.strip()
            if s:
                ref_seq = s
                break
    if ref_str.exists():
        for ln in ref_str.read_text().splitlines():
            s = ln.strip()
            if s:
                ref_2d = s
                break

    rmsd_tables = rmsd_table(rmsd_dir) if rmsd_dir.exists() else {}
    ss_rows = ss_table(ss_tsv) if ss_tsv.exists() else []
    cutoffs_str = ", ".join(f"{c:g}" for c in sorted(rmsd_tables)) + " Å" if rmsd_tables else "?"

    for lang in LANGS:
        build(run_dir, lang, total_pool, ref_seq, ref_2d, rmsd_tables, ss_rows, cutoffs_str)


if __name__ == "__main__":
    main()
