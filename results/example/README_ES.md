![](../../brand/assets/danslab-logo-report.png)

**DANSLAB · MMBAI** — Modelado Molecular, Bioinformática e IA · Universidad de la República, Salto, Uruguay

*[Read in English](README_EN.md) · [Leer en español](README_ES.md)*

# Subestados conformacionales — `example`

Análisis de clustering sobre el pool de frames de baja temperatura (niveles REMC 2–4, T≈1 en la escala relativa) de la corrida SimRNA de esta carpeta.

- **Frames en el pool low-T**: 30003
- **Secuencia**: `CCGCGAACGGUAACGUUCGCGG`
- **2D de referencia** (restricción SimRNA): `(((((((((....)))))))))`
- **Pipeline**: `scripts/run_cluster_populations.sh` → `extract_low_temp_frames.py` (--min-temp 2 --max-temp 4) → `clustering` (4, 8, 12, 16 Å) + `SimRNA_trafl2pdbs` + `ss_cluster.py`.

## RMSD clustering — top 10 por cutoff

### Cutoff 4.0 Å

| rank | cluster | frames | % |
|-----:|--------:|-------:|---:|
| 1 | 01 | 29832 | 99.44 |
| 2 | 02 | 69 | 0.23 |
| 3 | 03 | 68 | 0.23 |
| 4 | 04 | 10 | 0.03 |
| 5 | 05 | 8 | 0.03 |
| 6 | 06 | 3 | 0.01 |
| 7 | 07 | 2 | 0.01 |
| 8 | 08 | 2 | 0.01 |
| 9 | 09 | 1 | 0.00 |
| 10 | 10 | 1 | 0.00 |

_Top 10 cubren 29996 frames del pool (10 clusters reportados)._

### Cutoff 8.0 Å

| rank | cluster | frames | % |
|-----:|--------:|-------:|---:|
| 1 | 01 | 29999 | 100.00 |
| 2 | 02 | 1 | 0.00 |

_Top 10 cubren 30000 frames del pool (2 clusters reportados)._

### Cutoff 12.0 Å

| rank | cluster | frames | % |
|-----:|--------:|-------:|---:|
| 1 | 01 | 30000 | 100.00 |

_Top 10 cubren 30000 frames del pool (1 clusters reportados)._

### Cutoff 16.0 Å

| rank | cluster | frames | % |
|-----:|--------:|-------:|---:|
| 1 | 01 | 30000 | 100.00 |

_Top 10 cubren 30000 frames del pool (1 clusters reportados)._

## Clustering por patrón 2D (`.ss_detected`) — top 10

Cada subestado corresponde a un **patrón dot-bracket exacto**. La población es la fracción de frames del pool low-T que exhiben esa 2D.

| rank | 2D (dot-bracket) | frames | % |
|-----:|:-----------------|-------:|---:|
| 1 | `(((((((((....)))))))))` | 28968 | 96.55 |
| 2 | `((((((((......))))))))` | 234 | 0.78 |
| 3 | `.((((((((....)))))))).` | 211 | 0.70 |
| 4 | `(((((.(((....))).)))))` | 175 | 0.58 |
| 5 | `((((((.((....)).))))))` | 139 | 0.46 |
| 6 | `((((((.((....))).)))))` | 66 | 0.22 |
| 7 | `(((((((.(....).)))))))` | 61 | 0.20 |
| 8 | `((((.((((....)))).))))` | 37 | 0.12 |
| 9 | `(.(((((((....))))))).)` | 21 | 0.07 |
| 10 | `(((.(((((....))))).)))` | 14 | 0.05 |

_Archivos crudos: `low_temp.trafl`, `rmsd_clusters/`, `ss_clusters/`. Todos excluidos por `.gitignore` salvo estos reportes._

---

_DansLab · MMBAI — Universidad de la República, Salto, Uruguay. Generado por los builders de reportes de este repo; editar el generador, no el Markdown._
