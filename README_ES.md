# 2Dto3D-simRNA

Pipeline de principio a fin para pasar de una **secuencia de ARN + estructura secundaria en notación dot-bracket** a un **ensamble conformacional 3D** analizado, usando SimRNA con muestreo de Replica Exchange Monte Carlo a varias temperaturas (TREMD / REMC).

La versión en inglés de este documento está disponible en [README.md](README.md).

## Objetivo

Dado:
- una secuencia de ARN (una o más cadenas, con los códigos estándar A/C/G/U), y
- una estructura secundaria predicha en notación dot-bracket multilínea,

producir:
- una trayectoria REMC a múltiples temperaturas alrededor de la temperatura relativa estándar `T ≈ 1`,
- un conjunto de conformaciones 3D representativas (centroides de los clusters),
- una puntuación cuantitativa de cuánto cada ensamble muestreado reproduce la estructura 2D objetivo,
- gráficos diagnósticos para rankear variantes de diseño por estabilidad, diversidad y conservación de las regiones apareadas / no apareadas.

Los casos de uso objetivo son (a) **validación 3D** de diseños de ARN (¿las secuencias propuestas realmente se pliegan como se espera bajo el régimen de restricciones elegido?) y (b) **análisis de breathing** (¿cuánto fluctúa la estructura 3D alrededor del estado diseñado a temperaturas cercanas a las fisiológicas?).

## Por qué SimRNA

SimRNA es un motor de Monte Carlo coarse-grained desarrollado por GeneSilico que representa cada nucleótido con cinco pseudoátomos (dos de backbone: P, C4'; tres para la orientación de la base). Usa potenciales estadísticos derivados de estructuras de ARN resueltas y soporta ARN de cadena única y múltiple. El solvente es implícito — no es vacío, pero tampoco es agua explícita, lo cual puede introducir artefactos y motiva el protocolo de análisis de este repositorio.

El muestreo REMC corre varias réplicas en paralelo a distintas temperaturas relativas y las intercambia periódicamente para escapar de mínimos locales. La trayectoria a las temperaturas más bajas captura el fold diseñado; las réplicas a mayor temperatura exploran el paisaje de desplegamiento, lo que permite chequear tanto que el diseño es estable cerca de `T = 1` como que efectivamente se despliega bajo estrés térmico.

## Estructura del repositorio

```
.
├── inputs/
│   └── example/               # inputs de la corrida de referencia (versionado)
├── results/
│   └── example/               # resultados post-procesados de la corrida de referencia (versionado, sin .trafl)
├── CLAUDE.md                  # notas para sesiones de Claude Code sobre este repo
├── README.md                  # versión en inglés
└── README_ES.md               # este archivo
```

Solo la corrida de referencia bajo `inputs/example/` y `results/example/` está versionada. Cualquier otro directorio que crees dentro de `inputs/` o `results/` queda ignorado por git — podés usarlos libremente para tus propias corridas sin contaminar el repo. Los archivos de trayectoria pesados (`*.trafl`) se excluyen incluso para la corrida de referencia; el repo guarda solo los outputs post-procesados livianos (CSVs de scoring, PDBs de centroides, plots).

La **distribución académica de SimRNA** (binarios + potenciales estadísticos `data/`) **no** se distribuye desde este repo — la licencia no permite redistribución pública. Bajala aparte antes de correr el pipeline; ver [Instalación de SimRNA](#instalación-de-simrna) más abajo.

## Instalación de SimRNA

Descargá la distribución académica de SimRNA (este pipeline está calibrado contra v3.20 para Linux 64 bits) desde el sitio oficial de GeneSilico / IIMCB. La distribución incluye:

- `SimRNA` — motor principal MC / REMC
- `clustering` — clustering por RMSD de una trayectoria
- `SimRNA_trafl2pdbs` — extrae frames de un `.trafl` a PDB (opcionalmente con reconstrucción full-atom usando `AA`)
- `calc_rmsd_to_1st_frame` — RMSD por frame respecto al primero
- `traflView` — visualizador de trayectorias
- `trafl_extract_lowestE_frame.py` — helper **Python 2** que extrae el frame de menor energía
- `data/` — histogramas de los potenciales estadísticos, necesarios en runtime
- `SimRNA_academic_license.pdf` — leer antes de usar; la redistribución está restringida

Una vez descargada, agregá el directorio al `PATH` o invocá los binarios por path absoluto. Desde cualquier lugar donde lances una simulación, SimRNA necesita un directorio (o symlink) `data/` en el working directory — ver [Etapa 1](#etapa-1--preparar-el-directorio-de-corrida).

## Inputs

1. **Archivo de secuencia** — texto plano, nucleótidos `A/C/G/U`; múltiples cadenas separadas por espacios.
2. **Restricciones de estructura secundaria** (`*.str`, opcional pero usado en todo este pipeline) — string dot-bracket multilínea que coincide con la secuencia, una línea por cadena.
3. **Archivo de configuración** (`*.in` / `*.dat`) — parámetros de SimRNA. La distribución de SimRNA trae un `config.dat` / `config_n0.dat` mínimo al lado de los binarios. Valores típicos para TREMD:

   ```
   NUMBER_OF_ITERATIONS 16000000
   TRA_WRITE_IN_EVERY_N_ITERATIONS 16000
   INIT_TEMP 1.65
   FINAL_TEMP 0.90
   BONDS_WEIGHT 1.0
   ANGLES_WEIGHT 1.0
   TORS_ANGLES_WEIGHT 0.0      # deprecado, dejar en 0
   ETA_THETA_WEIGHT 0.4
   SECOND_STRC_RESTRAINTS_WEIGHT 0.3
   ```

4. **Simlink `data`** — SimRNA no arranca si no hay un directorio (o symlink) `data/` en el working directory. Es con diferencia la causa más común de fallos en las corridas. Apuntalo al `data/` que vino con tu descarga de SimRNA:

   ```bash
   ln -s /ruta/absoluta/a/tu_install_SimRNA/data data
   ```

5. **(Opcional) PDB de inicio** — con `-p` la secuencia se toma del PDB; con `-P` la ocupancia y el B-factor se interpretan como restricciones de posición átomo a átomo.

## Notas sobre los parámetros

- `NUMBER_OF_ITERATIONS` × `TRA_WRITE_IN_EVERY_N_ITERATIONS`: la razón determina cuántos frames (y en REMC, cuántos intentos de exchange) obtenés. Combinaciones típicas: 5M / 5k / 20 réplicas o 16M / 16k / 10 réplicas.
- `INIT_TEMP` / `FINAL_TEMP`: escala relativa — `T ≈ 1` es el punto fisiológico de operación. Para estructuras tipo TE el protocolo usa 16 réplicas cubriendo aproximadamente `0.9 → 1.65`. El extremo alto tiene que lograr desplegar la estructura, de lo contrario REMC no está haciendo su trabajo.
- `SECOND_STRC_RESTRAINTS_WEIGHT`: la perilla principal.
  - **Demasiado baja:** el ensamble cerca de `T = 1` se desvía de la estructura 2D objetivo.
  - **Demasiado alta:** la estructura sobrevive a las temperaturas más altas (REMC deja de funcionar) y la variabilidad en `T ≈ 1` se colapsa (sirve poco para análisis de breathing).

  Empezá alrededor de `0.05` para una corrida exploratoria, subí si el acuerdo con la estructura esperada es muy débil, bajá si las réplicas de alta T nunca despliegan. La variabilidad depende de este peso, así que corré **≥ 3 repeticiones con el peso elegido** para estimar el ruido.
- `TORS_ANGLES_WEIGHT` está ignorado en las versiones actuales de SimRNA y fue reemplazado por `ETA_THETA_WEIGHT`; dejalo en `0`.

## Pipeline

### Etapa 1 — Preparar el directorio de corrida

```bash
mkdir -p results/mi_corrida && cd results/mi_corrida
ln -s /ruta/abs/a/tu_install_SimRNA/data data
cp ../../inputs/mi_corrida/example            .     # secuencia
cp ../../inputs/mi_corrida/example.str        .     # restricciones dot-bracket
cp ../../inputs/mi_corrida/sim_config.in      .
```

> Solo `inputs/example/` y `results/example/` están versionados en git; todo lo demás dentro de esas dos carpetas queda local (ver [Estructura del repositorio](#estructura-del-repositorio)).

### Etapa 2 — Correr REMC (TREMD)

TREMD canónico de 16 réplicas a partir de secuencia + restricciones 2D:

```bash
nohup SimRNA -s example -S example.str -c sim_config.in -E 16 -o example >& example.log &
```

Alternativas:

```bash
# Partiendo de un PDB + restricciones 2D (MC a temperatura única):
nohup SimRNA -p example.pdb -S example.str -c sim_config.in -o example >& example.log &

# REMC desde secuencia sin restricciones (20 réplicas):
nohup SimRNA -s example -c sim_config.in -E 20 -o example >& example.log &
```

Outputs:

- `example_<NN>.trafl` — trayectoria por réplica (coarse-grained, un frame cada `TRA_WRITE_IN_EVERY_N_ITERATIONS` iteraciones)
- `example_<NN>-000001.pdb` — PDB del primer frame por réplica (se usa como template para `trafl2pdbs`)
- `example.ss_detected` — estructura secundaria detectada en la conformación inicial
- `example.log` — traza de energías y registro de intercambios entre réplicas

### Etapa 3 — Extraer frames de baja temperatura alrededor de `T ≈ 1`

Para la escalera estándar de 16 niveles `0.9 → 1.65`, las réplicas que quedan cerca de `T = 1` son típicamente los niveles `2`, `3` y `4` (temperaturas relativas `0.95`, `1.00`, `1.05`). Se juntan sus frames en una única trayectoria:

```bash
python simRNA/extract_low_temp_frames.py example.log \
  --base-name example --min-temp 2 --max-temp 4 --output around_1.trafl
```

### Etapa 4 — Clusterizar

Clusterizar los frames pooleados de baja temperatura con un cutoff de RMSD de 15 Å (ajustable según el sistema):

```bash
clustering around_1.trafl 1.0 15 >& clustering.log &
```

Esto genera archivos `around_1_thrs15.00A_clust<NN>.trafl`. Dentro de cada cluster los frames están ordenados por similitud al centroide; los RMSDs frame a centroide están en `clustering.log`.

### Etapa 5 — Extraer centroides / PDBs representativos

```bash
# Centroides de todos los clusters por encima de un cutoff de representación
python extract_major_clusters.py --cutoff <porcentaje>

# O el centroide (primer frame) de un cluster específico, reconstruido full-atom:
SimRNA_trafl2pdbs example_01-000001.pdb around_1_thrs15.00A_clustXX.trafl 1 AA
```

### Etapa 6 — Puntuar el acuerdo con la 2D

Comparar la estructura secundaria observada en el ensamble contra el dot-bracket objetivo:

```bash
# Scoring completo por temperatura + comparación entre niveles, agregando una fila al summary
python /ruta/a/Analyze_scorings.py \
  --name example --min-temp 1 --max-temp 16 \
  --structure example.str --label example \
  --summary ../summary.csv --threshold 0.8
```

Para cada temperatura emite tres scores: acuerdo de pares de bases, acuerdo de la secuencia completa (apareados + no apareados) y acuerdo solo de dots (solo no apareados). Los scores de dots y secuencia completa **no deben ir a cero a alta temperatura** — las regiones desapareadas son justamente lo que se espera en esas condiciones.

Comparar múltiples directorios (por ejemplo réplicas o variantes de diseño):

```bash
python /ruta/a/plot_scores_comparison.py summary.csv
```

Alternativa manual para el mismo paso:

```bash
python /ruta/a/extract_all_pdbs.py /ruta/a/dir --name example --suffixes 01,02,03,04,05
python /ruta/a/compare_ss.py example.str pdbs/
```

### Etapa 7 — Inspeccionar las trayectorias en ChimeraX

Combinar los PDBs frame-a-frame en un PDB multi-modelo que ChimeraX abre como trayectoria:

```bash
python multi-pdb.py --folder ./pdbs --output trajectory.pdb
```

## Flujo diagnóstico (qué mirar)

1. **Una corrida exploratoria** con un peso de restricciones bajo (`~0.05`), 16 réplicas, `0.9 → 1.65`.
2. Correr el scoring. Cerca de `T = 1` (niveles 2–4) el acuerdo de pares de bases tiene que ser alto pero no clavado en 1.0; en los dos niveles más altos (15–16) la estructura tiene que estar mayormente perdida.
3. Si la estructura objetivo no se recupera cerca de `T = 1` → **aumentar** `SECOND_STRC_RESTRAINTS_WEIGHT`.
4. Si la estructura sobrevive a las temperaturas más altas o la variabilidad en `T = 1` se colapsa → **bajar** el peso.
5. Una vez elegido el peso, repetir la corrida **≥ 3 veces** para cuantificar el ruido a ese peso.
6. Correr todas las variantes de diseño, hacer el análisis de ruptura y el whisker plot sobre pares de bases / dots / pares+dots. Pre-seleccionar devices por estabilidad, diversidad y conservación de las regiones desapareadas.
7. Clusterizar las variantes seleccionadas (cutoff de 15 Å como punto de partida) y tomar los centroides de los clusters por encima de un cutoff de representación elegido (por ej. 5%).
8. **(Opcional)** re-correr con un rango de temperatura ajustado para mayor sensibilidad.

## Aclaraciones

- **La solvatación es implícita.** La función de energía coarse-grained no ve el agua explícitamente; tratá las energías absolutas como relativas y validá las conclusiones con estadísticas a nivel de ensamble, no con frames individuales.
- **El largo de la secuencia** dispara el costo computacional de forma super-lineal; mantené las cadenas tan cortas como permita la pregunta biológica.
- **Los symlinks `data` rotos** son la causa principal de que este pipeline se reinicie en silencio — revisalos primero cuando una corrida no arranca o crashea de inmediato.
- **Los scripts Python de análisis downstream** (`extract_low_temp_frames.py`, `extract_major_clusters.py`, `Analyze_scorings.py`, `plot_scores_comparison.py`, `extract_all_pdbs.py`, `compare_ss.py`, `multi-pdb.py`) **todavía no están vendored en este repositorio** — viven en un árbol externo `simRNA/` y solo están referenciados por el protocolo documentado.

## Referencias

- **Manual de usuario de SimRNA v3.20** — viene con la distribución académica oficial como `SimRNA_UserManual_v3_20_20141002.pdf`.
- **Licencia académica de SimRNA** — viene como `SimRNA_academic_license.pdf` en la misma distribución; leerla antes de usar.
- **Reporte de método interno** (solo laboratorio, no está en este repo) — fuente de verdad para las decisiones de protocolo de arriba.
