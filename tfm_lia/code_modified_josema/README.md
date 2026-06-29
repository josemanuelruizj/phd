# Genome Window Analysis

Adaptacion de `tfm_lia/scripts_anabel` para analizar FASTA completos de genomas de referencia sin partirlos previamente por cromosomas.

El flujo genera tres tablas principales por genoma:

- `<species>_dinucleotides_by_window.tsv`: valores por ventana para los 16 dinucleotidos.
- `<species>_bhlh_CANNTG_by_window.tsv`: valores secuenciales por ventana para los 16 motivos bHLH `CANNTG`.
- `<species>_bhlh_CAN_CAN_functional_by_window.tsv`: valores funcionales por ventana para grupos no orientados `CAN-CAN`.

Tambien genera tres tablas resumen a partir de los grupos funcionales:

- `<species>_functional_CAN_CAN_top_<N>_over_under_by_group.tsv`: top N regiones sobre e infra representadas por cada grupo funcional en un unico archivo.
- `<species>_functional_CAN_CAN_log2_threshold_<X>.tsv`: regiones con `abs(log2_ratio_observed_expected) >= X`.
- `<species>_functional_CAN_CAN_top_<X>percent_over_under_by_group.tsv`: top X% de regiones mas positivas y mas negativas por cada grupo funcional.

## Uso en cluster

Desde esta carpeta:

```bash
sbatch run_one_genome.sbatch /ruta/al/genoma.fa.gz /ruta/a/resultados 50000
```

El nombre de especie se toma automaticamente del nombre del FASTA. Por ejemplo, `/ruta/hg38.fa.gz` genera salidas con prefijo `hg38`.

El job usa un solo nodo (`--nodes=1`, `--ntasks=1`) y paraleliza dentro del nodo con `--cpus-per-task`. Cambia los recursos al lanzar o editando el encabezado del `.sbatch`.

Por defecto se usa un solape de 4 bp entre ventanas, igual que el avance antiguo. Con ventanas de 50 kb, el script calcula automaticamente `STEP_SIZE=49996`.

Para cambiar el solape, pasalo como cuarto argumento:

```bash
sbatch run_one_genome.sbatch /ruta/al/genoma.fa.gz /ruta/a/resultados 50000 10
```

Tambien puedes forzar una etiqueta de especie o un `STEP_SIZE` concreto con variables de entorno:

```bash
SPECIES=human STEP_SIZE=49996 sbatch run_one_genome.sbatch /ruta/al/genoma.fa.gz /ruta/a/resultados 50000
```

Los resumenes funcionales usan por defecto `TOP_N=100`, `LOG2_THRESHOLD=1` y `TOP_PERCENT=1`. Se pueden cambiar al lanzar:

```bash
TOP_N=50 LOG2_THRESHOLD=2 TOP_PERCENT=5 sbatch run_one_genome.sbatch /ruta/al/genoma.fa.gz /ruta/a/resultados 50000
```

## Columnas principales

Ambas tablas incluyen `species`, `seq_id`, `window_index`, `start`, `end`, `window_length` y `acgt_fraction`.

La tabla de dinucleotidos incluye `dimer`, `observed_count`, `observed_frequency`, `expected_probability`, `expected_count` y `ratio_observed_expected`.

La tabla de bHLH incluye `motif`, `internal_dimer`, `observed_count`, `observed_frequency_all_hexamers`, `observed_frequency_bhlh`, `expected_probability_markov`, `expected_count_markov` y `ratio_observed_expected`.

La tabla funcional de bHLH agrupa secuencias equivalentes por el dimero funcional. Por ejemplo, `CATATG` es `CAT-CAT`, mientras que `CATCTG` y `CAGATG` se suman dentro de `CAT-CAG` porque representan el mismo heterodimero con orientacion opuesta. Incluye `functional_group`, `monomer_1`, `monomer_2`, `sequence_motifs`, `n_sequence_motifs`, `observed_count`, `observed_frequency_all_hexamers`, `expected_probability_markov`, `expected_count_markov` y `ratio_observed_expected`.

Todas las tablas principales incluyen tambien `log2_ratio_observed_expected`: valores negativos indican menos ocurrencias de las esperadas y valores positivos indican mas ocurrencias de las esperadas. Si `observed_count = 0` y `expected_count > 0`, el valor sera `-Inf`.

Por defecto se analiza tambien la cadena reversa complementaria, igual que hacia el script original al construir `fullwindow`.
