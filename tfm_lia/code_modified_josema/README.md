# Genome Window Analysis

Adaptacion de `tfm_lia/scripts_anabel` para analizar FASTA completos de genomas de referencia sin partirlos previamente por cromosomas.

El flujo se divide en dos pasos:

1. `genome_window_analysis.R`: calculo pesado por ventanas.
2. `summarize_functional_windows.R`: seleccion de regiones funcionales a partir de la tabla ya calculada.

El primer paso genera tres tablas principales por genoma:

- `<species>_dinucleotides_by_window.tsv`: valores por ventana para los 16 dinucleotidos.
- `<species>_bhlh_CANNTG_by_window.tsv`: valores secuenciales por ventana para los 16 motivos bHLH `CANNTG`.
- `<species>_bhlh_CAN_CAN_functional_by_window.tsv`: valores funcionales por ventana para grupos no orientados `CAN-CAN`.

El segundo paso genera tablas resumen a partir de los grupos funcionales sin recalcular el genoma:

- `<species>_functional_topN_by_log2_pseudocount.tsv`: top N regiones sobre e infra representadas por cada grupo funcional.
- `<species>_functional_topN_by_depletion_score.tsv`: top N regiones infrarepresentadas por diferencia absoluta esperado-observado.
- `<species>_functional_top<X>percent_by_log2_pseudocount.tsv`: top X% mas positivo y mas negativo por grupo funcional.
- `<species>_functional_log2_pseudocount_threshold_<X>.tsv`: regiones con `abs(log2_ratio_pseudocount) >= X`.
- `<species>_functional_padj_significant_<X>.tsv`: regiones significativas por p-value Poisson ajustado.
- `<species>_functional_zero_observed_high_expected.tsv`: regiones con `observed_count = 0` y esperado alto.

## Uso en cluster

### 1. Calcular probabilidades por ventanas

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

### 2. Crear summaries funcionales

Una vez generada la tabla funcional, los summaries se lanzan aparte:

```bash
sbatch run_functional_summaries.sbatch \
  /ruta/a/resultados/human_bhlh_CAN_CAN_functional_by_window.tsv \
  /ruta/a/resultados/summaries
```

Los summaries usan por defecto `TOP_N=100`, `TOP_PERCENT=1`, `LOG2_THRESHOLD=1`, `PVALUE_THRESHOLD=0.05`, `MIN_EXPECTED_ZERO=1` y `PSEUDOCOUNT=0.5`. Se pueden cambiar al lanzar:

```bash
TOP_N=50 LOG2_THRESHOLD=2 TOP_PERCENT=5 PVALUE_THRESHOLD=0.01 MIN_EXPECTED_ZERO=5 \
  sbatch run_functional_summaries.sbatch \
  /ruta/a/resultados/human_bhlh_CAN_CAN_functional_by_window.tsv \
  /ruta/a/resultados/summaries
```

## Que threshold afecta a cada summary

`TOP_N` afecta solo a:

- `<species>_functional_topN_by_log2_pseudocount.tsv`
- `<species>_functional_topN_by_depletion_score.tsv`

En `functional_topN_by_log2_pseudocount.tsv`, para cada `functional_group`, se escriben dos bloques en el mismo archivo:

- `representation = over`: las `TOP_N` ventanas con `log2_ratio_pseudocount` mas alto.
- `representation = under`: las `TOP_N` ventanas con `log2_ratio_pseudocount` mas bajo.

En `functional_topN_by_depletion_score.tsv`, se escriben solo regiones infrarepresentadas:

- `representation = under`
- ordenadas por `depletion_score = expected_count_markov - observed_count`

Esta tabla esta pensada para encontrar regiones donde faltan muchos sitios en terminos absolutos, aunque el log2 no sea extremo.

`TOP_PERCENT` afecta solo a:

- `<species>_functional_top<X>percent_by_log2_pseudocount.tsv`

Para cada `functional_group`, coge el X% mas positivo y el X% mas negativo segun `log2_ratio_pseudocount`. Igual que el top N, contiene ambos sentidos en el mismo archivo:

- `representation = over`
- `representation = under`

`LOG2_THRESHOLD` afecta solo a:

- `<species>_functional_log2_pseudocount_threshold_<X>.tsv`

Esta tabla usa el tamano del efecto, no el p-value. Incluye regiones que cumplen:

```text
abs(log2_ratio_pseudocount) >= LOG2_THRESHOLD
```

Por tanto incluye ambos sentidos:

- `representation = over` si `log2_ratio_pseudocount >= LOG2_THRESHOLD`
- `representation = under` si `log2_ratio_pseudocount <= -LOG2_THRESHOLD`

`PVALUE_THRESHOLD` afecta solo a:

- `<species>_functional_padj_significant_<X>.tsv`

Esta tabla usa p-values Poisson ajustados por BH. Incluye dos tipos de filas en el mismo archivo:

- `representation = over`: filas con `padj_over <= PVALUE_THRESHOLD`
- `representation = under`: filas con `padj_under <= PVALUE_THRESHOLD`

Es decir, no elige solo over ni solo under: mete los dos si pasan su p-value correspondiente.

`MIN_EXPECTED_ZERO` afecta solo a:

- `<species>_functional_zero_observed_high_expected.tsv`

Esta tabla esta pensada para rescatar ausencias biologicamente interesantes. Incluye solo:

```text
observed_count == 0
expected_count_markov >= MIN_EXPECTED_ZERO
```

Todas sus filas son:

```text
representation = under
```

Se ordena por `expected_count_markov` de mayor a menor, porque una ventana con `observed = 0` y `expected = 10` es mucho mas interesante que una con `observed = 0` y `expected = 0.01`.

`PSEUDOCOUNT` afecta a las tablas que usan `log2_ratio_pseudocount`:

- `<species>_functional_topN_by_log2_pseudocount.tsv`
- `<species>_functional_top<X>percent_by_log2_pseudocount.tsv`
- `<species>_functional_log2_pseudocount_threshold_<X>.tsv`

No afecta a los p-values Poisson ni a `depletion_score`.

## Por que hay p-value over y p-value under

Hay dos p-values porque son dos preguntas estadisticas distintas.

`p_over` pregunta:

```text
Si esperaba expected_count_markov sitios, cual es la probabilidad de observar observed_count o mas?
```

Sirve para detectar sobre-representacion. Ejemplo: esperabas 2 y observas 10.

`p_under` pregunta:

```text
Si esperaba expected_count_markov sitios, cual es la probabilidad de observar observed_count o menos?
```

Sirve para detectar infrarepresentacion. Ejemplo: esperabas 10 y observas 0.

Por eso el summary de p-value usa:

- `padj_over` para filas `representation = over`
- `padj_under` para filas `representation = under`

Ambos p-values se calculan para todas las filas, pero cada direccion usa el suyo.

## Columnas principales

Ambas tablas incluyen `species`, `seq_id`, `window_index`, `start`, `end`, `window_length` y `acgt_fraction`.

La tabla de dinucleotidos incluye `dimer`, `observed_count`, `observed_frequency`, `expected_probability`, `expected_count` y `ratio_observed_expected`.

La tabla de bHLH incluye `motif`, `internal_dimer`, `observed_count`, `observed_frequency_all_hexamers`, `observed_frequency_bhlh`, `expected_probability_markov`, `expected_count_markov` y `ratio_observed_expected`.

La tabla funcional de bHLH agrupa secuencias equivalentes por el dimero funcional. Por ejemplo, `CATATG` es `CAT-CAT`, mientras que `CATCTG` y `CAGATG` se suman dentro de `CAT-CAG` porque representan el mismo heterodimero con orientacion opuesta. Incluye `functional_group`, `monomer_1`, `monomer_2`, `sequence_motifs`, `n_sequence_motifs`, `observed_count`, `observed_frequency_all_hexamers`, `expected_probability_markov`, `expected_count_markov` y `ratio_observed_expected`.

Todas las tablas principales incluyen tambien:

- `log2_ratio_observed_expected`: log2 del ratio crudo. Si `observed_count = 0` y `expected_count > 0`, el valor sera `-Inf`.
- `log2_ratio_pseudocount`: log2 estable calculado como `log2((observed_count + 0.5) / (expected_count + 0.5))`.
- `depletion_score`: `expected_count - observed_count`.
- `enrichment_score`: `observed_count - expected_count`.
- `is_zero_observed`: marca regiones con `observed_count = 0` y esperado mayor que 0.
- `p_under`: p-value Poisson para infrarepresentacion.
- `p_over`: p-value Poisson para sobrerepresentacion.

El script de summaries calcula ademas `padj_under`, `padj_over`, `padj_under_by_group` y `padj_over_by_group` con correccion BH.

Para regiones con `observed_count = 0`, usa `functional_zero_observed_high_expected.tsv` y ordena por `expected_count_markov`: asi no se pierden regiones interesantes solo porque el ratio crudo sea 0 y el log2 crudo sea `-Inf`.

Por defecto se analiza tambien la cadena reversa complementaria, igual que hacia el script original al construir `fullwindow`.
