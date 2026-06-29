# Genome Window Analysis

Adaptacion de `tfm_lia/scripts_anabel` para analizar FASTA completos de genomas de referencia sin partirlos previamente por cromosomas.

El flujo genera tres tablas finales por genoma:

- `<species>_dinucleotides_by_window.tsv`: valores por ventana para los 16 dinucleotidos.
- `<species>_bhlh_CANNTG_by_window.tsv`: valores secuenciales por ventana para los 16 motivos bHLH `CANNTG`.
- `<species>_bhlh_CAN_CAN_functional_by_window.tsv`: valores funcionales por ventana para grupos no orientados `CAN-CAN`.

## Uso en cluster

Desde esta carpeta:

```bash
sbatch run_one_genome.sbatch /ruta/al/genoma.fa.gz human /ruta/a/resultados 50000
```

El job usa un solo nodo (`--nodes=1`, `--ntasks=1`) y paraleliza dentro del nodo con `--cpus-per-task`. Cambia los recursos al lanzar o editando el encabezado del `.sbatch`.

Para usar ventanas solapadas se puede definir `STEP_SIZE` al lanzar. Por ejemplo, para reproducir el avance antiguo de ventanas de 50 kb con solape de 4 bp:

```bash
STEP_SIZE=49996 sbatch run_one_genome.sbatch /ruta/al/genoma.fa.gz human /ruta/a/resultados 50000
```

## Columnas principales

Ambas tablas incluyen `species`, `seq_id`, `window_index`, `start`, `end`, `window_length` y `acgt_fraction`.

La tabla de dinucleotidos incluye `dimer`, `observed_count`, `observed_frequency`, `expected_probability`, `expected_count` y `ratio_observed_expected`.

La tabla de bHLH incluye `motif`, `internal_dimer`, `observed_count`, `observed_frequency_all_hexamers`, `observed_frequency_bhlh`, `expected_probability_markov`, `expected_count_markov` y `ratio_observed_expected`.

La tabla funcional de bHLH agrupa secuencias equivalentes por el dimero funcional. Por ejemplo, `CATATG` es `CAT-CAT`, mientras que `CATCTG` y `CAGATG` se suman dentro de `CAT-CAG` porque representan el mismo heterodimero con orientacion opuesta. Incluye `functional_group`, `monomer_1`, `monomer_2`, `sequence_motifs`, `n_sequence_motifs`, `observed_count`, `observed_frequency_all_hexamers`, `expected_probability_markov`, `expected_count_markov` y `ratio_observed_expected`.

Por defecto se analiza tambien la cadena reversa complementaria, igual que hacia el script original al construir `fullwindow`.
