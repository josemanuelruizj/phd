# Observed CG Frequency Plots

Estos plots estan pensados para responder a la pregunta:

```text
Las especies que metilan tienen menos CG real en el genoma?
```

Para esa pregunta se usa `observed_frequency` de las filas con `dimer == "CG"`, no `ratio_observed_expected`.

- `observed_frequency`: cantidad/frecuencia real de CG por ventana.
- `ratio_observed_expected`: deplecion relativa frente a lo esperado por el modelo local.

Si una especie ya tiene pocos CG, el modelo puede esperar pocos CG y el ratio puede acercarse a 1. Por eso, para demostrar menos CG en el genoma, la figura principal deberia ser `observed_frequency`.

## Uso en cluster

Desde `code_modified_josema`:

```bash
sbatch run_cg_observed_frequency_plots.sbatch \
  /users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb \
  /users/genomics/josema/Downloads/dnmt_methylation_info_download_species_UPDATED.tsv \
  /users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb/cg_observed_frequency_plots
```

El script lee cada `*_dinucleotides_by_window.tsv` por separado, se queda solo con `dimer == "CG"`, guarda una tabla ligera y luego pinta desde esa tabla.

## Salidas

```text
cg_observed_frequency_by_window.tsv
cg_observed_frequency_by_window_with_metadata.tsv
cg_observed_frequency_by_species.tsv
cg_observed_frequency_density_by_species.pdf/png
cg_observed_frequency_density_by_clade.pdf/png
cg_observed_frequency_density_by_methylation.pdf/png
cg_observed_frequency_violin_by_methylation.pdf/png
cg_observed_frequency_species_medians.pdf/png
```

## Reutilizar la tabla ligera

Si ya existe `cg_observed_frequency_by_window.tsv`, puedes rehacer plots sin volver a leer las tablas grandes:

```bash
Rscript plot_cg_observed_frequency.R \
  --cg_table /ruta/cg_observed_frequency_by_window.tsv \
  --methylation_table /ruta/dnmt_methylation_info_download_species_UPDATED.tsv \
  --output_dir /ruta/cg_observed_frequency_plots
```
