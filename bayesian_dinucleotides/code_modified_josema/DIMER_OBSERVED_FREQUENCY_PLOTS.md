# Observed Dinucleotide Frequency Plots

Estos plots muestran la frecuencia real observada de todos los grupos de dinucleotidos por ventana. Son distintos de los plots de `ratio_observed_expected`:

- `observed_frequency`: cuanto hay realmente de cada dinucleotido.
- `ratio_observed_expected`: si hay mas o menos de lo esperado por el modelo local de cada especie/ventana.

Para demostrar diferencias de composicion genómica absoluta, como que especies metiladoras tienen menos CG real, usa `observed_frequency`.

## Grupos de dimeros

Los dimeros se agrupan por reverso-complementario:

```text
CG, GC, CC/GG, AG/CT, GA/TC, CA/TG, AC/GT, AA/TT, TA, AT
```

## Uso en cluster

Desde `code_modified_josema`:

```bash
sbatch run_dimer_observed_frequency_plots.sbatch \
  /users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb \
  /users/genomics/josema/phd/bayesian_dinucleotides/code_modified_josema/dnmt_methylation_info_download_species_UPDATED.tsv \
  /users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb/dimer_observed_frequency_plots
```

El script lee cada `*_dinucleotides_by_window.tsv` por separado, guarda una tabla ligera y luego pinta desde esa tabla.

## Salidas

```text
dimer_observed_frequency_by_window.tsv
dimer_observed_frequency_by_window_with_metadata.tsv
dimer_observed_frequency_by_species.tsv
dimer_observed_frequency_density_by_species.pdf/png
dimer_observed_frequency_density_by_clade.pdf/png
dimer_observed_frequency_density_by_methylation.pdf/png
dimer_observed_frequency_violin_by_methylation.pdf/png
dimer_observed_frequency_species_medians.pdf/png
```

## Reutilizar la tabla ligera

```bash
Rscript plot_dimer_observed_frequency.R \
  --observed_table /ruta/dimer_observed_frequency_by_window.tsv \
  --methylation_table /users/genomics/josema/phd/bayesian_dinucleotides/code_modified_josema/dnmt_methylation_info_download_species_UPDATED.tsv \
  --output_dir /ruta/dimer_observed_frequency_plots
```

## Relacion con el plot solo de CG

`plot_cg_observed_frequency.R` es una version enfocada solo en CG. Este script hace lo mismo, pero para todos los grupos de dimeros.
