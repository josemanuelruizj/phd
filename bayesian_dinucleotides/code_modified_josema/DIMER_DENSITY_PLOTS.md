# Dimer Ratio Density Plots

Este flujo genera plots de densidad para las tablas `*_dinucleotides_by_window.tsv` sin cargar todas las tablas crudas a la vez.

El script hace dos pasos:

1. Lee cada tabla de dimeros por separado y calcula una densidad por especie y grupo de dinucleotidos.
2. Guarda una tabla ligera de densidades y genera tres plots conjuntos.

Los grupos de dinucleotidos se agrupan por reverso-complementario:

```text
CG, GC, CC/GG, AG/CT, GA/TC, CA/TG, AC/GT, AA/TT, TA, AT
```

## Uso en cluster

Desde `code_modified_josema`:

```bash
sbatch run_dimer_density_plots.sbatch \
  /users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb \
  /users/genomics/josema/phd/bayesian_dinucleotides/code_modified_josema/dnmt_methylation_info_download_species_UPDATED.tsv \
  /users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb/dimer_density_plots
```

Si la tabla de metilacion esta en otra ruta, cambia el segundo argumento. El script espera una tabla TSV con, como minimo:

```text
species
taxonomic_group_from_figure
methylation_CpG_OE
```

## Salidas

En la carpeta de salida se generan:

```text
dimer_ratio_density_values.tsv
dimer_ratio_density_values_with_metadata.tsv
dimer_ratio_density_by_species.pdf
dimer_ratio_density_by_species.png
dimer_ratio_density_by_clade.pdf
dimer_ratio_density_by_clade.png
dimer_ratio_density_by_methylation.pdf
dimer_ratio_density_by_methylation.png
```

`dimer_ratio_density_values.tsv` es la tabla ligera calculada desde los ficheros crudos. Si quieres cambiar solo el estilo del plot mas adelante, puedes reutilizarla con:

```bash
Rscript plot_dimer_ratio_densities.R \
  --density_table /ruta/dimer_ratio_density_values.tsv \
  --methylation_table /users/genomics/josema/phd/bayesian_dinucleotides/code_modified_josema/dnmt_methylation_info_download_species_UPDATED.tsv \
  --output_dir /ruta/dimer_density_plots
```

## Ajustes utiles

Cambiar rango del eje X:

```bash
X_MAX=3 sbatch run_dimer_density_plots.sbatch input_dir methylation.tsv output_dir
```

Cambiar la columna que se pinta:

```bash
RATIO_COL=log2_ratio_pseudocount sbatch run_dimer_density_plots.sbatch input_dir methylation.tsv output_dir
```

Por defecto usa `ratio_observed_expected` y pinta el rango `0..2.2`.
