library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(123)

# -----------------------------
# 1. Imaginary UMAP-like dataset
# -----------------------------

n <- 3000

umap_df <- tibble(
  UMAP_1 = c(
    rnorm(1000, mean = -3, sd = 0.8),
    rnorm(1000, mean =  0, sd = 0.9),
    rnorm(1000, mean =  3, sd = 0.7)
  ),
  UMAP_2 = c(
    rnorm(1000, mean =  1, sd = 0.8),
    rnorm(1000, mean = -1, sd = 0.9),
    rnorm(1000, mean =  1, sd = 0.7)
  ),
  celltype = rep(c("RG", "IPC", "Neuron"), each = 1000),
  subtype = sample(
    c("Subtype_A", "Subtype_B", "Subtype_C", "Subtype_D"),
    size = n,
    replace = TRUE
  )
)

# -----------------------------

plot_subset_umap <- function(df, subset_name) {
  
  subset_df <- df %>%
    filter(subtype == subset_name)
  
  ggplot() +
    # Perimeter / contour of the whole dataset
    geom_density_2d(
      data = df,
      aes(x = UMAP_1, y = UMAP_2),
      color = "grey60",
      linewidth = 0.4,
      bins = 6
    ) +
    
    # Points belonging only to this subset
    geom_point(
      data = subset_df,
      aes(x = UMAP_1, y = UMAP_2, color = celltype),
      size = 0.5,
      alpha = 0.8
    ) +
    
    coord_equal() +
    theme_classic() +
    labs(
      title = subset_name,
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Cell type"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}


subset_plots <- umap_df %>%
  pull(subtype) %>%
  unique() %>%
  sort() %>%
  lapply(function(x) plot_subset_umap(umap_df, x))

