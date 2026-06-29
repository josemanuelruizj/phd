library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)

# -----------------------------
# 1. Read files
# -----------------------------

expr <- read_excel("Downloads/summary_table_average_expr_cells.xlsx")
ages <- read_excel("Downloads/converted_ages_long.xlsx")

# -----------------------------
# 2. Define genes to plot
# -----------------------------

genes_all <- c(
  "CNR1",
  "CNR2",
  "CNRIP1",
  "DAGLA",
  "DAGLB",
  "FAAH",
  "MGLL",
  "NAPEPLD"
)

# If you want only the rest of genes, excluding CNR1:
genes_to_plot <- setdiff(genes_all, "CNR1")

# If you want all genes including CNR1, use this instead:
# genes_to_plot <- genes_all

# -----------------------------
# 3. Clean / harmonise expression table
# -----------------------------

expr_clean <- expr %>%
  mutate(
    DATASET = toupper(DATASET),
    Common_name = as.character(Common_name),
    age_join = AGE_combined %>%
      as.character() %>%
      str_remove("^.*_") %>%
      str_replace_all("\\.", ",")
  )

# -----------------------------
# 4. Clean / harmonise age table
# -----------------------------

ages_clean <- ages %>%
  mutate(
    DATASET = toupper(Dataset),
    age_join = `Real age` %>%
      as.character() %>%
      str_replace_all("\\.", ","),
    human_age = `Converted (heterochrony or days)`
  ) %>%
  select(DATASET, age_join, human_age)

# -----------------------------
# 5. Pivot expression table to long format
# -----------------------------

expr_long <- expr_clean %>%
  select(
    DATASET,
    Common_name,
    AGE_combined,
    age_join,
    n_cells,
    matches("_(avg_expression|percent_expressing)$")
  ) %>%
  pivot_longer(
    cols = matches("_(avg_expression|percent_expressing)$"),
    names_to = c("gene", "metric"),
    names_pattern = "^(.*)_(avg_expression|percent_expressing)$",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  ) %>%
  filter(gene %in% genes_all)

# -----------------------------
# 6. Join with converted ages
# -----------------------------

min_cells <- 25

plot_data_all <- expr_long %>%
  left_join(ages_clean, by = c("DATASET", "age_join")) %>%
  filter(!is.na(human_age)) %>%
  mutate(
    enough_cells = n_cells >= min_cells
  )

# -----------------------------
# 7. Keep only dorsal lineage
# -----------------------------

dorsal_lineage_order <- c(
  "Dorsal_NSC(SOX2+)",
  "Excit_IPC",
  "Glut_NEU"
)

plot_data_all <- plot_data_all %>%
  filter(Common_name %in% dorsal_lineage_order)

# -----------------------------
# 8. Normalize expression per dataset and gene
#    Only groups with enough cells are used for the max
# -----------------------------

plot_data_all <- plot_data_all %>%
  group_by(DATASET, gene) %>%
  mutate(
    max_expr_valid = max(avg_expression[enough_cells], na.rm = TRUE),
    max_expr_valid = ifelse(is.infinite(max_expr_valid), NA_real_, max_expr_valid),
    avg_expression_norm = ifelse(
      enough_cells & !is.na(max_expr_valid) & max_expr_valid > 0,
      avg_expression / max_expr_valid,
      NA_real_
    )
  ) %>%
  ungroup()

# -----------------------------
# 9. Order datasets and dorsal cell types
# -----------------------------

dataset_order <- c(
  "MANNO",
  "DIBELLA",
  "MICALI",
  "EZE",
  "BRAUN",
  "TREVINO"
)

plot_data_all <- plot_data_all %>%
  mutate(
    DATASET = factor(DATASET, levels = rev(dataset_order)),
    Common_name = factor(Common_name, levels = dorsal_lineage_order)
  ) %>%
  filter(
    !is.na(DATASET),
    !is.na(Common_name)
  )

# -----------------------------
# 10. Labels
# -----------------------------

celltype_labels <- c(
  "Dorsal_NSC(SOX2+)" = "Dorsal NSC (SOX2+)",
  "Excit_IPC" = "Excit IPC",
  "Glut_NEU" = "Glut NEU"
)

dataset_labels <- c(
  "MANNO"   = "La Manno et al. (2021)",
  "DIBELLA" = "Di Bella et al. (2021)",
  "MICALI"  = "Micali et al. (2023)",
  "EZE"     = "Eze et al. (2021)",
  "BRAUN"   = "Braun et al. (2023)",
  "TREVINO" = "Trevino et al. (2021)"
)

important_ages <- tibble::tibble(
  human_age = c(28, 56, 84, 112, 140, 168),
  label = c("4 PCW", "8 PCW", "12 PCW", "16 PCW", "20 PCW", "24 PCW")
)

# -----------------------------
# 11. Plot function
# -----------------------------

plot_gene_dorsal <- function(gene_name) {
  
  plot_data_gene <- plot_data_all %>%
    filter(gene == gene_name)
  
  plot_data_included <- plot_data_gene %>%
    filter(enough_cells)
  
  plot_data_excluded <- plot_data_gene %>%
    filter(!enough_cells)
  
  ggplot() +
    geom_vline(
      data = important_ages,
      aes(xintercept = human_age),
      linetype = "dashed",
      linewidth = 0.25,
      color = "grey75",
      alpha = 0.8
    ) +
    geom_point(
      data = plot_data_included,
      aes(
        x = human_age,
        y = DATASET,
        size = percent_expressing,
        fill = avg_expression_norm
      ),
      shape = 21,
      color = "grey15",
      stroke = 0.18,
      alpha = 0.95
    ) +
    geom_point(
      data = plot_data_excluded,
      aes(
        x = human_age,
        y = DATASET
      ),
      shape = 4,
      size = 2.4,
      stroke = 0.7,
      color = "grey45",
      alpha = 0.8
    ) +
    facet_wrap(
      ~ Common_name,
      ncol = 3,
      labeller = labeller(Common_name = celltype_labels)
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      direction = 1,
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "Mean expression\nscaled 0–1"
    ) +
    scale_size_continuous(
      range = c(1.2, 7),
      breaks = c(25, 50, 75, 100),
      limits = c(0, 100),
      name = paste0(gene_name, "+ cells (%)")
    ) +
    scale_x_continuous(
      breaks = seq(20, 170, by = 20),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_discrete(
      labels = dataset_labels
    ) +
    coord_cartesian(
      xlim = c(20, 170),
      clip = "off"
    ) +
    labs(
      x = "Human-equivalent developmental age (days)",
      y = NULL,
      title = paste0(gene_name, " expression dynamics across dorsal cortical lineages"),
      subtitle = paste0(
        "Crosses indicate groups with less than ", min_cells, " cells."
      )
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "left",
        title.hjust = 0.5,
        title.vjust = 0.5,
        barwidth = unit(4.5, "cm"),
        barheight = unit(0.35, "cm"),
        ticks = TRUE
      ),
      size = guide_legend(
        title.position = "left",
        title.hjust = 0.5,
        title.vjust = 0.5,
        direction = "horizontal",
        override.aes = list(
          fill = "grey70",
          color = "grey15",
          alpha = 0.95
        )
      )
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 9.5,
        color = "grey25",
        hjust = 0,
        margin = margin(t = 4, b = 8)
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 10,
        margin = margin(b = 5)
      ),
      axis.text.x = element_text(
        size = 9,
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(
        size = 9,
        color = "grey15"
      ),
      axis.title.x = element_text(
        size = 10,
        margin = margin(t = 8)
      ),
      axis.line = element_line(
        linewidth = 0.35,
        color = "grey20"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        color = "grey20"
      ),
      panel.spacing = unit(1.1, "lines"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.6, "cm"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

# -----------------------------
# 12. Generate plots
# -----------------------------

plots_dorsal <- map(
  genes_to_plot,
  plot_gene_dorsal
)

names(plots_dorsal) <- genes_to_plot

# Example: show one plot
plots_dorsal

