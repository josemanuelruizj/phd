library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)
library(tidyverse)
library(ggh4x)

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

genes_to_plot <- genes_all

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
    human_age = `Converted (heterochrony or days)`,
    human_age_pcw = human_age / 7
  ) %>%
  select(DATASET, age_join, human_age, human_age_pcw)

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
  filter(!is.na(human_age_pcw)) %>%
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
# 9. Order datasets, cell types and species groups
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
    Common_name = factor(Common_name, levels = dorsal_lineage_order),
    SPECIES_GROUP = case_when(
      DATASET %in% c("MANNO", "DIBELLA") ~ "Mouse",
      DATASET %in% c("MICALI") ~ "Macaque",
      DATASET %in% c("EZE", "BRAUN", "TREVINO") ~ "Human",
      TRUE ~ NA_character_
    ),
    SPECIES_GROUP = factor(
      SPECIES_GROUP,
      levels = c("Mouse", "Macaque", "Human")
    )
  ) %>%
  filter(
    !is.na(DATASET),
    !is.na(Common_name),
    !is.na(SPECIES_GROUP)
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
  human_age_pcw = seq(4, 24, by = 2),
  label = paste0(seq(4, 24, by = 2), " PCW")
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
      aes(xintercept = human_age_pcw),
      linetype = "dashed",
      linewidth = 0.22,
      color = "grey80",
      alpha = 0.75
    ) +
    geom_point(
      data = plot_data_included,
      aes(
        x = human_age_pcw,
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
        x = human_age_pcw,
        y = DATASET
      ),
      shape = 4,
      size = 2.4,
      stroke = 0.7,
      color = "grey45",
      alpha = 0.8
    ) +
    ggh4x::facet_grid2(
      rows = vars(SPECIES_GROUP),
      cols = vars(Common_name),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      labeller = labeller(
        Common_name = celltype_labels
      ),
      strip = ggh4x::strip_themed(
        background_x = list(
          ggplot2::element_rect(
            fill = "#CFE8F3",
            color = "grey55",
            linewidth = 0.4
          ),
          ggplot2::element_rect(
            fill = "#74A9CF",
            color = "grey55",
            linewidth = 0.4
          ),
          ggplot2::element_rect(
            fill = "#045A8D",
            color = "grey55",
            linewidth = 0.4
          )
        ),
        text_x = list(
          ggplot2::element_text(
            color = "black",
            face = "bold",
            size = 13
          ),
          ggplot2::element_text(
            color = "black",
            face = "bold",
            size = 13
          ),
          ggplot2::element_text(
            color = "white",
            face = "bold",
            size = 13
          )
        ),
        background_y = list(
          ggplot2::element_blank(),
          ggplot2::element_blank(),
          ggplot2::element_blank()
        ),
        text_y = list(
          ggplot2::element_text(
            color = "grey25",
            face = "bold",
            size = 8.5,
            angle = 0
          ),
          ggplot2::element_text(
            color = "grey25",
            face = "bold",
            size = 8.5,
            angle = 0
          ),
          ggplot2::element_text(
            color = "grey25",
            face = "bold",
            size = 8.5,
            angle = 0
          )
        )
      )
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
      breaks = seq(4, 24, by = 2),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_discrete(
      labels = dataset_labels,
      drop = TRUE
    ) +
    coord_cartesian(
      xlim = c(3, 25),
      clip = "off"
    ) +
    labs(
      x = "Human-equivalent developmental age (post-conception weeks)",
      y = NULL,
      title = paste0(gene_name, " expression dynamics across dorsal cortical lineages"),
      subtitle = paste0(
        "Datasets are grouped by species. Crosses indicate groups with fewer than ",
        min_cells, " cells."
      )
    ) +
    guides(
      fill = guide_colorbar(
        order = 1,
        title.position = "left",
        title.hjust = 0.5,
        title.vjust = 0.5,
        barwidth = unit(4.5, "cm"),
        barheight = unit(0.35, "cm"),
        ticks = TRUE
      ),
      size = guide_legend(
        order = 2,
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
        margin = ggplot2::margin(t = 4, b = 8)
      ),
      
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 13
      ),
      strip.placement = "outside",
      
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
        margin = ggplot2::margin(t = 8)
      ),
      axis.line = element_line(
        linewidth = 0.35,
        color = "grey20"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        color = "grey20"
      ),
      panel.spacing.x = unit(1.1, "lines"),
      panel.spacing.y = unit(1.4, "lines"),
      
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.6, "cm"),
      
      plot.margin = ggplot2::margin(8, 10, 8, 8)
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

# -----------------------------
# 13. Show example plots
# -----------------------------

plots_dorsal[["CNR1"]]
plots_dorsal[["CNR2"]]
plots_dorsal[["NAPEPLD"]]

# -----------------------------
# 14. Save all plots
# -----------------------------

dir.create("gene_dorsal_plots", showWarnings = FALSE)

walk(
  genes_to_plot,
  function(gene_name) {
    
    ggsave(
      filename = file.path('Downloads',
        "gene_dorsal_plots",
        paste0(gene_name, "_dorsal_lineage_dotplot_species_grouped.pdf")
      ),
      plot = plots_dorsal[[gene_name]],
      width = 18,
      height = 4,
      units = "in",
      device = cairo_pdf
    )
    
  }
)



library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)
library(tidyverse)
library(ggh4x)

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

genes_to_plot <- genes_all

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
    human_age = `Converted (heterochrony or days)`,
    human_age_pcw = human_age / 7
  ) %>%
  select(DATASET, age_join, human_age, human_age_pcw)

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
  filter(!is.na(human_age_pcw)) %>%
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
# 9. Order datasets, cell types and species groups
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
    Common_name = factor(Common_name, levels = dorsal_lineage_order),
    SPECIES_GROUP = case_when(
      DATASET %in% c("MANNO", "DIBELLA") ~ "Mouse",
      DATASET %in% c("MICALI") ~ "Macaque",
      DATASET %in% c("EZE", "BRAUN", "TREVINO") ~ "Human",
      TRUE ~ NA_character_
    ),
    SPECIES_GROUP = factor(
      SPECIES_GROUP,
      levels = c("Mouse", "Macaque", "Human")
    )
  ) %>%
  filter(
    !is.na(DATASET),
    !is.na(Common_name),
    !is.na(SPECIES_GROUP)
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
  human_age_pcw = seq(4, 24, by = 2),
  label = paste0(seq(4, 24, by = 2), " PCW")
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
      aes(xintercept = human_age_pcw),
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey80",
      alpha = 0.75
    ) +
    geom_point(
      data = plot_data_included,
      aes(
        x = human_age_pcw,
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
        x = human_age_pcw,
        y = DATASET
      ),
      shape = 4,
      size = 2.4,
      stroke = 0.7,
      color = "grey45",
      alpha = 0.8
    ) +
    ggh4x::facet_grid2(
      rows = vars(SPECIES_GROUP),
      cols = vars(Common_name),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      labeller = labeller(
        Common_name = celltype_labels
      ),
      strip = ggh4x::strip_themed(
        background_x = list(
          ggplot2::element_rect(
            fill = "#CFE8F3",
            color = "grey55",
            linewidth = 0.4
          ),
          ggplot2::element_rect(
            fill = "#74A9CF",
            color = "grey55",
            linewidth = 0.4
          ),
          ggplot2::element_rect(
            fill = "#045A8D",
            color = "grey55",
            linewidth = 0.4
          )
        ),
        text_x = list(
          ggplot2::element_text(
            color = "black",
            face = "bold",
            size = 13
          ),
          ggplot2::element_text(
            color = "black",
            face = "bold",
            size = 13
          ),
          ggplot2::element_text(
            color = "white",
            face = "bold",
            size = 13
          )
        ),
        background_y = list(
          ggplot2::element_blank(),
          ggplot2::element_blank(),
          ggplot2::element_blank()
        ),
        text_y = list(
          ggplot2::element_text(
            color = "grey25",
            face = "bold",
            size = 8.5,
            angle = 0
          ),
          ggplot2::element_text(
            color = "grey25",
            face = "bold",
            size = 8.5,
            angle = 0
          ),
          ggplot2::element_text(
            color = "grey25",
            face = "bold",
            size = 8.5,
            angle = 0
          )
        )
      )
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
      breaks = seq(4, 24, by = 2),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_discrete(
      labels = dataset_labels,
      drop = TRUE
    ) +
    coord_cartesian(
      xlim = c(3, 25),
      clip = "off"
    ) +
    labs(
      x = "Human-equivalent developmental age (post-conception weeks)",
      y = NULL,
      title = paste0(gene_name, " expression dynamics across dorsal cortical lineages"),
      subtitle = paste0(
        "Datasets are grouped by species. Crosses indicate groups with fewer than ",
        min_cells, " cells."
      )
    ) +
    guides(
      fill = guide_colorbar(
        order = 1,
        title.position = "left",
        title.hjust = 0.5,
        title.vjust = 0.5,
        barwidth = unit(4.5, "cm"),
        barheight = unit(0.35, "cm"),
        ticks = TRUE
      ),
      size = guide_legend(
        order = 2,
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
        margin = ggplot2::margin(t = 4, b = 8)
      ),
      
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 13
      ),
      strip.placement = "outside",
      
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
        margin = ggplot2::margin(t = 8)
      ),
      axis.line = element_line(
        linewidth = 0.35,
        color = "grey20"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        color = "grey20"
      ),
      panel.spacing.x = unit(1.1, "lines"),
      panel.spacing.y = unit(1.4, "lines"),
      
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.6, "cm"),
      
      plot.margin = ggplot2::margin(8, 10, 8, 8)
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

# -----------------------------
# 13. Show example plots
# -----------------------------

plots_dorsal[["CNR1"]]
plots_dorsal[["CNR2"]]
plots_dorsal[["NAPEPLD"]]

# -----------------------------
# 14. Save all plots
# -----------------------------

dir.create("gene_dorsal_plots", showWarnings = FALSE)

walk(
  genes_to_plot,
  function(gene_name) {
    
    ggsave(
      filename = file.path('Downloads',
        "gene_dorsal_plots",
        paste0(gene_name, "_dorsal_lineage_dotplot_species_grouped.pdf")
      ),
      plot = plots_dorsal[[gene_name]],
      width = 18,
      height = 4,
      units = "in",
      device = cairo_pdf
    )
    
    ggsave(
      filename = file.path('Downloads',
        "gene_dorsal_plots",
        paste0(gene_name, "_dorsal_lineage_dotplot_species_grouped.png")
      ),
      plot = plots_dorsal[[gene_name]],
      width = 18,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

