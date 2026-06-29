# -----------------------------
# Expression trajectory plot
# x = age
# y = expression level
# color = cell type
# size = % expressing
# -----------------------------

plot_expression_trajectory <- function(dataset_name = NULL) {
  
  plot_df <- plot_data_all %>%
    filter(gene %in% genes_to_plot) %>%
    filter(enough_cells) %>%
    mutate(
      Common_name_clean = recode(
        as.character(Common_name),
        "Dorsal_NSC(SOX2+)" = "Dorsal NSC (SOX2+)",
        "Excit_IPC" = "Excit IPC",
        "Glut_NEU" = "Glut NEU",
        "CR" = "CR"
      )
    )
  
  if (!is.null(dataset_name)) {
    plot_df <- plot_df %>%
      filter(DATASET == dataset_name)
  }
  
  ggplot(
    plot_df,
    aes(
      x = human_age,
      y = avg_expression_norm,
      color = Common_name_clean,
      group = Common_name_clean
    )
  ) +
    geom_vline(
      data = important_ages,
      aes(xintercept = human_age),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.25,
      color = "grey80",
      alpha = 0.8
    ) +
    geom_line(
      linewidth = 0.5,
      alpha = 0.65
    ) +
    geom_point(
      aes(size = percent_expressing),
      alpha = 0.9
    ) +
    facet_wrap(
      ~ gene,
      ncol = 4,
      scales = "free_y"
    ) +
    scale_size_continuous(
      range = c(1.2, 6),
      breaks = c(25, 50, 75, 100),
      limits = c(0, 100),
      name = "% expressing"
    ) +
    scale_x_continuous(
      breaks = seq(20, 170, by = 20),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    coord_cartesian(
      xlim = c(20, 170),
      clip = "off"
    ) +
    labs(
      x = "Human-equivalent developmental age (days)",
      y = "Mean expression scaled 0–1",
      color = "Cell type",
      title = ifelse(
        is.null(dataset_name),
        "Expression trajectories across dorsal cortical lineages",
        paste0("Expression trajectories in ", dataset_labels[[dataset_name]])
      ),
      subtitle = "Point size indicates the percentage of cells expressing each gene"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, color = "grey25"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
      axis.text.y = element_text(size = 8.5),
      axis.title = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.spacing = unit(1.1, "lines")
    )
}

plot_expression_trajectory()


# -----------------------------
# CNR1 trajectory plot
# x = age
# y = expression level
# color = cell type
# size = % expressing
# -----------------------------

plot_df_cnr1 <- plot_data_all %>%
  filter(
    gene == "CNR1",
    enough_cells
  ) %>%
  mutate(
    Common_name_clean = recode(
      as.character(Common_name),
      "Dorsal_NSC(SOX2+)" = "Dorsal NSC (SOX2+)",
      "Excit_IPC" = "Excit IPC",
      "Glut_NEU" = "Glut NEU",
      "CR" = "CR"
    ),
    DATASET_clean = recode(
      as.character(DATASET),
      "MANNO"   = "La Manno et al. (2021)",
      "DIBELLA" = "Di Bella et al. (2021)",
      "MICALI"  = "Micali et al. (2023)",
      "EZE"     = "Eze et al. (2021)",
      "BRAUN"   = "Braun et al. (2023)",
      "TREVINO" = "Trevino et al. (2021)"
    )
  )

p_cnr1_trajectory <- ggplot(
  plot_df_cnr1,
  aes(
    x = human_age,
    y = avg_expression_norm,
    color = Common_name_clean,
    group = Common_name_clean
  )
) +
  geom_vline(
    data = important_ages,
    aes(xintercept = human_age),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.25,
    color = "grey80",
    alpha = 0.8
  ) +
  geom_line(
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_point(
    aes(size = percent_expressing),
    alpha = 0.9
  ) +
  facet_wrap(
    ~ DATASET_clean,
    ncol = 3
  ) +
  scale_size_continuous(
    range = c(1.2, 6),
    breaks = c(25, 50, 75, 100),
    limits = c(0, 100),
    name = "CNR1+ cells (%)"
  ) +
  scale_x_continuous(
    breaks = seq(20, 170, by = 20),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  coord_cartesian(
    xlim = c(20, 170),
    clip = "off"
  ) +
  labs(
    x = "Human-equivalent developmental age (days)",
    y = "Mean CNR1 expression scaled 0–1",
    color = "Cell type",
    title = "CNR1 expression trajectories across dorsal cortical lineages",
    subtitle = paste0(
      "Color indicates cell type; point size indicates percentage of CNR1-expressing cells. ",
      "Groups with fewer than ", min_cells, " cells are excluded."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "grey25"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 8.5),
    axis.title = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.spacing = unit(1.1, "lines")
  )

p_cnr1_trajectory







# -----------------------------
# CNR1 trajectory plot
# x = age
# y = % cells expressing CNR1
# color = cell type
# size = mean CNR1 expression
# datasets ordered as requested
# -----------------------------

dataset_order_clean <- c(
  "La Manno et al. (2021)",
  "Di Bella et al. (2021)",
  "Micali et al. (2023)",
  "Eze et al. (2021)",
  "Braun et al. (2023)",
  "Trevino et al. (2021)"
)

important_ages <- tibble::tibble(
  human_age = c(28, 56, 84, 112, 140, 168),
  label = c("4 PCW", "8 PCW", "12 PCW", "16 PCW", "20 PCW", "24 PCW")
)

plot_df_cnr1 <- plot_data_all %>%
  filter(
    gene == "CNR1",
    enough_cells
  ) %>%
  mutate(
    Common_name_clean = recode(
      as.character(Common_name),
      "Dorsal_NSC(SOX2+)" = "Dorsal NSC (SOX2+)",
      "Excit_IPC" = "Excit IPC",
      "Glut_NEU" = "Glut NEU",
      "CR" = "CR"
    ),
    DATASET_clean = recode(
      as.character(DATASET),
      "MANNO"   = "La Manno et al. (2021)",
      "DIBELLA" = "Di Bella et al. (2021)",
      "MICALI"  = "Micali et al. (2023)",
      "EZE"     = "Eze et al. (2021)",
      "BRAUN"   = "Braun et al. (2023)",
      "TREVINO" = "Trevino et al. (2021)"
    ),
    DATASET_clean = factor(
      DATASET_clean,
      levels = dataset_order_clean
    ),
    Common_name_clean = factor(
      Common_name_clean,
      levels = c(
        "Dorsal NSC (SOX2+)",
        "Excit IPC",
        "Glut NEU",
        "CR"
      )
    )
  )

p_cnr1_trajectory <- ggplot(
  plot_df_cnr1,
  aes(
    x = human_age,
    y = percent_expressing,
    color = Common_name_clean,
    group = Common_name_clean
  )
) +
  geom_vline(
    data = important_ages,
    aes(xintercept = human_age),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.25,
    color = "grey80",
    alpha = 0.8
  ) +
  geom_line(
    linewidth = 0.6,
    alpha = 0.65
  ) +
  geom_point(
    aes(size = avg_expression_norm),
    alpha = 0.9
  ) +
  facet_wrap(
    ~ DATASET_clean,
    ncol = 3
  ) +
  scale_size_continuous(
    range = c(1.2, 6),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    name = "Mean CNR1 expression\nscaled 0–1"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_continuous(
    breaks = seq(20, 170, by = 20),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  coord_cartesian(
    xlim = c(20, 170),
    ylim = c(0, 100),
    clip = "off"
  ) +
  labs(
    x = "Human-equivalent developmental age (days)",
    y = "CNR1-expressing cells (%)",
    color = "Cell type",
    title = "CNR1 expression frequency across dorsal cortical lineages",
    subtitle = paste0(
      "Y-axis indicates the percentage of CNR1-expressing cells; ",
      "point size indicates mean CNR1 expression scaled per dataset. ",
      "Groups with fewer than ", min_cells, " cells are excluded."
    )
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      direction = "horizontal"
    ),
    size = guide_legend(
      title.position = "left",
      direction = "horizontal"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "grey25"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 8.5),
    axis.title = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.spacing = unit(1.1, "lines")
  )

p_cnr1_trajectory