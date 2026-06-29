# -----------------------------
# 9. Plot CNR1 across dorsal lineage
# -----------------------------

important_ages <- tibble::tibble(
  human_age = c(28, 56, 84, 112, 140, 168),
  label = c(
    "4 PCW",
    "8 PCW",
    "12 PCW",
    "16 PCW",
    "20 PCW",
    "24 PCW"
  )
)

p_cnr1_dorsal <- ggplot() +
  geom_vline(
    data = important_ages,
    aes(xintercept = human_age),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey45",
    alpha = 0.7
  ) +
  geom_point(
    data = plot_data_included,
    aes(
      x = human_age,
      y = DATASET,
      size = CNR1_percent_expressing,
      fill = CNR1_avg_expression_norm
    ),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  geom_point(
    data = plot_data_excluded,
    aes(
      x = human_age,
      y = DATASET
    ),
    shape = 4,
    size = 3,
    stroke = 0.8,
    color = "grey30"
  ) +
  facet_wrap(~ Common_name, ncol = 4) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = c(0, 1),
    name = "CNR1 avg expression\n0-max"
  ) +
  scale_size_continuous(
    range = c(1, 8),
    name = "% CNR1+"
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  labs(
    x = "Human equivalent age",
    y = "Dataset",
    title = "CNR1 expression across dorsal differentiation",
    subtitle = paste0(
      "X = group exists but excluded because n cells < ", min_cells,
      ". Dashed lines mark 4-week human-equivalent developmental intervals."
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5,
      title.vjust = 0.5,
      barwidth = unit(4, "cm"),
      barheight = unit(0.35, "cm")
    ),
    size = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      title.vjust = 0.5,
      direction = "horizontal"
    )
  )+
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p_cnr1_dorsal






################



# -----------------------------
# 9. Plot CNR1 across dorsal lineage
# -----------------------------

important_ages <- tibble::tibble(
  human_age = c(28, 56, 84, 112, 140, 168),
  label = c(
    "4 PCW",
    "8 PCW",
    "12 PCW",
    "16 PCW",
    "20 PCW",
    "24 PCW"
  )
)

p_cnr1_dorsal <- ggplot() +
  geom_vline(
    data = important_ages,
    aes(xintercept = human_age),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey45",
    alpha = 0.7
  ) +
  geom_point(
    data = plot_data_included,
    aes(
      x = human_age,
      y = DATASET,
      size = CNR1_avg_expression_norm,
      fill = CNR1_percent_expressing
    ),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  geom_point(
    data = plot_data_excluded,
    aes(
      x = human_age,
      y = DATASET
    ),
    shape = 4,
    size = 3,
    stroke = 0.8,
    color = "grey30"
  ) +
  facet_wrap(~ Common_name, ncol = 4) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = c(0, 100),
    name = "% CNR1+"
  ) +
  scale_size_continuous(
    range = c(1, 8),
    limits = c(0, 1),
    name = "CNR1 avg expression\n0-max"
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  labs(
    x = "Human equivalent age",
    y = "Dataset",
    title = "CNR1 expression across dorsal differentiation",
    subtitle = paste0(
      "X = group exists but excluded because n cells < ", min_cells,
      ". Dashed lines mark 4-week human-equivalent developmental intervals."
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5,
      title.vjust = 0.5,
      barwidth = unit(4, "cm"),
      barheight = unit(0.35, "cm")
    ),
    size = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      title.vjust = 0.5,
      direction = "horizontal"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(angle = 0)
  )

p_cnr1_dorsal




###################



# -----------------------------
# 9. Plot CNR1 across dorsal lineage
# -----------------------------

celltype_labels <- c(
  "Dorsal_NSC(SOX2+)" = "Dorsal NSC (SOX2+)",
  "Excit_IPC" = "Excit IPC",
  "Glut_NEU" = "Glut NEU",
  "CR" = "CR"
)


dataset_labels <- c(
  "MANNO"   = "La Manno et al. (2021)",
  "DIBELLA" = "Di Bella et al. (2021)",
  "MICALI"  = "Micali et al. (2023)",
  "EZE"     = "Eze et al. (2021)",
  "BRAUN"   = "Braun et al. (2023)",
  "TREVINO" = "Trevino et al. (2021)"
)

p_cnr1_dorsal <- ggplot() +
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
      size = CNR1_percent_expressing,
      fill = CNR1_avg_expression_norm
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
    y = NULL,
    title = "CNR1 expression dynamics across dorsal cortical lineages",
    subtitle = paste0(
      "Crosses indicate groups with fewer than ", min_cells, " cells."
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
  )+
  scale_y_discrete(
    labels = dataset_labels
  ) 




p_cnr1_dorsal



# -----------------------------
# 13. Save all plots
# -----------------------------

dir.create("gene_dorsal_plots", showWarnings = FALSE)

walk(
  genes_to_plot,
  function(gene_name) {
    ggsave(
      filename = file.path("gene_dorsal_plots", paste0(gene_name, "_dorsal_lineage_dotplot.pdf")),
      plot = plots_dorsal[[gene_name]],
      width = 9.73,
      height = 4.06,
      units = "in",
      device = cairo_pdf
    )
    
    ggsave(
      filename = file.path("gene_dorsal_plots", paste0(gene_name, "_dorsal_lineage_dotplot.png")),
      plot = plots_dorsal[[gene_name]],
      width = 9.73,
      height = 4.06,
      units = "in",
      dpi = 600
    )
  }
)
