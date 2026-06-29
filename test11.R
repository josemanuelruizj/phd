# -----------------------------
# CNR1 trajectory plot
# x = age
# y = % cells expressing CNR1
# color = dataset
# size = mean CNR1 expression
# facet = cell type
# custom facet title colors
# legends side by side
# -----------------------------

library(dplyr)
library(ggplot2)
library(ggh4x)

dataset_order_clean <- c(
  "La Manno et al. (2021)",
  "Di Bella et al. (2021)",
  "Micali et al. (2023)",
  "Eze et al. (2021)",
  "Braun et al. (2023)",
  "Trevino et al. (2021)"
)

dataset_colors <- c(
  "La Manno et al. (2021)" = "#1B7837",
  "Di Bella et al. (2021)" = "#A6DBA0",
  "Micali et al. (2023)"   = "#E66101",
  "Eze et al. (2021)"      = "#2166AC",
  "Braun et al. (2023)"    = "#67A9CF",
  "Trevino et al. (2021)"  = "#D1E5F0"
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
    Common_name_clean = factor(
      Common_name_clean,
      levels = c(
        "Dorsal NSC (SOX2+)",
        "Excit IPC",
        "Glut NEU",
        "CR"
      )
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
    )
  )

p_cnr1_trajectory_celltype <- ggplot(
  plot_df_cnr1,
  aes(
    x = human_age,
    y = percent_expressing,
    color = DATASET_clean,
    group = DATASET_clean
  )
) +
  geom_vline(
    data = important_ages,
    aes(xintercept = human_age),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey75",
    alpha = 0.85
  ) +
  geom_line(
    linewidth = 0.6,
    alpha = 0.65
  ) +
  geom_point(
    aes(size = avg_expression_norm),
    alpha = 0.9
  ) +
  ggh4x::facet_wrap2(
    ~ Common_name_clean,
    ncol = 3,
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
      )
    )
  ) +
  scale_color_manual(
    values = dataset_colors,
    drop = FALSE,
    name = "Dataset"
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
    title = "CNR1 expression frequency across dorsal cortical lineages",
    subtitle = paste0(
      "Each facet represents a dorsal cortical lineage. ",
      "Color indicates dataset; point size indicates mean CNR1 expression scaled per dataset. ",
      "Groups with fewer than ", min_cells, " cells are excluded."
    )
  ) +
  guides(
    color = guide_legend(
      order = 1,
      title.position = "left",
      direction = "horizontal",
      nrow = 2,
      override.aes = list(
        linewidth = 1,
        size = 3,
        alpha = 1
      )
    ),
    size = guide_legend(
      order = 2,
      title.position = "left",
      direction = "horizontal",
      nrow = 1,
      override.aes = list(
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
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8.5
    ),
    axis.text.y = element_text(
      size = 8.5
    ),
    axis.title = element_text(
      size = 10
    ),
    axis.line = element_line(
      linewidth = 0.35,
      color = "grey20"
    ),
    axis.ticks = element_line(
      linewidth = 0.3,
      color = "grey20"
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5),
    legend.spacing.x = unit(0.6, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.65, "cm"),
    panel.spacing = unit(1.1, "lines"),
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

p_cnr1_trajectory_celltype