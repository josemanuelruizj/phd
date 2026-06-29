library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(scales)

# -----------------------------
# 1. Read age file
# -----------------------------

ages <- read_excel("Downloads/converted_ages_long.xlsx")

# -----------------------------
# 2. Dataset order and labels
# -----------------------------

dataset_order <- c(
  "MANNO",
  "DIBELLA",
  "MICALI",
  "EZE",
  "BRAUN",
  "TREVINO"
)

dataset_labels <- c(
  "MANNO"   = "La Manno et al. (2021)",
  "DIBELLA" = "Di Bella et al. (2021)",
  "MICALI"  = "Micali et al. (2023)",
  "EZE"     = "Eze et al. (2021)",
  "BRAUN"   = "Braun et al. (2023)",
  "TREVINO" = "Trevino et al. (2021)"
)

species_labels <- c(
  "MANNO"   = "Mouse",
  "DIBELLA" = "Mouse",
  "MICALI"  = "Macaque",
  "EZE"     = "Human",
  "BRAUN"   = "Human",
  "TREVINO" = "Human"
)

species_colors <- c(
  "Mouse"   = "#4DAF4A",
  "Macaque" = "#E66101",
  "Human"   = "#2166AC"
)

# -----------------------------
# 3. Prepare data
# -----------------------------

age_plot_data <- ages %>%
  mutate(
    DATASET = toupper(Dataset),
    real_age_label = as.character(`Real age`),
    converted_age_days = `Converted (heterochrony or days)`,
    converted_age_pcw = converted_age_days / 7,
    
    DATASET = factor(DATASET, levels = rev(dataset_order)),
    
    DATASET_label = recode(
      as.character(DATASET),
      !!!dataset_labels
    ),
    
    SPECIES_GROUP = recode(
      as.character(DATASET),
      !!!species_labels
    ),
    
    SPECIES_GROUP = factor(
      SPECIES_GROUP,
      levels = c("Mouse", "Macaque", "Human")
    )
  ) %>%
  filter(
    !is.na(DATASET),
    !is.na(converted_age_pcw),
    !is.na(real_age_label),
    !is.na(SPECIES_GROUP)
  ) %>%
  group_by(DATASET) %>%
  arrange(converted_age_pcw, .by_group = TRUE) %>%
  mutate(
    # Alternate labels above and below the point
    label_offset = ifelse(row_number() %% 2 == 0, -0.28, 0.28),
    label_y = as.numeric(DATASET) + label_offset
  ) %>%
  ungroup()

# -----------------------------
# 4. Vertical reference lines every 2 PCW
# -----------------------------

important_ages <- tibble::tibble(
  human_age_pcw = seq(4, 24, by = 2),
  label = paste0(seq(4, 24, by = 2), " PCW")
)

# -----------------------------
# 5. Plot
# -----------------------------

p_age_timeline <- ggplot(
  age_plot_data,
  aes(
    x = converted_age_pcw,
    y = DATASET,
    color = SPECIES_GROUP
  )
) +
  geom_vline(
    data = important_ages,
    aes(xintercept = human_age_pcw),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey75",
    alpha = 0.85
  ) +
  geom_hline(
    aes(yintercept = as.numeric(DATASET)),
    color = "grey90",
    linewidth = 0.35
  ) +
  geom_point(
    size = 2.8,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = age_plot_data,
    aes(
      x = converted_age_pcw,
      y = DATASET,
      label = real_age_label
    ),
    inherit.aes = FALSE,
    size = 3,
    color = "grey15",
    
    # Alternate labels above and below the real point
    nudge_y = age_plot_data$label_offset,
    
    # Mostly vertical adjustment
    direction = "y",
    force = 1.2,
    force_pull = 2.5,
    
    box.padding = 0.15,
    point.padding = 0.12,
    min.segment.length = 0,
    
    segment.color = "grey55",
    segment.linewidth = 0.25,
    segment.alpha = 0.8,
    
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_color_manual(
    values = species_colors,
    name = "Species group"
  ) +
  scale_y_discrete(
    labels = dataset_labels
  ) +
  scale_x_continuous(
    breaks = seq(4, 24, by = 2),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  coord_cartesian(
    xlim = c(3, 25),
    clip = "off"
  ) +
  labs(
    x = "Human-equivalent developmental age (post-conception weeks)",
    y = NULL,
    title = "Dataset sampling ages projected onto human-equivalent developmental time",
    subtitle = "Point labels indicate the original age annotation reported for each dataset."
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
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5),
    legend.box = "horizontal",
    plot.margin = ggplot2::margin(8, 18, 8, 8)
  )

p_age_timeline

# -----------------------------
# 6. Save plot
# -----------------------------

ggsave(
  filename = "dataset_real_vs_human_equivalent_age_timeline.pdf",
  plot = p_age_timeline,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "dataset_real_vs_human_equivalent_age_timeline.png",
  plot = p_age_timeline,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)
