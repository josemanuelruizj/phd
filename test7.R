library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

# -----------------------------
# 1. Read files
# -----------------------------

expr <- read_excel("Downloads/summary_table_average_expr_cells.xlsx")
ages <- read_excel("Downloads/converted_ages_long.xlsx")

# -----------------------------
# 2. Keep only CNR1-relevant columns
# -----------------------------

expr_cnr1 <- expr %>%
  select(
    DATASET,
    Common_name,
    AGE_combined,
    n_cells,
    CNR1_avg_expression,
    CNR1_percent_expressing
  )

# -----------------------------
# 3. Clean / harmonise age labels
# -----------------------------

expr_clean <- expr_cnr1 %>%
  mutate(
    DATASET = toupper(DATASET),
    
    age_join = AGE_combined %>%
      as.character() %>%
      str_remove("^.*_") %>%
      str_replace_all("\\.", ",")
  )

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
# 4. Join with converted ages
# -----------------------------

plot_data <- expr_clean %>%
  left_join(ages_clean, by = c("DATASET", "age_join")) %>%
  filter(!is.na(human_age)) %>%
  group_by(DATASET) %>%
  mutate(
    CNR1_avg_expression_norm =
      CNR1_avg_expression / max(CNR1_avg_expression, na.rm = TRUE)
  ) %>%
  ungroup()

# -----------------------------
# 5. Keep only dorsal lineage
# -----------------------------

dorsal_lineage_order <- c(
  "Dorsal_NSC(SOX2+)",
  "Excit_IPC",
  "Glut_NEU",
  "CR"
)

plot_data <- plot_data %>%
  mutate(
    Common_name = as.character(Common_name)
  ) %>%
  filter(
    Common_name %in% dorsal_lineage_order
  )

# -----------------------------
# 6. Order datasets and dorsal cell types
# -----------------------------

dataset_order <- c(
  "MANNO",
  "DIBELLA",
  "MICALI",
  "EZE",
  "BRAUN",
  "TREVINO"
)

plot_data <- plot_data %>%
  mutate(
    DATASET = factor(DATASET, levels = rev(dataset_order)),
    Common_name = factor(Common_name, levels = dorsal_lineage_order)
  ) %>%
  filter(
    !is.na(DATASET),
    !is.na(Common_name)
  )

# -----------------------------
# 7. Plot CNR1 across dorsal lineage
# -----------------------------

p_cnr1_dorsal <- ggplot(
  plot_data,
  aes(
    x = human_age,
    y = DATASET,
    size = CNR1_percent_expressing,
    fill = CNR1_avg_expression_norm
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  facet_wrap(~ Common_name, ncol = 4) +
  scale_fill_gradient(
    low = "white",
    high = "blue",
    limits = c(0, 1),
    name = "CNR1 avg expression\n0-max per dataset"
  ) +
  scale_size_continuous(
    range = c(1, 8),
    name = "% CNR1+ cells"
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  labs(
    x = "Human equivalent age",
    y = "Dataset",
    title = "CNR1 expression across dorsal differentiation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_cnr1_dorsal





# -----------------------------
# 4. Join with converted ages
# -----------------------------

min_cells <- 25

plot_data <- expr_clean %>%
  left_join(ages_clean, by = c("DATASET", "age_join")) %>%
  filter(!is.na(human_age)) %>%
  mutate(
    Common_name = as.character(Common_name),
    enough_cells = n_cells >= min_cells
  )

# -----------------------------
# 5. Keep only dorsal lineage
# -----------------------------

dorsal_lineage_order <- c(
  "Dorsal_NSC(SOX2+)",
  "Excit_IPC",
  "Glut_NEU"
)

plot_data <- plot_data %>%
  filter(Common_name %in% dorsal_lineage_order)

# -----------------------------
# 6. Normalize expression only using groups with enough cells
# -----------------------------

plot_data <- plot_data %>%
  group_by(DATASET) %>%
  mutate(
    max_expr_valid = max(CNR1_avg_expression[enough_cells], na.rm = TRUE),
    CNR1_avg_expression_norm = ifelse(
      enough_cells,
      CNR1_avg_expression / max_expr_valid,
      NA_real_
    )
  ) %>%
  ungroup()

# -----------------------------
# 7. Order datasets and dorsal cell types
# -----------------------------

dataset_order <- c(
  "MANNO",
  "DIBELLA",
  "MICALI",
  "EZE",
  "BRAUN",
  "TREVINO"
)

plot_data <- plot_data %>%
  mutate(
    DATASET = factor(DATASET, levels = rev(dataset_order)),
    Common_name = factor(Common_name, levels = dorsal_lineage_order)
  ) %>%
  filter(
    !is.na(DATASET),
    !is.na(Common_name)
  )

# -----------------------------
# 8. Split included and excluded points
# -----------------------------

plot_data_included <- plot_data %>%
  filter(enough_cells)

plot_data_excluded <- plot_data %>%
  filter(!enough_cells)

# -----------------------------
# 9. Plot CNR1 across dorsal lineage
# -----------------------------

p_cnr1_dorsal <- ggplot() +
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
  scale_fill_gradient(
    low = "white",
    high = "blue",
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
    subtitle = paste0("X = group exists but excluded because n cells < ", min_cells)
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_cnr1_dorsal
