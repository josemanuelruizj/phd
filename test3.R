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
# 5. Use ALL datasets
# -----------------------------

dataset_order <- plot_data %>%
  distinct(DATASET) %>%
  arrange(DATASET) %>%
  pull(DATASET)

plot_data <- plot_data %>%
  mutate(
    DATASET = factor(DATASET, levels = dataset_order),
    Common_name = factor(Common_name)
  )

# -----------------------------
# 6. Plot CNR1 across all datasets
# -----------------------------

p_cnr1 <- ggplot(
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
  facet_wrap(~ Common_name) +
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
    title = "CNR1 expression across developmental time"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_cnr1