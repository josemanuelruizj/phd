library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

# Descarga manual el Excel desde Drive y pon la ruta aquí
library(googlesheets4)

sheet_url <- "https://docs.google.com/spreadsheets/d/1t6glRqrmuUh0w6Cu0xH6yF7Q4c8AedAO/edit?gid=2073198212#gid=2073198212"

raw_data <- read_sheet(
  sheet_url,
  sheet = "Copia de Original Josema",
  col_names = FALSE
)
colnames(raw_data) <- paste0("V", seq_len(ncol(raw_data)))

converted_ages <- raw_data %>%
  mutate(
    Dataset = ifelse(!V1 %in% c("heterochrony", "days"), V1, NA_character_)
  ) %>%
  fill(Dataset, .direction = "down") %>%
  mutate(row_type = case_when(
    V1 %in% c("heterochrony", "days") ~ "Converted",
    TRUE ~ "Real_age"
  )) %>%
  select(Dataset, row_type, V2:last_col()) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "col_id",
    values_to = "value"
  ) %>%
  filter(!is.na(value), value != "") %>%
  pivot_wider(
    names_from = row_type,
    values_from = value
  ) %>%
  filter(!is.na(Real_age), !is.na(Converted)) %>%
  rename(
    `Real age` = Real_age,
    `Converted (heterochrony or days)` = Converted
  ) %>%
  mutate(
    `Real age` = as.character(`Real age`),
    `Converted (heterochrony or days)` =
      as.numeric(str_replace_all(as.character(`Converted (heterochrony or days)`), ",", "."))
  )

write_xlsx(converted_ages, "converted_ages_long.xlsx")
write_tsv(converted_ages, "converted_ages_long.tsv")