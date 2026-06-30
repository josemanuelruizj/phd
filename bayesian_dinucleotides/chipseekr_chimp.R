BiocManager::install("TxDb.Ptroglodytes.UCSC.panTro6.refGene")
 BiocManager::install("org.Pt.eg.db", force = TRUE)
 
# Load necessary libraries
library(ChIPseeker)
library(TxDb.Ptroglodytes.UCSC.panTro6.refGene)  # Chimpanzee annotation
library(org.Pt.eg.db)  # Chimpanzee gene database
library(dplyr)
library(ggplot2)

# ============================================
# CREATE TxDb OBJECT FOR panTro6


txdb <- TxDb.Ptroglodytes.UCSC.panTro6.refGene

# ============================================
# SET WORKING DIRECTORY
# ============================================
setwd("/users/genomics/lia/significant_regions_log2_allspecies/chimpanzee")

# ============================================
# CREATE OUTPUT DIRECTORY
# ============================================
output_dir <- "functional_analysis_panTro6"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Directorio creado:", output_dir, "\n")
}

# List all .bed files
bed_files <- list.files(pattern = "\\.bed$", full.names = TRUE)

# View found files
cat("Archivos BED encontrados:\n")
print(basename(bed_files))
cat("Total de archivos:", length(bed_files), "\n\n")

# Create empty dataframe to combine all results
all_annotations_combined <- data.frame()

# Process each .bed file
for (bed_file in bed_files) {
  
  # Get full file name (without .bed)
  full_name <- gsub("\\.bed$", "", basename(bed_file))
  
  # Extract CAN-CAN pattern from file name
  can_pattern <- strsplit(full_name, "_")[[1]][1]
  
  # Also extract additional information from name if needed
  name_parts <- strsplit(full_name, "_")[[1]]
  region_type <- ifelse(length(name_parts) >= 2, name_parts[2], NA)
  region_number <- ifelse(length(name_parts) >= 3, name_parts[3], NA)
  
  cat("========================================\n")
  cat("Procesando:", full_name, "\n")
  cat("Patrón CAN-CAN:", can_pattern, "\n")
  cat("Tipo de región:", region_type, "\n")
  cat("Número:", region_number, "\n")
  
  # Read BED file
  peaks <- readPeakFile(bed_file)
  cat("Número de regiones:", length(peaks), "\n")
  
  # Annotate peaks using panTro6
  peak_anno <- annotatePeak(peaks, 
                            TxDb = txdb,
                            tssRegion = c(-3000, 3000),
                            annoDb = "org.Pt.eg.db",  # Changed to chimpanzee
                            verbose = FALSE)
  
  # Convert to dataframe
  anno_df <- as.data.frame(peak_anno)
  
  # Add columns with file information
  anno_df$CAN_pattern <- can_pattern
  anno_df$file_name <- full_name
  anno_df$region_type <- region_type
  anno_df$region_number <- region_number
  
  # Combine with main dataframe
  all_annotations_combined <- bind_rows(all_annotations_combined, anno_df)
  
  # Show summary
  cat("Resumen de anotación:\n")
  print(table(anno_df$annotation))
  cat("\n")
  
  # Save individual table in output directory
  output_file <- file.path(output_dir, paste0(full_name, "_annotation_panTro6.txt"))
  write.table(anno_df, 
              file = output_file,
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  cat("Tabla individual guardada en:", output_file, "\n\n")
}

# ============================================
# SAVE FINAL COMBINED TABLE
# ============================================

cat("========================================\n")
cat("Guardando tabla combinada final...\n")

# Save complete table in output directory
output_combined <- file.path(output_dir, "all_CAN_patterns_combined_annotation_panTro6.txt")
write.table(all_annotations_combined, 
            output_combined,
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("Tabla combinada guardada en:", output_combined, "\n")
cat("Dimensiones de la tabla:", nrow(all_annotations_combined), "filas x", 
    ncol(all_annotations_combined), "columnas\n")

# ============================================
# SUMMARY BY CAN-CAN PATTERN
# ============================================

cat("\n========================================\n")
cat("Generando resumen por patrón CAN-CAN...\n")

# Create summary by pattern
summary_by_pattern <- all_annotations_combined %>%
  group_by(CAN_pattern, annotation) %>%
  summarise(
    count = n(),
    .groups = 'drop'
  ) %>%
  group_by(CAN_pattern) %>%
  mutate(percentage = count/sum(count)*100)

output_summary_pattern <- file.path(output_dir, "summary_by_CAN_pattern_panTro6.txt")
write.table(summary_by_pattern, 
            output_summary_pattern,
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Show summary
cat("\nResumen de distribución genómica por patrón CAN-CAN:\n")
print(summary_by_pattern)

# ============================================
# SUMMARY BY FILE
# ============================================

cat("\n========================================\n")
cat("Generando resumen por archivo...\n")

summary_by_file <- all_annotations_combined %>%
  group_by(file_name, CAN_pattern, annotation) %>%
  summarise(
    count = n(),
    .groups = 'drop'
  ) %>%
  group_by(file_name) %>%
  mutate(percentage = count/sum(count)*100)

output_summary_file <- file.path(output_dir, "summary_by_file_panTro6.txt")
write.table(summary_by_file, 
            output_summary_file,
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# ============================================
# GENE LIST BY CAN-CAN PATTERN
# ============================================

cat("\n========================================\n")
cat("Generando lista de genes por patrón CAN-CAN...\n")

# Extract unique genes by pattern
genes_by_pattern <- all_annotations_combined %>%
  filter(!is.na(SYMBOL)) %>%
  select(CAN_pattern, SYMBOL, GENENAME, annotation, distanceToTSS) %>%
  distinct()

output_genes <- file.path(output_dir, "genes_by_CAN_pattern_panTro6.txt")
write.table(genes_by_pattern, 
            output_genes,
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Show number of unique genes by pattern
genes_summary <- genes_by_pattern %>%
  group_by(CAN_pattern) %>%
  summarise(
    unique_genes = n_distinct(SYMBOL),
    .groups = 'drop'
  )

cat("\nNúmero de genes únicos por patrón CAN-CAN:\n")
print(genes_summary)

# ============================================
# COMPARATIVE VISUALIZATIONS
# ============================================

cat("\n========================================\n")
cat("Generando visualizaciones comparativas...\n")

# Genomic distribution by CAN-CAN pattern
p1 <- ggplot(summary_by_pattern, aes(x = CAN_pattern, y = percentage, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Distribución genómica por patrón CAN-CAN (panTro6)",
       x = "Patrón CAN-CAN",
       y = "Porcentaje") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

output_p1 <- file.path(output_dir, "genomic_distribution_by_CAN_pattern_panTro6.pdf")
ggsave(output_p1, p1, width = 12, height = 8)

# Number of regions by pattern
regions_by_pattern <- all_annotations_combined %>%
  group_by(CAN_pattern, file_name) %>%
  summarise(
    total_regions = n(),
    .groups = 'drop'
  )

p2 <- ggplot(regions_by_pattern, aes(x = CAN_pattern, y = total_regions, fill = file_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Número de regiones por patrón CAN-CAN (panTro6)",
       x = "Patrón CAN-CAN",
       y = "Número de regiones") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

output_p2 <- file.path(output_dir, "regions_by_CAN_pattern_panTro6.pdf")
ggsave(output_p2, p2, width = 12, height = 8)

# ============================================
# COMPARATIVE DISTANCE TO TSS PLOTS
# ============================================

cat("\n========================================\n")
cat("Generando visualizaciones adicionales...\n")

# Distance to TSS distribution by pattern
p3 <- ggplot(all_annotations_combined, aes(x = CAN_pattern, y = distanceToTSS, fill = CAN_pattern)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distancia al TSS por patrón CAN-CAN (panTro6)",
       x = "Patrón CAN-CAN",
       y = "Distancia al TSS (bp)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10()

output_p3 <- file.path(output_dir, "distance_to_TSS_by_pattern_panTro6.pdf")
ggsave(output_p3, p3, width = 12, height = 8)

# ============================================
# FINAL SUMMARY
# ============================================

cat("\n========================================\n")
cat("ANÁLISIS COMPLETADO - GENOMA PAN TRO 6\n")
cat("========================================\n")
cat("Todos los archivos se guardaron en el directorio:", output_dir, "\n")
cat("Ruta completa:", file.path(getwd(), output_dir), "\n\n")
cat("Archivos generados:\n")
cat("  - all_CAN_patterns_combined_annotation_panTro6.txt: TABLA PRINCIPAL\n")
cat("  - summary_by_CAN_pattern_panTro6.txt: Resumen por patrón\n")
cat("  - summary_by_file_panTro6.txt: Resumen por archivo\n")
cat("  - genes_by_CAN_pattern_panTro6.txt: Lista de genes\n")
cat("  - genomic_distribution_by_CAN_pattern_panTro6.pdf: Gráfico distribución\n")
cat("  - regions_by_CAN_pattern_panTro6.pdf: Gráfico de regiones\n")
cat("  - distance_to_TSS_by_pattern_panTro6.pdf: Distancia al TSS\n")
cat("  - *_annotation_panTro6.txt: Tablas individuales\n")
cat("========================================\n")

# Save workspace
output_rdata <- file.path(output_dir, "chipseeker_CAN_analysis_panTro6.RData")
save.image(output_rdata)
cat("\nAmbiente guardado en:", output_rdata, "\n")

# Show first rows of combined table
cat("\nPrimeras 10 filas de la tabla combinada:\n")
print(head(all_annotations_combined, 10))

# Show complete location
cat("\n========================================\n")
cat("UBICACIÓN COMPLETA DE LOS RESULTADOS:\n")
cat(file.path(getwd(), output_dir), "\n")
cat("========================================\n")
