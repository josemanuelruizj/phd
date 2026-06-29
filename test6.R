# =============================================================================
# TWO-STEP UCELL ANNOTATION FOR DORSAL CORTICAL ORGANOID MULTIOME DATA
# 1) Broad major cell type classification
# 2) Fine subtype classification inside each major type
# =============================================================================

library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# =============================================================================
# 0. LOAD OBJECT AND BASIC SETUP
# =============================================================================

data <- readRDS("integratedmultiome_processed.rds")

DefaultAssay(data) <- "RNA"

# Join Seurat v5 RNA layers if present
if (length(Layers(data[["RNA"]])) > 1) {
  data[["RNA"]] <- JoinLayers(data[["RNA"]])
}

outdir <- "multiome_plts/subplot/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. FIRST UCELL: MAJOR CELL TYPE CLASSIFICATION
# =============================================================================

major_markers_raw <- list(
  
  Progenitor = c(
    "SOX2", "PAX6", "HES1", "HES5",
    "VIM", "NES", "PROM1", "HMGA2",
    "FABP7", "SLC1A3", "GLI3",
    "MKI67", "TOP2A", "ASPM", "CENPF"
  ),
  
  Intermediate_Progenitor = c(
    "EOMES", "INSM1", "PPP1R17",
    "NEUROG1", "NEUROG2", "NEUROD1",
    "HES6", "BTG2", "DLL3", "NHLH1"
  ),
  
  Neuronal = c(
    "DCX", "STMN2", "STMN1", "TUBB3",
    "NEUROD1", "NEUROD2", "NEUROD6",
    "ELAVL3", "ELAVL4", "SOX11", "SOX4",
    "GAP43", "MAP1B",
    "SLC17A6", "SLC17A7",
    "MAP2", "RBFOX3", "SNAP25", "SYT1"
  ),
  
  Astroglia_Like = c(
    "GFAP", "AQP4", "ALDH1L1",
    "SLC1A2", "SLC1A3", "GJA1",
    "CLU", "SPARCL1", "APOE",
    "SOX9", "FGFR3"
  ),
  
  Choroid_Hem_Like = c(
    "TTR", "AQP1", "CLIC6", "FOLR1",
    "KRT8", "KRT18", "KRT19",
    "TRPM3", "PCP4", "SLC4A10", "ENPP2",
    "OTX2", "RSPO2", "WNT3A", "WNT2B",
    "BMP4", "BMP7", "MSX1", "MSX2",
    "LHX5", "ZIC1", "ZIC2"
  ),
  
  Off_Target_Non_Neural = c(
    "DCN", "BGN", "COL1A1", "COL1A2",
    "COL3A1", "COL6A1", "LUM",
    "PDGFRA", "TAGLN", "ACTA2", "MMP2",
    "PECAM1", "VWF", "KDR", "FLT1",
    "CLDN5", "ESAM", "ENG", "RAMP2",
    "EMCN", "CDH5"
  )
)

# Keep only genes present in object
major_markers <- lapply(major_markers_raw, function(genes) {
  intersect(genes, rownames(data))
})

# Remove marker sets with too few detected genes
major_markers <- major_markers[lengths(major_markers) >= 2]

message("Major marker sets used:")
print(major_markers)

# Calculate UCell scores for major types
data <- AddModuleScore_UCell(
  obj = data,
  features = major_markers,
  name = '_major'
)

major_ucell_cols <- paste0(names(major_markers), "_major")
major_ucell_cols <- intersect(major_ucell_cols, colnames(data@meta.data))

# =============================================================================
# 2. ASSIGN MAJOR CELL TYPE
# =============================================================================

major_score_mat <- as.matrix(data@meta.data[, major_ucell_cols, drop = FALSE])

major_top_score <- apply(major_score_mat, 1, max, na.rm = TRUE)

major_second_score <- apply(major_score_mat, 1, function(x) {
  sx <- sort(x, decreasing = TRUE)
  if (length(sx) >= 2) sx[2] else NA_real_
})

major_delta <- major_top_score - major_second_score

major_top_label <- colnames(major_score_mat)[
  max.col(major_score_mat, ties.method = "first")
]

major_top_label <- gsub("_UCell$", "", major_top_label)
major_top_label <- gsub("_", " ", major_top_label)
major_top_label <- tools::toTitleCase(major_top_label)

data$celltype_major_raw <- major_top_label
data$celltype_major_score <- major_top_score
data$celltype_major_delta <- major_delta

# Thresholds for broad assignment
major_min_score <- 0.12
major_min_delta <- 0.015

data$celltype_major <- ifelse(
  data$celltype_major_score >= major_min_score &
    data$celltype_major_delta >= major_min_delta,
  data$celltype_major_raw,
  "Uncertain"
)

data$celltype_major <- dplyr::recode(
  data$celltype_major,
  "Intermediate Progenitor" = "Intermediate progenitor",
  "Astroglia Like" = "Astroglia-like",
  "Choroid Hem Like" = "Choroid/Hem-like",
  "Off Target Non Neural" = "Off-target/non-neural"
)

message("Major cell type counts:")
print(
  count(data@meta.data, celltype_major) %>%
    arrange(desc(n))
)

# =============================================================================
# 3. SECOND UCELL: SUBTYPE CLASSIFICATION
# =============================================================================

subtype_markers_raw <- list(
  
  # Progenitor subtypes
  Neuroepithelial_Radial_Glia = c(
    "SOX2", "PAX6", "HES1", "HES5",
    "VIM", "NES", "PROM1", "HMGA2",
    "FABP7", "SLC1A3", "GLI3"
  ),
  
  Cycling_Progenitor = c(
    "MKI67", "TOP2A", "ASPM", "CENPF",
    "MCM2", "MCM3", "MCM5", "PCNA",
    "UBE2C", "CDK1", "CCNB1", "CCNB2"
  ),
  
  Outer_Radial_Glia_Like = c(
    "HOPX", "TNC", "PTPRZ1", "MOXD1",
    "FAM107A", "LIFR", "ITGB5", "PDGFD",
    "CLU", "CRYAB"
  ),
  
  # IPC subtype
  Intermediate_Progenitor = c(
    "EOMES", "INSM1", "PPP1R17",
    "NEUROG1", "NEUROG2", "NEUROD1",
    "HES6", "BTG2", "DLL3", "NHLH1"
  ),
  
  # Neuronal subtypes
  Neuroblast_Immature_Neuron = c(
    "DCX", "STMN2", "STMN1", "TUBB3",
    "NEUROD1", "NEUROD2", "NEUROD6",
    "ELAVL3", "ELAVL4", "SOX11", "SOX4",
    "GAP43", "MAP1B"
  ),
  
  Excitatory_Neuron = c(
    "SLC17A6", "SLC17A7",
    "MAP2", "RBFOX3", "SNAP25", "SYT1",
    "CAMK2A", "GRIN1", "GRIN2B"
  ),
  
  Deep_Layer_Excitatory_Neuron = c(
    "TBR1", "BCL11B", "FEZF2",
    "SOX5", "FOXP2", "FOXP1",
    "ZFPM2", "CRYM", "LDB2"
  ),
  
  Upper_Layer_Excitatory_Neuron = c(
    "SATB2", "CUX1", "CUX2",
    "RORB", "POU3F2", "POU3F3",
    "LHX2", "MEF2C"
  ),
  
  # Astroglia-like
  Astroglia_Like = c(
    "GFAP", "AQP4", "ALDH1L1",
    "SLC1A2", "SLC1A3", "GJA1",
    "CLU", "SPARCL1", "APOE",
    "SOX9", "FGFR3"
  ),
  
  # Choroid / hem subtypes
  Choroid_Plexus_Like = c(
    "TTR", "AQP1", "CLIC6", "FOLR1",
    "KRT8", "KRT18", "KRT19",
    "TRPM3", "PCP4", "SLC4A10",
    "ENPP2"
  ),
  
  Hem_Choroid_Patterning_Like = c(
    "OTX2", "RSPO2", "WNT3A", "WNT2B",
    "BMP4", "BMP7", "MSX1", "MSX2",
    "LHX5", "ZIC1", "ZIC2"
  ),
  
  # Off-target subtypes
  Mesenchyme_Like = c(
    "DCN", "BGN", "COL1A1", "COL1A2",
    "COL3A1", "COL6A1", "LUM",
    "PDGFRA", "TAGLN", "ACTA2", "MMP2"
  ),
  
  Endothelial_Like = c(
    "PECAM1", "VWF", "KDR", "FLT1",
    "CLDN5", "ESAM", "ENG", "RAMP2",
    "EMCN", "CDH5"
  )
)

subtype_markers <- lapply(subtype_markers_raw, function(genes) {
  intersect(genes, rownames(data))
})

subtype_markers <- subtype_markers[lengths(subtype_markers) >= 2]

message("Subtype marker sets used:")
print(subtype_markers)

# Calculate UCell scores for subtypes
data <- AddModuleScore_UCell(
  obj = data,
  features = subtype_markers,
  name = '_minor'
)

subtype_ucell_cols <- paste0(names(subtype_markers), "_minor")
subtype_ucell_cols <- intersect(subtype_ucell_cols, colnames(data@meta.data))

# =============================================================================
# 4. FUNCTION TO ASSIGN SUBTYPE WITHIN A GIVEN MAJOR TYPE
# =============================================================================

assign_subtype_within_major <- function(
    metadata,
    allowed_cols,
    ambiguous_label,
    min_score = 0.13,
    min_delta = 0.02
) {
  
  allowed_cols <- intersect(allowed_cols, colnames(metadata))
  
  if (length(allowed_cols) == 0) {
    return(rep(ambiguous_label, nrow(metadata)))
  }
  
  score_mat <- as.matrix(metadata[, allowed_cols, drop = FALSE])
  
  top_score <- apply(score_mat, 1, max, na.rm = TRUE)
  
  second_score <- apply(score_mat, 1, function(x) {
    sx <- sort(x, decreasing = TRUE)
    if (length(sx) >= 2) sx[2] else NA_real_
  })
  
  delta_score <- top_score - second_score
  
  top_label <- colnames(score_mat)[
    max.col(score_mat, ties.method = "first")
  ]
  
  top_label <- gsub("_UCell$", "", top_label)
  top_label <- gsub("_", " ", top_label)
  top_label <- tools::toTitleCase(top_label)
  
  assigned <- ifelse(
    top_score >= min_score &
      (is.na(delta_score) | delta_score >= min_delta),
    top_label,
    ambiguous_label
  )
  
  return(assigned)
}

# =============================================================================
# 5. ASSIGN SUBTYPES WITHIN MAJOR TYPES
# =============================================================================

data$celltype_subtype <- "Uncertain"

# --------------------------
# Progenitor subtypes
# --------------------------

progenitor_cols <- c(
  "Neuroepithelial_Radial_Glia_UCell",
  "Cycling_Progenitor_UCell",
  "Outer_Radial_Glia_Like_UCell"
)

idx <- data$celltype_major == "Progenitor"

data$celltype_subtype[idx] <- assign_subtype_within_major(
  metadata = data@meta.data[idx, , drop = FALSE],
  allowed_cols = progenitor_cols,
  ambiguous_label = "Ambiguous progenitor",
  min_score = 0.13,
  min_delta = 0.02
)

# --------------------------
# Intermediate progenitor
# --------------------------

idx <- data$celltype_major == "Intermediate progenitor"

data$celltype_subtype[idx] <- "Intermediate progenitor"

# --------------------------
# Neuronal subtypes
# --------------------------

neuronal_cols <- c(
  "Neuroblast_Immature_Neuron_UCell",
  "Excitatory_Neuron_UCell",
  "Deep_Layer_Excitatory_Neuron_UCell",
  "Upper_Layer_Excitatory_Neuron_UCell"
)

idx <- data$celltype_major == "Neuronal"

data$celltype_subtype[idx] <- assign_subtype_within_major(
  metadata = data@meta.data[idx, , drop = FALSE],
  allowed_cols = neuronal_cols,
  ambiguous_label = "Ambiguous neuronal",
  min_score = 0.13,
  min_delta = 0.02
)

# --------------------------
# Astroglia-like
# --------------------------

idx <- data$celltype_major == "Astroglia-like"

data$celltype_subtype[idx] <- "Astroglia-like"

# --------------------------
# Choroid/Hem-like subtypes
# --------------------------

choroid_hem_cols <- c(
  "Choroid_Plexus_Like_UCell",
  "Hem_Choroid_Patterning_Like_UCell"
)

idx <- data$celltype_major == "Choroid/Hem-like"

data$celltype_subtype[idx] <- assign_subtype_within_major(
  metadata = data@meta.data[idx, , drop = FALSE],
  allowed_cols = choroid_hem_cols,
  ambiguous_label = "Ambiguous choroid/hem-like",
  min_score = 0.13,
  min_delta = 0.025
)

# --------------------------
# Off-target/non-neural subtypes
# --------------------------

offtarget_cols <- c(
  "Mesenchyme_Like_UCell",
  "Endothelial_Like_UCell"
)

idx <- data$celltype_major == "Off-target/non-neural"

data$celltype_subtype[idx] <- assign_subtype_within_major(
  metadata = data@meta.data[idx, , drop = FALSE],
  allowed_cols = offtarget_cols,
  ambiguous_label = "Ambiguous off-target",
  min_score = 0.13,
  min_delta = 0.02
)

# --------------------------
# Clean subtype names
# --------------------------

data$celltype_subtype <- dplyr::recode(
  data$celltype_subtype,
  "Neuroepithelial Radial Glia" = "Neuroepithelial radial glia",
  "Cycling Progenitor" = "Cycling progenitor",
  "Outer Radial Glia Like" = "Outer radial glia-like",
  "Intermediate Progenitor" = "Intermediate progenitor",
  "Neuroblast Immature Neuron" = "Neuroblast / immature neuron",
  "Excitatory Neuron" = "Excitatory neuron",
  "Deep Layer Excitatory Neuron" = "Deep-layer excitatory neuron",
  "Upper Layer Excitatory Neuron" = "Upper-layer excitatory neuron",
  "Astroglia Like" = "Astroglia-like",
  "Choroid Plexus Like" = "Choroid plexus-like",
  "Hem Choroid Patterning Like" = "Hem/choroid patterning-like",
  "Mesenchyme Like" = "Mesenchyme-like",
  "Endothelial Like" = "Endothelial-like"
)

# =============================================================================
# 6. COUNTS AND SUMMARY TABLES
# =============================================================================

major_counts <- count(data@meta.data, celltype_major) %>%
  arrange(desc(n))

subtype_counts <- count(data@meta.data, celltype_major, celltype_subtype) %>%
  arrange(celltype_major, desc(n))

message("Major counts:")
print(major_counts)

message("Subtype counts:")
print(subtype_counts)

write.csv(
  major_counts,
  file = file.path(outdir, "celltype_major_counts.csv"),
  row.names = FALSE
)

write.csv(
  subtype_counts,
  file = file.path(outdir, "celltype_major_subtype_counts.csv"),
  row.names = FALSE
)

# =============================================================================
# 7. UMAP PLOTS: MAJOR AND SUBTYPE
# =============================================================================

p_major <- DimPlot(
  object = data,
  reduction = "umap",
  group.by = "celltype_major",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Major cell types",
    subtitle = "First-pass UCell classification using broad dorsal cortical organoid signatures",
    colour = "Major type"
  ) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p_subtype <- DimPlot(
  object = data,
  reduction = "umap",
  group.by = "celltype_subtype",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Cell type subtypes",
    subtitle = "Second-pass UCell classification within each major type",
    colour = "Subtype"
  ) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p_annotation_panel <- p_major + p_subtype +
  plot_layout(ncol = 2)

print(p_major)
print(p_subtype)
print(p_annotation_panel)

ggsave(
  filename = file.path(outdir, "umap_celltype_major_ucell_first_pass.png"),
  plot = p_major,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "umap_celltype_major_ucell_first_pass.pdf"),
  plot = p_major,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(outdir, "umap_celltype_subtype_ucell_second_pass.png"),
  plot = p_subtype,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "umap_celltype_subtype_ucell_second_pass.pdf"),
  plot = p_subtype,
  width = 12,
  height = 9
)

ggsave(
  filename = file.path(outdir, "umap_celltype_major_and_subtype_ucell.png"),
  plot = p_annotation_panel,
  width = 22,
  height = 9,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(outdir, "umap_celltype_major_and_subtype_ucell.pdf"),
  plot = p_annotation_panel,
  width = 22,
  height = 9,
  limitsize = FALSE
)

# =============================================================================
# 8. VALIDATION HEATMAP: MEAN UCELL SCORES BY ASSIGNED SUBTYPE
# =============================================================================

all_ucell_cols <- unique(c(major_ucell_cols, subtype_ucell_cols))
all_ucell_cols <- intersect(all_ucell_cols, colnames(data@meta.data))

score_summary_subtype <- data@meta.data %>%
  select(celltype_major, celltype_subtype, all_of(all_ucell_cols)) %>%
  pivot_longer(
    cols = all_of(all_ucell_cols),
    names_to = "signature",
    values_to = "score"
  ) %>%
  mutate(
    signature = gsub("_UCell$", "", signature),
    signature = gsub("_", " ", signature),
    signature = tools::toTitleCase(signature)
  ) %>%
  group_by(celltype_major, celltype_subtype, signature) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

p_score_heatmap_subtype <- ggplot(score_summary_subtype, aes(
  x = celltype_subtype,
  y = signature,
  fill = mean_score
)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
    name = "Mean\nUCell score"
  ) +
  facet_grid(. ~ celltype_major, scales = "free_x", space = "free_x") +
  labs(
    title = "Mean UCell score by assigned subtype",
    subtitle = "Validation of two-step major/subtype classification",
    x = "Assigned subtype",
    y = "UCell signature"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

print(p_score_heatmap_subtype)

ggsave(
  filename = file.path(outdir, "heatmap_ucell_scores_by_subtype_two_step.png"),
  plot = p_score_heatmap_subtype,
  width = 18,
  height = 10,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(outdir, "heatmap_ucell_scores_by_subtype_two_step.pdf"),
  plot = p_score_heatmap_subtype,
  width = 18,
  height = 10,
  limitsize = FALSE
)

write.csv(
  score_summary_subtype,
  file = file.path(outdir, "ucell_score_summary_by_subtype.csv"),
  row.names = FALSE
)

# =============================================================================
# 9. CLUSTER-LEVEL SUMMARY
# =============================================================================

data$cluster_id <- Idents(data)

cluster_annotation_summary <- data@meta.data %>%
  count(cluster_id, celltype_major, celltype_subtype) %>%
  group_by(cluster_id) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(cluster_id, desc(freq)) %>%
  ungroup()

write.csv(
  cluster_annotation_summary,
  file = file.path(outdir, "cluster_annotation_summary_major_subtype.csv"),
  row.names = FALSE
)

message("Cluster annotation summary:")
print(cluster_annotation_summary)

# =============================================================================
# 10. PANEL: ONE SUBPLOT PER MAJOR TYPE WITH GLOBAL UMAP CONTOUR
# =============================================================================

umap_df <- Embeddings(data, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(
    cell_id = rownames(.),
    celltype_major = data@meta.data[cell_id, "celltype_major", drop = TRUE],
    celltype_subtype = data@meta.data[cell_id, "celltype_subtype", drop = TRUE]
  )

colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

umap_df <- umap_df %>%
  filter(!is.na(celltype_major), !is.na(celltype_subtype))

major_levels <- c(
  sort(setdiff(unique(umap_df$celltype_major), "Uncertain")),
  "Uncertain"
)

major_levels <- unique(major_levels[major_levels %in% unique(umap_df$celltype_major)])

umap_df$celltype_major <- factor(
  umap_df$celltype_major,
  levels = major_levels
)

p_all_major <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, colour = celltype_major)
) +
  geom_point(size = 0.25, alpha = 0.8) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Major cell types",
    subtitle = "All cells",
    colour = "Major type",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

major_subplots <- lapply(levels(umap_df$celltype_major), function(ct) {
  
  df_ct <- umap_df %>%
    filter(celltype_major == ct)
  
  ggplot() +
    geom_density_2d(
      data = umap_df,
      aes(x = UMAP_1, y = UMAP_2),
      colour = "grey60",
      linewidth = 0.25,
      alpha = 0.7,
      bins = 6
    ) +
    geom_point(
      data = df_ct,
      aes(x = UMAP_1, y = UMAP_2),
      size = 0.25,
      alpha = 0.9
    ) +
    coord_fixed(ratio = 1) +
    labs(
      title = ct,
      subtitle = paste0(nrow(df_ct), " cells"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6)
    )
})

p_major_subplots <- wrap_plots(major_subplots, ncol = 4)

final_major_panel <- p_all_major / p_major_subplots +
  plot_layout(heights = c(1.2, 3))

print(final_major_panel)

ggsave(
  filename = file.path(outdir, "umap_major_celltypes_global_and_subplots.png"),
  plot = final_major_panel,
  width = 24,
  height = 28,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(outdir, "umap_major_celltypes_global_and_subplots.pdf"),
  plot = final_major_panel,
  width = 24,
  height = 28,
  limitsize = FALSE
)


