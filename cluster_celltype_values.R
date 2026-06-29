library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(grid)
library(tidyverse)

data <- readRDS("integratedmultiome_processed.rds")

DefaultAssay(data) <- "RNA"
data[["RNA"]] <- JoinLayers(data[["RNA"]])

# -------------------------------------------------------------------------
# Marker genes per cell type


marker_features_raw <- list( Telencephalon = c( "FOXG1", "NFIA", "NFIB", "EMX1", "EMX2", "LHX2", "SIX3", "SIX6", "PAX6", "OTX1", "OTX2" ), Cell_Cycle = c( "MKI67", "ASPM", "TOP2A", "CENPF", "UBE2C", "CDK1", "CCNB1", "CCNB2", "AURKB", "BIRC5", "MCM2", "MCM5", "PCNA" ), Radial_Glia = c( "PAX6", "HES1", "GLI3", "SOX2", "SOX9", "VIM", "NES", "FABP7", "HES5", "HMGA2", "SLC1A3", "PROM1", "LIX1" ), Outer_Radial_Glia = c( "HOPX", "TNC", "PTPRZ1", "MOXD1", "FAM107A", "LIFR", "ITGB5", "PDGFD", "CRYAB", "CLU" ), Radial_Glia_Astrocyte = c( "TNC", "SLC1A3", "HOPX", "SOX9", "VIM", "FABP7", "CLU", "AQP4", "GFAP", "ALDH1L1", "SPARCL1", "GJA1" ), Astrocyte = c( "IL33", "AQP4", "GFAP", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1", "GJB6", "CLU", "SPARCL1", "APOE", "SOX9", "FGFR3", "EDNRB" ), Intermediate_Progenitors = c( "EOMES", "PPP1R17", "KCNQ3", "NHLH1", "TBR2", "NEUROG1", "NEUROG2", "INSM1", "ELAVL4", "BTG2", "HES6", "DLL3", "NEUROD1", "NEUROD4", "TOX3" ), Neuroblasts = c( "DCX", "STMN2", "STMN1", "TUBB3", "NEUROD1", "NEUROD2", "NEUROD6", "ELAVL3", "ELAVL4", "SOX11", "SOX4", "GAP43", "MAP1B" ), Excitatory_Neuron = c( "NEUROD6", "SLA", "NRP1", "SLC17A6", "SLC17A7", "TBR1", "BCL11B", "FEZF2", "SATB2", "CUX1", "CUX2", "RORB", "FOXP1", "FOXP2", "CAMK2A", "RBFOX3", "MAP2", "SNAP25", "SYT1" ), Maturing_Excitatory_Neuron = c( "RELN", "MEF2C", "SNAP25", "FGF12", "FEZF2", "BCL11B", "MAP2", "RBFOX3", "SYT1", "SYN1", "CAMK2A", "GRIN1", "GRIN2B", "DLG4", "CNTNAP2", "NRXN1", "NLGN1" ), Deep_Layer_Excitatory_Neuron = c( "TBR1", "BCL11B", "FEZF2", "SOX5", "FOXP2", "FOXP1", "ZFPM2", "NTNG2", "CRYM", "LDB2" ), Upper_Layer_Excitatory_Neuron = c( "SATB2", "CUX1", "CUX2", "RORB", "POU3F2", "POU3F3", "LHX2", "MEF2C", "CPNE4" ), Inhibitory_Neuron = c( "DLX1", "DLX2", "DLX6-AS1", "GAD1", "GAD2", "SLC32A1", "VGAT", "LHX6", "LHX8", "ARX", "SP8", "SP9", "ISL1", "ERBB4", "SCG2", "CALB2", "SST", "PVALB", "VIP" ), Interneuron_MGE = c( "LHX6", "LHX8", "NKX2-1", "SOX6", "MAF", "MAFB", "SST", "PVALB", "GAD1", "GAD2" ), Interneuron_CGE = c( "NR2F1", "NR2F2", "SP8", "SP9", "PROX1", "VIP", "CALB2", "RELN", "GAD1", "GAD2" ), Choroid = c( "OTX2", "TTR", "AQP1", "KRT8", "KRT18", "KRT19", "CLIC6", "FOLR1", "IGFBP2", "RSPO2" ), Hindbrain = c( "RSPO2", "HOXA3", "HOXB3", "GBX2", "HOXA2", "HOXB2", "HOXA4", "HOXB4", "KROX20", "EGR2", "PAX2", "PAX7" ), Retina = c( "VSX2", "RORB", "PAX6", "SIX3", "SIX6", "LHX2", "CRX", "OTX2", "NRL", "RHO", "ARR3", "POU4F2", "ATOH7" ), Mesenchyme = c( "DCN", "BGN", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "LUM", "PDGFRA", "VIM", "TAGLN", "ACTA2", "THY1", "MMP2" ), Neural_Crest = c( "SOX10", "FOXD3", "TFAP2A", "TFAP2B", "NGFR", "ERBB3", "MPZ", "PLP1", "SNAI1", "SNAI2", "EDNRB" ), Endothelial = c( "PECAM1", "VWF", "KDR", "FLT1", "CLDN5", "ESAM", "ENG", "RAMP2", "EMCN", "CDH5" ), Microglia = c( "CX3CR1", "P2RY12", "AIF1", "TYROBP", "C1QA", "C1QB", "C1QC", "CSF1R", "TMEM119", "SPI1", "LST1" ), Oligodendrocyte_Lineage = c( "OLIG1", "OLIG2", "SOX10", "PDGFRA", "CSPG4", "NKX2-2", "BCAS1", "MBP", "MOG", "PLP1", "MAG" ) )

# -------------------------------------------------------------------------

# Keep only genes present in the Seurat object
marker_features <- lapply(marker_features_raw, function(genes) {
  intersect(genes, rownames(data))
})

# Remove empty marker sets
marker_features <- marker_features[lengths(marker_features) > 0]

# Optional: print genes actually used per marker set
print(marker_features)

# Add module scores
data <- AddModuleScore(
  object = data,
  features = marker_features,
  name = "CellTypeScore"
)

# Rename score columns
score_cols <- paste0("CellTypeScore", seq_along(marker_features))
new_score_cols <- paste0(names(marker_features), "_Score")

colnames(data@meta.data)[
  match(score_cols, colnames(data@meta.data))
] <- new_score_cols

# -------------------------------------------------------------------------
# Create one FeaturePlot per cell type
# - no cluster labels
# - high-score cells plotted on top
# - square ratio
# - smaller colour scale
# - bigger titles
# - subtitles split into 2 lines
# -------------------------------------------------------------------------

plots <- lapply(seq_along(new_score_cols), function(i) {
  
  score_col <- new_score_cols[i]
  
  plot_title <- names(marker_features)[i] |>
    gsub("_", " ", x = _) |>
    tools::toTitleCase()
  
  genes_used <- marker_features[[i]]
  
  # Split marker genes into exactly 2 lines
  if (length(genes_used) > 1) {
    midpoint <- ceiling(length(genes_used) / 2)
    
    plot_subtitle <- paste(
      paste(genes_used[1:midpoint], collapse = ", "),
      paste(genes_used[(midpoint + 1):length(genes_used)], collapse = ", "),
      sep = "\n"
    )
  } else {
    plot_subtitle <- genes_used
  }
  
  FeaturePlot(
    object = data,
    features = score_col,
    label = FALSE,
    order = TRUE,
    raster = FALSE
  ) +
    scale_colour_gradientn(
      colours = rev(brewer.pal(n = 11, name = "RdBu"))
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      colour = "Module score"
    ) +
    coord_fixed(ratio = 1) +
    guides(
      colour = guide_colourbar(
        barheight = unit(2.2, "cm"),
        barwidth = unit(0.22, "cm"),
        title.position = "top"
      )
    ) +
    theme(
      aspect.ratio = 1,
      
      plot.title = element_text(
        size = 40,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 20,
        hjust = 0.5,
        lineheight = 0.9
      ),
      
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6)
    )
})

# Combine plots
combined_plot <- wrap_plots(plots, ncol = 6)

combined_plot

# -------------------------------------------------------------------------
# Save plot
# -------------------------------------------------------------------------

outdir <- "multiome_plts"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(outdir, "celltype_module_scores_pollen_markers.png"),
  plot = combined_plot,
  width = 50,
  height = 50,
  dpi = 300,
  limitsize = FALSE
)


# -------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(forcats)

# -------------------------------------------------------------------------
# Define cluster column
# -------------------------------------------------------------------------

# Option 1: use active identities
data$cluster_id <- Idents(data)

# Option 2: use a metadata column instead, for example:
# data$cluster_id <- data$seurat_clusters

# Check score columns
score_cols <- grep("_Score$", colnames(data@meta.data), value = TRUE)

score_cols


score_df <- FetchData(
  object = data,
  vars = c("cluster_id", score_cols)
)

score_long <- score_df %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "marker_set",
    values_to = "score"
  ) %>%
  mutate(
    marker_set = gsub("_Score$", "", marker_set),
    marker_set = gsub("_", " ", marker_set),
    marker_set = tools::toTitleCase(marker_set)
  )

avg_score_cluster <- score_long %>%
  group_by(cluster_id, marker_set) %>%
  summarise(
    avg_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    pct_score_positive = mean(score > 0, na.rm = TRUE) * 100,
    n_cells = n(),
    .groups = "drop"
  )

head(avg_score_cluster)


ggplot(avg_score_cluster, aes(
  x = cluster_id,
  y = marker_set
)) +
  geom_point(aes(
    size = pct_score_positive,
    color = avg_score
  )) +
  scale_color_gradientn(
    colours = rev(brewer.pal(n = 11, name = "RdBu")),
    name = "Average\nmodule score"
  ) +
  scale_size(
    range = c(0.5, 7),
    name = "% cells\nscore > 0"
  ) +
  labs(
    title = "Cell-type marker module scores per cluster",
    x = "Cluster",
    y = "Marker set"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )



ggsave(
  filename = file.path(outdir, "plot_score_markers_cluster.png"),
  width = 15,
  height = 10,
  dpi = 300,
  limitsize = FALSE
)





############## UCELL ###############


library(UCell)
library(Seurat)
library(dplyr)

DefaultAssay(data) <- "RNA"

marker_features <- lapply(marker_features_raw, function(genes) {
  intersect(genes, rownames(data))
})

marker_features <- marker_features[lengths(marker_features) >= 2]

data <- AddModuleScore_UCell(
  obj = data,
  features = marker_features,
  name =  "_UCell"
)

ucell_cols <- paste0(names(marker_features), "_UCell")



score_mat <- as.matrix(data@meta.data[, ucell_cols, drop = FALSE])

top_score <- apply(score_mat, 1, max, na.rm = TRUE)
second_score <- apply(score_mat, 1, function(x) sort(x, decreasing = TRUE)[2])
top_label <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]

data$celltype_ucell <- gsub("_UCell$", "", top_label)
data$celltype_ucell_score <- top_score
data$celltype_ucell_delta <- top_score - second_score

data$celltype_ucell_flag <- ifelse(
  data$celltype_ucell_score >= 0.15 &
    data$celltype_ucell_delta >= 0.03,
  data$celltype_ucell,
  "Uncertain"
)
count(data@meta.data,celltype_ucell_flag)


library(Seurat)
library(ggplot2)

p_ucell_celltypes <- DimPlot(
  object = data,
  reduction = "umap",
  group.by = "celltype_ucell_flag",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  coord_fixed(ratio = 1) +
  labs(
    title = "UCell-based cell type flags",
    subtitle = "Cell type assignment based on highest UCell marker signature",
    colour = "Cell type"
  ) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p_ucell_celltypes
outdir <- "multiome_plts"



ggsave(
  filename = file.path(outdir, "umap_ucell_celltype_flags.png"),
  plot = p_ucell_celltypes,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "umap_ucell_celltype_flags.pdf"),
  plot = p_ucell_celltypes,
  width = 10,
  height = 8
)


# Plot UCell scores for top cell type per cluster



library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# -------------------------------------------------------------------------
# Extract UMAP coordinates + celltype flags
# -------------------------------------------------------------------------

umap_df <- Embeddings(data, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(
    cell_id = rownames(.),
    celltype_ucell_flag = data@meta.data[cell_id, "celltype_ucell_flag", drop = TRUE]
  )

# Make sure column names are simple
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

# Remove cells without assignment
umap_df <- umap_df %>%
  filter(!is.na(celltype_ucell_flag))

# Optional: put Uncertain at the end
umap_df <- umap_df %>%
  mutate(
    celltype_ucell_flag = as.character(celltype_ucell_flag),
    celltype_ucell_flag = factor(
      celltype_ucell_flag,
      levels = c(
        sort(setdiff(unique(celltype_ucell_flag), "Uncertain")),
        "Uncertain"
      )
    )
  )

# -------------------------------------------------------------------------
# Global UMAP with all detected celltypes
# -------------------------------------------------------------------------

p_all_celltypes <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, colour = celltype_ucell_flag)
) +
  geom_point(size = 0.25, alpha = 0.8) +
  coord_fixed(ratio = 1) +
  labs(
    title = "UCell-based cell type flags",
    subtitle = "All detected cell types",
    colour = "Cell type",
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

p_all_celltypes


# -------------------------------------------------------------------------
# One subplot per celltype
# Global density contour + points from only that celltype
# -------------------------------------------------------------------------

celltypes_to_plot <- levels(umap_df$celltype_ucell_flag)

celltype_plots <- lapply(celltypes_to_plot, function(ct) {
  
  df_ct <- umap_df %>%
    filter(celltype_ucell_flag == ct)
  
  ggplot() +
    # Global contour: all cells
    geom_density_2d(
      data = umap_df,
      aes(x = UMAP_1, y = UMAP_2),
      colour = "grey60",
      linewidth = 0.25,
      alpha = 0.7,
      bins = 6
    ) +
    # Points from selected celltype only
    geom_point(
      data = df_ct,
      aes(x = UMAP_1, y = UMAP_2),
      size = 0.25,
      alpha = 0.9
    ) +
    coord_fixed(ratio = 1) +
    labs(
      title = tools::toTitleCase(gsub("_", " ", ct)),
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

p_celltype_subplots <- wrap_plots(celltype_plots, ncol = 5)

p_celltype_subplots# -------------------------------------------------------------------------
# One subplot per celltype
# Global density contour + points from only that celltype
# -------------------------------------------------------------------------

celltypes_to_plot <- levels(umap_df$celltype_ucell_flag)

celltype_plots <- lapply(celltypes_to_plot, function(ct) {
  
  df_ct <- umap_df %>%
    filter(celltype_ucell_flag == ct)
  
  ggplot() +
    # Global contour: all cells
    geom_density_2d(
      data = umap_df,
      aes(x = UMAP_1, y = UMAP_2),
      colour = "grey60",
      linewidth = 0.25,
      alpha = 0.7,
      bins = 6
    ) +
    # Points from selected celltype only
    geom_point(
      data = df_ct,
      aes(x = UMAP_1, y = UMAP_2),
      size = 0.25,
      alpha = 0.9
    ) +
    coord_fixed(ratio = 1) +
    labs(
      title = tools::toTitleCase(gsub("_", " ", ct)),
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

p_celltype_subplots <- wrap_plots(celltype_plots, ncol = 5)

p_celltype_subplots
##########################################333




# -------------------------------------------------------------------------
# One subplot per celltype
# Global density contour + points from only that celltype
# -------------------------------------------------------------------------

celltypes_to_plot <- levels(umap_df$celltype_ucell_flag)

celltype_plots <- lapply(celltypes_to_plot, function(ct) {
  
  df_ct <- umap_df %>%
    filter(celltype_ucell_flag == ct)
  
  ggplot() +
    # Global contour: all cells
    geom_density_2d(
      data = umap_df,
      aes(x = UMAP_1, y = UMAP_2),
      colour = "grey60",
      linewidth = 0.25,
      alpha = 0.7,
      bins = 6
    ) +
    # Points from selected celltype only
    geom_point(
      data = df_ct,
      aes(x = UMAP_1, y = UMAP_2),
      size = 0.25,
      alpha = 0.9
    ) +
    coord_fixed(ratio = 1) +
    labs(
      title = tools::toTitleCase(gsub("_", " ", ct)),
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

p_celltype_subplots <- wrap_plots(celltype_plots, ncol = 5)

p_celltype_subplots







final_celltype_panel <- p_all_celltypes / p_celltype_subplots +
  plot_layout(heights = c(1.2, 4))

final_celltype_panel