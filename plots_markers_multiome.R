library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

# -----------------------------
# 0. Define working environment
# -----------------------------

general_env <- "integrated_multiome/"
plots_env <- "plots/"

setwd(general_env)

if (!dir.exists(plots_env)) {
  dir.create(plots_env, recursive = TRUE)
}

getwd()

# -----------------------------
# 1. Read data
# -----------------------------

integrated_multiome <- readRDS("integratedmultiome_processed.rds")

DefaultAssay(integrated_multiome) <- "RNA"

# Seurat v5: join RNA layers if needed
integrated_multiome[["RNA"]] <- JoinLayers(integrated_multiome[["RNA"]])

# Scale RNA data
integrated_multiome <- ScaleData(
  integrated_multiome,
  assay = "RNA",
  features = rownames(integrated_multiome)
)

# -----------------------------
# 2. Set cluster identity
[# -----------------------------

Idents(integrated_multiome) <- "seurat_clusters"

# Keep cluster names as character to avoid factor/level problems later
integrated_multiome$seurat_clusters <- as.character(integrated_multiome$seurat_clusters)

# -----------------------------
# 3. Major cell type markers
# Logical order:
# progenitors → neurogenic progenitors → neurons → glia → other
# -----------------------------

major_markers <- list(
  "Radial glia / progenitors" = c(
    "SOX2", "PAX6", "VIM", "NES", "HES1", "FABP7"
  ),
  "Outer / truncated radial glia" = c(
    "HOPX", "MOXD1", "PTN", "ID4", "SAMD4A", "CCN2"
  ),
  "Intermediate progenitors" = c(
    "EOMES", "NEUROG1", "NEUROG2", "SOX4", "SOX11"
  ),
  "Neurons - general" = c(
    "DCX", "STMN2", "MAP2", "GAP43", "SNAP25", "RBFOX3"
  ),
  "Excitatory neurons" = c(
    "SLC17A7", "SLC17A6", "TBR1", "BCL11B", "SATB2"
  ),
  "Inhibitory neurons" = c(
    "GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "ASCL1"
  ),
  "Ependymal / ciliated" = c(
    "FOXJ1", "GMNC", "TPPP3", "PIFO"
  ),
  "Astro / glial" = c(
    "GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A3"
  ),
  "Oligodendrocyte lineage" = c(
    "OLIG1", "OLIG2", "SOX10", "PDGFRA", "MBP", "PLP1"
  ),
  "Microglia / immune" = c(
    "P2RY12", "CX3CR1", "AIF1", "C1QA", "TYROBP"
  ),
  "Cycling" = c(
    "MKI67", "TOP2A", "PCNA", "HMGB2"
  )
)

major_markers_vec <- unlist(major_markers, use.names = FALSE)
major_markers_vec <- unique(major_markers_vec)
major_markers_vec <- major_markers_vec[major_markers_vec %in% rownames(integrated_multiome)]

# -----------------------------
# 4. Subtype markers
# Logical order:
# RG states → IPC/neurogenesis → neurons → interneurons → glia/other
# -----------------------------

subtype_markers <- list(
  "Ventricular radial glia / neuroepithelial" = c(
    "SOX2", "PAX6", "VIM", "NES", "HES1", "HES5", "FABP7"
  ),
  "Outer radial glia" = c(
    "HOPX", "MOXD1", "PTN", "ID4", "FAM107A"
  ),
  "Truncated / gliogenic radial glia" = c(
    "SAMD4A", "CCN2", "TNC", "CLU"
  ),
  "IPC / neurogenic progenitors" = c(
    "EOMES", "NEUROG1", "NEUROG2", "INSM1", "SOX4", "SOX11"
  ),
  "Newborn neurons" = c(
    "DCX", "STMN2", "GAP43", "TUBB3", "NEUROD1", "NEUROD2"
  ),
  "Excitatory neurons - deep layer" = c(
    "TBR1", "BCL11B", "FEZF2", "SOX5"
  ),
  "Excitatory neurons - upper layer" = c(
    "SATB2", "CUX1", "CUX2", "POU3F2", "RORB"
  ),
  "Excitatory neurons - general / maturing" = c(
    "SLC17A7", "SLC17A6", "MAP2", "SNAP25", "RBFOX3", "SYT1"
  ),
  "Interneurons - general" = c(
    "GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "DLX5", "ASCL1"
  ),
  "Interneurons - MGE-like" = c(
    "LHX6", "NKX2-1", "SST", "PVALB"
  ),
  "Interneurons - CGE-like" = c(
    "NR2F2", "SP8", "PROX1", "RELN", "VIP"
  ),
  "Astrocyte-like" = c(
    "GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A3", "SLC1A2"
  ),
  "Oligodendrocyte lineage" = c(
    "OLIG1", "OLIG2", "SOX10", "PDGFRA", "CSPG4", "MBP", "PLP1"
  ),
  "Ependymal / ciliated" = c(
    "FOXJ1", "GMNC", "TPPP3", "PIFO", "DNAH5"
  ),
  "Microglia / immune" = c(
    "P2RY12", "CX3CR1", "AIF1", "C1QA", "C1QB", "TYROBP"
  ),
  "Cycling / mitotic" = c(
    "MKI67", "TOP2A", "PCNA", "HMGB2", "UBE2C", "CENPF"
  )
)

subtype_markers_vec <- unlist(subtype_markers, use.names = FALSE)
subtype_markers_vec <- unique(subtype_markers_vec)
subtype_markers_vec <- subtype_markers_vec[subtype_markers_vec %in% rownames(integrated_multiome)]

# -----------------------------
# 5. All markers used for cluster similarity ordering
# -----------------------------

all_marker_genes <- unique(c(major_markers_vec, subtype_markers_vec))
all_marker_genes <- all_marker_genes[all_marker_genes %in% rownames(integrated_multiome)]

if (length(all_marker_genes) < 2) {
  stop("Too few marker genes found in the object. Check gene names and rownames(integrated_multiome).")
}

# -----------------------------
# 6. Reorder clusters by similarity
# Important fix:
# AverageExpression can rename numeric clusters 0, 1, 2 as g0, g1, g2.
# We convert them back before applying factor levels.
# -----------------------------

Idents(integrated_multiome) <- "seurat_clusters"

avg_expr <- AverageExpression(
  object = integrated_multiome,
  assays = "RNA",
  features = all_marker_genes,
  group.by = "seurat_clusters",
  slot = "data"
)$RNA

# Fix Seurat renaming of numeric identities: g0 -> 0, g1 -> 1, etc.
colnames(avg_expr) <- sub("^g(?=[0-9])", "", colnames(avg_expr), perl = TRUE)

# Remove genes with zero variance across clusters
avg_expr <- avg_expr[apply(avg_expr, 1, var) > 0, ]

if (nrow(avg_expr) < 2) {
  stop("Too few variable marker genes across clusters to compute similarity ordering.")
}

# Scale each gene across clusters
avg_expr_scaled <- t(scale(t(avg_expr)))

# Replace possible NA values after scaling
avg_expr_scaled[is.na(avg_expr_scaled)] <- 0

# Hierarchical clustering of clusters
cluster_dist <- dist(t(avg_expr_scaled), method = "euclidean")
cluster_hclust <- hclust(cluster_dist, method = "ward.D2")

cluster_order_similarity <- colnames(avg_expr_scaled)[cluster_hclust$order]

print("Cluster order by similarity:")
print(cluster_order_similarity)

print("Clusters present in metadata:")
print(sort(unique(integrated_multiome$seurat_clusters)))

# Safety check before applying factor levels
missing_clusters <- setdiff(unique(integrated_multiome$seurat_clusters), cluster_order_similarity)

if (length(missing_clusters) > 0) {
  warning(
    "Some clusters in metadata are not present in cluster_order_similarity: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# Apply new cluster order to metadata
integrated_multiome$seurat_clusters <- factor(
  integrated_multiome$seurat_clusters,
  levels = cluster_order_similarity
)

Idents(integrated_multiome) <- "seurat_clusters"

print("Cluster table after applying similarity order:")
print(table(integrated_multiome$seurat_clusters, useNA = "ifany"))

if (any(is.na(integrated_multiome$seurat_clusters))) {
  stop("Some cells became NA after applying cluster order. Check cluster labels before and after AverageExpression.")
}

# -----------------------------
# 7. Helper function to add marker group labels
# -----------------------------

make_marker_annotation <- function(marker_list, present_markers) {
  data.frame(
    marker = unlist(marker_list, use.names = FALSE),
    group = rep(names(marker_list), lengths(marker_list))
  ) %>%
    filter(marker %in% present_markers) %>%
    distinct(marker, .keep_all = TRUE)
}

major_annot <- make_marker_annotation(major_markers, major_markers_vec)
subtype_annot <- make_marker_annotation(subtype_markers, subtype_markers_vec)

# -----------------------------
# 8. Major cell type DotPlot
# -----------------------------

p_major <- DotPlot(
  integrated_multiome,
  features = major_markers_vec,
  group.by = "seurat_clusters",
  assay = "RNA"
) +
  scale_color_gradientn(
    colors = c("grey90", "lightskyblue", "dodgerblue3", "navy"),
    name = "Average\nexpression"
  ) +
  scale_size(
    range = c(0.5, 6),
    name = "% cells\nexpressing"
  ) +
  labs(
    title = "Major cell type markers by cluster",
    subtitle = "Clusters ordered by similarity based on marker expression",
    x = "Markers",
    y = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0),
    panel.grid = element_line(color = "grey92"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

# -----------------------------
# 9. Subtype DotPlot
# -----------------------------

p_subtype <- DotPlot(
  integrated_multiome,
  features = subtype_markers_vec,
  group.by = "seurat_clusters",
  assay = "RNA"
) +
  scale_color_gradientn(
    colors = c("grey90", "lightskyblue", "dodgerblue3", "navy"),
    name = "Average\nexpression"
  ) +
  scale_size(
    range = c(0.5, 6),
    name = "% cells\nexpressing"
  ) +
  labs(
    title = "Subtype markers by cluster",
    subtitle = "Clusters ordered by similarity based on marker expression",
    x = "Markers",
    y = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0),
    panel.grid = element_line(color = "grey92"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

# -----------------------------
# 10. Split major markers visually by biological groups
# -----------------------------

major_data <- p_major$data %>%
  left_join(major_annot, by = c("features.plot" = "marker")) %>%
  mutate(
    group = factor(group, levels = names(major_markers)),
    features.plot = factor(features.plot, levels = major_markers_vec),
    id = factor(as.character(id), levels = cluster_order_similarity)
  )

p_major_grouped <- ggplot(
  major_data,
  aes(
    x = features.plot,
    y = id,
    size = pct.exp,
    color = avg.exp.scaled
  )
) +
  geom_point() +
  facet_grid(
    . ~ group,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_gradientn(
    colors = c("grey90", "lightskyblue", "dodgerblue3", "navy"),
    name = "Scaled average\nexpression"
  ) +
  scale_size(
    range = c(0.5, 6),
    name = "% cells\nexpressing"
  ) +
  labs(
    title = "Major cell type markers by cluster",
    subtitle = "Clusters ordered by similarity based on marker expression",
    x = "Markers",
    y = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 9),
    strip.text.x = element_text(
      face = "bold",
      size = 8,
      angle = 0
    ),
    plot.title = element_text(face = "bold", hjust = 0),
    panel.grid = element_line(color = "grey92"),
    legend.position = "right"
  )

# -----------------------------
# 11. Split subtype markers visually by biological groups
# -----------------------------

subtype_data <- p_subtype$data %>%
  left_join(subtype_annot, by = c("features.plot" = "marker")) %>%
  mutate(
    group = factor(group, levels = names(subtype_markers)),
    features.plot = factor(features.plot, levels = subtype_markers_vec),
    id = factor(as.character(id), levels = cluster_order_similarity)
  )

p_subtype_grouped <- ggplot(
  subtype_data,
  aes(
    x = features.plot,
    y = id,
    size = pct.exp,
    color = avg.exp.scaled
  )
) +
  geom_point() +
  facet_grid(
    . ~ group,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_gradientn(
    colors = c("grey90", "lightskyblue", "dodgerblue3", "navy"),
    name = "Scaled average\nexpression"
  ) +
  scale_size(
    range = c(0.5, 6),
    name = "% cells\nexpressing"
  ) +
  labs(
    title = "Subtype markers by cluster",
    subtitle = "Clusters ordered by similarity based on marker expression",
    x = "Markers",
    y = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 9),
    strip.text.x = element_text(
      face = "bold",
      size = 8,
      angle = 0
    ),
    plot.title = element_text(face = "bold", hjust = 0),
    panel.grid = element_line(color = "grey92"),
    legend.position = "right"
  )

# -----------------------------
# 12. Print plots
# -----------------------------

p_major_grouped
p_subtype_grouped

# -----------------------------
# 13. Save plots
# -----------------------------

ggsave(
  filename = file.path(plots_env, "major_celltype_markers_by_cluster_similarity_ordered.pdf"),
  plot = p_major_grouped,
  width = 18,
  height = 7,
  useDingbats = FALSE
)

ggsave(
  filename = file.path(plots_env, "subtype_markers_by_cluster_similarity_ordered.pdf"),
  plot = p_subtype_grouped,
  width = 24,
  height = 7,
  useDingbats = FALSE
)

ggsave(
  filename = file.path(plots_env, "major_celltype_markers_by_cluster_similarity_ordered.png"),
  plot = p_major_grouped,
  width = 18,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(plots_env, "subtype_markers_by_cluster_similarity_ordered.png"),
  plot = p_subtype_grouped,
  width = 24,
  height = 7,
  dpi = 300
)

# -----------------------------
# 14. Also save basic Seurat-style DotPlots, non-faceted
# -----------------------------

p_major_basic <- p_major
p_subtype_basic <- p_subtype

ggsave(
  filename = file.path(plots_env, "basic_major_celltype_markers_by_cluster_similarity_ordered.pdf"),
  plot = p_major_basic,
  width = 18,
  height = 7,
  useDingbats = FALSE
)

ggsave(
  filename = file.path(plots_env, "basic_subtype_markers_by_cluster_similarity_ordered.pdf"),
  plot = p_subtype_basic,
  width = 24,
  height = 7,
  useDingbats = FALSE
)

ggsave(
  filename = file.path(plots_env, "basic_major_celltype_markers_by_cluster_similarity_ordered.png"),
  plot = p_major_basic,
  width = 18,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(plots_env, "basic_subtype_markers_by_cluster_similarity_ordered.png"),
  plot = p_subtype_basic,
  width = 24,
  height = 7,
  dpi = 300
)