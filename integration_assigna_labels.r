# kanton analysis 

library(Seurat)
library(Signac)
library(tidyverse)
library(Matrix)

counts <- readMM('/users/genomics/josema/scdatasets/kanton/human_cell_counts_consensus.mtx')
metadata <- read_delim('/users/genomics/josema/scdatasets/kanton/metadata_human_cells.tsv') %>%
  as.data.frame() %>% 
  column_to_rownames('...1')
genes <-  read_delim('/users/genomics/josema/scdatasets/kanton/genes_consensus.txt', col_names = F)

rownames(counts) <- genes$X1
colnames(counts) <- rownames(metadata)

kanton <- CreateSeuratObject(counts = counts,meta.data = metadata)

# 1. Cargar la tabla de correspondencia (sin colnames)
colnames(genes) <- c("Ensembl", "GeneName")  # Renombramos las columnas correctamente

# 2. Extraer la matriz de expresión del objeto Seurat
expr_matrix <- GetAssayData(kanton, slot = "counts")  # Usa "data" si prefieres normalizados

# 3. Asegurar correspondencia entre genes en Seurat y la tabla
common_genes <- intersect(rownames(expr_matrix), genes$Ensembl)
conversion_filtered <- genes[genes$Ensembl %in% common_genes, ]

# 4. Identificar los genes que tienen nombres repetidos
gene_counts <- table(conversion_filtered$GeneName)
duplicated_genes <- names(gene_counts[gene_counts > 1])

# 5. Identificar los Ensembl IDs de los genes duplicados
duplicated_ensembls <- conversion_filtered$Ensembl[conversion_filtered$GeneName %in% duplicated_genes]

# 6. Mantener genes únicos y quitar los duplicados
expr_matrix_corrected <- expr_matrix[!(rownames(expr_matrix) %in% duplicated_ensembls), ]

# 7. Para cada gen duplicado, sumar la expresión de sus distintos Ensembl IDs
for (gene in duplicated_genes) {
  ensembl_ids <- conversion_filtered$Ensembl[conversion_filtered$GeneName == gene]
  
  # Si hay más de un Ensembl ID, sumamos sus expresiones
  summed_expression <- colSums(expr_matrix[ensembl_ids, , drop = FALSE])
  
  # Agregamos la nueva fila con el nombre del gen
  expr_matrix_corrected <- rbind(expr_matrix_corrected, summed_expression)
  rownames(expr_matrix_corrected)[nrow(expr_matrix_corrected)] <- gene
}

name_map <- setNames(genes$GeneName, genes$Ensembl)

gene_vector_updated <- ifelse(rownames(expr_matrix_corrected) %in% names(name_map), name_map[rownames(expr_matrix_corrected)], rownames(expr_matrix_corrected))

rownames(expr_matrix_corrected) <-  gene_vector_updated

# 8. Crear un nuevo objeto Seurat con la matriz corregida
seurat_obj_merged <- CreateSeuratObject(counts = expr_matrix_corrected)

# 9. Transferir metadata (opcional)
seurat_obj_merged@meta.data <- kanton@meta.data


# Porcesamiento del kanton

kanton <- NormalizeData(seurat_obj_merged)
kanton <- FindVariableFeatures(kanton, selection.method = "vst", nfeatures = 2000)
kanton <- ScaleData(kanton)
kanton <- RunPCA(kanton, npcs = 30)


# 10. Guardar el objeto corregido
saveRDS(kanton, "/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")




############################## 

options(future.globals.maxSize = 100000000000)  # Adjust max memory if needed
future::plan("multicore", workers = 35)

###########predict annotations batch 1##########
human <- readRDS('/users/genomics/josema/multiome/rmd_files/human_batch1_filtered.rmd')
kanton <- readRDS("/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")

anchors <- FindTransferAnchors(
  reference = kanton,
  query = human,
  dims = 1:30
)

# Transfer cell type labels from kanton to human_batch1
human_batch1 <- TransferData(
  anchorset = anchors,
  refdata = kanton$PredCellType   , 
  dims = 1:30
)

# Store the predicted cell types
human$predicted_cell_type <- human_batch1$predicted.id

# Save the updated Seurat object
saveRDS(human_batch1, "/users/genomics/josema/multiome/rmd_files/predicted_annotations_human_batch_1.rds")

# Check transferred labels
head(human_batch1$predicted_cell_type)

pdf('distribution_celltype_plot.pdf')
DimPlot(human, group.by  = 'predicted_cell_type', label = T)
DimPlot(human, group.by  = 'seurat_clusters', label = T)

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(seurat_clusters) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=percent, fill=predicted_cell_type), position = 'stack')

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(predicted_cell_type) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=predicted_cell_type, y=percent, fill=seurat_clusters), position = 'stack')
dev.off()




###########predict human_chimp d30#########



human <- readRDS('/users/genomics/josema/multiome/CH507_FIPS1_D30_unfiltered.rds')
#kanton <- readRDS("/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")

DefaultAssay(human) <-  'RNA'

anchors <- FindTransferAnchors(
  reference = kanton,
  query = human,
  dims = 1:30
)

# Transfer cell type labels from kanton to human_batch1
human_predictions <- TransferData(
  anchorset = anchors,
  refdata = kanton$PredCellType   , 
  dims = 1:30
)

# Store the predicted cell types
human$predicted_cell_type <- human_predictions$predicted.id

# Save the updated Seurat object
saveRDS(human_predictions, "/users/genomics/josema/multiome/rmd_files/predicted_annotations_CH507_FIPS1_D30_unfiltered.rds")

# Check transferred labels
head(human_batch1$predicted_cell_type)

pdf('distribution_celltype_plot_CH507_FIPS1_D30_unfiltered.pdf')
DimPlot(human, group.by  = 'predicted_cell_type', label = T)
DimPlot(human, group.by  = 'seurat_clusters', label = T)

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(seurat_clusters) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=percent, fill=predicted_cell_type), position = 'stack')

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(predicted_cell_type) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=predicted_cell_type, y=percent, fill=seurat_clusters), position = 'stack')
dev.off()



###########predict human_chimp d90#####


human <- readRDS('/users/genomics/josema/multiome/CH507_FIPS1_D90_unfiltered.rds')
#kanton <- readRDS("/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")
DefaultAssay(human) <-  'RNA'

anchors <- FindTransferAnchors(
  reference = kanton,
  query = human,
  dims = 1:30
)

# Transfer cell type labels from kanton to human_batch1
human_predictions <- TransferData(
  anchorset = anchors,
  refdata = kanton$PredCellType   , 
  dims = 1:30
)

# Store the predicted cell types
human$predicted_cell_type <- human_predictions$predicted.id

# Save the updated Seurat object
saveRDS(human_predictions, "/users/genomics/josema/multiome/rmd_files/predicted_annotations_CH507_FIPS1_D90_unfiltered.rds")

# Check transferred labels
head(human_batch1$predicted_cell_type)

pdf('distribution_celltype_plot_CH507_FIPS1_D90_unfiltered.pdf')
DimPlot(human, group.by  = 'predicted_cell_type', label = T)
DimPlot(human, group.by  = 'seurat_clusters', label = T)

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(seurat_clusters) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=percent, fill=predicted_cell_type), position = 'stack')

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(predicted_cell_type) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=predicted_cell_type, y=percent, fill=seurat_clusters), position = 'stack')
dev.off()




###########predict human_pongo d30#####


human <- readRDS('/users/genomics/josema/multiome/Orango_FIPS2_D30_unfiltered.rds')
#kanton <- readRDS("/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")

DefaultAssay(human) <-  'RNA'

anchors <- FindTransferAnchors(
  reference = kanton,
  query = human,
  dims = 1:30
)

# Transfer cell type labels from kanton to human_batch1
human_predictions <- TransferData(
  anchorset = anchors,
  refdata = kanton$PredCellType   , 
  dims = 1:30
)

# Store the predicted cell types
human$predicted_cell_type <- human_predictions$predicted.id

# Save the updated Seurat object
saveRDS(human_predictions, "/users/genomics/josema/multiome/rmd_files/predicted_annotations_Orango_FIPS2_D30_unfiltered.rds")

# Check transferred labels
head(human_batch1$predicted_cell_type)

pdf('distribution_celltype_plot_Orango_FIPS2_D30_unfiltered.pdf')
DimPlot(human, group.by  = 'predicted_cell_type', label = T)
DimPlot(human, group.by  = 'seurat_clusters', label = T)

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(seurat_clusters) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=percent, fill=predicted_cell_type), position = 'stack')

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(predicted_cell_type) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=predicted_cell_type, y=percent, fill=seurat_clusters), position = 'stack')
dev.off()



###########predict human_pongo d90#####


human <- readRDS('/users/genomics/josema/multiome/Orango_FIPS2_D90_unfiltered.rds')
#kanton <- readRDS("/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")

DefaultAssay(human) <-  'RNA'

anchors <- FindTransferAnchors(
  reference = kanton,
  query = human,
  dims = 1:30
)

# Transfer cell type labels from kanton to human_batch1
human_predictions <- TransferData(
  anchorset = anchors,
  refdata = kanton$PredCellType   , 
  dims = 1:30
)

# Store the predicted cell types
human$predicted_cell_type <- human_predictions$predicted.id

# Save the updated Seurat object
saveRDS(human_predictions, "/users/genomics/josema/multiome/rmd_files/predicted_annotations_Orango_FIPS2_D90_unfiltered.rds")

# Check transferred labels
head(human_batch1$predicted_cell_type)

pdf('distribution_celltype_plot_Orango_FIPS2_D90_unfiltered.pdf')
DimPlot(human, group.by  = 'predicted_cell_type', label = T)
DimPlot(human, group.by  = 'seurat_clusters', label = T)

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(seurat_clusters) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=percent, fill=predicted_cell_type), position = 'stack')

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(predicted_cell_type) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=predicted_cell_type, y=percent, fill=seurat_clusters), position = 'stack')
dev.off()




###########predict gorilla_chimp d30#####


human <- readRDS('/users/genomics/josema/multiome/Gorilla_CHA322_D30_unfiltered.rds')
kanton <- readRDS("/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")

DefaultAssay(human) <-  'RNA'

anchors <- FindTransferAnchors(
  reference = kanton,
  query = human,
  dims = 1:30
)

# Transfer cell type labels from kanton to human_batch1
human_predictions <- TransferData(
  anchorset = anchors,
  refdata = kanton$PredCellType   , 
  dims = 1:30
)

# Store the predicted cell types
human$predicted_cell_type <- human_predictions$predicted.id

# Save the updated Seurat object
saveRDS(human_predictions, "/users/genomics/josema/multiome/rmd_files/predicted_annotations_Gorilla_CHA322_D30_unfiltered.rds")

# Check transferred labels
head(human_batch1$predicted_cell_type)

pdf('distribution_celltype_plot_Gorilla_CHA322_D30_unfiltered.pdf')
DimPlot(human, group.by  = 'predicted_cell_type', label = T)
DimPlot(human, group.by  = 'seurat_clusters', label = T)

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(seurat_clusters) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=percent, fill=predicted_cell_type), position = 'stack')

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(predicted_cell_type) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=predicted_cell_type, y=percent, fill=seurat_clusters), position = 'stack')
dev.off()



###########predict gorilla_chimp d90#####


human <- readRDS('/users/genomics/josema/multiome/Gorilla_CHA322_D90_unfiltered.rds')
#kanton <- readRDS("/users/genomics/josema/scdatasets/kanton/seurat_corrected.rds")

DefaultAssay(human) <-  'RNA'

anchors <- FindTransferAnchors(
  reference = kanton,
  query = human,
  dims = 1:30
)

# Transfer cell type labels from kanton to human_batch1
human_predictions <- TransferData(
  anchorset = anchors,
  refdata = kanton$PredCellType   , 
  dims = 1:30
)

# Store the predicted cell types
human$predicted_cell_type <- human_predictions$predicted.id

# Save the updated Seurat object
saveRDS(human_predictions, "/users/genomics/josema/multiome/rmd_files/predicted_annotations_Gorilla_CHA322_D90_unfiltered.rds")

# Check transferred labels
head(human_batch1$predicted_cell_type)

pdf('distribution_celltype_plot_Gorilla_CHA322_D90_unfiltered.pdf')
DimPlot(human, group.by  = 'predicted_cell_type', label = T)
DimPlot(human, group.by  = 'seurat_clusters', label = T)

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(seurat_clusters) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=percent, fill=predicted_cell_type), position = 'stack')

human@meta.data %>% 
  as.data.frame() %>% 
  count(seurat_clusters, predicted_cell_type) %>% 
  group_by(predicted_cell_type) %>% 
  mutate(percent= n/sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=predicted_cell_type, y=percent, fill=seurat_clusters), position = 'stack')
dev.off()

