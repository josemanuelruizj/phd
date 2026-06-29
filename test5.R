ssh shiva
library(Seurat)
multi<-readRDS("integratedmultiome.rds")
multi[["RNA"]] <- JoinLayers(multi[["RNA"]])


Idents(multi) <- "run_name"   # o la columna de metadata que quieras usar
multi_5k_per_group <- subset(
  x = multi,
  downsample = 1000
)

outdir <- "datosWorkshop"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- multi_5k_per_group[[]]
expr_mat <- LayerData(multi_5k_per_group, assay = "RNA", layer = "data")
counts_mat <- LayerData(multi_5k_per_group, assay = "RNA", layer = "counts")
scale_mat <- LayerData(multi_5k_per_group, assay = "RNA", layer = "scale.data")

meta <- meta[colnames(expr_mat), , drop = FALSE]

saveRDS(meta, file.path(outdir, "multi_metadata.rds"))
saveRDS(expr_mat, file.path(outdir, "multi_expr_mat.rds"))
saveRDS(counts_mat, file.path(outdir, "multi_counts_mat.rds"))
saveRDS(scale_mat, file.path(outdir, "multi_scale_mat.rds"))

library(data.table)
library(Matrix)

fwrite(
  as.data.table(meta, keep.rownames = "cell"),
  file.path(outdir, "multi_metadata.tsv.gz"),
  sep = "\t"
)

expr_dt <- as.data.table(as.matrix(expr_mat), keep.rownames = "gene")
counts_dt <- as.data.table(as.matrix(counts_mat), keep.rownames = "gene")
scale_dt <- as.data.table(as.matrix(scale_mat), keep.rownames = "gene")

writeMM(expr_mat, file.path(outdir, "multi_expr_mat.mtx"))
writeLines(rownames(expr_mat), file.path(outdir, "genes.txt"))
writeLines(colnames(expr_mat), file.path(outdir, "cells.txt"))


writeMM(counts_mat, file.path(outdir, "multi_counts_mat.mtx"))
writeLines(rownames(counts_mat), file.path(outdir, "genes_count.txt"))
writeLines(colnames(counts_mat), file.path(outdir, "cells_count.txt"))


saveRDS(multi_5k_per_group, file.path(outdir, "multi_5k_per_group.rds")
)
