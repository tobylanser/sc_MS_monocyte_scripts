library(zellkonverter)
library(dplyr)
library(Seurat)
library(scCustomize)
library(reticulate)
library(SeuratData)
library(SeuratDisk)

# devtools::install_github('satijalab/seurat-data')
# reticulate::install_python(version = '3.10')
use_condaenv("omicverse")
# BiocManager::install("zellkonverter")

setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")

merged <- readRDS("Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds")

table(merged$celltypes)
# saveRDS(merged, "objs/myeloid_BP2_unintegrated_annotated_v3.rds", compress = F)
# obj <- Convert_Assay(seurat_object = obj, convert_to = "V3")
as.anndata(x = merged, file_path = "h5ad/", file_name = "myeloid_BP2_unintegrated_annotated2.h5ad")
sce_obj <- as.SingleCellExperiment(merged, assay = c("RNA"))
rm(merged)
gc()
writeH5AD(sce_obj, "h5ad/myeloid_BP2_unintegrated_annotated2.h5ad", verbose = T)
?writeH5AD

as.anndata(x = merged, file_path = "h5ad/", file_name = "v3_allGenes_myeloid_unintegrated_annotated.h5ad")

SaveH5Seurat(merged, filename = "Rds/myeloid_BP2_unintegrated_annotated2.v3.h5Seurat")
gc()
Convert("Rds/myeloid_BP2_unintegrated_annotated2.v3.h5Seurat", dest = "h5ad")

