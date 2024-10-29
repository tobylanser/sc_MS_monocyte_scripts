setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")

myeloid <- readRDS("Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds")
# saveRDS(myeloid, "Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds", compress = F)


myeloid <- myeloid[, sample(colnames(myeloid), size =100000, replace=F)]

meta_data <- data.frame(myeloid@meta.data)
gc()
library(data.table)
library(Seurat)
path <- ("compass/")
fwrite(meta_data, 
       paste0(path, 'compass_input/myeloidmeta_data.csv'), quote=F, sep=",",
       row.names = T, col.names = T,  nThread = 28)
fwrite(meta_data, 
       paste0(path, 'compass_input/myeloidmeta_data.txt'), quote=F, sep="\t",
       row.names = T, col.names = T,  nThread = 28)
gc()
fwrite(as.data.frame(GetAssayData(object = myeloid, layer = "data")), 
       paste0(path, 'compass_input/myeloid_counts.csv'), quote=F, sep=",",
       row.names = T, col.names = T,  nThread = 28)
gc()
fwrite(as.data.frame(GetAssayData(object = myeloid, layer = "data")), 
       paste0(path, 'compass_input/myeloid_counts.txt'), quote=F, sep="\t",
       row.names = T, col.names = T,  nThread = 28)


