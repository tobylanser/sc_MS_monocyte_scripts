dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC",
    "CXX14=g++", # or clang++ but you may need a version postfix
    file = M, sep = "\n", append = TRUE)

# devtools::install_github("jaleesr/BITFAM")
library(Seurat)
library(BITFAM)
# remotes::install_github("stan-dev/loo")

setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")

myeloid <- readRDS("Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds")
myeloid <- myeloid[, sample(colnames(myeloid), size =150000, replace=F)]
gc()
data_matrix_normalized <- BITFAM_preprocess(raw_data = myeloid@assays$RNA@counts)
rm(myeloid)
gc()
# data_matrix_normalized
BITFAM_res <- BITFAM(data = data_matrix_normalized, species = "human", scATAC_obj = NA, ncores = parallel::detectCores())

Z <- BITFAM_activities(BITFAM_res)

library(data.table)

fwrite(Z, 
       paste0('Z_BITFAM_out.csv'), quote=F, sep=",",
       row.names = T, col.names = T,  nThread = 28)



