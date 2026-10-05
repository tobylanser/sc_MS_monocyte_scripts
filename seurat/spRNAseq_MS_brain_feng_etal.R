# Seurat v5
library(openxlsx)
library(MAST)
library(Seurat)
library(harmony)
library(reticulate)
library(BPCells)
library(dplyr)
library(RColorBrewer)
library(SeuratWrappers)
library(Azimuth)
library(Matrix)
library(sctransform)
library(car)
library(scater)
library(ggplot2)
library(patchwork)
library(ggrepel)
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")
library(BiocParallel)
register(MulticoreParam(14))
library(paletteer)
library(nichenetr)
library(tidyverse)
library(circlize)
library(ggpubr)
library(VennDiagram)
library(pheatmap)
library(gplots)
library(GOplot)
library(readr)
library(robustbase)
library(EnhancedVolcano)
library(alluvial)
library(reshape2)
library(scico)
library(ggsci)
library(rcartocolor)
library(ggside)
library(viridis)
library(ggstatsplot)
library(grid)
library(shadowtext)
library(tidyr)



#########################         BPCells           ###################

setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNaseq_spatial_feng_et_al/spatial/")

prot.combined <- readRDS("./GSE284005_merfish_all.rds") # pre-done






