library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(sctransform)
library(openxlsx)
library(BiocParallel)
register(MulticoreParam(14))
options(future.globals.maxSize = 8e+09)
library(MAST)
library(harmony)
library(reticulate)
library(BPCells)
library(RColorBrewer)
library(SeuratWrappers)
# library(Azimuth)
library(Matrix)
library(car)
library(scater)
library(ggrepel)
# options(Seurat.object.assay.version = "v5")
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
library(DropletUtils)


setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNaseq_spatial_feng_et_al/spatial/")

prot.combined <- readRDS("./GSE284005_merfish_all.rds")

# DimPlot(prot.combined, reduction = "umap") # no reduction in this object


my_comparisons <- list(c("healthyWM","NAWM"),c("healthyWM","DMWM"),c("healthyWM","Rim"),
                       c("healthyWM","GM"),c("NAWM","DMWM"),c("NAWM","Rim"),
                       c("NAWM","GM"),c("DMWM","Rim"),c("DMWM","GM"),c("Rim","GM"))
colors <- brewer.pal(n=8,name = "Dark2")
colors <- colors[c(1,3,2,4,5)]

gene <- "IL1B"
p1 <- VlnPlot(prot.combined, features = gene,
              pt.size = 0.05, raster = F, group.by = "Region_banksy_major", cols = colors,
              # idents = c("16","13","19","2")
              # idents = c("pos")
              # idents = c("SDC2+LPCAT1+")
              # idents = c("16","13","19","2","26","21")
              # idents = c("Astrocytes")
              # idents = c("CD8+ T cells")
              # idents = c("mNK")
              # idents = c("Classical")
) + theme(legend.position = "none") + xlab("")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 20, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000000, max(p1[[1]][["data"]][[gene]])+0.8*max(p1[[1]][["data"]][[gene]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.01
p1













ImageDimPlot(prot.combined, fov = "s2r1", cols = "polychrome", axes = TRUE)

p1 <- ImageDimPlot(prot.combined, fov = "s2r1", cols = "red", cells = WhichCells(prot.combined, idents = 1))
p2 <- ImageDimPlot(prot.combined, fov = "s2r1", cols = "red", cells = WhichCells(prot.combined, idents = 15))
p1 + p2

p1 <- ImageFeaturePlot(prot.combined, features = "Slc17a7")
p2 <- ImageDimPlot(prot.combined, molecules = "Slc17a7", nmols = 10000, alpha = 0.3, mols.cols = "red")
p1 + p2


p1 <- ImageDimPlot(prot.combined, fov = "s2r1", alpha = 0.3, molecules = c("Slc17a7", "Olig1"), nmols = 10000)
markers.14 <- FindMarkers(prot.combined, ident.1 = "14")
p2 <- ImageDimPlot(prot.combined, fov = "s2r1", alpha = 0.3, molecules = rownames(markers.14)[1:4],
                   nmols = 10000)
p1 + p2


# create a Crop
cropped.coords <- Crop(prot.combined[["s2r1"]], x = c(1750, 3000), y = c(3750, 5250), coords = "plot")
# set a new field of view (fov)
prot.combined[["hippo"]] <- cropped.coords


# visualize FOV using default settings (no cell boundaries)
p1 <- ImageDimPlot(prot.combined, fov = "hippo", axes = TRUE, size = 0.7, border.color = "white", cols = "polychrome",
                   coord.fixed = FALSE)

# visualize FOV with full cell segmentations
DefaultBoundary(prot.combined[["hippo"]]) <- "segmentation"
p2 <- ImageDimPlot(prot.combined, fov = "hippo", axes = TRUE, border.color = "white", border.size = 0.1,
                   cols = "polychrome", coord.fixed = FALSE)

# simplify cell segmentations
prot.combined[["hippo"]][["simplified.segmentations"]] <- Simplify(coords = prot.combined[["hippo"]][["segmentation"]],
                                                                tol = 3)
DefaultBoundary(prot.combined[["hippo"]]) <- "simplified.segmentations"

# visualize FOV with simplified cell segmentations
DefaultBoundary(prot.combined[["hippo"]]) <- "simplified.segmentations"
p3 <- ImageDimPlot(prot.combined, fov = "hippo", axes = TRUE, border.color = "white", border.size = 0.1,
                   cols = "polychrome", coord.fixed = FALSE)

p1 + p2 + p3

# Since there is nothing behind the segmentations, alpha will slightly mute colors
ImageDimPlot(prot.combined, fov = "hippo", molecules = rownames(markers.14)[1:4], cols = "polychrome",
             mols.size = 1, alpha = 0.5, mols.cols = c("red", "blue", "yellow", "green"))






