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


setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/spRNaseq_alsema_et_al/spatial/")

##################               Visium 
####      V1

# setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNAseq_fagiani_el.al_2025/spaceranger/GSE277435_RAW/GSM8522424/")
# m <- Read10X(data.dir = "/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNAseq_fagiani_el.al_2025/spaceranger/GSE277435_RAW/GSM8522424/filtered_feature_bc_matrix/")
# 
# # 2. Write out as 10x HDF5
# DropletUtils::write10xCounts(
#   # file       = "filtered_feature_bc_matrix.h5",
#   path = "/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNAseq_fagiani_el.al_2025/spaceranger/GSE277435_RAW/GSM8522424/filtered_feature_bc_matrix.h5",
#   x          = m,
#   type       = "HDF5",
#   version    = "3",          # 3 for feature-barcode format
#   overwrite  = TRUE,
#   gene.id    = rownames(m),
#   gene.symbol = rownames(m)  # or your gene symbols if different
# )
# 
# V1 <- Load10X_Spatial(data.dir = "/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNAseq_fagiani_el.al_2025/spaceranger/GSE277435_RAW/GSM8522424/")
# gc()


######################   merging all slides

file.dir <- "../spaceranger/"
data.list <- c()
meta1 <- read.delim("meta.txt", header = T)

files.set <- c(
  "A1", #1
  "A2",#2
  "A3",#3
  "A4",#4
  "C1",#5
  "C2",#6
  "CG1",#7
  "CG2",#8
  "CG3",#9
  "M1",#10
  "M2",#11
  "M3",#12
  "M4",#13
  "M5",#14
  "M6",#15
  "N1",#16
  "N2",#17
  "N3",#18
  "NG1",#19
  "NG2",#20
  "NG3",#21
  "SG1",#22
  "SG2",#23
  "SG3",#24
  "SG4",#25
  "SG5")#26

for (i in 1:length(files.set)) {
  # path1 <- paste0(file.dir, files.set[i],"/filtered_feature_bc_matrix/")
  # m <- Read10X(data.dir = path1)
  # path2 <- paste0(file.dir, files.set[i],"/filtered_feature_bc_matrix.h5")
  # DropletUtils::write10xCounts(
  #   path = path2,
  #   x          = m,
  #   type       = "HDF5",
  #   version    = "3",          # 3 for feature-barcode format
  #   overwrite  = TRUE,
  #   gene.id    = rownames(m),
  #   gene.symbol = rownames(m)  # or your gene symbols if different
  # )
  dataset_name <- files.set[i]
  path3 <- paste0(file.dir, files.set[i])
  mat <- Load10X_Spatial(data.dir = path3) %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  mat$sample <- dataset_name
  data.list[[i]] <- mat
  rm(mat)
  rm(m)
  gc()
}
# Name layers
names(data.list) <- files.set

# Merge layers and create seurat obj during merging 
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
prot.combined <- merge(data.list[[1]], y = data.list[2:length(data.list)], 
                       add.cell.ids = files.set, merge.data = T
)

rm(data.list)
gc()

VariableFeatures(prot.combined) <- features

### Normalize and scale merged obj
# prot.combined <- NormalizeData(prot.combined, normalization.method = "LogNormalize")
prot.combined <- FindVariableFeatures(prot.combined, selection.method = "vst", nfeatures = 3000)
prot.combined <- ScaleData(prot.combined #, vars.to.regress = c("percent.mt","nCount_RNA"), model.use = "linear" #latent.data = "nFeature_RNA", 
)#, features = all.genes)
gc()
### Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(prot.combined), 10)

### plot variable features with and without labels
#plot1 <- VariableFeaturePlot(prot.combined, raster = F)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, raster = F)
#plot1 + plot2

### Dimensionality reduction and integration
prot.combined <- RunPCA(prot.combined, npcs = 50)
gc()
ElbowPlot(prot.combined, ndims = 50)
prot.combined <- FindNeighbors(prot.combined, dims = 1:30, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = 2, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
prot.combined <- RunUMAP(prot.combined, #umap.method = "umap-learn", 
                         dims = 1:30, reduction = "pca", reduction.name = "umap")
gc()

DimPlot(prot.combined, reduction = "umap", label = T, raster = F) # + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

FeaturePlot(prot.combined, features = "SDC2", raster = F, split.by = "pathology")

colors <- brewer.pal(n=8,name = "Dark2")
colors <- colors[c(1,2)]
SpatialDimPlot(prot.combined, label = F, repel = T, group.by = "SDC2_status",
               # cols = colors, label.size = 4,
               
               images = c("slice1.5","slice1.6", ## CWM
                          "slice1.16","slice1.17","slice1.18", ## NAWM
                          "slice1","slice1.2","slice1.3","slice1.4", ## act WM lesion
                          "slice1.10","slice1.11","slice1.12","slice1.14","slice1.15"  ## act/inact "slice1.13",
                          ),
               ncol = 5) + theme(legend.position = "none")

Idents(prot.combined) <- "SDC2_status"
cells <- CellsByIdentities(prot.combined, idents = "SDC2+LPCAT1+")

p1 <- SpatialDimPlot(prot.combined, #image.scale = "hires",
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                    combine = T, images = c("slice1.5") ## CWM
) + NoLegend()
p2 <- SpatialDimPlot(prot.combined, #image.scale = "hires",
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.6") ## CWM
) + NoLegend()

p3 <- SpatialDimPlot(prot.combined,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                    combine = T, images = c("slice1.16") ## NAWM
) + NoLegend()
p4 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("grey50","#FFFF00"), facet.highlight = T, 
                     combine = T, images = c("slice1.17") ## NAWM
) + NoLegend()
p5 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("grey50","#FFFF00"), facet.highlight = T, 
                     combine = T, images = c("slice1.18") ## NAWM
) + NoLegend()

p6 <- SpatialDimPlot(prot.combined,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                    combine = T, images = c("slice1")  ## act WM lesion
) + NoLegend()
p7 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.2")  ## act WM lesion
) + NoLegend()
p8 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.3")  ## act WM lesion
) + NoLegend()
p9 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.4")  ## act WM lesion
) + NoLegend()

p10 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.10") ## act/inact
) + NoLegend()
p11 <- SpatialDimPlot(prot.combined,
                      cells.highlight = cells[setdiff(names(cells), "NA")],
                      cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                      combine = T, images = c("slice1.11") ## act/inact
) + NoLegend()
p12 <- SpatialDimPlot(prot.combined,
                      cells.highlight = cells[setdiff(names(cells), "NA")],
                      cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                      combine = T, images = c("slice1.12") ## act/inact
) + NoLegend()
# p13 <- SpatialDimPlot(prot.combined,
#                       cells.highlight = cells[setdiff(names(cells), "NA")],
#                       cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
#                       combine = T, images = c("slice1.13") ## act/inact
# ) + NoLegend()
p14 <- SpatialDimPlot(prot.combined,
                      cells.highlight = cells[setdiff(names(cells), "NA")],
                      cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                      combine = T, images = c("slice1.14") ## act/inact
) + NoLegend()
p15 <- SpatialDimPlot(prot.combined,
                      cells.highlight = cells[setdiff(names(cells), "NA")],
                      cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                      combine = T, images = c("slice1.15") ## act/inact
) + NoLegend()

(p1 + p2 + 
    p3 + p4 + p5 + 
    p6 + p7 + p8 + p9 + 
    p10 + p11 + p12 + p14 + p15) + 
  plot_layout(nrow = 3, ncol = 5, guides = 'collect')





p1 <- SpatialDimPlot(prot.combined, #image.scale = "hires",
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.7") ## CGM
) + NoLegend()
p2 <- SpatialDimPlot(prot.combined, #image.scale = "hires",
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.8") ## CGM
) + NoLegend()
p3 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.9") ## CGM
) + NoLegend()

p4 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.19") ## NAGM
) + NoLegend()
p5 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.20") ## NAGM
) + NoLegend()

p6 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.21")  ## NAGM
) + NoLegend()
p7 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.22")  ## subpial
) + NoLegend()
p8 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.23")  ## subpial
) + NoLegend()
p9 <- SpatialDimPlot(prot.combined,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                     combine = T, images = c("slice1.24")  ## subpial
) + NoLegend()

p10 <- SpatialDimPlot(prot.combined,
                      cells.highlight = cells[setdiff(names(cells), "NA")],
                      cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                      combine = T, images = c("slice1.25") ## subpial
) + NoLegend()
p11 <- SpatialDimPlot(prot.combined,
                      cells.highlight = cells[setdiff(names(cells), "NA")],
                      cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                      combine = T, images = c("slice1.26") ## subpial
) + NoLegend()


(p1 + p2 + p3 + p5 + 
    p6 + p7 + p9 + p10) + 
  plot_layout(nrow = 2, ncol = 4, guides = 'collect')


SpatialFeaturePlot(prot.combined, features = c("LPCAT1"), max.cutoff = 3,
                   slot = "data", images = c("slice1.5","slice1.6", ## CWM
                                             "slice1.16","slice1.17","slice1.18", ## NAWM
                                             "slice1","slice1.2","slice1.3","slice1.4", ## act WM lesion
                                             "slice1.10","slice1.11","slice1.12","slice1.14","slice1.15"  ## act/inact "slice1.13",
                                             ), ncol = 5
)

SpatialFeaturePlot(prot.combined, features = c("LPCAT1"), max.cutoff = 3,
                   slot = "data", images = c("slice1.7","slice1.8","slice1.9", ## CGM
                                             "slice1.20","slice1.21", ## NAGM "slice1.19",
                                              "slice1.22","slice1.23","slice1.24",  ## subpial
                                              "slice1.25","slice1.26" ## subpial
                                             ), ncol = 5
)

patho <- c("CWM","NAWM","active","act/inact","CGM","NAGM","subpial")
my_comparisons <- list(c("CGM","NAGM"),c("CGM","subpial"),c("NAGM","subpial")#,
                       # c("CWM","NAWM"),c("CWM","active"),c("CWM","act/inact"),
                       # c("NAWM","active"),c("NAWM","act/inact"),c("active","act/inact")
                       )
colors <- brewer.pal(n=8,name = "Dark2")
colors <- colors[c(1,3,2)]

Idents(prot.combined) <- "pathology"

gene <- "LPCAT1"
p1 <- VlnPlot(prot.combined, features = gene,
              pt.size = 0.05, raster = F, group.by = "pathology", cols = colors,
              idents = c("CGM","NAGM","subpial")
              # idents = c("CWM","NAWM","active","act/inact")
) + theme(legend.position = "none") + xlab("")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 20, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000005, max(p1[[1]][["data"]][[gene]])+.3*max(p1[[1]][["data"]][[gene]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1

clus4.markers <- c("SDC2","PTGS2","HBEGF","IL1B","FOSL2","PPIF","GPR183","RGCC", "TRIB1", "NFKB1","EGR1","NLRP3","EGR2","EGR3","CXCL8",  
                   "G0S2")

DotPlot(prot.combined, features = clus4.markers,
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        idents = c("CWM","NAWM","active","act/inact"
                   # "CGM","NAGM","subpial"
                   ),
        cluster.idents = F, 
        group.by = "pathology",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip() + ylab("") + xlab("")



## Add metadata
meta <- read.delim("meta.txt", header = T)

# disease <- c("CTRL","MS","MS","MS")
# prog <- c("CTRL","RRMS","PMS","PMS")
# diagnosis <- c("CTR","RRMS","PPMS","SPMS")
prot.combined$pathology <- prot.combined$sample
prot.combined$disease <- prot.combined$sample
prot.combined$age <- prot.combined$sample
prot.combined$sex <- prot.combined$sample
prot.combined$donor <- prot.combined$sample

for (i in 1:length(meta$GEM)) {
  prot.combined$pathology <- recode(prot.combined$pathology, "meta$GEM[i] = meta$Tissue_type[i]")
  prot.combined$age <- recode(prot.combined$age, "meta$GEM[i] = meta$Age[i]")
  prot.combined$sex <- recode(prot.combined$sex, "meta$GEM[i] = meta$Sex[i]")
  prot.combined$disease <- recode(prot.combined$disease, "meta$GEM[i] = meta$disease[i]")
  prot.combined$donor <- recode(prot.combined$donor, "meta$GEM[i] = meta$Donor..[i]")
}

patho <- c("CWM","NAWM","active","act/inact","CGM","NAGM","subpial")
prot.combined$pathology <- factor(prot.combined$pathology, levels = patho)
samples <- c(
  "C1",
  "C2",
  "N1",
  "N2",
  "N3",
  "A1",
  "A2",
  "A3",
  "A4",
  "M1",
  "M2",
  "M3",
  "M4",
  "M5",
  "M6",
  "CG1",
  "CG2",
  "CG3",
  "NG1",
  "NG2",
  "NG3",
  "SG1",
  "SG2",
  "SG3",
  "SG4",
  "SG5")
prot.combined$sample <- factor(prot.combined$sample, levels = samples)

## Save and Load V5 data
saveRDS(object = prot.combined, file = "obj_unintegrated_merged.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_merged.Rds")



## cell types annotation

DotPlot(prot.combined, features = c( "GFAP","SLC1A3","AQP4","LCN2", "GJA1", "SLC1A2","FGFR3","NKAIN4",   #Astrocytes
                                     "SDC2",
                                     "FBLN1","FBLN5", # fibroblasts
                                     "CHRM3", # cholinergic neurons
                                     "TH","SLC18A2", # dopaminergic neurons
                                     "TAGLN","MYH11", # vascular smooth muscle cells
                                     "CFAP44","CFAP43", # ependymal cells
                                     "SLC17A6","SLC17A7","NRGN","CAMK2A", "SATB2", "COL5A1","SDK2","NEFM","HTR2C",    #Excitatory_neurons
                                     "SLC32A1","GAD1","GAD2","TAC1","PENK","SST","NPY","MYBPC1","PVALB","GABBR2",   #Inhibitory_neurons
                                     "OLIG2", "MBP","MOBP","PLP1","MOG","CLDN11","MYRF","GALC","ERMN","MAG",   #Oligodendrocytes
                                     "VCAN","CSPG4","PDGFRA", "SOX10","NEU4", "PCDH15","GPR37L1","C1QL1","CDO1","EPN2",   #Oligodendrocyte_precursor_cells
                                     "AMBP","HIGD1B","COX4I2", "AOC3","PDE5A","PTH1R","P2RY14","ABCC9","KCNJ8","CD248",  #Pericytes
                                     "FLT1","CLDN5", "VTN","ITM2A", "VWF", "FAM167B","BMX","CLEC1B",    #Endothelial_cells
                                     "P2RY12","CSF1R","C3","APOE","CD74","CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5","AIF1","IL1B","IRF8",   #Microglia
                                     "MS4A4A","CD163","MAFB","TNFAIP2","IL15","ASAH1","PLA2G7", "PLXND1","EMILIN2","SIGLEC1","F13A1","MARCO","GAS7","GDA", # monocytes
                                     "CD68", #"CD14","FCGR3A","FCGR1A","TFRC","CCR5","ITGAM","CCR2","HP","SELL", #Macrophages
                                     "CD8A","CD4","CD3E","CD3D","CD19","CD22","IGKC","PTPRC","CD27"
),
col.max = 20, 
dot.scale = 10, 
cluster.idents = T, 
# group.by = "type_broad",
#scale = F,
#split.by = "cohort"
) + RotatedAxis() + ylab("") + xlab("")










prot.combined[["Spatial"]] <- JoinLayers(prot.combined[["Spatial"]])

SDC2.pos <- WhichCells(prot.combined, expression = SDC2 > 0 &  AQP4 == 0, slot = "counts")
prot.combined$SDC2.pos <- colnames(prot.combined) %in% SDC2.pos
table(prot.combined$SDC2.pos)


SDC2_expr <- GetAssayData(prot.combined, layer = "counts")["SDC2", ]
GFAP_expr <- GetAssayData(prot.combined, layer = "counts")["GFAP", ]
prot.combined <- AddMetaData(prot.combined,
                             metadata = ifelse(SDC2_expr > 0 & GFAP_expr == 0, "pos", "neg"),
                             col.name = "SDC2_status")


SDC2.pos <- WhichCells(prot.combined, expression = SDC2 > 0 &  LPCAT1 > 0, slot = "counts")
prot.combined$SDC2.pos <- colnames(prot.combined) %in% SDC2.pos
table(prot.combined$SDC2.pos)


SDC2_expr <- GetAssayData(prot.combined, layer = "counts")["SDC2", ]
LPCAT1_expr <- GetAssayData(prot.combined, layer = "counts")["LPCAT1", ]
AQP4_expr <- GetAssayData(prot.combined, layer = "counts")["AQP4", ]
prot.combined <- AddMetaData(prot.combined,
                             metadata = ifelse(SDC2_expr > 0 & LPCAT1_expr > 0, "SDC2+LPCAT1+", "SDC2-LPCAT1-"),
                             col.name = "SDC2_status")
table(prot.combined$SDC2_status)
# 
# prot.combined2 <- subset(prot.combined, subset = SDC2 > 0.1)

Idents(prot.combined) <- "SDC2_status"

zk.response0 <- FindMarkers(prot.combined, ident.1 = "pos",
                            ident.2 = "neg",
                            slot = "data",
                            assay = "Spatial",
                            features = NULL,
                            logfc.threshold = 0,
                            test.use = "wilcox",
                            min.pct = 0.0,
                            min.diff.pct = -Inf,
                            verbose = TRUE,
                            only.pos = FALSE,
                            max.cells.per.ident = Inf,
                            random.seed = 1,
                            latent.vars = NULL,
                            min.cells.feature = 3,
                            min.cells.group = 3,
                            pseudocount.use = 1,
                            mean.fxn = NULL,
                            fc.name = NULL,
                            base = 2,
                            densify = FALSE,
                            recorrect_umi = TRUE
)
# zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_SDC2_pos_x_neg_DEGs.xlsx")

















#########################################################################
# DefaultAssay(V1) <- "Spatial.008um"
vln.plot <- VlnPlot(V1, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 10)) + NoLegend()
count.plot <- SpatialFeaturePlot(V1, features = "nCount_Spatial") + theme(legend.position = "right")
# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot

# DefaultAssay(V1) <- "Spatial.008um"
V1 <- NormalizeData(V1)
gc()

# switch spatial resolution to 2um from 8um
DefaultAssay(CT1_A1) <- "Spatial.002um"
p1 <- SpatialFeaturePlot(CT1_A1, features = "Apoe", min.cutoff = 5, max.cutoff = 8) + ggtitle("Apoe expression (2um)")
p1
# switch back to 8um
DefaultAssay(CT1_A1) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(CT1_A1, features = "Apoe") + ggtitle("Apoe expression (8um)")
p1 | p2

## Unsupervised clustering
DefaultAssay(CT1_A1) <- "Spatial.008um"
CT1_A1 <- FindVariableFeatures(CT1_A1)
CT1_A1 <- ScaleData(CT1_A1)
gc()
# we select 50,0000 cells and create a new 'sketch' assay
CT1_A1 <- SketchData(
  object = CT1_A1,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
gc()

# switch analysis to sketched cells
DefaultAssay(CT1_A1) <- "sketch"

# perform clustering workflow
CT1_A1 <- FindVariableFeatures(CT1_A1)
CT1_A1 <- ScaleData(CT1_A1)
CT1_A1 <- RunPCA(CT1_A1, assay = "sketch", reduction.name = "pca.sketch")
ElbowPlot(CT1_A1, ndims = 50, reduction = "pca.sketch")
gc()
CT1_A1 <- FindNeighbors(CT1_A1, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
CT1_A1 <- FindClusters(CT1_A1, cluster.name = "seurat_cluster.sketched", resolution = 3)
CT1_A1 <- RunUMAP(CT1_A1, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
gc()

CT1_A1 <- ProjectData(
  object = CT1_A1,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
gc()

DefaultAssay(CT1_A1) <- "sketch"
Idents(CT1_A1) <- "seurat_cluster.sketched"
p1 <- DimPlot(CT1_A1, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(CT1_A1) <- "Spatial.008um"
Idents(CT1_A1) <- "seurat_cluster.projected"
p2 <- DimPlot(CT1_A1, reduction = "full.umap.sketch", label = T, raster = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "right")
p2
FeaturePlot(CT1_A1, features = "Apoe", raster = F)
p1 | p2

SpatialDimPlot(CT1_A1, label = T, repel = T, label.size = 4)

# Save and Load data
saveRDS(object = CT1_A1, file = "obj_CT1_A1.Rds")
#rm(prot.combined)
CT1_A1 <- readRDS("./obj_CT1_A1.Rds")

DotPlot(CT1_A1, features = c( "Apoe",
                              "Gfap", "Aqp4", "Lcn2", "Gja1", "Slc1a2", "Fgfr3", "Nkain4",   #astrocytes
                              "Flt1", "Cldn5", "Vtn", "Itm2a", "Vwf", "Fam167b", "Bmx", "Clec1b",    #endothelial_cells
                              "Slc17a6",  "Slc17a7",  "Nrgn", "Camk2a", "Satb2", "Col5a1", "Sdk2", "Nefm",    #excitatory_neurons
                              "Slc32a1",  "Gad1", "Gad2", "Tac1", "Penk", "Sst",  "Npy",  "Mybpc1", "Pvalb", "Gabbr2",   #inhibitory_neurons
                              "P2ry12", "Csf1r",  "Cd74", "C3", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia
                              "Olig2",  "Mbp",  "Mobp", "Plp1", "Mog",  "Cldn11", "Myrf", "Galc", "Ermn", "Mag",   #oligodendrocytes
                              "Vcan", "Cspg4", "Pdgfra", "Sox10", "Neu4", "Pcdg15", "Gpr37l1", "C1ql1", "Cdo1", "Epn2",   #oligodendrocyte_precursor_cells
                              "Ambp",  "Higd1b", "Cox4i2", "Aoc3", "Pde5a",  "Pth1r",  "P2ry14", "Abcc9", "Kcnj8", "Cd248", #Pericytes
                              "Ccr2","Cd68","Cd11b","Cd14","Fcgr3", #monocytes
                              "Ccr5","Itgam","Trfc","Fcgr1" #macrophage
                              
),
col.max = 20,
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()

Idents(CT1_A1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(CT1_A1, idents = c(0:8))
p <- SpatialDimPlot(CT1_A1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

### astrocytes
Idents(CT1_A1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(CT1_A1, idents = c(34,7,18,32))
p <- SpatialDimPlot(CT1_A1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p
### myeloids
Idents(CT1_A1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(CT1_A1, idents = c(17,29,49))
p <- SpatialDimPlot(CT1_A1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p
# ggsave(filename = "CT1_A1_spatial_feat_clus0-8.png",
#        plot = p,
#        width = 10,
#        height = 10,
#        dpi = 600,
#        device = "png")

## find and visualize the top gene expression markers for each cluster
# Create downsampled object to make visualization either
DefaultAssay(CT1_A1) <- "Spatial.008um"
Idents(CT1_A1) <- "seurat_cluster.projected"
object_subset <- subset(CT1_A1, cells = Cells(CT1_A1[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE) %>%
  group_by(cluster) #%>%
# dplyr::filter(avg_log2FC > 1)%>%
# slice_head(n = 50)
write.xlsx(as.data.frame(markers), rowNames = T, file="wilcox_clus_all_markers_IPSI_50pcs_res3.xlsx")
markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) #+ NoLegend()
p

ggsave(filename = "CT1_A1_heat_top10.png",
       plot = p,
       width = 20,
       height = 20,
       dpi = 600,
       device = "png")





### IPSI separately

IPSI_D1 <- Load10X_Spatial(data.dir = "/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/izzy/visium_cancer_joe/spaceranger/IPSI_D1/outs/", bin.size = c(2, 8, 16))
Assays(IPSI_D1)
gc()

DefaultAssay(IPSI_D1) <- "Spatial.008um"
vln.plot <- VlnPlot(IPSI_D1, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 10)) + NoLegend()
count.plot <- SpatialFeaturePlot(IPSI_D1, features = "nCount_Spatial.008um") + theme(legend.position = "right")
# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot

DefaultAssay(IPSI_D1) <- "Spatial.002um"
IPSI_D1 <- NormalizeData(IPSI_D1)
DefaultAssay(IPSI_D1) <- "Spatial.008um"
IPSI_D1 <- NormalizeData(IPSI_D1)
DefaultAssay(IPSI_D1) <- "Spatial.016um"
IPSI_D1 <- NormalizeData(IPSI_D1)
gc()

# switch spatial resolution to 2um from 8um
DefaultAssay(IPSI_D1) <- "Spatial.002um"
p1 <- SpatialFeaturePlot(IPSI_D1, features = "Apoe", slot = "counts") + ggtitle("Apoe expression (2um)")
p1
# switch back to 8um
DefaultAssay(IPSI_D1) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(IPSI_D1, features = "Apoe", slot = "counts") + ggtitle("Apoe expression (8um)")
p2
p1 | p2

## Unsupervised clustering
DefaultAssay(IPSI_D1) <- "Spatial.008um"
IPSI_D1 <- FindVariableFeatures(IPSI_D1)
gc()
IPSI_D1 <- ScaleData(IPSI_D1)
gc()
# we select 50,0000 cells and create a new 'sketch' assay
IPSI_D1 <- SketchData(
  object = IPSI_D1,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
gc()

# switch analysis to sketched cells
DefaultAssay(IPSI_D1) <- "sketch"

# perform clustering workflow
IPSI_D1 <- FindVariableFeatures(IPSI_D1)
IPSI_D1 <- ScaleData(IPSI_D1)
IPSI_D1 <- RunPCA(IPSI_D1, assay = "sketch", reduction.name = "pca.sketch")
ElbowPlot(IPSI_D1, ndims = 50, reduction = "pca.sketch")
gc()
IPSI_D1 <- FindNeighbors(IPSI_D1, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
IPSI_D1 <- FindClusters(IPSI_D1, cluster.name = "seurat_cluster.sketched", resolution = 3)
IPSI_D1 <- RunUMAP(IPSI_D1, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
gc()

IPSI_D1 <- ProjectData(
  object = IPSI_D1,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
gc()

DefaultAssay(IPSI_D1) <- "sketch"
Idents(IPSI_D1) <- "seurat_cluster.sketched"
p1 <- DimPlot(IPSI_D1, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(IPSI_D1) <- "Spatial.008um"
Idents(IPSI_D1) <- "seurat_cluster.projected"
p2 <- DimPlot(IPSI_D1, reduction = "full.umap.sketch", label = T, raster = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "right")
p2
FeaturePlot(IPSI_D1, features = "Apoe", raster = F)
p1 | p2

SpatialDimPlot(IPSI_D1, label = T, repel = T, label.size = 4)

DefaultAssay(IPSI_D1) <- "Spatial.008um"
SpatialFeaturePlot(IPSI_D1, features = c("Apoe")
                   ,slot = "data"
)

# Save and Load data
saveRDS(object = IPSI_D1, file = "obj_IPSI_D1.Rds")
#rm(prot.combined)
IPSI_D1 <- readRDS("./obj_IPSI_D1.Rds")



Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(0:8))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p
# ggsave(filename = "IPSI_D1_spatial_feat_clus0-8.png",
#        plot = p,
#        width = 10,
#        height = 10,
#        dpi = 600,
#        device = "png")

Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(9:17))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p


Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(18:26))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p


Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(27:35))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p


Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(36:44))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p


Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(45:52))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(12,25,13,37))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

cells <- CellsByIdentities(IPSI_D1, idents = c(7,22,16,13,20,15,23,11,37))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

cells <- CellsByIdentities(IPSI_D1, idents = c(12,25,13,37))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p


### astrocytes
Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(7,22,42,16))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p
### myeloids
Idents(IPSI_D1) <- "seurat_cluster.projected"
cells <- CellsByIdentities(IPSI_D1, idents = c(12,25,13,36))
p <- SpatialDimPlot(IPSI_D1,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

## find and visualize the top gene expression markers for each cluster
# Create downsampled object to make visualization either
DefaultAssay(IPSI_D1) <- "Spatial.008um"
Idents(IPSI_D1) <- "seurat_cluster.projected"
object_subset <- subset(IPSI_D1, cells = Cells(IPSI_D1[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE) %>%
  group_by(cluster) #%>%
# dplyr::filter(avg_log2FC > 1)%>%
# slice_head(n = 50)
write.xlsx(as.data.frame(markers), rowNames = T, file="wilcox_clus_all_markers_IPSI_50pcs_res3.xlsx")
markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) #+ NoLegend()
ggsave(filename = "IPSI_D1_heat_top10_HD.png",
       plot = p,
       width = 30,
       height = 30,
       #dpi = 600,
       device = "png")


# Astrocytes <- c("Gfap", "EAAT1", "AQP4", "LCN2", "GJA1", "SLC1A2", "FGFR3", "NKAIN4")
# Excitatory_neurons <- c("SLC17A6",  "SLC17A7",  "NRGN", "CAMK2A", "SATB2", "COL5A1", "SDK2", "NEFM")
# Inhibitory_neurons <- c("SLC32A1",  "GAD1", "GAD2", "TAC1", "PENK", "SST",  "NPY",  "MYBPC1", "PVALB", "GABBR2")
# Microglia <- c("IBA-1", "P2RY12", "CSF1R",  "CD74", "C3", "CST3", "HEXB", "C1QA", "CX3CR1", "AIF-1")
# Oligodendrocytes <- c("OLIG2",  "MBP",  "MOBP", "PLP1", "MOG",  "CLDN11", "MYRF", "GALC", "ERMN", "MAG")
# Oligodendrocyte_precursor_cells <- c("VCAN", "CSPG4", "PDGFRA", "SOX10", "NEU4", "PCDG15", "GPR37L1", "C1QL1", "CDO1", "EPN2")
# Macrophages <- c("CD14", "CD16", "CD64","CD71", "CCR5","Cc68", "Ccr2", "Cd11b")
#Endothelial_cells <- c("FLT1", "CLDN5", "VTN", "ITM2A", "VWF", "FAM167B", "BMX", "CLEC1B")
#Pericytes <- c("AMBP",  "HIGD1B", "COX4I2", "AOC3", "PDE5A",  "PTH1R",  "P2RY14", "ABCC9", "KCNJ8", "CD248")

DotPlot(IPSI_D1, features = c( "Apoe",
                               "Gfap", "Aqp4", "Lcn2", "Gja1", "Slc1a2", "Fgfr3", "Nkain4",   #astrocytes
                               "Flt1", "Cldn5", "Vtn", "Itm2a", "Vwf", "Fam167b", "Bmx", "Clec1b",    #endothelial_cells
                               "Slc17a6",  "Slc17a7",  "Nrgn", "Camk2a", "Satb2", "Col5a1", "Sdk2", "Nefm",    #excitatory_neurons
                               "Slc32a1",  "Gad1", "Gad2", "Tac1", "Penk", "Sst",  "Npy",  "Mybpc1", "Pvalb", "Gabbr2",   #inhibitory_neurons
                               "P2ry12", "Csf1r",  "Cd74", "C3", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia
                               "Olig2",  "Mbp",  "Mobp", "Plp1", "Mog",  "Cldn11", "Myrf", "Galc", "Ermn", "Mag",   #oligodendrocytes
                               "Vcan", "Cspg4", "Pdgfra", "Sox10", "Neu4", "Pcdg15", "Gpr37l1", "C1ql1", "Cdo1", "Epn2",   #oligodendrocyte_precursor_cells
                               "Ambp",  "Higd1b", "Cox4i2", "Aoc3", "Pde5a",  "Pth1r",  "P2ry14", "Abcc9", "Kcnj8", "Cd248", #Pericytes
                               "Ccr2","Cd68","Cd11b","Cd14","Fcgr3", #monocytes
                               "Ccr5","Itgam","Trfc","Fcgr1" #macrophage
                               
),
#cols = c("blue","blue","blue"),#"green","yellow","gray","pink","brown","lightblue"), 
col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),#"Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"),
#c("Classical Mono_1","Classical Mono_2","Intermediate Mono","Nonclassical Mono", 
#"pDC_AD","pDC_C", 
#"mo-DC_AD","mo-DC_C"
#),
#idents = "CD8+ TEM",
# c("NK_4_C","NK_4_AD","NK_8_C","NK_8_AD","NK_21_C","NK_21_AD",
#  "CD8+ NKT-like_C", "CD8+ NKT-like_AD"#, "NK_AD", "NK_C"
# c("NK_AD","NK_C","Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C","Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"
#c("62","67","49","47"),
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()


p1 <- VlnPlot(IPSI_D1, features = c("Apoe"),#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
              #split.by = "disease",
              pt.size = 0.05,
              raster = F,
              #ncol = 1,
              #group.by = "disease",
              #slot = "counts",
              #add.noise = F,
              #log = T,
              #sort = "increasing",
              #idents = c("0","1","6","13","14","18") #endothelial
              #c("24","26","30")#microglia
) + scale_y_continuous(limits = c(0.000,8.5)) + #geom_boxplot(width=0.1, color="black", alpha=0.2) +
  stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = 95)
p1$layers[[2]]$aes_params$alpha <- 0.05
p1

###############################
## astrocytes

DotPlot(IPSI_D1, features = c( "Apoe",
                               "Gfap", "Aqp4", "Lcn2", "Gja1", "Slc1a2", "Fgfr3", "Nkain4"   #astrocytes
),
col.max = 20, 
idents = c("7","22","42","16"),
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
scale = F,
#split.by = "cohort"
) + RotatedAxis()


## myeloids

DotPlot(IPSI_D1, features = c( "Apoe",
                               "P2ry12", "Csf1r",  "Cd74", "C3", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia
                               "Ccr2","Cd68","Cd11b","Cd14","Fcgr3", #monocytes
                               "Ccr5","Itgam","Trfc","Fcgr1" #macrophage                    
),
col.max = 20, 
idents = c("12","25","13","36"),
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
scale = F,
#split.by = "cohort"
) + RotatedAxis()

# DoHeatmap(IPSI_D1, features = c( "Apoe",
#                                   "P2ry12","Csf1r","Cd74","C3","Cst3","Hexb","C1qa","Aif1", "Tmem119",  #microglia
#                                   "Ccr2","Cd68","Cd14","Fcgr3", #monocytes
#                                   "Ccr5","Itgam","Fcgr1" #macrophage                    
# ), 
#           slot = "data",
#           assay = "Spatial.008um",
#           cells = 1:1000, 
#           size = 2.5,
#           #disp.max = 2.5,
#           #disp.min = 0
# )

#### CT1
# astrocytes

DotPlot(CT1_A1, features = c( "Apoe",
                              "Gfap", "Aqp4", "Lcn2", "Gja1", "Slc1a2", "Fgfr3", "Nkain4"  #astrocytes
),
col.max = 20, 
idents = c("34","7","18","32"),
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
scale = F,
#split.by = "cohort"
) + RotatedAxis()


## myeloids

DotPlot(CT1_A1, features = c( "Apoe",
                              "P2ry12", "Csf1r",  "Cd74", "C3", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia
                              "Ccr2","Cd68","Cd11b","Cd14","Fcgr3", #monocytes
                              "Ccr5","Itgam","Trfc","Fcgr1" #macrophage                    
),
col.max = 20, 
idents = c("17","29","49"),
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
scale = F,
#split.by = "cohort"
) + RotatedAxis()
















############# ##########   merging CT1 and IPSI    #####################################

##################               Visium HD    ########################

## Controlateral
CT1_A1 <- Load10X_Spatial(data.dir = "/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/izzy/visium_cancer_joe/spaceranger/CT1_A1/outs/", bin.size = c(8))
Assays(CT1_A1)
gc()
# DefaultAssay(CT1_A1) <- "Spatial.002um"
# CT1_A1 <- NormalizeData(CT1_A1)
DefaultAssay(CT1_A1) <- "Spatial.008um"
CT1_A1 <- NormalizeData(CT1_A1)
# DefaultAssay(CT1_A1) <- "Spatial.016um"
# CT1_A1 <- NormalizeData(CT1_A1)
gc()
# DefaultAssay(CT1_A1) <- "Spatial.008um"
CT1_A1 <- FindVariableFeatures(CT1_A1)
CT1_A1 <- ScaleData(CT1_A1)
gc()
CT1_A1$orig.ident <- "Contralateral"
CT1_A1$ID <- "Contralateral"

## Ipsilateral
IPSI_D1 <- Load10X_Spatial(data.dir = "/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/izzy/visium_cancer_joe/spaceranger/IPSI_D1/outs/", bin.size = c(8))
Assays(IPSI_D1)
gc()
# DefaultAssay(IPSI_D1) <- "Spatial.002um"
# IPSI_D1 <- NormalizeData(IPSI_D1)
DefaultAssay(IPSI_D1) <- "Spatial.008um"
IPSI_D1 <- NormalizeData(IPSI_D1)
# DefaultAssay(IPSI_D1) <- "Spatial.016um"
# IPSI_D1 <- NormalizeData(IPSI_D1)
gc()
# DefaultAssay(IPSI_D1) <- "Spatial.008um"
IPSI_D1 <- FindVariableFeatures(IPSI_D1)
IPSI_D1 <- ScaleData(IPSI_D1)
gc()
IPSI_D1$orig.ident <- "Ipsilateral"
IPSI_D1$ID <- "Ipsilateral"

files.set <- c("CT1","IPSI")

brain.merge <- merge(CT1_A1, y = IPSI_D1, 
                     add.cell.ids = files.set, 
                     #merge.data = T
)
DefaultAssay(brain.merge) <- "Spatial.008um"
gc()
VariableFeatures(brain.merge) <- c(VariableFeatures(CT1_A1), VariableFeatures(IPSI_D1))
# rm(CT1_A1)
# rm(IPSI_D1)
gc()

# we select 50,0000 cells and create a new 'sketch' assay
brain.merge <- SketchData(
  object = brain.merge,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
gc()

# switch analysis to sketched cells
DefaultAssay(brain.merge) <- "sketch"

# perform clustering workflow
brain.merge <- FindVariableFeatures(brain.merge)
brain.merge <- ScaleData(brain.merge)
brain.merge <- RunPCA(brain.merge, assay = "sketch", reduction.name = "pca.sketch")
ElbowPlot(brain.merge, ndims = 50, reduction = "pca.sketch")
gc()
brain.merge <- FindNeighbors(brain.merge, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
brain.merge <- FindClusters(brain.merge, cluster.name = "seurat_cluster.sketched", resolution = 3)
brain.merge <- RunUMAP(brain.merge, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
gc()

brain.merge <- ProjectData(
  object = brain.merge,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
gc()

DefaultAssay(brain.merge) <- "sketch"
Idents(brain.merge) <- "seurat_cluster.sketched"
p1 <- DimPlot(brain.merge, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(brain.merge) <- "Spatial.008um"
Idents(brain.merge) <- "seurat_cluster.projected"
p2 <- DimPlot(brain.merge, reduction = "full.umap.sketch", label = T, raster = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")
p2
p1 | p2
FeaturePlot(brain.merge, features = "Apoe", raster = F)

SpatialDimPlot(brain.merge, label = T, repel = T, label.size = 4)

# Save and Load data
saveRDS(object = brain.merge, file = "obj_brain.merge.Rds")
#rm(prot.combined)
brain.merge <- readRDS("./obj_brain.merge.Rds")


SpatialFeaturePlot(brain.merge, features = c("Apoe")
                   ,slot = "data"
)

p2 <- SpatialFeaturePlot(brain.merge, features = "Apoe", slot = "counts") #+ ggtitle("Apoe expression (8um)")
p2

#### plots

DotPlot(brain.merge, features = c( "Apoe",
                                   "Gfap", "Aqp4", "Lcn2", "Gja1", "Slc1a2", "Fgfr3", "Nkain4",   #astrocytes
                                   "Flt1", "Cldn5", "Vtn", "Itm2a", "Vwf", "Fam167b", "Bmx", "Clec1b",    #endothelial_cells
                                   "Slc17a6",  "Slc17a7",  "Nrgn", "Camk2a", "Satb2", "Col5a1", "Sdk2", "Nefm",    #excitatory_neurons
                                   "Slc32a1",  "Gad1", "Gad2", "Tac1", "Penk", "Sst",  "Npy",  "Mybpc1", "Pvalb", "Gabbr2",   #inhibitory_neurons
                                   "P2ry12", "Csf1r",  "Cd74", "C3", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia
                                   "Olig2",  "Mbp",  "Mobp", "Plp1", "Mog",  "Cldn11", "Myrf", "Galc", "Ermn", "Mag",   #oligodendrocytes
                                   "Vcan", "Cspg4", "Pdgfra", "Sox10", "Neu4", "Pcdg15", "Gpr37l1", "C1ql1", "Cdo1", "Epn2",   #oligodendrocyte_precursor_cells
                                   "Ambp",  "Higd1b", "Cox4i2", "Aoc3", "Pde5a",  "Pth1r",  "P2ry14", "Abcc9", "Kcnj8", "Cd248", #Pericytes
                                   "Ccr2","Cd68","Cd11b","Cd14","Fcgr3", #monocytes
                                   "Ccr5","Itgam","Trfc","Fcgr1" #macrophage
                                   
),
col.max = 20, #idents = c("62","67","49","47"),
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()


p1 <- VlnPlot(brain.merge, features = c("Apoe"),#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
              #split.by = "disease",
              pt.size = 0.05,
              raster = F,
              #ncol = 1,
              group.by = "ID",
              #slot = "counts",
              #add.noise = F,
              #log = T,
              #sort = "increasing",
              idents = c("17","56","16","29","50") #endothelial
              #c("24","26","30")#microglia
) + scale_y_continuous(limits = c(0.000,8.5)) + #geom_boxplot(width=0.1, color="black", alpha=0.2) +
  stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95)
p1$layers[[2]]$aes_params$alpha <- 0.1
p1

###############################
## astrocytes

DotPlot(brain.merge, features = c( "Apoe",
                                   "Gfap", "Aqp4", "Lcn2", "Gja1", "Slc1a2", "Fgfr3", "Nkain4"   #astrocytes
),
col.max = 20, 
idents = c("7","22","42","16"),
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
scale = F,
#split.by = "cohort"
) + RotatedAxis()


## myeloids

DotPlot(brain.merge, features = c( "Apoe",
                                   "P2ry12", "Csf1r",  "Cd74", "C3", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia
                                   "Ccr2","Cd68","Cd11b","Cd14","Fcgr3", #monocytes
                                   "Ccr5","Itgam","Trfc","Fcgr1" #macrophage                    
),
cols = c("blue","blue"),
col.max = 20, 
idents = c("17","56","16","29","50"),
dot.scale = 10, 
cluster.idents = F, #group.by = "patho",
#scale = F,
split.by = "ID"
) + RotatedAxis()


## find and visualize the top gene expression markers for each cluster
# Create downsampled object to make visualization either
DefaultAssay(brain.merge) <- "Spatial.008um"
Idents(brain.merge) <- "seurat_cluster.projected"
object_subset <- subset(brain.merge, cells = Cells(brain.merge[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE) %>%
  group_by(cluster) #%>%
# dplyr::filter(avg_log2FC > 1)%>%
# slice_head(n = 50)
write.xlsx(as.data.frame(markers), rowNames = T, file="wilcox_clus_all_markers_IPSI_50pcs_res3.xlsx")
markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) #+ NoLegend()
ggsave(filename = "brain.merge_heat_top10_HD.png",
       plot = p,
       width = 35,
       height = 40,
       #dpi = 600,
       device = "png")



Idents(brain.merge) <- "seurat_cluster.projected"
cells <- CellsByIdentities(brain.merge, idents = c(54:60))
p <- SpatialDimPlot(brain.merge,
                    images = "slice1.008um",
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

Idents(brain.merge) <- "seurat_cluster.projected"
cells <- CellsByIdentities(brain.merge, idents = c(54:60))
p <- SpatialDimPlot(brain.merge,
                    images = "slice1.008um.2",
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p










######################################      First attempt (SCTransform --> did not work, too heavy)
############      control in A1
# CT1_A1 <- Load10X_Spatial(data.dir = "/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/izzy/visium_cancer_joe/spaceranger/CT1_A1/outs/binned_outputs/square_002um/", slice = "slice1")#, slice = "slice1") # dir should contain filtered_feature_bc_matrix.h5
# 
# plot1 <- VlnPlot(CT1_A1, features = "nCount_Spatial", pt.size = 0.0, raster = F) + NoLegend()
# plot2 <- SpatialFeaturePlot(CT1_A1, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)
# # P1 <- plot1 + plot2
# # P1
# 
# ### Normalization  -> standard approaches (such as the LogNormalize() function), which force each data point to have the same underlying ‘size’ after normalization, can be problematic.
# CT1_A1 <- SCTransform(CT1_A1, assay = "Spatial")
# 
# ### Gene expression visualization
# SpatialFeaturePlot(CT1_A1, features = c("Apoe", "Cd36"), slot = "counts")
# 
# # plot <- SpatialFeaturePlot(CT1_A1, features = c("Ttr")) + theme(legend.text = element_text(size = 0),
# #                                                                legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
# # jpeg(filename = "../output/images/spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
# # print(plot)
# # dev.off()
# # p1 <- SpatialFeaturePlot(CT1_A1, features = "Ttr", pt.size.factor = 1) ## This will scale the size of the spots. Default is 1.6
# # p2 <- SpatialFeaturePlot(CT1_A1, features = "Ttr", alpha = c(0.1, 1)) ## minimum and maximum transparency, setting to alpha c(0.1, 1) to downweight the transparency of points with lower expression
# # p1 + p2
# 
# 
# #### Dimensionality reduction, clustering, and visualization
# CT1_A1 <- RunPCA(CT1_A1, assay = "SCT")
# Elbowplot(CT1_A1, npcs = 50)
# CT1_A1 <- FindNeighbors(CT1_A1, reduction = "pca", dims = 1:30)
# CT1_A1 <- FindClusters(CT1_A1)
# CT1_A1 <- RunUMAP(CT1_A1, reduction = "pca", dims = 1:30)
# 
# p1 <- DimPlot(CT1_A1, reduction = "umap", label = TRUE)
# p2 <- SpatialDimPlot(CT1_A1, label = TRUE, label.size = 3)
# p1 + p2
# SpatialDimPlot(CT1_A1, cells.highlight = CellsByIdentities(
#   object = CT1_A1, idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3)
# 
# 
# ####### Identification of Spatially Variable Features
# de_markers <- FindMarkers(CT1_A1, ident.1 = 5, ident.2 = 6)
# SpatialFeaturePlot(object = CT1_A1, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
# ## An alternative approach, implemented in FindSpatiallyVariables(), is to search for features exhibiting spatial patterning in the absence of pre-annotation
# ## The default method (method = 'markvariogram), is inspired by the Trendsceek, which models spatial transcriptomics data as a mark point process and computes a ‘variogram’, which identifies genes whose expression level is dependent on their spatial location
# CT1_A1 <- FindSpatiallyVariableFeatures(CT1_A1, assay = "SCT", features = VariableFeatures(CT1_A1)[1:1000],
#                                        selection.method = "moransi")
# # Now we visualize the expression of the top 6 features identified by this measure
# top.features <- head(SpatiallyVariableFeatures(CT1_A1, selection.method = "moransi"), 6)
# SpatialFeaturePlot(CT1_A1, features = top.features, ncol = 3, alpha = c(0.1, 1))
# 
# 
# 
# 
# 
# ###########        control in A1 and IPSI in D1
# CT1_A1 <- Load10X_Spatial(data.dir = "/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/izzy/visium_cancer_joe/spaceranger/CT1_A1/outs/binned_outputs/square_008um/", slice = "slice1") #%>% SCTransform(assay = "Spatial")
# gc()
# CT1_A1 <- SCTransform(CT1_A1, assay = "Spatial")
# gc()
# IPSI_D1 <- Load10X_Spatial(data.dir = "/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/izzy/visium_cancer_joe/spaceranger/IPSI_D1/outs/binned_outputs/square_002um/", slice = "slice1") #%>% SCTransform(assay = "Spatial")
# gc()
# IPSI_D1 <- SCTransform(IPSI_D1, assay = "Spatial")
# gc()
# 
# plot1 <- SpatialFeaturePlot(CT1_A1, features = "Apoe", slot = "counts") #+ theme(legend.position = "right")
# plot2 <- SpatialFeaturePlot(IPSI_D1, features = "Apoe", slot = "counts") #+ theme(legend.position = "right")
# wrap_plots(plot1, plot2)
# 
# 
# brain.merge <- merge(CT1_A1, IPSI_D1)
# 
# DefaultAssay(brain.merge) <- "SCT"
# VariableFeatures(brain.merge) <- c(VariableFeatures(CT1_A1), VariableFeatures(IPSI_D1))
# brain.merge <- RunPCA(brain.merge, verbose = FALSE)
# Elbowplot(brain.merge, npcs = 50)
# brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
# brain.merge <- FindClusters(brain.merge)
# brain.merge <- RunUMAP(brain.merge, dims = 1:30)
# 
# DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
# 
# SpatialDimPlot(brain.merge)
# 
# SpatialFeaturePlot(brain.merge, features = c("Apoe", "Cd36"))
# 
# 
# # Save and Load data
# saveRDS(object = prot.combined, file = "obj_unintegrated_CSF.Rds")
# #rm(prot.combined)
# prot.combined <- readRDS("./obj_unintegrated_CSF.Rds")
# 
# 
# 




















