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
# library(Azimuth)
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

setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snATACseq_RNAseq_spatial_elkjaer/Seurat/")

file.dir <- "../cellranger_GEX/"
data.list <- c()

files.set <- c(
  "C1_WM",
  "C3_WM",
  "C4_WM",
  "MS1_AL",
  "MS1_NAWM",
  "MS3_CA",
  "MS3_RL",
  "MS6_AL",
  "MS6_CA",
  "MS6_NAWM",
  "MS6_RL",
  "MS7_CA",
  "MS7_RL",
  "MS9_NAWM"
)

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i],"/outs/raw_feature_bc_matrix.h5")
  data <- Read10X_h5(filename = path)
  dataset_name <- files.set[i]
  mat <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 200) %>% 
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5) %>%
    NormalizeData(normalization.method = "LogNormalize") %>% FindVariableFeatures(selection.method = "vst")
  mat$sample <- dataset_name
  mat$lesion <- sub("^.*[0-9]_(.*$)","\\1", dataset_name)
  mat$disease <- sub("(^.*)[0-9]_.*$","\\1", dataset_name)
  mat$patient <- sub("(^.*[0-9])_.*$","\\1", dataset_name)
  data.list[[i]] <- mat
  rm(mat)
  rm(data)
}
# Name layers
names(data.list) <- files.set

# Merge layers and create seurat obj during merging 
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
prot.combined <- merge(data.list[[1]], y = data.list[2:length(data.list)], 
                       add.cell.ids = files.set, merge.data = T)


rm(data.list)
gc()
VariableFeatures(prot.combined) <- features

# Visualize QC metrics as a violin plot
#VlnPlot(prot.combined, pt.size = 0.0,group.by = "cohort",raster = F,
#        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(prot.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",raster = F)
#plot2 <- FeatureScatter(prot.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)
#plot1 + plot2

### Normalize and scale merged obj
prot.combined <- NormalizeData(prot.combined, normalization.method = "LogNormalize")
prot.combined <- FindVariableFeatures(prot.combined, selection.method = "vst", nfeatures = 3000)
prot.combined <- ScaleData(prot.combined, vars.to.regress = c("percent.mt","nCount_RNA"), #latent.data = "nFeature_RNA", 
                           model.use = "linear")#, features = all.genes)
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
prot.combined <- FindNeighbors(prot.combined, dims = 1:40, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = .7, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
prot.combined <- RunUMAP(prot.combined, #umap.method = "umap-learn", 
                         dims = 1:40, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
        label = T,
        repel = T,
        # group.by = "seurat_clusters"
        # ncol = 3,
        # split.by = "patient"
)#, combine = F)

## Dimensionality reduction of integrated data
prot.combined <- RunHarmony(prot.combined, group.by.vars = "sample")
gc()
ElbowPlot(prot.combined, ndims = 50, reduction = "harmony")
prot.combined <- RunUMAP(prot.combined, dims = 1:40, reduction = "harmony", reduction.name = "umap")
prot.combined <- FindNeighbors(prot.combined, reduction = "harmony", dims = 1:40)
prot.combined <- FindClusters(prot.combined, resolution = 1.0, cluster.name = "harmony_clusters")
gc()

DimPlot(prot.combined, reduction = "umap", raster = F,
        ncol = 3,
        label = T,
        repel = T,
        # group.by = "seurat_clusters",
        split.by = "lesion"
) #, combine = F)

## Save and Load V5 data
saveRDS(object = prot.combined, file = "obj_harmony_scaled_by.MT.counts.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_harmony_scaled_by.MT.counts.Rds")


prot.combined$lesion <- factor(prot.combined$lesion, levels = c("WM","AL","NAWM","CA","RL"))

# disease <- c("CTRL","RRMS","PMS","PMS")
# diagnosis <- c("CTR","RRMS","PPMS","SPMS")
# prot.combined$disease <- prot.combined$diagnosis
# 
# for (i in 1:length(disease)) {
#   prot.combined$disease <- recode(prot.combined$disease, "diagnosis[i] = disease[i]")
# }
# prot.combined$disease_lesion <- paste(prot.combined$disease, prot.combined$lesion_type, sep = "_")
# prot.combined$disease_lesion <- factor(prot.combined$disease_lesion, levels = c("CTRL_GM","CTRL_WM",
#                                                                                 "RRMS_AL","RRMS_NAWM",
#                                                                                 "PMS_NAGM","PMS_GML",
#                                                                                 "PMS_AL","PMS_NAWM",
#                                                                                 "PMS_CAL","PMS_CIL","PMS_RL"))




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



prot.combined[["RNA"]] <- JoinLayers(prot.combined[["RNA"]])
gc()
prot.combined <- BuildClusterTree(prot.combined, assay = "RNA", reduction = "umap", reorder = T)
gc()
prot.markers <- FindAllMarkers(prot.combined, assay = "RNA", only.pos = T, min.pct = 0.3, logfc.threshold = 0.33) %>% group_by(cluster)
write.xlsx(as.data.frame(prot.markers), rowNames = T, file="all.clusters.markers.xlsx")
prot.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5


DotPlot(prot.combined, features = unique(top5$gene),
        col.max = 20, 
        dot.scale = 10, 
        #cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()

Idents(prot.combined) <- "harmony_clusters"

prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Oligodendrocytes",
                              `1`= "Oligodendrocytes",                     # 
                              `2`= "Oligodendrocytes",  
                              `3`= "Microglia",
                              `4`= "Oligodendrocytes", 
                              `5`= "Oligodendrocytes",
                              `6`= "Oligodendrocytes",
                              `7`= "Oligodendrocytes",
                              `8`= "Astrocytes",
                              `9`= "Oligodendrocytes",
                              `10`= "Oligodendrocytes",
                              `11`= "Oligodendrocytes",
                              `12`= "Oligodendrocytes",
                              `13`= "Oligodendrocyte precursor cells",
                              `14`= "Oligodendrocytes",                   # 
                              `15`= "Oligodendrocytes",
                              `16`= "Excitatory neurons",
                              `17`= "Astrocytes",
                              `18`= "Oligodendrocytes",
                              `19`= "Microglia",
                              `20`= "Inhibitory neurons",
                              `21`= "Oligodendrocytes",
                              `22`= "MgND microglia",  ##  *************** MgND
                              `23`= "Inhibitory neurons",
                              `24`= "Monocytes",     ## monocytes
                              `25`= "Oligodendrocytes",
                              `26`= "Microglia",          # 
                              `27`= "T cells",              #
                              `28`= "B cells",
                              `29`= "Excitatory neurons",
                              `30`= "Endothelial cells",           # 
                              `31`= "Neuroendocrine cells",
                              `32`= "Fibroblasts",           # 
                              `33`= "Oligodendrocyte precursor cells",         
                              `34`= "Pericytes",
                              `35`= "Astrocytes",
                              `36`= "Astrocytes",
                              `37`= "Oligodendrocyte precursor cells",
                              `38`= "Ependymal cells",             
                              `39`= "Excitatory neurons",                      # 
                              `40`= "Oligodendrocyte precursor cells",
                              `41`= "Excitatory neurons",                  # 
                              `42`= "Oligodendrocytes",
                              `43`= "Oligodendrocytes",
                              `44`= "Oligodendrocytes"
)

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
# Idents(prot.combined) <- "seurat_clusters"

## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats.annotated.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_3000feats.annotated.Rds")





clus4.markers <- c("SDC2","PPIF","PTGS2","NLRP3","HBEGF","NFKB1","FOSL2","GPR183","EGR1","EGR3","CXCL8","EGR2","IL1B",
                   "G0S2","RGCC","TRIB1")

prot.combined <- AddModuleScore(
  prot.combined,
  features = clus4.markers,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "clus4.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)










gene <- "SDC2"
FeaturePlot(prot.combined, features = gene, 
            raster = F,
            # pt.size = 2.,
            # ncol = 3,
            reduction = "umap", 
            # min.cutoff = 0.7,
            # max.cutoff = 1,
            cols = c("gray90","#9E1021"),
            # split.by = "lesion"
) + theme(legend.position = "right")

# Idents(prot.combined) <- "celltypes"
Idents(prot.combined) <- "type_broad"

my_comparisons <- list(c("CTRL WM","PLWM"),c("CTRL WM","DMWM"),c("DMWM","PLWM"))
colors <- brewer.pal(n=8,name = "Dark2")
colors <- colors[c(1,3,2)]
gene <- "clus4.markers1"
p1 <- VlnPlot(prot.combined, features = gene,
              pt.size = 0.05, raster = F, group.by = "lesion", #cols = colors,
              # idents = c("Microglia")
              idents = c("8")
              # idents = c("Astrocytes")
              # idents = c("CD8+ T cells")
              # idents = c("mNK")
              # idents = c("Classical")
) + theme(legend.position = "none") + xlab("")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene]])+0.3*max(p1[[1]][["data"]][[gene]]))) #+
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
  xlab("")
p1$layers[[2]]$aes_params$alpha <- 0.5
p1



eicosanoids <- c("PTGS2","FADS1","FADS2","FADS3","EGR2","EGR3","PLA2G10","PLA2G4A","DAGLB","PLA2G4E",
                 "PLA2G3", "PLBD1","ACSL1","ACSL2","ACSL3","ACSL4",
                 "PLA2G4F","PLA2G6","PLA2G7","PLA2G5","PLA2G4C","ABCC4","PTGES2","SLCO2A1")

plcs <- c("ABCC4","SLCO2A1",
          "PTGS1","PTGS2","PTGES2",
          "PTGES","PTGES3",
          "ALOX5","ALOX5AP","ALOX12","ALOX15",
          "ELOVL2","FADS1","ELOVL5","FADS2","FADS3",
          "DAGLA","DAGLB",
          "PLA2G5","PLA2G2F","PLA2G12A","PLA2G10","PLA2G12B","PLA2G4A","PLA2G4E","PLA2G7","PLA2G4C","PLA2G6",
          "PLA2G2D","PLA2G2E","PLA2G2C","PLA2G3","PLA2G1B","PLA2G2A","PLA2G4F","PLBD1",
          "ACSL3","ACSL4","ACSL1", "ACSL2"
)

plcs2 <- c("ACSL3","ACSL4","ACSL1",
           "PLA2G7",
           "FADS1","FADS2","FADS3",
           "ALOX5AP","LTA4H","PTGES3","TBXAS1"
)

ptgs <- c("ABCC4","SLCO2A1","PTGES2","PTGS2")

eicosanoids2 <- c("ABCC4","PTGES2","PTGS2","FADS1","ACSL1","PLA2G10","PLA2G4A","PLBD1")

genes1 <- c("SDC2","HBEGF","CD14","FCGR3A","ABCC4","PLA2G7","LTA4H","FADS3","ACSL4",  #"PTGS2",
            "PTGES3","TBXAS1")

plcs2 <- c("ACSL3","ACSL4","ACSL1",  "ACAT1","ACAT2",
           "PLA2G7",
           "FADS3","PTGS2",  "FADS1","FADS2","ALOX5AP","LTA4H",
           "PTGES3","TBXAS1"
)

ptgs <- c("ABCC4","SLCO2A1","PTGES2","PTGS2")

genes1 <- c("SDC2","HBEGF","CD14","FCGR3A","ABCC4","PLA2G7","LTA4H","FADS3","ACSL4",  #"PTGS2",
            "PTGES3","TBXAS1")
clus4.markers <- c("SDC2","PTGS2","EGR2","EGR3","CXCL8", "G0S2"#  "RGCC", ,"TRIB1", "PPIF","NLRP3","NFKB1","GPR183","EGR1","HBEGF","IL1B","FOSL2"#,
) ## myeloids
clus4.markers <- c("SDC2","EGR3","PPIF","TRIB1", #"PTGS2","HBEGF","IL1B","FOSL2","EGR2","CXCL8",  "NLRP3","NFKB1","GPR183","EGR1","RGCC",
                   "G0S2") ## monocytes
clus4.markers <- c("SDC2","EGR3","PPIF","PTGS2","EGR2","CXCL8",  "RGCC" #, #"NLRP3","NFKB1","GPR183","EGR1","TRIB1", "HBEGF","IL1B","FOSL2","G0S2"
) ## monocytes microglia


DotPlot(prot.combined, features = plcs2,
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        idents = c("Microglia",
                   #"MgND microglia",
                   "Monocytes"
        ),
        cluster.idents = F, 
        group.by = "lesion",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip() + xlab("") + ylab("")




############### DEG analysis

cell.types <- c(
  "Microglia"
)
# cell.types <- as.character(1:7)

prot.combined$cluster.disease <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")
prot.combined$celltypes.disease <- paste(prot.combined$type_broad, prot.combined$disease, sep = "_")
prot.combined$celltypes.lesion <- paste(prot.combined$type_broad, prot.combined$lesion_type, sep = "_")
prot.combined$cluster.lesion <- paste(prot.combined$seurat_clusters, prot.combined$lesion_type, sep = "_")
# Idents(prot.combined) <- "cluster.disease"
Idents(prot.combined) <- "celltypes.lesion"


# ## wilcox
lesions <- c("AL","GM","GML","NAGM","NAWM","WM"
             # ,"CAL" #,"CIL","RL",
)
up.list <- c()
down.list <- c()

for (i in 1:length(cell.types)) {
  for ( e in 1:length(lesions)) {
    zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_CIL"),
                                ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                slot = "data",
                                assay = "RNA",
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
    zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_CIL_x_",lesions[e],"_", cell.types[i], "_by.age.sex_DEGs.xlsx"))
    zk.response1 <- zk.response0[zk.response0$avg_log2FC > 0.322,]
    zk.response2 <- zk.response0[zk.response0$avg_log2FC < -0.322,]
    gc()
    zk.response0b <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_CAL"),
                                 ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                 slot = "data",
                                 assay = "RNA",
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
    zk.response0b <- zk.response0b[zk.response0b$p_val_adj < 0.1,]
    # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_CAL_x_",lesions[e],"_", cell.types[i], "_DEGs.xlsx"))
    zk.response1b <- zk.response0b[zk.response0b$avg_log2FC > 0.322,]
    zk.response2b <- zk.response0b[zk.response0b$avg_log2FC < -0.322,]
    gc()
    zk.response0c <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_RL"),
                                 ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                 slot = "data",
                                 assay = "RNA",
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
    zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_RL_x_",lesions[e],"_", cell.types[i], "_DEGs.xlsx"))
    zk.response1c <- zk.response0c[zk.response0c$avg_log2FC > 0.322,]
    zk.response2c <- zk.response0c[zk.response0c$avg_log2FC < -0.322,]
    gc()
    zk.response1e <- zk.response1b[row.names(zk.response1b) %in% row.names(zk.response1),]
    zk.response1e <- zk.response1e[row.names(zk.response1e) %in% row.names(zk.response1c),]
    zk.response2e <- zk.response2b[row.names(zk.response2b) %in% row.names(zk.response2),]
    zk.response2e <- zk.response2e[row.names(zk.response2e) %in% row.names(zk.response2c),]
    up.list[[e]] <- row.names(zk.response1e)
    down.list[[e]] <- row.names(zk.response2e)
    rm(zk.response1)
    rm(zk.response2)
    rm(zk.response0)
    rm(zk.response1b)
    rm(zk.response2b)
    rm(zk.response0b)
    rm(zk.response1c)
    rm(zk.response2c)
    rm(zk.response0c)
    rm(zk.response1e)
    rm(zk.response2e)
    gc()
  }
}

names(up.list) <- lesions
up.list2 <- t(plyr::ldply(up.list, rbind))
colnames(up.list2) <- up.list2[1,]
up.list2 <- up.list2[-c(1), ] 
write.xlsx(as.data.frame(up.list2), rowNames = F,file="wilcox_microglia_chronic_lesions_converged_up_per.non.chronic_DEGs.xlsx")
names(down.list) <- lesions
down.list2 <- t(plyr::ldply(down.list, rbind))
colnames(down.list2) <- down.list2[1,]
down.list2 <- down.list2[-c(1), ] 
write.xlsx(as.data.frame(down.list2), rowNames = F,file="wilcox_microglia_chronic_lesions_converged_down_per.non.chronic_DEGs.xlsx")

lesions <- c("AL","GM","GML","NAGM","NAWM","WM"
             # ,"CAL" #,"CIL","RL",
)
common_up_genes <- Reduce(intersect, list(up.list2[,1],up.list2[,3],up.list2[,4],up.list2[,5],up.list2[,6]))



for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(myeloids, ident.1 = c(paste0(cell.types[i], "_CAL"),
                                                    paste0(cell.types[i], "_CIL"),
                                                    paste0(cell.types[i], "_RL")),
                              ident.2 = c(paste0(cell.types[i], "_AL"),
                                          paste0(cell.types[i], "_GM"),
                                          paste0(cell.types[i], "_GML"),
                                          paste0(cell.types[i], "_NAGM"),
                                          paste0(cell.types[i], "_NAWM"),
                                          paste0(cell.types[i], "_WM")),
                              slot = "data",
                              assay = "RNA",
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
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_chronic_lesions_x_all_others_", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}


## MAST

lesions <- c("AL","GM","GML","NAGM","NAWM","WM"
             #,"CAL" #,"CIL","RL",
)
for (i in 1:length(cell.types)) {
  for ( e in 1:length(lesions)) {
    zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_RL"),
                                ident.2 = paste0(cell.types[i], "_",lesions[e]),
                                slot = "data",
                                assay = "RNA",
                                features = NULL,
                                logfc.threshold = 0,
                                test.use = "MAST",
                                min.pct = 0.0,
                                min.diff.pct = -Inf,
                                verbose = TRUE,
                                only.pos = FALSE,
                                max.cells.per.ident = Inf,
                                random.seed = 1,
                                latent.vars = c("sex","age_cat"),
                                min.cells.feature = 3,
                                min.cells.group = 3,
                                pseudocount.use = 0.001,
                                mean.fxn = NULL,
                                fc.name = NULL,
                                base = 2,
                                densify = FALSE,
                                #recorrect_umi = TRUE
    )
    zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_RL_x_",lesions[e],"_", cell.types[i], "_by.age.sex_DEGs.xlsx"))
    rm(zk.response0)
    gc()
  }
}


for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(cell.types[i], "_CAL"),
                              ident.2 = paste0(cell.types[i], "_AL"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.0,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("sex","age_cat"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 0.001,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              #recorrect_umi = TRUE
  )
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_CAL_x_AL_", cell.types[i], "_minpct0_DEGs.xlsx"))
  rm(zk.response0)
  gc()
  # zk.response0 <- FindMarkers(myeloids, ident.1 = c(paste0(cell.types[i], "_CAL"),paste0(cell.types[i], "_CIL"),paste0(cell.types[i], "_RL")),
  #                             ident.2 = paste0(cell.types[i], "_HC"),
  #                             slot = "data",
  #                             assay = "RNA",
  #                             features = NULL,
  #                             logfc.threshold = 0,
  #                             test.use = "MAST",
  #                             min.pct = 0.0,
  #                             min.diff.pct = -Inf,
  #                             verbose = TRUE,
  #                             only.pos = FALSE,
  #                             max.cells.per.ident = Inf,
  #                             random.seed = 1,
  #                             latent.vars = c("sex","age_cat"),
  #                             min.cells.feature = 3,
  #                             min.cells.group = 3,
  #                             pseudocount.use = 0.001,
  #                             mean.fxn = NULL,
  #                             fc.name = NULL,
  #                             base = 2,
  #                             densify = FALSE,
  #                             #recorrect_umi = TRUE
  # )
  # zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_PMS_x_HC_clus", cell.types[i], "_minpct0.1_DEGs.xlsx"))
  # rm(zk.response0)
  # gc()
  # zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(cell.types[i], "_RRMS"),
  #                             ident.2 = paste0(cell.types[i], "_HC"),
  #                             slot = "data",
  #                             assay = "RNA",
  #                             features = NULL,
  #                             logfc.threshold = 0,
  #                             test.use = "MAST",
  #                             min.pct = 0.1,
  #                             min.diff.pct = -Inf,
  #                             verbose = TRUE,
  #                             only.pos = FALSE,
  #                             max.cells.per.ident = Inf,
  #                             random.seed = 1,
  #                             latent.vars = c("sex"),
  #                             min.cells.feature = 3,
  #                             min.cells.group = 3,
  #                             pseudocount.use = 0.001,
  #                             mean.fxn = NULL,
  #                             fc.name = NULL,
  #                             base = 2,
  #                             densify = FALSE,
  #                             #recorrect_umi = TRUE
  # )
  # zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  # write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_RRMS_x_HC_clus", cell.types[i], "_minpct0.1_DEGs.xlsx"))
  # rm(zk.response0)
  # gc()
}


















