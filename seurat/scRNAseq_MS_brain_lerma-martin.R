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

setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/scRNaseq_spatial_lerma-martin_et_al/Seurat/")

# prot.combined <- readRDS("./GSE301908_sn_all.rds") # pre-done

file.dir <- "../cellranger/"
data.list <- c()
# meta1 <- read.delim("meta.txt", header = T)

files.set <- c(
  "GSM8563681_CO37",
  "GSM8563682_CO40",
  "GSM8563683_CO41",
  "GSM8563684_CO45",
  "GSM8563685_CO74",
  "GSM8563686_CO85",
  "GSM8563687_MS197D",
  "GSM8563688_MS229",
  "GSM8563689_MS377N",
  "GSM8563690_MS377T",
  "GSM8563691_MS377I",
  "GSM8563692_MS411",
  "GSM8563693_MS497I",
  "GSM8563694_MS497T",
  "GSM8563695_MS549H",
  "GSM8563696_MS549T")

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i], "/filtered_feature_bc_matrix/")
  data <- Read10X(data.dir = path)
  dataset_name <- files.set[i]
  # condition <- sub("(.*)-.*", "\\1", files.set[i])
  # treatment1 <- treatment[i]
  mat <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 500) %>% 
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5) %>%
    NormalizeData(normalization.method = "LogNormalize") %>% FindVariableFeatures(selection.method = "vst")
  #SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)# %>%
  #RunPCA(npcs = 50)
  # mat$treatment <- treatment1
  mat$sample <- dataset_name
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
prot.combined <- FindClusters(prot.combined, resolution = 1.5, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
prot.combined <- RunUMAP(prot.combined, #umap.method = "umap-learn", 
                         dims = 1:40, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
        #ncol = 2,
        label = T,
        repel = T,
        #group.by = "seurat_clusters"
        # split.by = "disease"
)#, combine = F)


## Save and Load V5 data
saveRDS(object = prot.combined, file = "obj_unintegrated_scaled_by.MT.counts.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_scaled_by.MT.counts.Rds")

## Add metadata
meta <- read.delim("meta.txt", header = T)

prot.combined$GEM <- sub("(.*)_.*", "\\1", prot.combined$sample)
prot.combined$pathology <- prot.combined$GEM
prot.combined$disease <- prot.combined$GEM
prot.combined$age <- prot.combined$GEM
prot.combined$sex <- prot.combined$GEM
prot.combined$donor <- prot.combined$GEM

for (i in 1:length(meta$GEM)) {
  prot.combined$pathology <- recode(prot.combined$pathology, "meta$GEM[i] = meta$Lesion.type[i]")
  prot.combined$age <- recode(prot.combined$age, "meta$GEM[i] = meta$Age[i]")
  prot.combined$sex <- recode(prot.combined$sex, "meta$GEM[i] = meta$Sex[i]")
  prot.combined$disease <- recode(prot.combined$disease, "meta$GEM[i] = meta$Condition[i]")
  prot.combined$donor <- recode(prot.combined$donor, "meta$GEM[i] = meta$patient[i]")
}
prot.combined$pathology <- factor(prot.combined$pathology, levels = c("CTRL","CA","CI"))

############################ cell types annotation
################ label transfer from previously sequenced scRNAseq
allen <- readRDS("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/public_data/MS/snRNAseq_fagiani_el.al_2025/Seurat/obj_unintegrated_3000feats.annotated.Rds")


transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = prot.combined,
  dims = 1:40,
  reduction = 'pcaproject'
)
gc()
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$celltypes,
  weight.reduction = prot.combined[['pca']],
  dims = 1:40
)
gc()
prot.combined <- AddMetaData(object = prot.combined, metadata = predicted.labels)

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
        #ncol = 2,
        label = T,
        repel = T,
        # group.by = "predicted.id"
        group.by = "celltypes"
        # split.by = "treatment"
)#, combine = F)


Idents(prot.combined) <- "predicted.id"

############ regular cell type annotation

DotPlot(prot.combined, features = c( "GFAP","SLC1A3","AQP4","LCN2", "GJA1", "SLC1A2","FGFR3","NKAIN4",   #Astrocytes
                                     "SDC2",
                                     "FBLN1","FBLN5", # fibroblasts
                                     "CHRM3", # cholinergic neurons
                                     "TH","SLC18A2", # dopaminergic neurons
                                     "TAGLN","MYH11", # vascular smooth muscle cells
                                     "CFAP44","CFAP43", # ependymal cells
                                     "SLC17A6","SLC17A7","NRGN","CAMK2A", "SATB2", "COL5A1","SDK2","NEFM","HTR2C",   #Excitatory_neurons
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
prot.combined <- BuildClusterTree(prot.combined, assay = "RNA", reduction = "umap.unintegrated", reorder = T)
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

Idents(prot.combined) <- "unintegrated_clusters"

prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Oligodendrocytes",
                              `1`= "Oligodendrocytes",                     # 
                              `2`= "Oligodendrocytes",  
                              `3`= "Oligodendrocytes",
                              `4`= "Oligodendrocytes", 
                              `5`= "Oligodendrocytes",
                              `6`= "Oligodendrocytes",
                              `7`= "Oligodendrocytes",
                              `8`= "Oligodendrocytes",
                              `9`= "Astrocytes",
                              `10`= "Oligodendrocytes",
                              `11`= "Oligodendrocytes",
                              `12`= "Astrocytes",
                              `13`= "Glutamatergic neurons",
                              `14`= "Oligodendrocytes",                   # 
                              `15`= "Myeloids",
                              `16`= "Oligodendrocytes precursor cells",
                              `17`= "Oligodendrocytes",
                              `18`= "Myeloids",
                              `19`= "Myeloids",
                              `20`= "Oligodendrocytes",
                              `21`= "Oligodendrocytes",
                              `22`= "Myeloids",  ##  *************** 
                              `23`= "Astrocytes",
                              `24`= "Oligodendrocytes",     
                              `25`= "Oligodendrocytes precursor cells",
                              `26`= "Excitatory cholinergic neurons",          # 
                              `27`= "Astrocytes",              #
                              `28`= "Astrocytes",
                              `29`= "Oligodendrocytes",
                              `30`= "Endothelial cells",           # 
                              `31`= "Astrocytes",
                              `32`= "Inhibitory cholinergic neurons",           # 
                              `33`= "Pericytes",         
                              `34`= "T cells",
                              `35`= "Myeloids",
                              `36`= "Astrocytes",
                              `37`= "Ependymal cells",
                              `38`= "Myeloids",             
                              `39`= "Astrocytes",                      # 
                              `40`= "Fibroblasts",
                              `41`= "Oligodendrocytes",                  # 
                              `42`= "Apoptotic cells",
                              `43`= "B cells",
                              `44`= "Astrocytes",
                              `45`= "GABAergic neurons",
                              `46`= "Excitatory cholinergic neurons",
                              `47`= "Inhibitory cholinergic neurons",
                              `48`= "Oligodendrocytes precursor cells",
                              `49`= "Vascular smooth muscle cells",
                              `50`= "Astrocytes"
)

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
# Idents(prot.combined) <- "seurat_clusters"

## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats.annotated.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_3000feats.annotated.Rds")



############### analysis

clus4.markers <- c("SDC2","PPIF","PTGS2","NLRP3","HBEGF","IL1B","NFKB1","FOSL2","GPR183","EGR1","EGR2","EGR3","CXCL8",
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
            reduction = "umap.unintegrated", 
            # min.cutoff = 0.7,
            # max.cutoff = 1,
            cols = c("gray90","#9E1021"),
            split.by = "pathology"
) + theme(legend.position = "right") + xlab("UMAP_1") + ylab("UMAP_2")

Idents(prot.combined) <- "celltypes"

my_comparisons <- list(c("CTRL","CA"),c("CTRL","CI"),c("CA","CI"))
colors <- brewer.pal(n=8,name = "Dark2")
colors <- colors[c(1,3,2)]
gene <- "SDC2"
p1 <- VlnPlot(prot.combined, features = gene,
              pt.size = 0.05, raster = F, group.by = "pathology", cols = colors,
              idents = c("Myeloids")
              # idents = c("Microglia","MgND microglia","Monocytes")
              # idents = c("Astrocytes")
              # idents = c("CD8+ T cells")
              # idents = c("mNK")
              # idents = c("Classical")
) + theme(legend.position = "none") + xlab("")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 20, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.00001, max(p1[[1]][["data"]][[gene]])+0.3*max(p1[[1]][["data"]][[gene]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
  xlab("")
p1$layers[[2]]$aes_params$alpha <- 0.1
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
clus4.markers <- c("SDC2","PTGS2","EGR2","EGR3","CXCL8", "G0S2",  "RGCC", "TRIB1", "PPIF","NLRP3","NFKB1","GPR183","EGR1","HBEGF","IL1B","FOSL2"#,
) ## myeloids
clus4.markers <- c("SDC2","EGR3","PPIF","TRIB1", #"PTGS2","HBEGF","IL1B","FOSL2","EGR2","CXCL8",  "NLRP3","NFKB1","GPR183","EGR1","RGCC",
                   "G0S2") ## monocytes
clus4.markers <- c("SDC2","EGR3","PPIF","PTGS2","EGR2","CXCL8",  "RGCC" #, #"NLRP3","NFKB1","GPR183","EGR1","TRIB1", "HBEGF","IL1B","FOSL2","G0S2"
) ## monocytes microglia

clus4.markers <- c("SDC2","SDC4","G0S2","EGR2","EGR3","CCL3L3","IRAK2","DUSP2","B3GNT5","TSPOAP1","HBEGF","PTGS2","TRIB1","PPIF","RGCC", # 4
                   "NLRP3","NFKB1","IL1B","CXCL8", "GPR183")


plcs2 <- c("ACSL3","ACSL4","ACSL1",  "ACAT1","ACAT2",
           "PLA2G7",#"CES1","ABCA1","ABCG1",
           # "FADS3","PTGS2",  "FADS1","FADS2","ALOX5AP",
           "LTA4H",
           "PTGES3","TBXAS1"
)
clus4.markers <- c("SDC2","HBEGF","SDC4","G0S2","TSPO","PPIF","RGCC" )


DotPlot(prot.combined, features = plcs2,
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        idents = c(#"Microglia",
                   #"MgND microglia",
                   "Myeloids"
        ),
        cluster.idents = F, 
        group.by = "pathology",
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

########## subset microglia

myeloids <- subset(prot.combined, subset = celltypes %in% c("Myeloids"))

### Normalize and scale merged obj
myeloids <- NormalizeData(myeloids, normalization.method = "LogNormalize")
myeloids <- FindVariableFeatures(myeloids, selection.method = "vst", nfeatures = 1000)
myeloids <- ScaleData(myeloids, vars.to.regress = c("percent.mt","nCount_RNA"), #latent.data = "nFeature_RNA", 
                      model.use = "linear")#, features = all.genes)
gc()

### Dimensionality reduction and integration
myeloids <- RunPCA(myeloids, npcs = 50)
gc()
ElbowPlot(myeloids, ndims = 50)
myeloids <- FindNeighbors(myeloids, dims = 1:30, reduction = "pca")
myeloids <- FindClusters(myeloids, resolution = .5, cluster.name = "unintegrated_clusters")
myeloids <- RunUMAP(myeloids, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

colors1 <- paletteer_d("RColorBrewer::Set3")[1:10]
colors1 <- colors1[c(8,4,2,5)]
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

DimPlot(myeloids, reduction = "umap.unintegrated", 
        raster = F,# pt.size = 2,
        ncol = 3,
        label = T,
        repel = T,
        cols = colors1,
        # group.by = "unintegrated_clusters",
        group.by = "celltypes",
        # group.by = "predicted.id"
        split.by = "pathology"
) + ggtitle("White matter") + xlab("UMAP_1") + ylab("UMAP_2") + #xlim(-8,9) + ylim(-7,6) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow(type = "closed")))

myeloids <- BuildClusterTree(myeloids, assay = "RNA", reduction = "umap.unintegrated", reorder = T)
prot.markers <- FindAllMarkers(myeloids, assay = "RNA", only.pos = T, min.pct = 0.2, logfc.threshold = 0.33) %>% group_by(cluster)
write.xlsx(as.data.frame(prot.markers), rowNames = T, file="all.clusters.myeloids.markers.xlsx")

markers1 <- c("P2RY12","P2RY13","CX3CR1","SELPLG","TMEM119","CSF1R","SLC2A5","C3",  # "SLC2A3",
              # "IL1B","TLR2","ALOX5", #"HIF1A","IL4R",
              "CD74","CST3","APOE", #"IRF8","CST7","CLEC7A","AIF1","HEXB","PTGES3","NFKB1","C1QA", 
              "MS4A4A","CD163","SIGLEC1","F13A1","MARCO","IL12RB2","IL15",
              "PLXND1","TNFAIP2","EMILIN2","MAFB","GAS7",#"CD4","GDA",#"TBXAS1",
              "GPR183","PLA2G7","ASAH1","FOSL2","ADIPOR1","LPCAT1","SDC2" #, "PPIF","PTGS2","FLVCR2","HBEGF",
              # "ATP6V1B2","ATP6V1E1","ATP6V1H", ## "ATP6V1C1","ATP6V0A1","ATP6V1A",
              # "ATP6V0A3","ATP6V0C","ATP6V0D1","ATP6V0E1"   #"ATP6V0E2","ATP6V0D2","ATP6V0A4",    "ATP6V0A2"
              #,"EGR1","EGR2","EGR3"#,"G0S2","NLRP3","RGCC","CXCL8","TRIB1",
              # "CD8A","CD8B","CD3E","IL2RB","KLRB1","IL7R","PRF1","GZMH","GZMA","GZMB","NKG7",
              # "CD19","CD22","IGHM"  #,"TRAC","TRBC1"#,"TRDC","TRGC1","GZMK","CD3D",
)
# markers1 <- markers1[37:1]

DotPlot(myeloids, features = markers1,
        col.max = 20,
        dot.scale = 10, 
        cols = "RdBu",
        # cluster.idents = T, 
        group.by = "celltypes",
        # group.by = "unintegrated_clusters",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip() + xlab("") + ylab("")

Idents(myeloids) <- "unintegrated_clusters"
myeloids <- RenameIdents(myeloids, 
                         `0`= "Homeostatic microglia",
                         `1`= "Homeostatic microglia",
                         `2`= "Homeostatic microglia",
                         `3`= "MgND microglia", # Foam cells
                         `4`= "Homeostatic microglia",
                         `5`= "MgND microglia",
                         `6`= "MgND microglia",
                         `7`= "Homeostatic microglia", # # Foam cells
                         `8`= "MgND microglia", # Inflammatory microglia
                         `9`= "Homeostatic microglia",
                         `10`= "Foam cells",
                         `11`= "Homeostatic microglia",
                         `12`= "Foam cells",
                         `13`= "Monocytes",
                         `14`= "Foam cells" #,
                         # `15`= "Monocytes",
                         # `16`= "Homeostatic microglia" # Inflammatory microglia
)
myeloids$celltypes <- sub("(.*)_.*", "\\1", Idents(myeloids))
Idents(myeloids) <- "celltypes"
myeloids$celltypes <- factor(myeloids$celltypes, levels = c("Homeostatic microglia",
                                                            # "Inflammatory microglia",
                                                            "MgND microglia",
                                                            "Monocytes",
                                                            "Foam cells"
))


## Save and Load data
saveRDS(object = myeloids, file = "myeloids_reg.out.mito.ncounts.patients.annotated.v3.Rds")
# #rm(myeloids)
myeloids <- readRDS("./myeloids_reg.out.mito.ncounts.patients.annotated.v3.Rds")

##### freq
num.cells <- as.data.frame(table(myeloids$pathology, myeloids$pathology))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(myeloids$pathology, myeloids$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.cells.myeloids_disease_lesion.csv")

# paletteLength <- 40
# colors <- colorRampPalette( rev(brewer.pal(11, "Set3")))(paletteLength)
# n <- 60
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
colors1 <- paletteer_d("RColorBrewer::Set3")[1:10]
colors1 <- colors1[c(8,4,2,5)]
barplot(tfreq.num.cells.celltype, col = colors1, #col = sample(col_vector, n), 
        legend.text = rownames(tfreq.num.cells.celltype),
        xlim = c(0,7.5), main = "Myeloid cell types distribution", xlab = "", ylab = "Cells percentage (%)", cex.names = 0.9)

myeloids$lesion.patient <- paste(myeloids$pathology, myeloids$sample, sep = "_")
num.cells <- as.data.frame(table(myeloids$lesion.patient, myeloids$lesion.patient))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(myeloids$lesion.patient, myeloids$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.cells.myeloids_lesion_patient.csv")

mydata <- read.csv("freq.num.cells.myeloids_lesion_patient2.csv", header = T)

coluna1 <- c("CTRL WM","PLWM","DMWM")

my_comparisons <- list(c("CTRL WM","PLWM"),c("CTRL WM","DMWM"),c("PLWM","DMWM")
)

p1 <- ggplot(mydata, aes_string(x="factor(lesion,levels=coluna1)", y="Foam.cells", color="lesion")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  # scale_color_brewer(palette="Set1") + 
  scale_color_brewer(palette="Dark2") +
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("cells percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD8+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  # ggtitle("COX2+ Classical monocytes") +
  # ggtitle("Nonclassical monocytes") +
  ggtitle("Foam cells") +
  # ggtitle("Intermediate monocytes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1


clus4.markers <- c("SDC2","PPIF","PTGS2","NLRP3","HBEGF","IL1B","NFKB1","FOSL2","GPR183","EGR1","EGR2","EGR3","CXCL8",
                   "G0S2","RGCC","TRIB1")

clus4.markers <- c("PPIF","PTGS2","PTGES3","NLRP3","HBEGF","FOSL2","GPR183",
                   "CXCL8","G0S2","RGCC","SDC2")

myeloids <- AddModuleScore(
  myeloids,
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
FeaturePlot(myeloids, features = gene, 
            raster = F,
            # pt.size = 2.,
            # ncol = 3,
            reduction = "umap.unintegrated", 
            # min.cutoff = 0.7,
            # max.cutoff = 1,
            cols = c("gray90","#9E1021"),
            # split.by = "lesion_type"
) + theme(legend.position = "right")

# Idents(myeloids) <- "celltypes"
Idents(myeloids) <- "type_broad"

my_comparisons <- list(c("CTRL","CA"),c("CTRL","CI"),c("CA","CI"))
colors <- brewer.pal(n=8,name = "Dark2")
colors <- colors[c(1,3,2)]

my_comparisons <- list(c("Homeostatic microglia","MgND microglia"),c("Homeostatic microglia","Monocytes"),c("Homeostatic microglia","Foam cells"),
                       c("MgND microglia","Monocytes"),c("MgND microglia","Foam cells"),c("Monocytes","Foam cells"))
colors <- paletteer_d("RColorBrewer::Set3")[1:10]
colors <- colors[c(8,4,2,5)]
gene <- "SDC2"
p1 <- VlnPlot(myeloids, features = gene,
              pt.size = 0.05, raster = F, 
              # group.by = "celltypes", 
              cols = colors,
              group.by = "pathology",
              # idents = c("Microglia")
              # idents = c("Astrocytes")
              # idents = c("CD8+ T cells")
              # idents = c("mNK")
              # idents = c("Classical")
) + theme(legend.position = "none") + xlab("") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 20, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000001, max(p1[[1]][["data"]][[gene]])+0.3*max(p1[[1]][["data"]][[gene]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
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
           "FADS3", #"FADS2","FADS1",
           "ALOX5AP","LTA4H","PTGES3","TBXAS1"
)

ptgs <- c("ABCC4","SLCO2A1","PTGES2","PTGS2")

eicosanoids2 <- c("ABCC4","PTGES2","PTGS2","FADS1","ACSL1","PLA2G10","PLA2G4A","PLBD1")

genes1 <- c("SDC2","HBEGF","CD14","FCGR3A","ABCC4","PLA2G7","LTA4H","FADS3","ACSL4",  #"PTGS2",
            "PTGES3","TBXAS1")

myeloids$celltypes.lesion <- paste(myeloids$celltypes, myeloids$lesion_type)

DotPlot(myeloids, features = ferroptosis,
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        # idents = "Microglia",
        cluster.idents = F, 
        group.by = "celltypes",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip() + xlab("") + ylab("")

choline  <- c("SCARB1","MFSD2A","FABP3","LPCAT1","LDLR","APOC1","RAB38","ACSL3","ABCA3","CAPN2")
lipids <- c("ACSL1","ACSL3","ACSL5","PLIN2")

genes1 <- c("CPT2","CPT1A","SLC25A20","LCLAT1","HADHA","ALDH9A1", # carnitine
            "ADIPOR1", # lipid metabolism and storage
            "PLSCR1", # phosphatidyl serine
            "GPX4","TFRC","SLC25A37","STEAP4","SLC11A1","SLC7A11",  # ferroptosis
            # "UGCG",  # glycosphingolipids
            "PGS1", # phosphatidyl glycerol
            # "TSPO","TSPOAP1",#"STARD4",#"SREBF2", # cholesterol
            "FLVCR2","LPCAT1","CHKA","CHKB","ETNK1","ETNK2","PCYT1A","PCYT1B","PCYT2","CEPT1","SELENOI","GPAT4","PLA2G4A","PTDSS1"  # choline "GPCPD1",,"LCAT","PLB1","PLA2G7","PLCE1","PLCH1","ATP8A2"
)  
genes1 <- c("MBOAT1","MBOAT2","HHATL","MBOAT4","MBOAT7","LPCAT4","LPCAT3","LPCAT2","LPCAT1","PLA2G4A") ## lands cycle

genes1 <- c("CHPT1","CHKB","ETNK1","ETNK2","PCYT1A","PCYT1B","PCYT2","CEPT1","SELENOI") ## kennedy pathway

genes1 <- c("SLC25A48","CHDH","CHAT","SLC18A3") # "FLVCR2","ALDH9A1" 

genes1 <- c("MFSD2A","MFSD2B","SPNS1","SPNS2","ATP8B1","SGMS1","SGMS2","SMPD1")  ## sphingomyelin metabolism

genes1 <- c("ACAA1","CYBA","PARK7","ECH1","GSTK1","DHRS4","RPS27A", ## peroxisome related associated with PIRA2
            "UBA52","UBB","UBE2D1","FPR2","GNAI2","PEX16","PMVK",  # "ACSL3","RAC1","ACSL1","RAC2",
            "GSTP1","ITGB2","TYROBP","CFL1","CYP1B1","DUSP1",  # "RHOB","CLEC7A","ITGAM",
            "ERN1","RPS3","TXN","KLF4","PPIF","KDM6B","PRDX5","PRDX4","PRDX3","PRDX2","PRDX1","SOD1"
) # "SOD2","SCP2","HSD17B4","CYBB","IDH1","UBE2D3","FYN","TRPM2","OSER1","EEF2","KLF2","SIRPA",

genes1 <- c("ACAA1","PARK7","RAC2","RPS27A","UBA52","FPR2","ACSL1","ACSL3", ## peroxisome related associated with foam cells cluster
            ## "SOD1","PRDX5","PRDX4","PRDX3","RPS3","ITGAM","ITGB2","TYROBP","PEX16","PMVK","GNAI2","UBB","UBE2D1","CYBA","ECH1","GSTK1","DHRS4","RAC1",
            "GSTP1","CLEC7A","RHOB","CFL1","CYP1B1","DUSP1","ERN1","TXN","KLF4","PPIF","KDM6B","SIRPA","PRDX2","PRDX1"
) # "SOD2","SCP2","HSD17B4","CYBB","IDH1","UBE2D3","FYN","TRPM2","OSER1","EEF2","KLF2",

genes1 <- c("AMD1","OAZ1","ODC1","SAT1","SRM","SMS","SMOX","DHPS","DOHH","EIF5A")  ## polyamines

genes1 <- c("LYST","NCF2","BIN2","RAB14","CLCN3","PRKCD","SPG11","ITGAL","P2RX7","PECAM1", ## phagocytosis "CD93","CD36","IRF8",
            "TICAM2","VAV1","CCR2","ICAM3","PAK1","ARHGAP25","ITGB2","ABL1", ### phagocytosis  "FCN1","ANXA1","ITGAM","CORO1A","RAC1","FCGR1A","NCF4","CD14","VAV2","HCK",
            "ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1E1","ATP6V1H","ATP6V0A2",  ## ROS "ATP6V0A1",
            "LIPA","DNAJC13","HEATR5A","DENND1A","FKBP15","LRRK2","RIN2","SNX10"  ## endocytosis  "ASGR1","PYCARD","SNX17","CAP1",
)

endocytosis <- c("CD93","LYST","SNX17","CLCN3","CORO1A","ITGAL","P2RX7","PECAM1",
                 "LIPA","DNAJC13","FCGR1A","CAP1","CCR2","HEATR5A","ASGR1","PYCARD",
                 "DENND1A","FKBP15","ITGB2","HCK","LRRK2","CD36","RIN2","SNX10")
genes1 <- c(
  #"CYBB","RAC1","RAC2","HVCN1","NCF1",
  "ATP6V1D","ATP6V1F","ATP6V1G1","ATP6V1G2","ATP6V1B1","ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1E1","ATP6V1H", ## ROS
  "ATP6V0A3","ATP6V0C","ATP6V0D1","ATP6V0E1","ATP6V0D2","ATP6V0A4","ATP6V0A1"  ## "ATP6V1C2","ATP6V0A2","ATP6V1G3","ATP6V1E2","ATP6V0E2",
)

genes1 <- c("FLVCR2","ATG7",  "ATG5","SLC33A1","FTL","FTH1","PCBP2", ## "SLC39A8","NCOA4","TFRC","SLC40A1",
            "SLC11A1"  ,"GSS","GCLM","STEAP4","GPX4" ## "TMEM164","SLC7A11","SLC25A37","STEAP3", # ferroptosis #
)

genes1 <- c("FAR1","FAR2","GNPAT","AGPS","AGPAT3","PEDS1") ## plasmalogen peroxisome biosynthesis

genes1 <- c("PTGS2","HBEGF","IL1B","FOSL2","EGR2","EGR3","CXCL8",  "PPIF","NLRP3","NFKB1","EGR1","TRIB1",  # clus4.markers  "GPR183","RGCC",
            "G0S2","SDC2")

prot.combined$celltypes.lesion <- paste(prot.combined$pathology, prot.combined$celltypes, sep = "_")
Idents(prot.combined) <- "celltypes.lesion"
Idents(prot.combined) <- "celltypes"

prot.combined$celltypes <- factor(prot.combined$celltypes, levels = c("Microglia","Foam cells"))
prot.combined$celltypes.lesion <- factor(prot.combined$celltypes.lesion, levels = c("control WM_Microglia","control WM_Foam cells",
                                                                                    "control cortex_Microglia","control cortex_Foam cells",
                                                                                    "myelinated cortex_Microglia","myelinated cortex_Foam cells",
                                                                                    "demyelinated cortex_Microglia","demyelinated cortex_Foam cells"))

DotPlot(myeloids, features = genes1,
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        # idents = c(#"Microglia","Foam cells"
        #   "control WM_Microglia","control cortex_Microglia",#"control cortex_Foam cells",
        #   "myelinated cortex_Microglia",#"myelinated cortex_Foam cells",
        #   "demyelinated cortex_Microglia","demyelinated cortex_Foam cells"
        # ),
        cluster.idents = F, 
        group.by = "celltypes",
        # scale = F,
        #split.by = "disease"
) + xlab("") + ylab("") + coord_flip() + theme(axis.text.x = element_text(angle = 50, hjust = 1))#+ RotatedAxis() 

DotPlot(myeloids, features = genes1,
        # cols = c("gray90","#9E1021"),
        cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        # idents = c("Microglia","Foam cells"
        # ),
        cluster.idents = F, 
        group.by = "pathology",
        # group.by = "disease",
        # scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip() + xlab("") + ylab("")




















