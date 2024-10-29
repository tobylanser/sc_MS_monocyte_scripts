library(Seurat)
library(scCustomize)
library(ggplot2)
library(SCP)
# BiocManager::install("dittoSeq")
library(dittoSeq)


setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
gc()
# merged_females <- readRDS("data/merged_seurat/females/females_merged_object.rds")

merged <- readRDS("Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds")

table(merged$seurat_clusters)
Idents(merged) <- "seurat_clusters"
merged <- RenameIdents(merged, `17` = "10")
merged$seurat_clusters <- Idents(merged)
merged <- RenameIdents(merged, `0` = "0", `1` = "1", `2` = "2", 
                    `3` = "3", `4` = "4", `5` = "5", `6` = "6",
                    `8` = "7", `9` = "8", `10` = "9")
merged$seurat_clusters <- Idents(merged)
merged$seurat_clusters <- factor(merged$seurat_clusters, 
                                 levels=c("0", 
                                          "1",
                                          "2",
                                          "3", 
                                          "4",
                                          "5",
                                          "6",
                                          "7",
                                          "8",
                                          "9"))


levels(x = merged) <- c("0", 
                          "1",
                          "2",
                          "3", 
                          "4",
                          "5",
                          "6",
                          "7",
                          "8",
                          "9",
                          "10")


table(merged$seurat_clusters)
saveRDS(merged, "Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds", compress = F)

clust_colors <- c(
  "#F4CCCC",
  "#E6E6FA",
  "#FFFACD",
  "#DDA0DD",
  "#ADD8E6",
  "#b2df8a",
  "#FFF0F5", 
  "#D3D3D3",
  "#F0F8FF",
  "#FFDAB9",
  "#87CEEB",
  "#FFE5B4",
  "#C7E5C9"
)

clust_colors <- c(
  "#E6B8B8",  # Darkened #F4CCCC
  "#CFCFEF",  # Darkened #E6E6FA
  "#EFE49A",  # Darkened #FFFACD
  "#C185C1",  # Darkened #DDA0DD
  "#8EB5D6",  # Darkened #ADD8E6
  "#9ACB6F",  # Darkened #b2df8a
  "#E6C4C9",  # Darkened #FFF0F5
  "#B0B0B0",  # Darkened #D3D3D3
  "#D4E1F1",  # Darkened #F0F8FF
  "#E6C099",  # Darkened #FFDAB9
  "#6EA8CC",  # Darkened #87CEEB
  "#E6C092",  # Darkened #FFE5B4
  "#A0CFA1"   # Darkened #C7E5C9
)



sccust_celltype_cols <- c("#7ba39d", "#cab2d6","#fdbf6f", "#a6cee3")

SC <- DimPlot(merged, group.by = "seurat_clusters", cols = clust_colors, label = T,
                       reduction = "umap.unintegrated", raster = T)
SC
ggsave("figures_paper2/fig1/cluster_UMAP.pdf", device = "pdf", 
       width = 12, height = 12, units = "in", dpi = 300, plot = SC)

CT <- DimPlot(merged, group.by = "celltypes", cols = sccust_celltype_cols, label = T,
              reduction = "umap.unintegrated", raster = T)
CT
ggsave("figures_paper2/fig1/celltype_UMAP.pdf", device = "pdf", 
       width = 12, height = 12, units = "in", dpi = 300, plot = CT)

table(Idents(merged))
merged$celltype_disease <- factor(merged$celltype_disease, 
                               levels=c("Classical_HC", 
                                        "Classical_MS",
                                        "Intermediate_HC",
                                        "Intermediate_MS", 
                                        "Nonclassical_HC",
                                        "Nonclassical_MS",
                                        "mo_DC_HC",
                                        "mo_DC_MS"))
library(dplyr)

disease_state_cols <- c("#7ba39d", "#fdbf6f", "#5f75ae", "#fb9a99", "#1f78b4")

sccust_celltype_cols <- c("#7ba39d","#fdbf6f", "#cab2d6","#a6cee3")
ditto_celltype_cols <- c("#7ba39d","#fdbf6f","#a6cee3",  "#cab2d6")


DimPlot_scCustom(merged, group.by = "celltype", colors_use = sccust_celltype_cols, raster = F)

barplot <- dittoBarPlot(merged, color.panel = ditto_celltype_cols, group.by = "GEM", var = "celltypes")
barplot
ggsave("figures_paper2/fig1/GEM_celltype_barplot.pdf", device = "pdf", 
       width = 15, height = 10, units = "in", dpi = 600, plot = barplot)

table(merged$disease)
barplot <- dittoBarPlot(merged, color.panel = ditto_celltype_cols, group.by = "prog", var = "celltypes")
barplot
ggsave("figures_paper2/fig1/disease_celltype_barplot.pdf", device = "pdf", 
       width = 8, height = 10, units = "in", dpi = 600, plot = barplot)

barplot <- dittoBarPlot(merged, color.panel = ditto_celltype_cols, group.by = "celltypes", var = "prog")
barplot

table(merged$seurat_clusters)

barplot <- dittoBarPlot(merged, color.panel = clust_colors,
                        group.by = "disease", var = "seurat_clusters",
                        var.labels.rename = c("0", "1", "2", "3", "4", "5",
                                              "6", "7", "8", "9"))
barplot
ggsave("figures_paper2/fig1/disease_cluster_barplot.pdf", device = "pdf", 
       width = 8, height = 10, units = "in", dpi = 600, plot = barplot)

