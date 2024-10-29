library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)

setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper")
gc()
# merged_females <- readRDS("data/merged_seurat/females/females_merged_object.rds")

myeloids <- readRDS("Rds/myeloid_BP2_unintegrated_annotated2.v3.Rds")

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


Idents(myeloids) <- "seurat_clusters"
genes1 <- c("PTGS2","ABCC4","PLA2G4A","FOS","FOSL2","EGR2")                                                                               # violins for cluster 4
genes1 <- c("PTGS2","ABCC4","PLA2G4A","FOS","FOSL2","EGR2","PTGES3","TBXAS1")                   # violins for nonclassical and intermediate monocytes (clusters 2 and 6)
genes1 <- c("PTGER4","CD83","PTGER2","NFKB1","NFKBIZ","IL1B","RIPK2","PLAUR","NAMPT","BTG2","ATF3","TBXAS1",
            "HLA-DRB5","LGALS2","LGALS3","MAT2A","RELA","TNFRSF10B")                                                # violins for mo-DC (cluster 8)
my_comparisons <- list( c("progressor","RR"), c("progressor","non-progressor"), c("progressor","HC"),
                        c("RR","HC"),c("RR","non-progressor"),c("non-progressor","HC"))
for (i in 1:length(genes1)){
  p1 <- VlnPlot(myeloids, features = genes1[i],
                pt.size = 0.05,
                raster = FALSE,
                group.by = "prog",
                idents = "8"
  ) 
  p1 <- p1 + 
    stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+2.8)) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + 
    scale_fill_manual(values = clust_colors)  # Apply custom colors
  
  p1$layers[[2]]$aes_params$alpha <- 0.1
  
  ggsave(filename = paste0("../../figures_paper2/violins/vln_",genes1[i],"_mean_mo.DC_by.prog.pdf"),
         plot = p1,
         width = 7.5,
         height = 12,
         device = "pdf")
  
  rm(p1)
}
