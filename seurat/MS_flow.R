# libraries
library(car)
library(pheatmap)
library(openxlsx)
library(RColorBrewer)
library(ggplot2)
library(matrixStats)
library(tidyverse)
library(zoo)
library(paletteer)
library(EnhancedVolcano)
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
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(ggtext)

setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/MS/scRNAseq_monocytes/flow/")

mydata <- read.xlsx("percentage_03-Sep-2025_patrick_MFI_all.xlsx")
mydata0 <- mydata[!is.na(mydata$Statistic),]
mydata0$sample <- sub("(.*) Well_0[0-9][0-9]_Plate_001.fcs.*","\\1",mydata0$Name) 
mydata0$group <- sub(".* Well_0[0-9][0-9]_Plate_001.fcs/(.*$)","\\1",mydata0$Name)
mydata0 <- mydata0[!(mydata0$sample %in% c("B9")),]


# meta.data <- read.xlsx("meta.xlsx")
meta.data <- read.csv("final_meta_no.out.txt")
meta.data <- meta.data[meta.data$FLOW_CYTOMETRY == "yes",]
mydata0$disease <- mydata0$sample
mydata0$sex <- mydata0$sample
mydata0$EDSS <- mydata0$sample
mydata0$PIRA <- mydata0$sample

for (i in 1:length(meta.data$Flow_ID)) {
  mydata0$disease <- car::recode(mydata0$disease, "meta.data$Flow_ID[i] = meta.data$MS3[i]")
  mydata0$sex <- car::recode(mydata0$sex, "meta.data$Flow_ID[i] = meta.data$sex[i]")
  mydata0$EDSS <- car::recode(mydata0$EDSS, "meta.data$Flow_ID[i] = meta.data$edss[i]")
  mydata0$PIRA <- car::recode(mydata0$PIRA, "meta.data$Flow_ID[i] = meta.data$PIRA3[i]")
}

# cellgroups <- factor(mydata0$group)

target.cell.groups <- c(
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/intermediate mo",
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/intermediate mo/COX2+",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/non-classical mo",
                        # "Cells/Single Cells/Live/CD45+/CD3- B/Lin-/non-classical mo B" ,
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/CD362+"#,
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/COX2+",
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/EGR2+",
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/GPR183+",
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/HBEGF+"#,
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+/CD25+CD27+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD27+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+/FOXP3+ CD27+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD8+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD8+/CD27+"
                        #,"","","","",""
                        )

cell.groups <- c(
                        "Intermediate monocytes",
                        # "Intermediate monocytes COX2+",
                        "Nonclassical monocytes",
                        # "Nonclassical monocytes B" ,
                        "Classical monocytes",
                        "Classical monocytes SDC2+"#, # "Classical monocytes CD362+",
                        # "Classical monocytes COX2+",
                        # "Classical monocytes EGR2+",
                        # "Classical monocytes GPR183+",
                        # "Classical monocytes HBEGF+"#,
                        # "CD4+CD25+ T cells",
                        # "CD4+CD25+CD27+ T cells",
                        # "CD4+CD27+ T cells",
                        # "CD4+FOXP3+ T cells",
                        # "CD4+FOXP3+CD27+ T cells",
                        # "CD25+ T cells",
                        # "CD25+CD27+ T cells",
                        # "CD27+ T cells"#,
                        # "FOXP3+ T cells",
                        # "FOXP3+CD27+ T cells"#,
                        # "CD4+ T cells",
                        # "CD8+ T cells",
                        # "CD8+CD27+ T cells"
                        #,"","","","",""
)

# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo", ]
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+", ]
mydataRef1 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3- A/Lin-", ]
mydataRef2 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3- B/Lin-", ]
mydataRef <- mydataRef1
mydataRef$`#Cells` <- mydataRef1$`#Cells` + mydataRef2$`#Cells`
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+", ]
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD8+", ]
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD4+", ]
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+", ]
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+", ]
mydata1 <- mydata0[mydata0$group == target.cell.groups[4], ]
mydata2 <- data.frame(percent = 100*(mydata1$`#Cells` / mydataRef$`#Cells`), 
                      # percent = 100-(100*(mydata1$`#Cells` / mydataRef$`#Cells`)), # negative percentage
                      disease = mydata1$disease,
                      PIRA = mydata1$PIRA,
                      sample = mydata1$sample,
                      sex = mydata1$sex,
                      ncells = mydata1$`#Cells`)
mydata2$disease.sex <- paste(mydata2$disease, mydata2$sex, sep = "_")
mydata2$disease <- factor(mydata2$disease, levels = c("HC","RRMS","ProgMS"))

write.xlsx(as.data.frame(mydata2), rowNames = T,file="/media/patrick/GERVAZIO/Bioinfo/weiner_lab/MS/scRNAseq_monocytes/Seurat/figures_paper4/flow/SDC2_on.CD3.CD20.CD56.neg.CD45.pos.xlsx")
write.xlsx(as.data.frame(mydata2), rowNames = T,file="/media/patrick/GERVAZIO/Bioinfo/weiner_lab/MS/scRNAseq_monocytes/Seurat/figures_paper4/flow/SDC2_on.classic.mono.xlsx")
# mydata2$disease.sex <- factor(mydata2$disease.sex, levels = c("HC","RRMS","ProgMS"))
colors <- brewer.pal(n=4,name = "Dark2")
colors <- colors[c(1,3,2,4)]
plt <- ggbetweenstats(data = mydata2, x = PIRA, y = percent,
                      # plot.type = "violin",
                      p.adjust.method = "none", 
                      xlab = "", ylab = "percentage (%)",
                      ggplot.component = list(theme(axis.text.x = element_text(size = 14, face = "bold"),
                                                    axis.text.y = element_text(size = 14, face = "bold"),
                                                    axis.title.y = element_text(size = 14, face = "bold"))),
                      ggsignif.args = list(textsize = 4),
                      boxplot.args = list(width = 0),
                      point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 1, size = 3, stroke = 0),
                      violin.args = list(width = 0.5, alpha = 0.2),
                      package = "RColorBrewer", palette = "Dark2",
                      results.subtitle = FALSE,
                      mean.plotting = FALSE,
                      mean.point.args = list(size = 0, color = "darkred"),
                      mean.label.args = list(size = 0),
                      messages = FALSE,
                      # package = "ggsci", palette = "nrc_npg", 
                      # package = "ggsci", palette = "uniform_startrek",
                      # package = "yarrr", palette = "info2",
                      type = "r"
)
plt

coluna1 <- c("HC","RRMS","ProgMS")
my_comparisons <- list(c("HC","RRMS"),
                       c("HC","ProgMS"),
                       c("ProgMS","RRMS"))
coluna1 <- c("HC","RRMS","ProgMS","PIRA")
my_comparisons <- list(c("HC","RRMS"),c("HC","PIRA"),
                       c("HC","ProgMS"),c("PIRA","RRMS"),
                       c("ProgMS","RRMS"),c("PIRA","ProgMS"))


p1 <- ggplot(mydata2, aes_string(x="factor(PIRA,levels=coluna1)", y="percent", color="PIRA")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+CD3-CD20-CD56- cells percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD8+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  # ggtitle("COX2+ Classical monocytes") +
  # ggtitle("Nonclassical monocytes") +
  ggtitle("COX2+ Intermediate monocytes") +
  # ggtitle("Intermediate monocytes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1

for (i in 1:length(target.cell.groups)) {
  mydata1 <- mydata0[mydata0$group == target.cell.groups[i], ]
  mydata2 <- data.frame(percent = 100*(mydata1$`#Cells` / mydataRef$`#Cells`), 
                        PIRA = mydata1$PIRA,
                        disease = mydata1$disease,
                        sample = mydata1$sample,
                        sex = mydata1$sex,
                        ncells = mydata1$`#Cells`)
  # mydata2$disease <- factor(mydata2$disease, levels = c("HC","RRMS","PMS"))
  # p1 <- ggplot(mydata2, aes_string(x="factor(disease,levels=coluna1)", y="percent", color="disease")) +
  p1 <- ggplot(mydata2, aes_string(x="factor(PIRA,levels=coluna1)", y="percent", color="PIRA")) +
    geom_violin(trim=T) + 
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
    geom_jitter(position=position_jitter(0.2)) +
    geom_boxplot(width=0.1) + 
    scale_color_brewer(palette="Dark2") + 
    #scale_fill_brewer(palette="Dark2") + 
    stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
    #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
    # stat_compare_means() +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
    # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
    # theme_minimal() + rotate_x_text(angle = 45)+ ylab("T cells percentage (%)") + xlab("") +
    theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+CD3-CD20-CD56- cells percentage (%)") + xlab("") +
    # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD8+ T cells percentage (%)") + xlab("") +
    # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ T cells percentage (%)") + xlab("") +
    # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
    # theme_minimal() + rotate_x_text(angle = 45)+ ylab("Treg percentage (%)") + xlab("") +
    # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
    # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
    ggtitle(cell.groups[i]) +
    # ggtitle("COX2+ Classical monocytes") +
    # ggtitle("Nonclassical monocytes") +
    # ggtitle("COX2+ Intermediate monocytes") +
    # ggtitle("Intermediate monocytes") +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          axis.text.x = element_text(size = 14, face = "bold"))
  ggsave(#filename = paste0("Lin.neg.A.B/wilcox_percentage_",cell.groups[i],"_on.CD45.pos_by.disease.png"),
         # filename = paste0("final_figs/pdf/percentage/wilcox_percentage_",cell.groups[i],"_on.CD3.CD20.CD56.neg.CD45.pos_by.disease.pdf"),
         # filename = paste0("final_figs/pdf/percentage/wilcox_percentage_",cell.groups[i],"_on.classic.mono_by.disease.pdf"),
         # filename = paste0("final_figs/pdf/percentage/wilcox_percentage_",cell.groups[i],"_on.classic.mono_by.PIRA.pdf"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_Tregs_by.disease.png"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_T.cells_by.disease.png"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_classic.mono_by.disease.png"),
         # filename = paste0("wilcox/percentage/vln_percentage_",cell.groups[i],"_on_classic.mono_by.PIRA.png"),
         # filename = paste0("wilcox/percentage/wilcox_percentage_",cell.groups[i],"_on.CD45.pos_by.PIRA.png"),
         filename = paste0("final_figs/pdf/percentage/wilcox_percentage_",cell.groups[i],"_on.CD3.CD20.CD56.neg.CD45.pos_by.PIRA.pdf"),
         plot = p1,
         # width = 3.,
         width = 3.5,
         height = 5,
         dpi = 600,
         device = "pdf")
  rm(p1)
  rm(mydata1)
  rm(mydata2)
  gc()
}

for (i in 1:length(target.cell.groups)) {
  mydata1 <- mydata0[mydata0$group == target.cell.groups[i], ]
  mydata2 <- data.frame(percent = 100*(mydata1$`#Cells` / mydataRef$`#Cells`), 
                        disease = mydata1$disease,
                        sample = mydata1$sample,
                        sex = mydata1$sex,
                        ncells = mydata1$`#Cells`)
  mydata2$disease <- factor(mydata2$disease, levels = c("HC","RRMS","ProgMS"))
  p1 <- ggbetweenstats(data = mydata2, x = disease, y = percent,
                       # plot.type = "violin",
                       p.adjust.method = "none", 
                       xlab = "", 
                       # ylab = "CD45+CD3-CD20-CD56- cells percentage (%)",
                       ylab = "classical monocytes percentage (%)",
                       # ylab = "CD45+ cells percentage (%)",
                       ggplot.component = list(theme(axis.text.x = element_text(size = 14, face = "bold"),
                                                     axis.text.y = element_text(size = 14, face = "bold"),
                                                     axis.title.y = element_text(size = 10, face = "bold"))),
                       ggsignif.args = list(textsize = 4),
                       boxplot.args = list(width = 0),
                       point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 1, size = 3, stroke = 0),
                       violin.args = list(width = 0.5, alpha = 0.2),
                       package = "RColorBrewer", palette = "Dark2",
                       results.subtitle = FALSE,
                       title = cell.groups[i],
                       # mean.plotting = FALSE,
                       # mean.point.args = list(size = 0, color = "darkred"),
                       # mean.label.args = list(size = 0),
                       # messages = FALSE,
                       # package = "ggsci", palette = "nrc_npg", 
                       # package = "ggsci", palette = "uniform_startrek",
                       # package = "yarrr", palette = "info2",
                       type = "r"
  )
  ggsave(#filename = paste0("ggbetweenstats/robust_holm_percentage_",cell.groups[i],"_on.CD45.pos_by.disease.png"),
         filename = paste0("ggbetweenstats/robust_none_percentage_",cell.groups[i],"_on.classic.mono_by.disease.png"),
         # filename = paste0("ggbetweenstats/robust_none_percentage_",cell.groups[i],"_on.CD3.CD20.CD56.neg.CD45.pos_by.disease.png"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_Tregs_by.disease.png"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_T.cells_by.disease.png"),
         plot = p1,
         width = 4,
         height = 7.5,
         dpi = 600,
         device = "png")
  rm(p1)
  rm(mydata1)
  rm(mydata2)
  gc()
}












coluna1 <- c("HC_M","HC_F","RRMS_M","RRMS_F","PMS_M","PMS_F")

my_comparisons <- list(c("HC_M","HC_F"),
                       c("RRMS_M","RRMS_F"),
                       c("PMS_M","PMS_F"))

my_comparisons <- list(c("HC_M","RRMS_M"),
                       c("HC_M","PMS_M"),
                       c("PMS_M","RRMS_M"),
                       c("HC_F","RRMS_F"),
                       c("HC_F","PMS_F"),
                       c("PMS_F","RRMS_F"))

p1 <- ggplot(mydata2, aes_string(x="factor(disease.sex,levels=coluna1)", y="percent", color="disease.sex")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD8+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  ggtitle("Classical monocytes CD362+") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1

mydata1 <- mydata0[mydata0$group == target.cell.groups[11], ]
df_freqs <- mydata1
p1 <- ggplot(df_freqs, aes_string(x="factor(disease,levels=coluna1)", y="`#Cells`", color="disease")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("Number of cells") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  ggtitle("FOXP3+ CD27+") +
  theme(plot.title = element_text(hjust = 0.5))
p1

p1 <- ggplot(df_freqs, aes_string(x="factor(disease,levels=coluna1)", y="`Statistic`", color="disease")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  ggtitle("FOXP3+ CD27+") +
  theme(plot.title = element_text(hjust = 0.5))
p1


##################   MFI   ######################

target.cell.groups <- c(
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/intermediate mo/COX2+",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/CD362+"#,
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/COX2+",
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/EGR2+",
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/GPR183+",
                        # "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/HBEGF+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+/CD25+CD27+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD27+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+/FOXP3+ CD27+",
                        # "Cells/Single Cells/Live/CD45+/CD3+/CD8+/CD27+"
                        #,"","","","",""
)

cell.groups <- c(
  # "COX2 - Intermediate monocytes",
  "SDC2 - Classical monocytes"#,  # "CD362 - Classical monocytes",
  # "COX2 - Classical monocytes",
  # "EGR2 - Classical monocytes",
  # "GPR183 - Classical monocytes",
  # "HBEGF - Classical monocytes",
  # "CD25 - CD4+ T cells",
  # "CD27 - CD4+CD25+ T cells",
  # "CD27 - CD4+ T cells",
  # "FOXP3 - CD4+ T cells",
  # "CD27 - CD4+FOXP3+ T cells",
  # "CD27 - CD8+ T cells"
  #,"","","","",""
)


group_idx <- which(mydata0$group == target.cell.groups[2])
group_mean <- group_idx + 1
group_median <- group_idx + 2

mydata2 <- data.frame(mean = mydata0$Statistic[group_mean],
                      median = mydata0$Statistic[group_median],
                      disease = mydata0$disease[group_idx],
                      PIRA = mydata0$PIRA[group_idx],
                      sample = mydata0$sample[group_idx],
                      sex = mydata0$sex[group_idx])
mydata2$disease.sex <- paste(mydata2$disease, mydata2$sex, sep = "_")


coluna1 <- c("HC","RRMS","ProgMS")
my_comparisons <- list(c("HC","RRMS"),
                       c("HC","ProgMS"),
                       c("ProgMS","RRMS"))
coluna1 <- c("HC","RRMS","ProgMS","PIRA")
my_comparisons <- list(c("HC","RRMS"),c("HC","PIRA"),
                       c("HC","ProgMS"),c("PIRA","RRMS"),
                       c("ProgMS","RRMS"),c("PIRA","ProgMS"))

p1 <- ggplot(mydata2, aes_string(x="factor(disease,levels=coluna1)", y="mean", color="disease")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("MFI mean") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  ggtitle("HBEGF - Classical monocytes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1

p1 <- ggplot(mydata2, aes_string(x="factor(PIRA,levels=coluna1)", y="mean", color="PIRA")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("MFI mean") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  ggtitle("CD362 - Classical monocytes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1

p1 <- ggplot(mydata2, aes_string(x="factor(disease,levels=coluna1)", y="median", color="disease")) + 
  geom_violin(trim=T) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
  geom_jitter(position=position_jitter(0.2)) +
  geom_boxplot(width=0.1) + 
  scale_color_brewer(palette="Dark2") + 
  #scale_fill_brewer(palette="Dark2") + 
  stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
  #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  # stat_compare_means() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + # Add pairwise comparisons p-value
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD45+ cells percentage (%)") + xlab("") +
  theme_minimal() + rotate_x_text(angle = 45)+ ylab("MFI median") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("CD25+CD4+ percentage (%)") + xlab("") +
  # theme_minimal() + rotate_x_text(angle = 45)+ ylab("classical monocytes percentage (%)") + xlab("") +
  # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
  ggtitle("CD27 - FOXP3+CD4+ T cells") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"))
p1


for (i in 1:length(target.cell.groups)) {
  group_idx <- which(mydata0$group == target.cell.groups[i])
  group_mean <- group_idx + 1
  group_median <- group_idx + 2
  mydata2 <- data.frame(mean = mydata0$Statistic[group_mean],
                        median = mydata0$Statistic[group_median],
                        disease = mydata0$disease[group_idx],
                        PIRA = mydata0$PIRA[group_idx],
                        sample = mydata0$sample[group_idx],
                        sex = mydata0$sex[group_idx])
  mydata2$disease.sex <- paste(mydata2$disease, mydata2$sex, sep = "_")
  p1 <- ggplot(mydata2, aes_string(x="factor(disease,levels=coluna1)", y="mean", color="disease")) +
  # p1 <- ggplot(mydata2, aes_string(x="factor(PIRA,levels=coluna1)", y="mean", color="PIRA")) +
    geom_violin(trim=T) + 
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
    geom_jitter(position=position_jitter(0.2)) +
    geom_boxplot(width=0.1) + 
    scale_color_brewer(palette="Dark2") + 
    #scale_fill_brewer(palette="Dark2") + 
    stat_summary(fun = mean, geom='point', size = 3, colour = "darkred") +
    #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
    # stat_compare_means() +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") + # Add pairwise comparisons p-value
    theme_minimal() + rotate_x_text(angle = 45)+ ylab("MFI mean") + xlab("") +
    # ggtitle(expression("PTGE"[2]*" / "*"TXB"[2])) +
    ggtitle(cell.groups[i]) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 14, face = "bold"))
  ggsave(#filename = paste0("final_figs/pdf/MFI/vln_mean_MFI_",cell.groups[i],"_by.PIRA.pdf"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_Tregs_by.PIRA.png"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_T.cells_by.PIRA.png"),
         filename = paste0("final_figs/pdf/MFI/vln_mean_MFI_",cell.groups[i],"_by.disease.pdf"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_Tregs_by.disease.png"),
         # filename = paste0("vln_percentage_",cell.groups[i],"_on_T.cells_by.disease.png"),
         plot = p1,
         width = 3.5,
         # width = 3.,
         height = 5,
         dpi = 600,
         device = "pdf")
  rm(p1)
  rm(mydata1)
  rm(mydata2)
  gc()
}




######################        Correlation      #####################################



target.cell.groups <- c("Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/CD362+",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/COX2+",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/EGR2+",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/GPR183+",
                        "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/HBEGF+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+/CD25+CD27+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD27+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+/FOXP3+ CD27+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD4+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD8+",
                        "Cells/Single Cells/Live/CD45+/CD3+/CD8+/CD27+"
                        #,"","","","",""
)

# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo", ]
mydataRef1 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3- A/Lin-", ]
mydataRef2 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3- B/Lin-", ]
mydataRef <- mydataRef1
mydataRef$`#Cells` <- mydataRef1$`#Cells` + mydataRef2$`#Cells`
# mydataRef2 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+", ]
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+", ]
# mydataRef <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD8+", ]
# mydataRef1 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD4+", ]
mydataRef1 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+", ]
# mydataRef1 <- mydata0[mydata0$group == "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+", ]
mydataCD27 <- mydata0[mydata0$group == target.cell.groups[11], ]
mydataCD362 <- mydata0[mydata0$group == target.cell.groups[2], ]
mydata2 <- data.frame(percentCD27CD4 = 100*(mydataCD27$`#Cells` / mydataRef1$`#Cells`), 
                      percentCD362 = 100*(mydataCD362$`#Cells` / mydataRef$`#Cells`),
                      disease = mydataCD27$disease,
                      sample = mydataCD27$sample,
                      EDSS = mydataCD27$EDSS,
                      sex = mydataCD27$sex
                      )
mydata2$disease.sex <- paste(mydata2$disease, mydata2$sex, sep = "_")

# mydata3 <- mydata2[mydata2$disease == "PMS",]
mydata3 <- mydata2[!(mydata2$EDSS %in% c("0")),]
mydata3$EDSS <- as.numeric(mydata3$EDSS)
mydata4 <- mydata3[mydata3$disease == "RRMS",]


ggscatter(mydata3, x = "percentCD27CD4", y = "EDSS", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = 90, label.y = 10) +
  annotate("text", x = 92, y = 10.5, label = (paste0("slope==", coef(lm(mydata3$EDSS~mydata3$percentCD27CD4))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("EDSS")) +
  xlab(expression("CD27+ percentage on CD4+ T cells (%)"))


ggscatter(mydata3, x = "percentCD27CD4", y = "EDSS", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = 90, label.y = 10) +
  annotate("text", x = 92, y = 10.5, label = (paste0("slope==", coef(lm(mydata3$EDSS~mydata3$percentCD27CD4))[2])), parse = T) + 
  ggtitle(expression("All MS samples (RRMS and PMS)")) + ylab(expression("EDSS")) +
  xlab(expression("CD27+ percentage on CD4+ T cells (%)"))

ggscatter(mydata4, x = "percentCD27CD4", y = "EDSS", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = 90, label.y = 10) +
  annotate("text", x = 92, y = 10.5, label = (paste0("slope==", coef(lm(mydata4$EDSS~mydata4$percentCD27CD4))[2])), parse = T) + 
  ggtitle(expression("PMS patients")) + ylab(expression("EDSS")) +
  xlab(expression("CD27+ percentage on CD4+ T cells (%)"))



ggscatterstats(mydata3, percentCD27CD4, EDSS#, type = "nonparametric"
               )

ggscatter(mydata3, x = "percentCD362", y = "EDSS", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = 0.5, label.y = 10) +
  annotate("text", x = 3, y = 11, label = (paste0("slope==", coef(lm(mydata3$EDSS~mydata3$percentCD362))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("EDSS")) +
  # xlab(expression("SCD2+ percentage on CD14+CD3-CD20-CD56-CD45+ (%)"))
  xlab(expression("SDC2+ cells percentage on classical monocytes (%)"))

ggscatter(mydata4, x = "percentCD362", y = "EDSS", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = 0.1, label.y = 10) +
  annotate("text", x = 2.5, y = 11, label = (paste0("slope==", coef(lm(mydata4$EDSS~mydata4$percentCD362))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("EDSS")) +
  # xlab(expression("CD362+ percentage on classical monocytes (%)"))+
  xlab(expression("CD362+ classical monocytes percentage on CD3-CD20-CD56-CD45+ (%)"))

### MFI corr


target.cell.groups <- c(#"Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo",
  "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/CD362+",
  "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/COX2+",
  "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/EGR2+",
  "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/GPR183+",
  "Cells/Single Cells/Live/CD45+/CD3- A/Lin-/classical mo/HBEGF+",
  "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+",
  "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD25+/CD25+CD27+",
  "Cells/Single Cells/Live/CD45+/CD3+/CD4+/CD27+",
  "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+",
  "Cells/Single Cells/Live/CD45+/CD3+/CD4+/FOXP3+/FOXP3+ CD27+",
  #"Cells/Single Cells/Live/CD45+/CD3+/CD4+",
  #"Cells/Single Cells/Live/CD45+/CD3+/CD8+",
  "Cells/Single Cells/Live/CD45+/CD3+/CD8+/CD27+"
  #,"","","","",""
)

group_idx <- which(mydata0$group == target.cell.groups[1])
group_mean <- group_idx + 1
group_median <- group_idx + 2

mydata2 <- data.frame(mean = mydata0$Statistic[group_mean],
                      median = mydata0$Statistic[group_median],
                      disease = mydata0$disease[group_idx],
                      sample = mydata0$sample[group_idx],
                      sex = mydata0$sex[group_idx],
                      EDSS = mydata0$EDSS[group_idx])
mydata2$disease.sex <- paste(mydata2$disease, mydata2$sex, sep = "_")
mydata3 <- mydata2[!(mydata2$EDSS %in% c("0")),]
mydata3 <- mydata3[!is.na(mydata3$EDSS),]
mydata3 <- mydata3[!(mydata3$median > 5000),]
mydata3$EDSS <- as.numeric(mydata3$EDSS)



ggscatter(mydata3, x = "mean", y = "EDSS", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = 3900, label.y = 10) +
  annotate("text", x = 4500, y = 10.5, label = (paste0("slope==", coef(lm(mydata3$EDSS~mydata3$mean))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("EDSS")) +
  xlab(expression("SDC2 MFI mean on classical monocytes"))


ggscatter(mydata3, x = "median", y = "EDSS", #label = rownames(mydata3),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "spearman", label.x = 2900, label.y = 10) +
  annotate("text", x = 3500, y = 10.5, label = (paste0("slope==", coef(lm(mydata3$EDSS~mydata3$median))[2])), parse = T) + 
  ggtitle(expression("")) + ylab(expression("EDSS")) +
  xlab(expression("SDC2 MFI median on classical monocytes"))






