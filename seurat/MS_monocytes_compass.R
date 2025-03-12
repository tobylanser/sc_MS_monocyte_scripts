
# Install compassR from YosefLab.
# devtools::install_github("YosefLab/compassR")
setwd("/home/overflow/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper/compass")
setwd("/home/toby/Desktop/Partners HealthCare Dropbox/Toby Lanser/MS_monocytes_paper/compass")
library(compassR)

# micro <- read.table("../COMPASS_outs/single_cell_benari_outs/micropools.tsv", sep = "\t", header = T)
# meta <- read.csv("ilc3_metadata.csv")
# meta$cluster <- micro$microcluster
# write.csv(meta, "ilc3_metadata.csv")
compass_settings <- CompassSettings$new(
  user_data_directory = system.file(cell_metadata_path = "compass_input/myeloidmeta_data.csv",
                                    gene_metadata = "extdata/RECON2/gene_metadata.csv",
                                    reaction_metadata_path = "extdata/RECON2/reaction_metadata.csv",
                                    package = "compassR"),
  gene_id_col_name = "HGNC.symbol",
  cell_id_col_name = "Name"
)



compass_settings$metabolite_metadata_path
compass_settings$cell_metadata_path <- "compass_input/myeloidmeta_data.csv"
compass_settings$compass_reaction_scores_path <- "COMPASS_outs/reactions.tsv"
compass_settings$reaction_metadata_path <- "extdata/RECON2/reaction_metadata.csv"
compass_settings$linear_gene_expression_matrix_path <- "compass_input/myeloid_counts.csv"
compass_data <- CompassData$new(compass_settings)

compass_analyzer <- CompassAnalyzer$new(compass_settings)

library(tidyverse)
gc()
group_A_cell_ids <-
  compass_data$cell_metadata %>%
  filter(Condition == "non-progressor_mo-DC") %>%
  pull(microcluster)
group_B_cell_ids <-
  compass_data$cell_metadata %>%
  filter(Condition == "RR_mo-DC") %>%
  pull(microcluster)
wilcoxon_results <- compass_analyzer$conduct_wilcoxon_test(
  compass_data$reaction_consistencies,
  group_A_cell_ids,
  group_B_cell_ids,
  for_metareactions = F
)

library(ggrepel)
facets <- c( "Glycolysis", "TCA cycle",
             "Fatty acid oxidation", "Amino acid metabolism")

compass_scores_by_cell_type <-
  wilcoxon_results %>%
  left_join(
    dplyr::select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"),
    by = "reaction_id"
  ) %>%
  left_join(
    compass_data$reaction_metadata,
    by = "reaction_no_direction"
  ) %>%
  # Keep only "confident reactions", as defined in our paper.
  filter(!is.na(EC_number)) %>%
  filter(confidence == "0" | confidence == "4") %>%
  # Exclude non-mitochondrially localized reactions from TCA.
  mutate(subsystem = case_when(
    reaction_id == "SPMDOX_pos" ~ "Arginine and Proline Metabolism",
    subsystem == "Citric acid cycle" & !grepl("[m]", formula, fixed = TRUE) ~ "Other",
    TRUE ~ subsystem
  )) %>%
  # Assign reactions to the appropriate subsystem.
  mutate(
    subsystem_priority = factor(subsystem) %>%
      fct_recode(
        "Glycolysis" = "Glycolysis/gluconeogenesis",
        "TCA cycle" = "Citric acid cycle"
      ) %>%
      fct_collapse("Amino acid metabolism" = c(
        "Alanine and aspartate metabolism",
        "Arginine and Proline Metabolism",
        "beta-Alanine metabolism",
        "Cysteine Metabolism",
        "D-alanine metabolism",
        "Folate metabolism",
        "Glutamate metabolism",
        "Glycine, serine, alanine and threonine metabolism",
        "Histidine metabolism",
        "Lysine metabolism",
        "Methionine and cysteine metabolism",
        "Taurine and hypotaurine metabolism",
        "Tryptophan metabolism",
        "Tyrosine metabolism",
        "Urea cycle",
        "Valine, leucine, and isoleucine metabolism"
      )) %>%
      fct_other(keep = facets) %>%
      fct_relevel(facets)
  ) %>%
  # Keep only the subsystems for which we want to plot a facet.
  filter(subsystem_priority != "Other") %>%
  # Lower-bound the adjusted p-value.
  mutate(adjusted_p_value = if_else(
    subsystem_priority == "Amino acid metabolism" & adjusted_p_value <= 1e-12,
    1e-12,
    adjusted_p_value
  )) %>%
  # Assign descriptive labels to various reactions.
  mutate(label = case_when(
    reaction_id == "PGM_neg" ~ "phosphoglycerate mutase (PGAM)",
    reaction_id == "LDH_L_neg" ~ "lactate dehydrogenase",
    reaction_id == "PDHm_pos" ~ "pyruvate dehydrogenase (PDH)",
    reaction_id == "TPI_neg" ~ "triosephosphate isomerase (DHAP forming)",
    reaction_id == "FACOAL1821_neg" ~ "long-chain fatty-acid-CoA ligase",
    reaction_id == "r1257_pos" ~ "long-chain fatty-acid-CoA ligase",
    reaction_id == "FACOAL1831_neg" ~ "long-chain fatty-acid-CoA ligase",
    reaction_id == "CSNATr_neg" ~ "carnitine O-acetyltransferase",
    reaction_id == "C160CPT1_pos" ~ "carnitine O-palmitoyltransferase",
    reaction_id == "ACONTm_pos" ~ "aconitate hydratase",
    reaction_id == "SUCOASm_pos" ~ "succinate-CoA ligase",
    reaction_id == "AKGDm_pos" ~ "alpha-ketoglutarate dehydrogenase",
    reaction_id == "SUCD1m_pos" ~ "succinate dehydrogenase",
    reaction_id == "ICDHyrm_pos" ~ "isocitrate dehydrogenase",
    reaction_id == "CK_pos" ~ "creatine\nkinase",
    reaction_id == "PGCD_pos" ~ "phosphoglycerate dehydrogenase",
    reaction_id == "ARGSS_pos" ~ "arginosuccinate synthase",
    reaction_id == "r0281_neg" ~ "putrescine diamine oxidase",
    reaction_id == "SPMDOX_pos" ~ "spermidine dehydrogenase (spermidine -> GABA)",
    reaction_id == "ARGDCm_pos" ~ "arginine decarboxylase",
    reaction_id == "AGMTm_pos" ~ "agmatinase",
    reaction_id == "GHMT2r_pos" ~ "serine hydroxymethyltransferase",
    reaction_id == "APPMS_pos" ~ "adenosylhomocysteinase",
    reaction_id == "METAT_pos" ~ "methionine adenosyltransferase",
    reaction_id == "METS_pos" ~ "methionine\nsynthase",
    reaction_id == "ARGN_pos" ~ "arginase",
    TRUE ~ ""
  ))

ggplot(
  compass_scores_by_cell_type,
  aes(
    x = cohens_d,
    y = -log10(adjusted_p_value),
    color = subsystem_priority
  )
) +
  ggtitle("Differential COMPASS Scores for PPMS vs. HC Cells") +
  xlab("Cohen's d") + ylab("-log(BH-adjusted p-value)") +
  xlim(-0.5, 0.5) +
  facet_wrap(vars(subsystem_priority), scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(
    "Glycolysis" = "#662D8C",
    "TCA cycle" = "#B87013",
    "Fatty acid oxidation" = "#0B0D9D",
    "Amino acid metabolism" = "#B82130"
  )) +
  guides(color = FALSE) +
  geom_point(size = 1, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype="dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype="dashed", color = "blue") +
  geom_text_repel(
    aes(label = label),
    min.segment.length = 0.1,
    point.padding = 0.5,
    size = 2,
    seed = 7,
    max.overlaps = Inf
  ) +
  theme_bw()
ggsave("plots/mo-DC_nonprog_vs_RR_Compass_volcanos.pdf", device = "pdf", 
       width = 15, height = 10, units = "in", dpi = 600)

cohens_d_by_subsystem <-
  wilcoxon_results %>%
  left_join(
    dplyr::select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"),
    by = "reaction_id"
  ) %>%
  left_join(
    compass_data$reaction_metadata,
    by = "reaction_no_direction"
  ) %>%
  # Keep only "confident reactions", as defined in our paper.
  filter(!is.na(EC_number)) %>%
  filter(confidence == "0" | confidence == "4") %>%
  # Keep only "interesting subsystems", as defined in our paper.
  filter(!(subsystem == "Miscellaneous" | subsystem == "Unassigned")) %>%
  filter(!(startsWith(subsystem, "Transport") | startsWith(subsystem, "Exchange"))) %>%
  # Keep only subsystems of non-negligible size.
  group_by(subsystem) %>%
  filter(n() > 5) %>%
  ungroup() %>%
  # Order subsystems in a manner that will lend itself to a visually aesthetic plot.
  mutate(
    subsystem_priority = factor(subsystem) %>%
      fct_reorder2(
        cohens_d,
        adjusted_p_value,
        .fun = function(cohens_d, adjusted_p_value) {
          abs(median(cohens_d[adjusted_p_value < 0.1]))
        },
        .desc = FALSE
      )
  )

ggplot(
  cohens_d_by_subsystem,
  aes(
    x = subsystem_priority,
    y = cohens_d,
    color = if_else(cohens_d > 0, "up_regulated", "down_regulated"),
    alpha = if_else(adjusted_p_value < 0.1, "significant", "insignificant")
  )
) +
  ggtitle("Up- and Down-Regulated Reactions Cross Pathway Boundaries") +
  xlab("") + ylab("Cohen's d") +
  scale_color_manual(
    values = c(up_regulated = "#ca0020", down_regulated = "#0571b0"),
    guide = FALSE
  ) +
  scale_alpha_manual(
    name = "",
    values = c(significant = 1, insignificant = 0.25),
    labels = c(significant = "BH-adjusted p-value < 0.1", insignificant = "insignificant")
  ) +
  coord_flip() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave("plots/mo-DC_nonprog_vs_RR_Compass_CrossPathways.pdf", device = "pdf", 
       width = 15, height = 10, units = "in", dpi = 600)


write.csv(wilcoxon_results, "post_compass_stats/mo-DC_nonprog_vs_RR_wilxocon_results.csv")
write.csv(compass_scores_by_cell_type, "post_compass_stats/mo-DC_nonprog_vs_RR_compass_scores_by_treatment.csv")
write.csv(cohens_d_by_subsystem, "post_compass_stats/mo-DC_nonprog_vs_RR_compass_scores_by_subsystem.csv")

