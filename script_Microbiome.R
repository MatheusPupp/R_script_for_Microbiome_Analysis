# LITTER AND WATER SAMPLE ANALYSIS SCRIPT
setwd("C:\\Users\\mathe\\Documents")

library(qiime2R)
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(vegan)
library(picante)
library(ALDEx2)
library(zCompositions)
library(MASS)
library(metagenomeSeq)
library(dendextend)
library(selbal)
library(HPM)
library(rms)
library(breakaway)
library(dplyr)
library(microViz)
library(ggplot2)
library(ANCOMBC)
library(FSA)

# boxlink!
#"/Users/lacona/Library/CloudStorage/Box-Box/"

# Creating the original phyloseq object 
# metadata has to be .tsv

Jolie_Analysis<-qza_to_phyloseq(
  features="table-Jolie3.qza",
  tree="rooted-tree-Jolie3.qza",
  taxonomy="taxonomy-Jolie3.qza",
  metadata = "Jolie_metadata_Andrea_cleaned.tsv"
)

###################
# FILTERING AND SORTING

# Sort the samples by number of reads
sort(phyloseq::sample_sums(Jolie_Analysis))
     
# Filter out everything not related to bacteria

psJ <- Jolie_Analysis %>%
  subset_taxa(Kingdom == "d__Bacteria" &  Family  != "Mitochondria" &             #filter out mitochondria
                Class   != "Chloroplast")

# Filter out everything with less than 5000 reads and empty otus
psJ <- prune_samples(sample_sums(psJ) > 5000, psJ)

# Get unique sampling days
unique_days <- unique(sample_data(psJ)$Day)
print(unique_days)
# Assuming your sample metadata has a column named "Day"
physeq_day20 <- subset_samples(psJ, Day == "20")
physeq_day28 <- subset_samples(psJ, Day == "28")

# Get the total counts per taxa
taxa_counts <- taxa_sums(physeq_day20)
taxa_counts <- taxa_sums(physeq_day28)

# Identify taxa with structural zeros across groups
taxa_with_zeros <- taxa_counts[taxa_counts == 0]
print(taxa_with_zeros)

# arange the treatments for day 20
physeq_day20@sam_data$Water.Treatment <- factor(physeq_day20@sam_data$Water.Treatment, 
                                    levels=c("Nothing", "Acetic", "Citric"))

# ANCOMBC STEP IF analyzing day 20 and 28 separately
# Day 20 analysis (only water treatments)
result_day20 <- ancombc(
  data = physeq_day20,
  assay.type = NULL,
  assay_name = "counts",
  rank = NULL,
  tax_level = "Family",
  formula = "Water.Treatment",
  p_adj_method = "holm",
  prv_cut = 0.1,
  lib_cut = 0,
  group = "Water.Treatment",
  struc_zero = FALSE,
  neg_lb = FALSE,
  tol = 1e-05,
  max_iter = 100,
  conserve = FALSE,
  alpha = 0.05,
  global = FALSE,
  n_cl = 1,
  verbose = TRUE
)

# arange the treatments for day 28
physeq_day28@sam_data$Water.Treatment <- factor(physeq_day28@sam_data$Water.Treatment, 
                                                levels=c("Nothing", "Acetic", "Citric"))
physeq_day28@sam_data$Litter.Treatment <- factor(physeq_day28@sam_data$Litter.Treatment, 
                                                levels=c("Nothing", "Compost", "FormicAcidSalt"))

# Day 28 analysis (both water and litter treatments)
result_day28 <- ancombc(
  data = physeq_day28,
  assay.type = NULL,
  assay_name = "counts",
  rank = NULL,
  tax_level = "Family",
  formula = "Water.Treatment + Litter.Treatment",
  p_adj_method = "holm",
  prv_cut = 0.1,
  lib_cut = 0,
  group = "Water.Treatment",
  struc_zero = FALSE,
  neg_lb = FALSE,
  tol = 1e-05,
  max_iter = 100,
  conserve = FALSE,
  alpha = 0.05,
  global = FALSE,
  n_cl = 1,
  verbose = TRUE
)

# Name the columns for the p-values tables
col_name20 = c("taxon", "(Intercept)", "Water.TreatmentAcetic", "Water.TreatmentCitric")
col_name28 = c("taxon", "(Intercept)", "Water.TreatmentAcetic", "Water.TreatmentCitric", "Litter.TreatmentCompost", "Litter.TreatmentFormicAcidSalt")

# ANCOMBC p-values
tab_p20 = result_day20$res$p_val
colnames(tab_p20) = col_name20
tab_p20 %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name20[-1], digits = 3)

tab_p28 = result_day28$res$p_val
colnames(tab_p28) = col_name28
tab_p28 %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name28[-1], digits = 3)

# ANCOMBC Adjusted p-values (q_val = adjusted p-values)
tab_q20 = result_day20$res$q_val
colnames(tab_q20) = col_name20
tab_q20 %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name20[-1], digits = 3)
tab_q_day20_filtered <- tab_q20 %>%
  filter(`(Intercept)` < 0.05 | `Water.TreatmentAcetic` < 0.05 | `Water.TreatmentCitric` < 0.05)

tab_q28 = result_day28$res$q_val
colnames(tab_q28) = col_name28
tab_q28 %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name28[-1], digits = 3)
tab_q_day28_filtered <- tab_q28 %>%
  filter(`(Intercept)` < 0.05 | `Water.TreatmentAcetic` < 0.05 | `Water.TreatmentCitric` < 0.05 | `Litter.TreatmentCompost` < 0.05 | `Litter.TreatmentFormicAcidSalt` < 0.05)

# Filtering only lines with p-values < 0.05 for result_day20
tab_q_day20_filtered <- tab_q20 %>%
  filter(`(Intercept)` < 0.05 | `Water.TreatmentAcetic` < 0.05 | `Water.TreatmentCitric` < 0.05)
# Filtering only lines with p-values < 0.05 for result_day28
tab_q_day28_filtered <- tab_q28 %>%
  filter(`(Intercept)` < 0.05 | `Water.TreatmentAcetic` < 0.05 | `Water.TreatmentCitric` < 0.05 | 
           `Litter.TreatmentCompost` < 0.05 | `Litter.TreatmentFormicAcidSalt` < 0.05)

# Displaying filtered tables with taxon names
tab_q_day20_filtered %>%
  datatable(caption = "Adjusted P-values < 0.05 for Day 20")

tab_q_day28_filtered %>%
  datatable(caption = "Adjusted P-values < 0.05 for Day 28")

# Save the adjusted p-values table of day 20 in TSV format
write.table(tab_q_day20_filtered, file = "adjusted_p_values_day20.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Save the adjusted p-values table of day 28 in TSV format
write.table(tab_q_day28_filtered, file = "adjusted_p_values_day28.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Looking for Log Fold Changes (LFC)

# Day 20
# Extracting LFC and preparing table structure
df_lfc = data.frame(result_day20$res$lfc[, -1] * result_day20$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = result_day20$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Renaming columns with treatments
colnames(df_lfc) <- c("taxon_id", "Nothing", "Acetic", "Citric")

# Filtering to keep only Acetic and Citric treatments
df_fig_water = df_lfc %>%
  dplyr::select(taxon_id, Acetic, Citric) %>%
  tidyr::pivot_longer(cols = c("Acetic", "Citric"), names_to = "treatment", values_to = "LFC") %>%
  dplyr::filter(!is.na(LFC) & LFC != 0) %>%  # Removendo taxons com LFC igual a 0
  dplyr::mutate(direct = ifelse(LFC > 0, "Positive LFC", "Negative LFC"))

# Adjusting factors for the graph
df_fig_water$taxon_id = factor(df_fig_water$taxon_id, levels = unique(df_fig_water$taxon_id))
df_fig_water$direct = factor(df_fig_water$direct, levels = c("Positive LFC", "Negative LFC"))
df_fig_water$treatment = factor(df_fig_water$treatment, levels = c("Citric", "Acetic"))  # Citric primeiro

# Creating the graph with ggplot
p_water_d20 = ggplot(data = df_fig_water, aes(y = taxon_id, x = LFC, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  facet_wrap(~treatment, ncol = 2, scales = "free_y") + # Plots lado a lado
  labs(y = "Taxon", x = "Log Fold Change (LFC)", title = "Log Fold Changes for Acetic and Citric Treatments on day 20") + 
  scale_fill_discrete(name = "Direction") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size = 14), # Tamanho da fonte do título
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 10, color = "black"),  # Tamanho e cor dos taxons (letra preta)
        axis.text.x = element_text(size = 11, color = "black"),  # Ajustando o tamanho do eixo x e cor
        strip.text = element_text(size = 12, color = "black"),  # Tamanho e cor do título de cada gráfico
        strip.text.x = element_blank(), # Removendo as labels das facetas (taxon label)
        plot.margin = unit(c(1,1,1,1), "cm"), 
        aspect.ratio = 1.5,  # Ajustando a proporção para maior altura
        panel.spacing = unit(1, "lines"))  # Aumentando o espaçamento entre os gráficos

# Displaying the graph
p_water_d20

# DAY 28
# WATER 
# Extracting CFL and preparing the table structure for the Water Treatments on day 28
df_lfc_water_d28 = data.frame(result_day28$res$lfc[, -1] * result_day28$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = result_day28$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Renaming the columns with the water treatments
colnames(df_lfc_water_d28) <- c("taxon_id", "Water.TreatmentsAcetic", "Water.TreatmentsCitric")

# Filtering to keep only Acetic and Citric treatments
df_fig_water_d28 = df_lfc_water_d28 %>%
  dplyr::select(taxon_id, Water.TreatmentsAcetic, Water.TreatmentsCitric) %>%
  tidyr::pivot_longer(cols = c("Water.TreatmentsAcetic", "Water.TreatmentsCitric"), names_to = "treatment", values_to = "LFC") %>%
  dplyr::filter(!is.na(LFC) & LFC != 0) %>%  # Removendo taxons com LFC igual a 0
  dplyr::mutate(direct = ifelse(LFC > 0, "Positive LFC", "Negative LFC"))

# Adjusting factors for the graph
df_fig_water_d28$taxon_id = factor(df_fig_water_d28$taxon_id, levels = unique(df_fig_water$taxon_id))
df_fig_water_d28$direct = factor(df_fig_water_d28$direct, levels = c("Positive LFC", "Negative LFC"))
df_fig_water_d28$treatment = factor(df_fig_water_d28$treatment, levels = c("Water.TreatmentsCitric", "Water.TreatmentsAcetic"))  # Citric primeiro

# Creating the Water Treatments graph with ggplot
p_water_d28 = ggplot(data = df_fig_water, aes(y = taxon_id, x = LFC, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  facet_wrap(~treatment, ncol = 2, scales = "free_y") + # Plots side by side
  labs(y = "Taxon", x = "Log Fold Change (LFC)", title = "Log Fold Changes for Citric and Acetic Water Treatments on day 28") + 
  scale_fill_discrete(name = "Direction") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size = 14), # Title font size
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 7, color = "black"),  # Size and color of taxons 
        axis.text.x = element_text(size = 8, color = "black"),  # Adjusting the x-axis size and color
        strip.text = element_text(size = 12, color = "black"),  # Size and color of the title of each graphic
        strip.text.x = element_blank(), # Removing facet labels (taxon label)
        plot.margin = unit(c(1,1,1,1), "cm"), 
        aspect.ratio = 1.5,  # Adjusting the aspect ratio for greater height
        panel.spacing = unit(1, "lines"))  # Increasing the spacing between graphics

# Displaying the Water Treatments graph
p_water_d28

# LITTER 
# Extracting LFC and preparing the table structure for the Litter Treatments on day 28
df_lfc_litter = data.frame(result_day28$res$lfc[, -1] * result_day28$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = result_day28$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Renaming columns with bed treatments
colnames(df_lfc_litter) <- c("taxon_id", "Litter.TreatmentsCompost", "Litter.TreatmentsFormicAcidSalt")

# Filtering to keep only Compost and Formic Acid Salt treatments
df_fig_litter = df_lfc_litter %>%
  dplyr::select(taxon_id, "Litter.TreatmentsCompost", "Litter.TreatmentsFormicAcidSalt") %>%
  tidyr::pivot_longer(cols = c("Litter.TreatmentsCompost", "Litter.TreatmentsFormicAcidSalt"), names_to = "treatment", values_to = "LFC") %>%
  dplyr::filter(!is.na(LFC) & LFC != 0) %>%  # Removing taxons with LFC equal to 0
  dplyr::mutate(direct = ifelse(LFC > 0, "Positive LFC", "Negative LFC"))

# Adjusting factors for the graph
df_fig_litter$taxon_id = factor(df_fig_litter$taxon_id, levels = unique(df_fig_litter$taxon_id))
df_fig_litter$direct = factor(df_fig_litter$direct, levels = c("Positive LFC", "Negative LFC"))
df_fig_litter$treatment = factor(df_fig_litter$treatment, levels = c("Litter.TreatmentsCompost", "Litter.TreatmentsFormicAcidSalt"))  # Compost first

# Creating the graph for the Litter Treatments with ggplot
p_litter_d28 = ggplot(data = df_fig_litter, aes(y = taxon_id, x = LFC, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  facet_wrap(~treatment, ncol = 2, scales = "free_y") + # Plots side by side
  labs(y = "Taxon", x = "Log Fold Change (LFC)", title = "Log Fold Changes for Compost and Formic Acid Salt Litter Treatments") + 
  scale_fill_discrete(name = "Direction") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size = 14), # Title font size
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 7, color = "black"),  # Taxons size and color
        axis.text.x = element_text(size = 8, color = "black"),  # Adjusting the x-axis size and color
        strip.text = element_text(size = 12, color = "black"),  # Size and color of the title of each graphic
        strip.text.x = element_blank(), # Removing facet labels (taxon label)
        plot.margin = unit(c(1,1,1,1), "cm"), 
        aspect.ratio = 1.5,  # Adjusting the aspect ratio for greater height
        panel.spacing = unit(1, "lines"))  # Increasing the spacing between graphics

# Displaying the Litter Treatments graphic
p_litter_d28

##########################
# ALPHA DIVERSITY

# DAY 20
#Generate a data.frame with adiv measures
adiv20 <- data.frame(
  "Observed" = phyloseq::estimate_richness(physeq_day20, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(physeq_day20, measures = "Shannon"),
  "Group" = phyloseq::sample_data(physeq_day20)$Water.Treatment)
#Generate Pielou
"Pielou" = adiv20$Pielou <- adiv20$Shannon / log(adiv20$Observed)
#Plot adiv measures
alpha_plot20 <- adiv20 %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Pielou")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Pielou"))) %>%
  ggplot(aes(group=Group, x = factor(Group, level=c("Nothing", "Acetic", "Citric")), y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Group), height = 0, width = .2, cex = 2.5) +
  scale_color_manual(values=c( "#53B74C", "#B3A033", "#E77D72",
                               "#55BCC2", "#6F9BF8",  "#E46EDD" )) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  ggtitle("Day 20") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(color = "black", size = 20, hjust=1, angle=45),
        axis.text.y = element_text(color = "black", size = 14),
        strip.text = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, color = "black"),
        panel.spacing.x = unit(2, "lines"))
ggsave("alpha diversity OSS 20.png", alpha_plot20, height = 6, width = 13, units = "in", dpi = 300 )

#Kruskal test of location
kruskal.test(Observed ~ Group, data=adiv20)
kruskal.test(Shannon ~ Group, data = adiv20)            
kruskal.test(Pielou ~ Group, data = adiv20)

# Pairwise comparisons (with Bonferroni correction)
dunnTest(Shannon ~ Group,
         data = adiv20,
         method= "bonferroni")

# Pairwise comparisons (with Bonferroni correction)
dunnTest(Observed ~ Group,
         data = adiv20,
         method= "bonferroni")

# DAY 28 WATER

#Generate a data.frame with adiv measures
adiv28w <- data.frame(
  "Observed" = phyloseq::estimate_richness(physeq_day28, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(physeq_day28, measures = "Shannon"),
  "Group" = phyloseq::sample_data(physeq_day28)$Water.Treatment)
#Generate Pielou
"Pielou" = adiv28w$Pielou <- adiv28w$Shannon / log(adiv28w$Observed)

#Plot adiv measures
alpha_plot28w <- adiv28w %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Pielou")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Pielou"))) %>%
  ggplot(aes(group=Group, x = factor(Group, level=c("Nothing", "Acetic", "Citric")), y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Group), height = 0, width = .2, cex = 2.5) +
  scale_color_manual(values=c( "#53B74C", "#B3A033", "#E77D72",
                               "#55BCC2", "#6F9BF8",  "#E46EDD" )) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  ggtitle("Day 28") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(color = "black", size = 20, hjust=1, angle=45),
        axis.text.y = element_text(color = "black", size = 14),
        strip.text = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, color = "black"),
        panel.spacing.x = unit(2, "lines"))
ggsave("alpha diversity OSS 28w.png", alpha_plot28w, height = 6, width = 13, units = "in", dpi = 300 )

#Kruskal test of location
kruskal.test(Observed ~ Group, data=adiv28w)
kruskal.test(Shannon ~ Group, data = adiv28w)            
kruskal.test(Pielou ~ Group, data = adiv28w)

# DAY 28 LITTER

#Generate a data.frame with adiv measures
adiv28L <- data.frame(
  "Observed" = phyloseq::estimate_richness(physeq_day28, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(physeq_day28, measures = "Shannon"),
  "Group" = phyloseq::sample_data(physeq_day28)$Litter.Treatment)
#Generate Pielou
"Pielou" = adiv28L$Pielou <- adiv28L$Shannon / log(adiv28L$Observed)

#Plot adiv measures
alpha_plot28L <- adiv28L %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Pielou")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Pielou"))) %>%
  ggplot(aes(group=Group, x = factor(Group, level=c("Nothing", "Compost", "Formic Acid Salt")), y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Group), height = 0, width = .2, cex = 2.5) +
  scale_color_manual(values=c( "#53B74C", "#B3A033", "#E77D72",
                               "#55BCC2", "#6F9BF8",  "#E46EDD" )) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  ggtitle("Day 28") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(color = "black", size = 20, hjust=1, angle=45),
        axis.text.y = element_text(color = "black", size = 14),
        strip.text = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, color = "black"),
        panel.spacing.x = unit(2, "lines"))
ggsave("alpha diversity OSS 28L.png", alpha_plot28L, height = 6, width = 13, units = "in", dpi = 300 )


#Kruskal test of location
kruskal.test(Observed ~ Group, data=adiv28L)
kruskal.test(Shannon ~ Group, data = adiv28L)            
kruskal.test(Pielou ~ Group, data = adiv28L)


########################

# BETA DIVERSITY

# DAY 20
ps20_clr <- microbiome::transform(physeq_day20, "clr")

#PCA via phyloseq
ord_clr20 <- phyloseq::ordinate(ps20_clr, "RDA")

#Plot scree plot
pcoa20 <-phyloseq::plot_scree(ord_clr20) + 
  geom_bar(stat="identity", fill = "royalblue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle("Day 20") +
  theme_bw() 
theme(axis.text.x = element_text(size = 12, angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size = 12))
ggsave("pcoa 20.png", pcoa20, height = 8, width = 13, units = "in", dpi = 300 )

#Scale axes and plot ordination with first two PCs
clr1 <- ord_clr20$CA$eig[1] / sum(ord_clr20$CA$eig)
clr2 <- ord_clr20$CA$eig[2] / sum(ord_clr20$CA$eig)


ord_plot20 <- phyloseq::plot_ordination(ps20_clr, ord_clr20, type="samples", color="Water.Treatment") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Water.Treatment), linetype = 2) +
  ggtitle("Day 20") +
  theme_bw()
ggsave("ordination plot 20.png", ord_plot20, height = 8, width = 10, units = "in", dpi = 300 )

#Extract the OTU table:
otu_table20 <- otu_table(physeq_day20)

#Convert the OTU table to a matrix:
otu_matrix20 <- as.matrix(otu_table20)

# Plot the ordination
plot_ordination(physeq_day20, ord_clr20, type = "samples", color = "Water.Treatment") +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Water.Treatment), linetype = 2) +
  ggtitle("Day 20") +
  theme_bw()

#Perform PERMANOVA with Bray-Curtis dissimilarity
bray_curtis_distance20 <- phyloseq::distance(physeq_day20, method ="bray") 
adonis2_resultbraycurtis20 <- adonis2(bray_curtis_distance20 ~ phyloseq::sample_data(ps20_clr)$Water.Treatment)
pairwise.adonis(bray_curtis_distance20, phyloseq::sample_data(physeq_day20)$Water.Treatment)


# Perform PERMANOVA with Unweighted UniFrac distance
unweighted_unifrac20 <- UniFrac(physeq=physeq_day20, weighted=FALSE, normalized=T, parallel=T, fast=T)
adonis2_resultunweighted20 <- adonis2(unweighted_unifrac20 ~ phyloseq::sample_data(ps20_clr)$Water.Treatment)
pairwise.adonis(unweighted_unifrac20, phyloseq::sample_data(physeq_day20)$Water.Treatment)


# Perform PERMANOVA with weighted UniFrac distance
weighted_unifrac20 = UniFrac(physeq=physeq_day20, weighted=TRUE, normalized=T, parallel=T, fast=T)
adonis2_resultweighted20 <- adonis2(weighted_unifrac20 ~ phyloseq::sample_data(ps20_clr)$Water.Treatment)
pairwise.adonis(weighted_unifrac20, phyloseq::sample_data(physeq_day20)$Water.Treatment)

# DAY 28 WATER
ps28w_clr <- microbiome::transform(physeq_day28, "clr")

#PCA via phyloseq
ord_clr28w <- phyloseq::ordinate(ps28w_clr, "RDA")

#Plot scree plot
pcoa28w <-phyloseq::plot_scree(ord_clr28w) + 
  geom_bar(stat="identity", fill = "royalblue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle("Day 28") +
  theme_bw() 
theme(axis.text.x = element_text(size = 12, angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size = 12))
ggsave("pcoa 28w.png", pcoa20, height = 8, width = 13, units = "in", dpi = 300 )

#Scale axes and plot ordination with first two PCs
clr128w <- ord_clr28w$CA$eig[1] / sum(ord_clr28w$CA$eig)
clr228w <- ord_clr28w$CA$eig[2] / sum(ord_clr28w$CA$eig)


ord_plot28w <- phyloseq::plot_ordination(ps28w_clr, ord_clr28w, type="samples", color="Water.Treatment") + 
  geom_point(size = 2) +
  coord_fixed(clr228w / clr128w) +
  stat_ellipse(aes(group = Water.Treatment), linetype = 2) +
  ggtitle("Day 28") +
  theme_bw()
ggsave("ordination plot 28w.png", ord_plot28w, height = 8, width = 10, units = "in", dpi = 300 )

#Extract the OTU table:
otu_table28w <- otu_table(physeq_day28)

#Convert the OTU table to a matrix:
otu_matrix28w <- as.matrix(otu_table28w)

# Plot the ordination
plot_ordination(physeq_day28, ord_clr28w, type = "samples", color = "Water.Treatment") +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Water.Treatment), linetype = 2) +
  ggtitle("Day 28") +
  theme_bw()

#Perform PERMANOVA with Bray-Curtis dissimilarity
bray_curtis_distance28w <- phyloseq::distance(physeq_day28, method= "bray") 
adonis2_resultbraycurtis28w <- adonis2(bray_curtis_distance28w ~ phyloseq::sample_data(ps28w_clr)$Water.Treatment)
pairwise.adonis(bray_curtis_distance28w, phyloseq::sample_data(physeq_day28)$Water.Treatment)

# Perform PERMANOVA with Unweighted UniFrac distance
unweighted_unifrac28w <- UniFrac(physeq=physeq_day28, weighted=FALSE, normalized=T, parallel=T, fast=T)
adonis2_resultunweighted28w <- adonis2(unweighted_unifrac28w ~ phyloseq::sample_data(ps28w_clr)$Water.Treatment)
pairwise.adonis(unweighted_unifrac28w, phyloseq::sample_data(physeq_day28)$Water.Treatment)


# Perform PERMANOVA with weighted UniFrac distance
weighted_unifrac28w = UniFrac(physeq=physeq_day28, weighted=TRUE, normalized=T, parallel=T, fast=T)
adonis2_resultweighted28w <- adonis2(weighted_unifrac28w ~ phyloseq::sample_data(ps28w_clr)$Water.Treatment)
pairwise.adonis(weighted_unifrac28w, phyloseq::sample_data(physeq_day28)$Water.Treatment)

# DAY 28 LITTER
ps28L_clr <- microbiome::transform(physeq_day28, "clr")

#PCA via phyloseq
ord_clr28L <- phyloseq::ordinate(ps28L_clr, "RDA")

#Plot scree plot
pcoa28L <-phyloseq::plot_scree(ord_clr28L) + 
  geom_bar(stat="identity", fill = "royalblue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle("Day 28") +
  theme_bw() 
theme(axis.text.x = element_text(size = 12, angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size = 12))
ggsave("pcoa 28L.png", pcoa20, height = 8, width = 13, units = "in", dpi = 300 )

#Scale axes and plot ordination with first two PCs
clr128L <- ord_clr28L$CA$eig[1] / sum(ord_clr28L$CA$eig)
clr228L <- ord_clr28L$CA$eig[2] / sum(ord_clr28L$CA$eig)


ord_plot28L <- phyloseq::plot_ordination(ps28L_clr, ord_clr28L, type="samples", color="Litter.Treatment") + 
  geom_point(size = 2) +
  coord_fixed(clr228L / clr128L) +
  stat_ellipse(aes(group = Litter.Treatment), linetype = 2) +
  ggtitle("Day 28") +
  theme_bw()
ggsave("ordination plot 28L.png", ord_plot28w, height = 8, width = 10, units = "in", dpi = 300 )

#Extract the OTU table:
otu_table28L <- otu_table(physeq_day28)

#Convert the OTU table to a matrix:
otu_matrix28L <- as.matrix(otu_table28L)

# Plot the ordination
plot_ordination(physeq_day28, ord_clr28L, type = "samples", color = "Litter.Treatment") +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Litter.Treatment), linetype = 2) +
  ggtitle("Day 28") +
  theme_bw()

#Perform PERMANOVA with Bray-Curtis dissimilarity
bray_curtis_distance28L <- phyloseq::distance(physeq_day28, method ="bray") 
adonis2_resultbraycurtis28L <- adonis2(bray_curtis_distance28L ~ phyloseq::sample_data(ps28L_clr)$Litter.Treatment)
pairwise.adonis(bray_curtis_distance28L, phyloseq::sample_data(physeq_day28)$Litter.Treatment)

# Perform PERMANOVA with Unweighted UniFrac distance
unweighted_unifrac28L <- UniFrac(physeq=physeq_day28, weighted=FALSE, normalized=T, parallel=T, fast=T)
adonis2_resultunweighted28L <- adonis2(unweighted_unifrac28L ~ phyloseq::sample_data(ps28L_clr)$Litter.Treatment)
pairwise.adonis(unweighted_unifrac28L, phyloseq::sample_data(physeq_day28)$Litter.Treatment)

# Perform PERMANOVA with weighted UniFrac distance
weighted_unifrac28L = UniFrac(physeq=physeq_day28, weighted=TRUE, normalized=T, parallel=T, fast=T)
adonis2_resultweighted28L <- adonis2(weighted_unifrac28L ~ phyloseq::sample_data(ps28L_clr)$Litter.Treatment)
pairwise.adonis(weighted_unifrac28L, phyloseq::sample_data(physeq_day28)$Litter.Treatment)

#################################################
