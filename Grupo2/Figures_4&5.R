# Analysis of bacterial genomes evolved in the absence of phage
# Reena Debray
# Feb 1, 2022
# Modified: Amairani Cancino Bello (04/Sep/2024) 
# This script runs the ReproHack_data_figures4_5.RData variable to create the graphs.


#--- Packages ----
library(ggplot2)
library(viridis)
# Load objects from the .RData file
# The variable contains : 
# Saving objects needed for the figures from 2_Genome_data_analysis.R script.
# save(p_day6, p_day36, annotation_day6, test_gene_df_day6, annotation_day6, breseq_gene_props_day6, test_gene_df_day36, breseq_gene_props_day36, 
#    file = "ReproHack_data_figures4_5.RData")
getwd()
load("ReproHack_data_figures4_5.RData")

# Figure 4A y 4B
#--- Fig. 4A ---- 
# Creation and display of the graph
p_day6 <- ggplot(test_gene_df_day6, aes(label, sim_coef)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
  theme_classic(base_size = 18) +
  xlab("") +
  ylab("Overlap in acquired mutations\n(Pairwise Sørensen-Dice similarity)") +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 18, face = "bold")) +
  ylim(0.1, 0.55) +
  annotate("text", x = 1.5, y = 0.55, label = annotation_day6, size = 6, hjust = 0.5, vjust = 1.5, color = "black")

print(p_day6)

#--- Plot Fig 4B---- 
ggplot(breseq_gene_props_day6, aes(top_gene, res_gene, fill = percent_pops * 100)) +
  geom_tile(color = "grey80") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")  # Define all adjustments in one call
  ) +
  scale_fill_viridis(
    name = "Populations\nwith mutation (%)",
    limits = c(0, 100),  # Define scale limits
    breaks = c(0, 25, 50, 75, 100)  # Define specific breaks to display
  ) +
  ylab("") +
  xlab("") +
  ggtitle("Day 6 of experimental evolution")


#--- Fig. 5A----
p_day36 <- ggplot(test_gene_df_day36, aes(label, sim_coef)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
  theme_classic(base_size = 18) +
  xlab("") +
  ylab("Overlap in acquired mutations\n(Paired Sørensen-Dice similarity)") +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 18, face = "bold")) +
  ylim(0.1, 0.55)

print(p_day36)

#--- Fig. 5B----
ggplot(breseq_gene_props_day36, aes(top_gene, res_gene, fill = percent_pops * 100)) +
  geom_tile(color = "grey80") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")) +
  scale_fill_viridis(name = "Populations\nwith mutation (%)", limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  ylab("") +
  xlab("") +
  ggtitle("Day 36 of experimental evolution")

