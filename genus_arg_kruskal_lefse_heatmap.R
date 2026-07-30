# ==============================================================================
# Genus-ARG Kruskal-Wallis Testing, LEfSe, and Heatmap
# Peshawar Wastewater Metagenomics -- underlies Figure 4A and the genus-ARG
# significance testing reported in the Results (Kruskal-Wallis / FDR section).
#
# Input:
#   CARD_abundance_mixed.csv  Genus-stratified CARD abundance table, columns:
#                             CARD (format "ARO:ID|k__...;g__Genus.s__..."),
#                             PHW1, PHW2, PCW1, PCW2, PSLW1, PSLW2
#
# Output:
#   Genus_ARG_mean_matrix_all.csv-equivalent summary tables (printed), a
#   Kruskal-Wallis/FDR table, and a heatmap of the top 30 significant
#   genus-ARG pairs.
#
# Required packages:
#   install.packages(c("tidyverse", "reshape2", "broom", "pheatmap"))
#   # lefser is a Bioconductor package:
#   # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   # BiocManager::install("lefser")
# ==============================================================================

library(tidyverse)
library(reshape2)
library(broom)
library(pheatmap)
library(lefser)

CARD_abundance_mixed <- read.csv("data/CARD_abundance_mixed.csv", check.names = FALSE)

# ------------------------------------------------------------------------------
# Part 1 -- Preliminary exploration: a handful of named clinical pathogens
# ------------------------------------------------------------------------------
stratified_data <- CARD_abundance_mixed %>%
  separate(CARD, into = c("ARO", "Taxonomy"), sep = "\\|") %>%
  filter(Taxonomy != "unclassified")

target_pathogens <- stratified_data %>%
  filter(grepl("Klebsiella|Acinetobacter|Escherichia|Salmonella|Citrobacter", Taxonomy)) %>%
  mutate(Pathogen = str_extract(Taxonomy, "g__[A-Za-z]+")) %>%
  mutate(Pathogen = str_replace(Pathogen, "g__", ""))

groups <- list(
  Hospital = c("PHW1", "PHW2"),
  Community = c("PCW1", "PCW2"),
  Slaughterhouse = c("PSLW1", "PSLW2")
)

sector_summary <- target_pathogens %>%
  pivot_longer(cols = starts_with("P"), names_to = "Sample", values_to = "Abundance") %>%
  mutate(Sector = case_when(
    Sample %in% groups$Hospital ~ "Hospital",
    Sample %in% groups$Community ~ "Community",
    Sample %in% groups$Slaughterhouse ~ "Slaughterhouse"
  )) %>%
  group_by(Sector, Pathogen, ARO) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop") %>%
  filter(Mean_Abundance > 0)

stats_results <- sector_summary %>%
  group_by(Pathogen) %>%
  do(tidy(kruskal.test(Mean_Abundance ~ Sector, data = .)))

ggplot(sector_summary, aes(x = Pathogen, y = Mean_Abundance, fill = Sector)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean Abundance of Pathogen-Associated ARGs by Sector",
       y = "Mean Abundance (CPC)", x = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------------------------------------------------------
# Part 2 -- Full genus-ARO analysis: Kruskal-Wallis + FDR, LEfSe, heatmap
# ------------------------------------------------------------------------------
peshawar_cols <- c("CARD", "PHW1", "PHW2", "PCW1", "PCW2", "PSLW1", "PSLW2")
CARD_peshawar <- CARD_abundance_mixed %>% select(all_of(peshawar_cols))

tax_data <- CARD_peshawar %>%
  separate(CARD, into = c("ARO", "Taxonomy"), sep = "\\|", extra = "merge", fill = "right") %>%
  filter(!is.na(Taxonomy), !grepl("unclassified", Taxonomy, ignore.case = TRUE))

extract_genus <- function(tax) {
  if (is.na(tax)) return(NA)
  match <- str_extract(tax, "g__[A-Za-z]+")
  if (is.na(match)) return(NA)
  str_replace(match, "g__", "")
}
tax_data$Genus <- sapply(tax_data$Taxonomy, extract_genus)
tax_data <- tax_data %>% filter(!is.na(Genus))

long_data <- tax_data %>%
  pivot_longer(cols = c("PHW1", "PHW2", "PCW1", "PCW2", "PSLW1", "PSLW2"),
               names_to = "Sample", values_to = "Abundance") %>%
  mutate(Sector = case_when(
    Sample %in% c("PHW1", "PHW2") ~ "Hospital",
    Sample %in% c("PCW1", "PCW2") ~ "Community",
    Sample %in% c("PSLW1", "PSLW2") ~ "Slaughterhouse"
  ))

sector_means <- long_data %>%
  group_by(Sector, Genus, ARO) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(Mean_Abundance > 0)

test_data <- long_data %>%
  group_by(Genus, ARO, Sample, Sector) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

kruskal_results <- test_data %>%
  group_by(Genus, ARO) %>%
  summarise(
    p_value = ifelse(n_distinct(Sector) >= 3 & n() >= 3,
                     kruskal.test(Abundance ~ Sector, data = pick(everything()))$p.value,
                     NA),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"),
         significant = ifelse(p_adjusted < 0.05, "Yes", "No"))

cat("Total Genus-ARO pairs tested:", nrow(kruskal_results), "\n")
cat("Significant after FDR (p_adj < 0.05):", sum(kruskal_results$significant == "Yes"), "\n")
print(head(kruskal_results %>% arrange(p_adjusted), 20))

# Prepare a feature x sample matrix and run LEfSe (genus|ARO as feature)
long_data$Feature <- paste0(long_data$Genus, "|", long_data$ARO)
lefse_mat <- long_data %>%
  select(Feature, Sample, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)

metadata <- data.frame(
  sample = names(lefse_mat)[-1],
  Sector = case_when(
    names(lefse_mat)[-1] %in% c("PHW1", "PHW2") ~ "Hospital",
    names(lefse_mat)[-1] %in% c("PCW1", "PCW2") ~ "Community",
    names(lefse_mat)[-1] %in% c("PSLW1", "PSLW2") ~ "Slaughterhouse"
  )
)

lefse_mat_matrix <- as.matrix(lefse_mat[, -1])
rownames(lefse_mat_matrix) <- lefse_mat$Feature

lefse_result <- lefser(lefse_mat_matrix, groups = metadata$Sector, cols = metadata$Sector)
print(head(lefse_result[order(lefse_result$lda, decreasing = TRUE), ], 10))

# Heatmap of the top 30 most significant genus-ARO pairs
top_pairs <- kruskal_results %>%
  arrange(p_adjusted) %>%
  head(30) %>%
  mutate(Feature = paste0(Genus, "|", ARO))

heatmap_data <- sector_means %>%
  filter(paste0(Genus, "|", ARO) %in% top_pairs$Feature) %>%
  select(Sector, Feature, Mean_Abundance) %>%
  pivot_wider(names_from = Sector, values_from = Mean_Abundance, values_fill = 0)

heatmap_mat <- as.matrix(heatmap_data[, c("Hospital", "Community", "Slaughterhouse")])
rownames(heatmap_mat) <- heatmap_data$Feature
heatmap_mat_log <- log1p(heatmap_mat)

pheatmap(heatmap_mat_log,
         main = "Top 30 significant Genus-ARG pairs (FDR < 0.05)",
         cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         fontsize_row = 8)
