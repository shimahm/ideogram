require(ggideogram)
library(rlang)
library(ggplot2)
library(dplyr)

# Load data
karyotype_a_bna <- read.table("updated_data_with_arms.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Check original data
cat("Original unique chromosomes:", unique(karyotype_a_bna$chr), "\n")
cat("Original rows per chromosome and arm:\n")
print(table(karyotype_a_bna$chr, karyotype_a_bna$arm))

# Validate data for missing values
cat("Any missing values in key columns:\n")
cat("chr:", any(is.na(karyotype_a_bna$chr)), "\n")
cat("arm:", any(is.na(karyotype_a_bna$arm)), "\n")
cat("start:", any(is.na(karyotype_a_bna$start)), "\n")
cat("end:", any(is.na(karyotype_a_bna$end)), "\n")

# Aggregate data to remove binning (one row per chr + arm)
karyotype_a_bna_unbinned <- karyotype_a_bna %>%
  group_by(chr, arm) %>%
  summarize(
    start = min(start, na.rm = TRUE),
    end = max(end, na.rm = TRUE),
    .groups = "drop"
  )

# Check aggregated data
cat("Unbinned unique chromosomes:", unique(karyotype_a_bna_unbinned$chr), "\n")
cat("Unbinned rows per chromosome and arm:\n")
print(table(karyotype_a_bna_unbinned$chr, karyotype_a_bna_unbinned$arm))

# Split into A and C chromosomes
karyotype_A <- karyotype_a_bna_unbinned %>% filter(grepl("^A", chr))
karyotype_C <- karyotype_a_bna_unbinned %>% filter(grepl("^C", chr))

# Check split data
cat("A chromosomes:", unique(karyotype_A$chr), "\n")
cat("C chromosomes:", unique(karyotype_C$chr), "\n")

# Define cytoband_colors for white fill
cytoband_colors <- c("p" = "white", "q" = "white")

# Plot for A chromosomes
kardy_A <- ggplot(karyotype_A) +
  # First copy (left)
  geom_ideogram(aes(x = chr, ymin = start, ymax = end, 
                    chrom = chr, fill = arm, arm = arm), 
                radius = unit(9, 'pt'), width = 0.4, linewidth = 0.8, 
                colour = 'black', show.legend = FALSE, just = -0.05) +
  # Second copy (right)
  geom_ideogram(aes(x = chr, ymin = start, ymax = end, 
                    chrom = chr, fill = arm, arm = arm), 
                radius = unit(9, 'pt'), width = 0.4, linewidth = 0.8, 
                colour = 'black', show.legend = FALSE, just = 1.05) +
  scale_fill_manual(values = cytoband_colors) +
  scale_alpha(range = c(0.1, 1)) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for C chromosomes
kardy_C <- ggplot(karyotype_C) +
  # First copy (left)
  geom_ideogram(aes(x = chr, ymin = start, ymax = end, 
                    chrom = chr, fill = arm, arm = arm), 
                radius = unit(9, 'pt'), width = 0.4, linewidth = 0.8, 
                colour = 'black', show.legend = FALSE, just = -0.05) +
  # Second copy (right)
  geom_ideogram(aes(x = chr, ymin = start, ymax = end, 
                    chrom = chr, fill = arm, arm = arm), 
                radius = unit(9, 'pt'), width = 0.4, linewidth = 0.8, 
                colour = 'black', show.legend = FALSE, just = 1.05) +
  scale_fill_manual(values = cytoband_colors) +
  scale_alpha(range = c(0.1, 1)) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
print(kardy_A)
print(kardy_C)

# Save plots to files
ggsave("karyotype_A_diploid.png", plot = kardy_A, width = 12, height = 8, dpi = 300)
ggsave("karyotype_C_diploid.png", plot = kardy_C, width = 12, height = 8, dpi = 300)