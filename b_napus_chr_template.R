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

# Define cytoband_colors for white fill
cytoband_colors <- c("p" = "white", "q" = "white")

# Plotting the unbinned data
kardy <- ggplot(karyotype_a_bna_unbinned) +
  geom_ideogram(aes(x = chr, ymin = start, ymax = end, 
                    chrom = chr, fill = arm, 
                    arm = arm), 
                radius = unit(9, 'pt'), width = 0.5, linewidth = 0.8, 
                colour = 'black', show.legend = FALSE) + # Hide legend since both arms are white
  scale_fill_manual(values = cytoband_colors) +
  scale_alpha(range = c(0.1, 1)) +
  theme_void()

# Display plot
print(kardy)

# Save the plot to a file
ggsave("karyotype_plot1_unbinned_white.png", plot = kardy, width = 10, height = 8, dpi = 300)
