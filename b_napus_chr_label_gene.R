require(ggideogram)
library(rlang)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load karyotype data
karyotype_a_bna <- read.table("updated_data_with_arms.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Check original karyotype data
cat("Original unique chromosomes:", unique(karyotype_a_bna$chr), "\n")
cat("Original rows per chromosome and arm:\n")
print(table(karyotype_a_bna$chr, karyotype_a_bna$arm))

# Validate karyotype data for missing values
cat("Any missing values in key columns (karyotype):\n")
cat("chr:", any(is.na(karyotype_a_bna$chr)), "\n")
cat("arm:", any(is.na(karyotype_a_bna$arm)), "\n")
cat("start:", any(is.na(karyotype_a_bna$start)), "\n")
cat("end:", any(is.na(karyotype_a_bna$end)), "\n")

# Aggregate karyotype data to remove binning (one row per chr + arm)
karyotype_a_bna_unbinned <- karyotype_a_bna %>%
  group_by(chr, arm) %>%
  summarize(
    start = min(start, na.rm = TRUE),
    end = max(end, na.rm = TRUE),
    .groups = "drop"
  )

# Check aggregated karyotype data
cat("Unbinned unique chromosomes:", unique(karyotype_a_bna_unbinned$chr), "\n")
cat("Unbinned rows per chromosome and arm:\n")
print(table(karyotype_a_bna_unbinned$chr, karyotype_a_bna_unbinned$arm))

# Load gene data (no header, take only first four columns)
genes <- read.table("annotated_genes_plot.txt", sep = "\t", header = FALSE, 
                    stringsAsFactors = FALSE, 
                    col.names = c("chr", "start", "end", "gene_name", 
                                  paste0("extra_", seq_len(10))), 
                    fill = TRUE, 
                    colClasses = c("character", "numeric", "numeric", "character", 
                                   rep("NULL", 10))) %>% 
  mutate(Mid = (start + end) / 2)

# Check gene data
cat("Unique chromosomes in gene data:", unique(genes$chr), "\n")
cat("Any missing values in gene data:\n")
cat("chr:", any(is.na(genes$chr)), "\n")
cat("start:", any(is.na(genes$start)), "\n")
cat("end:", any(is.na(genes$end)), "\n")
cat("gene_name:", any(is.na(genes$gene_name)), "\n")

# Validate that gene chromosomes match karyotype chromosomes
cat("Gene chromosomes not in karyotype:", 
    setdiff(unique(genes$chr), unique(karyotype_a_bna_unbinned$chr)), "\n")

# Define cytoband_colors for yellow fill
cytoband_colors <- c("p" = "#e8e87d", "q" = "#e8e87d")

# Plotting the unbinned data with genes
kardy <- ggplot(karyotype_a_bna_unbinned) +
  # Chromosomes
  geom_ideogram(aes(x = 1, ymin = start, ymax = end, 
                    chrom = chr, fill = arm, arm = arm), 
                radius = unit(4, 'pt'), width = 0.1, 
                linewidth = 0.2, colour = 'black', 
                show.legend = FALSE) +
  scale_fill_manual(values = cytoband_colors) +
  # Genes as segments
  geom_segment(data = genes, 
               aes(x = 1, xend = 1, y = start, yend = end), 
               color = "red", size = 0.5, 
               position = position_nudge(x = 0.2)) +
  # Gene labels
  geom_text_repel(data = genes, 
                  aes(x = 1, y = Mid, label = gene_name), 
                  size = 2, hjust = 0, nudge_x = 0.5, 
                  color = "black", 
                  arrow = arrow(length = unit(0.015, "npc")),
                  direction = "y",
                  force_pull = 0,
                  segment.size = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  max.overlaps = Inf, 
                  box.padding = 0.5, 
                  point.padding = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), 
                     limits = c(0.9, 3)) +
  facet_wrap(vars(chr), scales = 'free_y') +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank()
  )

# Display plot
print(kardy)

# Save the plot to a file
ggsave("karyotype_plot1_unbinned_yellow_genes.png", plot = kardy, width = 12, height = 10, dpi = 300)