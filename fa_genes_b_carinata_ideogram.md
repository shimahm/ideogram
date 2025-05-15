# BLAST Gene Analysis Pipeline (15 May 2025)

## Step 1: Create the BLAST Database

Run the following command to generate a BLAST nucleotide database from your reference genome:

```bash
makeblastdb -in /home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/GCA_040584065.1_ASM4058406v1_genomic.fna -dbtype nucl -out Brassica_car_DB
```

---

## Step 2: Run BLAST — Separate Result Files per Gene

```bash
for gene in /home/mahmoudi/Bigdata/computing/Shima/gene_sequences/*.fa; do
  blastn \
    -query "$gene" \
    -db Brassica_car_DB \
    -out "/home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/result/$(basename "${gene%.fa}")_blast_results.txt" \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
done
```

---

## Step 2: Run BLAST — Combine Results into One File

```bash
#!/bin/bash

# Set paths
gene_dir="/home/mahmoudi/Bigdata/computing/Shima/gene_sequences/"
blast_db="Brassica_car_DB"
combined_out="/home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/all_blast_results_combined.txt"

# Create or empty the combined result file
echo -e "gene_name\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "$combined_out"

# Loop through each gene fasta file
for gene in "$gene_dir"*.fa; do
  gene_name=$(basename "${gene%.fa}")
  out_file="${gene%.fa}_blast_results.txt"
  
  # Run BLAST and save individual result
  blastn \
    -query "$gene" \
    -db "$blast_db" \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$out_file"

  # Add gene name as first column and append to combined file
  awk -v gn="$gene_name" '{print gn"\t"$0}' "$out_file" >> "$combined_out"
done
```

---

## Combine Already-Generated BLAST Files into One

```bash
#!/bin/bash

# Set paths
input_dir="/home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/result/"
combined_out="/home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/all_results_combined.txt"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$combined_out")"

# Create header for the combined file
header="gene_name\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
echo -e "$header" > "$combined_out"

# Process each BLAST result file
for result_file in "$input_dir"*.txt; do
    [ -f "$result_file" ] || continue
    gene_name=$(basename "$result_file" | sed 's/_blast_results\.txt$//; s/\.txt$//')
    awk -v gn="$gene_name" 'BEGIN {OFS="\t"} {print gn, $0}' "$result_file" >> "$combined_out"
done

echo "Success! Combined results saved to:"
echo "$combined_out"
```

---

## Combine and Annotate with Gene Common Names

```bash
#!/bin/bash

# Set paths
input_dir="/home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/result/"
combined_out="/home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/all_results_combined_with_names.txt"

# Create gene name mapping (GeneID -> Common Name)
declare -A gene_names=(
    [AT3G10370]="SDP6"
    [AT3G16950]="LPD1"
    [AT3G18850]="LPAT5"
    [AT3G25110]="FATA"
    [AT3G51520]="DGAT2"
    [AT4G30950]="FAD6"
    [AT1G74960]="FAB1"
    [AT1G75020]="LPAT4"
    [AT1G78690]="LPLAT"
    [AT5G13640]="PDAT"
    [AT1G01610]="GPAT4"
    [AT2G29980]="FAD3"
    [AT2G30200]="EMB3147"
    [AT2G38040]="CAC3"
    [AT2G41540]="GPDHc1"
    [AT2G43710]="Fab2"
    [AT2G45670]="LPEAT2"
    [AT3G02600]="LPP3"
    [AT3G11170]="FAD7"
    [AT3G15820]="PDCT"
    [AT4G01950]="GPAT3"
    [AT5G05580]="FAD8"
    [AT5G15530]="BCCP2"
    [AT1G32200]="ATS1"
    [AT1G51260]="LPAT3"
    [AT2G38110]="GPAT6"
    [AT3G54340]="WRI1"
    [AT3G57650]="LPAT2"
    [AT1G08510]="FATB"
    [AT1G21970]="LEC1"
    [AT1G28300]="LEC2"
    [AT1G48300]="DGAT3"
    [AT1G62640]="KAS III"
    [AT1G12640]="LPCAT1"
    [AT1G63050]="LPCAT2"
    [AT2G19450]="TAG1, DGAT1"
    [AT3G55360]="CER10"
    [AT1G06520]="GPAT1"
    [AT1G02390]="GPAT2"
    [AT3G11430]="GPAT5"
    [AT5G06090]="GPAT7"
    [AT4G34520]="FAE1"
    [AT5G35360]="CAC2"
    [AT1G80950]="LPEAT1"
    [AT3G12120]="FAD2"
)

# Header
echo -e "gene_name\tgene_id\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "$combined_out"

# Process
for result_file in "$input_dir"*.txt; do
    gene_id=$(basename "$result_file" | grep -oE 'AT[1-5]G[0-9]{5}' | head -1)
    gene_name="${gene_names[$gene_id]:-$gene_id}"
    awk -v gn="$gene_name" -v id="$gene_id" 'BEGIN {OFS="\t"} {print gn, id, $0}' "$result_file" >> "$combined_out"
done

echo "Success! Combined results with gene names saved to:"
echo "$combined_out"
```

---

## R Script for Ideogram Plot

```r
setwd("Y:/Bigdata/computing/Shima/b_carinata/blast_gene/plot")
require(ggideogram)
library(rlang)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load karyotype data
karyotype_a_bna_unbinned <- read.table("brassica_c_plot.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load gene data
genes <- read.table("annotated_genes_plot.txt", sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE,
                    col.names = c("chr", "start", "end", "gene_name",
                                  paste0("extra_", seq_len(10))),
                    fill = TRUE,
                    colClasses = c("character", "numeric", "numeric", "character",
                                   rep("NULL", 10))) %>%
  mutate(Mid = (start + end) / 2)

# Mapping
genbank_to_chromosome <- c(
  "CM081008.1" = "B1", "CM081009.1" = "B2", "CM081010.1" = "B3",
  "CM081011.1" = "B4", "CM081012.1" = "B5", "CM081013.1" = "B6",
  "CM081014.1" = "B7", "CM081015.1" = "B8", "CM081016.1" = "C1",
  "CM081017.1" = "C2", "CM081018.1" = "C3", "CM081019.1" = "C4",
  "CM081020.1" = "C5", "CM081021.1" = "C6", "CM081022.1" = "C7",
  "CM081023.1" = "C8", "CM081024.1" = "C9"
)

genes$chr <- genbank_to_chromosome[genes$chr]

# Check
cat("Unique chromosomes in gene data:", unique(genes$chr), "\n")

# Colors
cytoband_colors <- c("B" = "#5aae61", "C" = "#fdae6b")

# Plot
kardy <- ggplot(karyotype_a_bna_unbinned) +
  geom_ideogram(aes(x = 1, ymin = start, ymax = end,
                    chrom = chr, fill = color_group, arm = arm),
                radius = unit(4, 'pt'), width = 0.1,
                linewidth = 0.2, colour = 'black', show.legend = FALSE) +
  scale_fill_manual(values = cytoband_colors) +
  geom_segment(data = genes,
               aes(x = 1, xend = 1, y = start, yend = end),
               color = "red", size = 0.5,
               position = position_nudge(x = 0.2)) +
  geom_text_repel(data = genes,
                  aes(x = 1, y = Mid, label = gene_name),
                  size = 2, hjust = 0, nudge_x = 0.5,
                  color = "black",
                  arrow = arrow(length = unit(0.015, "npc")),
                  direction = "y", force_pull = 0,
                  segment.size = 0.2, segment.curvature = -0.1,
                  segment.ncp = 3, segment.angle = 20,
                  max.overlaps = Inf, box.padding = 0.5,
                  point.padding = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)),
                     limits = c(0.9, 3)) +
  facet_wrap(vars(chr), scales = 'free_y') +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

# Display and save
print(kardy)
ggsave("karyotype_plot1_unbinned_green_orange_genes_2.png", plot = kardy, width = 12, height = 10, dpi = 300)
```

---

## Karyotype Format Example

```
chr     arm   start       end         color_group
A01     p     1           17000000    A
A01     q     17000001    32958928    A
A02     p     1           20000000    A
A02     q     20000001    33432960    A
A03     p     1           35000000    A
A03     q     35000001    39685748    A
```

