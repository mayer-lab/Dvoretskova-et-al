# Load required libraries
library(readxl)
library(dplyr)
library(VennDiagram)

# Read the data from Excel file
data <- read_excel("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/Table_S4_E16_clusters_DE.xlsx")

# Get column titles of the data
column_titles <- colnames(data)

# Assuming your data frame is named "data"

# Group the data by cell_type and apply filtering
filtered_data <- data %>%
  group_by(cell_type) %>%
  filter(p_val_adj < 0.05 & (avg_logFC < -1.0 | avg_logFC > 1.0))

# Number of rows in the filtered data
nrow(filtered_data)

# Extract upregulated and downregulated genes
upregulated_genes <- filtered_data %>%
  filter(avg_logFC > 1.0) %>%
  select(gene)

downregulated_genes <- filtered_data %>%
  filter(avg_logFC < -1.0) %>%
  select(gene)

# Number of unique upregulated and downregulated genes
num_upregulated <- length(unique(upregulated_genes$gene))
num_downregulated <- length(unique(downregulated_genes$gene))
total_intersection <- num_upregulated + num_downregulated

# Load ChipSeq data
chipseq_data <- read.table("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/Dvoretskova_et_al_2023-DATASETS/enhancer_promoter_target_genes_MEIS2_gorkin_071122.txt", header = TRUE)

# Create a Venn diagram
venn_data <- list(
  Upregulated = unique(upregulated_genes$gene),
  Downregulated = unique(downregulated_genes$gene),
  ChipSeq = unique(chipseq_data$gene_symbol)
)

venn_diagram_file <- "venn_diagram.png"
venn.diagram(
  x = venn_data,
  category.names = c("Upregulated", "Downregulated", "ChipSeq"),
  filename = venn_diagram_file  # Specify the filename and extension as a PNG
)

# Calculate the lengths of the intersections
intersection_Up_Down_Chip <- length(intersect(intersect(venn_data$Upregulated, venn_data$Downregulated), venn_data$ChipSeq))
intersection_Up_Down <- length(setdiff(intersect(venn_data$Upregulated, venn_data$Downregulated), venn_data$ChipSeq))
intersection_Up_Chip <- length(setdiff(intersect(venn_data$Upregulated, venn_data$ChipSeq), venn_data$Downregulated))
intersection_Down_Chip <- length(setdiff(intersect(venn_data$Downregulated, venn_data$ChipSeq), venn_data$Upregulated))

# Print the numbers of the intersections
cat("Intersection of Upregulated, Downregulated, and ChipSeq:", intersection_Up_Down_Chip, "\n")
cat("Intersection of Upregulated and Downregulated:", intersection_Up_Down, "\n")
cat("Intersection of Upregulated and ChipSeq:", intersection_Up_Chip, "\n")
cat("Intersection of Downregulated and ChipSeq:", intersection_Down_Chip, "\n")

# Final message
cat("Total number of genes in Upregulated and Downregulated categories:", total_intersection, "\n")
cat("Venn diagram saved as", venn_diagram_file, "\n")
