# Load required libraries
library(Seurat)
library(Libra)
library(readr)
library(ggplot2)
library(forcats)
library(EnhancedVolcano)
library(grid)

# Load Seurat objects and perform cell type assignment
E16_CCA_inhib <- readRDS("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/E16_CCA_inhib.rds")
E16_CCA_inhib@meta.data[, "CCA_Assignment"] <- plyr::mapvalues(E16_CCA_inhib@meta.data$CCA_Assignment,
                                                               from = c("IN:Tshz1/Pbx1"),
                                                               to = c("PN:Tshz1/Pbx1"))

# Find DE genes using libra
DE <- run_de(E16_CCA_inhib, cell_type_col = "CCA_Assignment", replicate_col = "orig.ident", label_col = "gRNA", de_family = "pseudobulk", de_method = "edgeR", de_type = "LRT")
DE = run_de(E16_CCA_inhib, cell_type_col = "neuron_class", replicate_col = "orig.ident", label_col = "gRNA", de_family = "pseudobulk", de_method = "edgeR", de_type = "LRT")

# Calculate the total number of DEGs for each inhibitory cluster.
DEGs_clusters <- DE %>%
  count(cell_type) %>%
  arrange(desc(n))

celltype_num <- E16_CCA_inhib@meta.data %>%
  count(CCA_Assignment) %>%
  rename(cell_type = CCA_Assignment, num_cells = n)

DEGs_clusters <- left_join(DEGs_clusters, celltype_num)
DEGs_clusters$cell_type <- forcats::fct_reorder(DEGs_clusters$cell_type, DEGs_clusters$n, .desc = TRUE)

# Plot lollipop chart for DEGs in each inhibitory cluster
p <- ggplot(DEGs_clusters, aes(cell_type, n)) +
  geom_segment(aes(xend=cell_type, y=0, yend=n)) +
  geom_point(mapping = aes(size = num_cells, colour = cell_type)) + 
  scale_colour_viridis_d() +
  scale_size_continuous("Number of Cells", range = c(2,10), breaks = c(500,1000,1500,2000,2500), labels = c("500", "1000", "1500", "2000", "2500")) +
  labs(y = "# of DEGs", x = "Cell Type") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_classic()

# Function to create EnhancedVolcano plot
createVolcanoPlot <- function(data, title, interest_genes) {
  EnhancedVolcano(data, 
                  lab = data$gene,
                  x = 'avg_logFC', 
                  y = 'p_val_adj',
                  FCcutoff = 0.25,
                  pCutoff = 0.05,
                  title = title,
                  #xlim = c(-4,4),
                  ylim = c(0,28),
                  legendPosition = 'none',  # Remove the legend
                  drawConnectors = TRUE,
                  boxedLabels = TRUE,
                  legendLabels = c('NS.','LogFC','p-value',"p-value and LogFC"),
                  selectLab = interest_genes,
                  max.overlaps = Inf, 
                  border = "full",
                  labSize = 3,
                  xlab = "log2(Fold Change)",
                  ylab = expression('FDR adjusted -log'["10"]*italic('P')))
}

# E16 Volcano Plots for different cell-classes 
E16_proj_DE_filter <- read_csv("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/E16_proj_DE_filter.csv") # Load DE genes from a CSV file
proj_interest_genes <- unique(c("Adora2a", "Drd1", "Gucy1a3", "Gucy1b3", "Six3", "Nxph1","Tle4",
                                "Htr3a", "Maf", "Sp8", "Tcf4", "Ngf", "Prox1", "Arx", "Ibsp", "Sox14", "Npas1","Cck","Nr2f1"))
createVolcanoPlot(E16_proj_DE_filter, "DEGs in Projection Neurons", proj_interest_genes)
ggsave(filename = "PN_Volcano_4x6_ylim.pdf", width = 4, height = 6, units = "in")

E16_IN_DE_filter <- read_csv("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/E16_interneurons_DE.csv")
IN_interest_genes <- unique(c("Ibsp","Cck","Ngf","Scgn","Slc30a3","Cacna2d3","Neb","Wnt5a","Nr3c2","Pdzrn3","Arrdc4","Rwdd3"))
createVolcanoPlot(E16_IN_DE_filter, "DEGs in Interneurons", IN_interest_genes)
ggsave(filename = "IN_Volcano_4x6.pdf", width = 4, height = 6, units = "in")

E16_M_DE_filter <- read_csv("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/E16_Meis2_mitotic.csv")
M_interest_genes <- unique(c("Wnt5a","Pdzrn3","Npas1","Zic3","Nefm","Esrrg","Tiam2","Serpini1","Erbb4","Rbfox3","Npas1","Gpd1","Sorcs3","Vstm2b"))
createVolcanoPlot(E16_M_DE_filter, "DEGs in Mitotic Cells", M_interest_genes)
ggsave(filename = "Mitotic_Volcano_4x6.pdf", width = 4, height = 6, units = "in")
