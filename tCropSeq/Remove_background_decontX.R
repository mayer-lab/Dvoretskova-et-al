if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("celda")



library(celda)
library(SingleCellExperiment)

setwd("/data/mayerlab/cmayer/meis2_paper/results")

P7GE_sub <- readRDS("/datastore_share/Users/mayho/meis2_paper/seurat_objects/P7GE_sub.rds")
DimPlot(P7GE_sub)
DimPlot(P7GE_sub, group.by = "level4")
Idents(object = P7GE_sub) <- "level4"
DimPlot(P7GE_sub)

P7GE_sub@meta.data[, "level4"] <- plyr::mapvalues(P7GE_sub@meta.data$level4,
                                                               from = c("Cck int"),
                                                               to = c("Interneurons"))


## DecontX II
markers <- c("Dlx6os1","Stmn2","Igfbpl1","Tmsb10","Zfp704","Pdgfra","Lhfpl3","Cntn1","C1ql1","Ostf1","Top2a","2810417H13Rik","Ube2c","Hist1h1b","Hmgb2","Gfap","Igfbp5","Ntrk2","Prdx6","Clu","Zbtb18","Uchl1","Nrgn","Eif1b","Epha5","Sparcl1","Apoe","Slc4a4","Slc1a2","Pla2g7","Cpne4","Rprml","Synpr","Grin2b","Isoc1","Bcas1","S100a13","Sirt2","Mbp","Plp1","Foxp2","Sst","Gap43","Pcp4","Npy","Syt1","Gucy1a3","Tshz1","Serpini1","Reln","Cxcl12","Snhg11","Ndnf","Meg3","Cck",
  "Maf","Cnr1","Cxcl14","Trh","Nrsn1","Ndrg4","Gad2")

idents.order <- c("OB prec","OPCs","Mitotic Ascl1","Astro Gfap","Excitatory","Astro Aqp4","OB Cpne4","Oligo","DRD1/DRD2 MSN","ITC","Cajal Retzius","Interneurons","OB Th")


counts <- GetAssayData(object = P7GE_sub, slot = "counts")
sce <- SingleCellExperiment(list(counts = counts))
sce <- decontX(sce)
P7GE_sub[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce))



DefaultAssay(P7GE_sub) <- "RNA"
Idents(object = P7GE_sub) <- factor(Idents(object = P7GE_sub), levels = idents.order)
DotPlot(P7GE_sub, features = markers, assay = "RNA")

DefaultAssay(P7GE_sub) <- "decontXcounts"
Idents(object = P7GE_sub) <- factor(Idents(object = P7GE_sub), levels = idents.order)
DotPlot(P7GE_sub, features = markers, assay = "decontXcounts")



# Load required libraries
library(ggplot2)
library(gridExtra)


# Set the figure size (adjust as needed)
options(repr.plot.width=14, repr.plot.height=8)


p1 <- DotPlot(P7GE_sub, features = markers, assay = "RNA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  ggtitle("RNA DotPlot")

# Create the second plot with rotated x-axis labels (45 degrees)
p2 <- DotPlot(P7GE_sub, features = markers, assay = "decontXcounts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  ggtitle("decontXcounts DotPlot")

# Combine the two plots vertically using grid.arrange
final_plot <- grid.arrange(p1, p2, ncol = 1)

# Save the final plot as a PDF
ggsave("high_profile_figure.pdf", final_plot, width = 14, height = 8)

# Display the final plot (optional, comment this out if you only want to save the PDF)
print(final_plot)





