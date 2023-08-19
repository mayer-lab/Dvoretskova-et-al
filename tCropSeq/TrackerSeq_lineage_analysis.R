# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(philentropy)
library(dendextend)
library(UpSetR)
library(ComplexUpset)

# Lineage Barcode distribution between classes.  
# Author: May Ho 

# Load Seurat object
E16_lin_subset <- readRDS("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/E16_lin_subset.rds")

# Data Preprocessing: Calculate clone sizes and filter out single-cell clones
E16_lin_subset_clones <- na.omit(E16_lin_subset@meta.data[, c("gRNA", "cloneID")])
clone_size_tot <- E16_lin_subset_clones %>%
  group_by(gRNA, cloneID) %>%
  summarise(clone_size = n())
clone_size_remove <- clone_size_tot[clone_size_tot$clone_size == 1,]
clone_size_tot <- clone_size_tot[!(clone_size_tot$cloneID %in% clone_size_remove$cloneID),]
clone_size_tot <- clone_size_tot %>%
  ungroup() %>%
  group_by(gRNA, clone_size) %>%
  summarise(n = n())

# Calculate mean and standard deviation for clone sizes
clone_size_lacZ <- clone_size_tot[clone_size_tot$gRNA == "glacZ",]
clone_size_meis2 <- clone_size_tot[clone_size_tot$gRNA == "gMeis2",]
clone_size_data <- data.frame(
  gRNA = c("gMeis2", "glacZ"),
  clone_size = c(mean(clone_size_meis2$clone_size), mean(clone_size_lacZ$clone_size)),
  sd = c(sd(clone_size_meis2$clone_size), sd(clone_size_lacZ$clone_size))
)

# Visualize the result using a bar plot with error bars
ggplot(clone_size_data) +
  geom_bar(aes(x = gRNA, y = clone_size), stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(x = gRNA, ymin = clone_size - sd, ymax = clone_size + sd), width = 0.1) +
  theme_classic()

# Clonal Intersection UpSet diagram
merge_cloneseq <- data.frame(class = character(),
                             dataset = character(),
                             value = numeric(),
                             gRNA = character(),
                             stringsAsFactors = FALSE)

gRNAs <- unique(E16_lin_subset@meta.data$gRNA)
for (i in gRNAs) {
  seur_obj <- subset(E16_lin_subset, subset = gRNA == i)
  pool <- FetchData(seur_obj, c("cloneID", "neuron_class"))
  pool <- pool[complete.cases(pool), ]
  pool <- pool %>%
    ungroup() %>%
    group_by(cloneID) %>%
    filter(n() != 1)
  pool <- as.data.frame(pool)
  cloneseq <- fromList(split(pool$cloneID, pool$neuron_class))
  cloneseq$gRNA <- i
  merge_cloneseq <- rbind(merge_cloneseq, cloneseq)
}

# Visualize clonal intersection between mitotic, interneurons, and projection neurons
ComplexUpset::upset(
  merge_cloneseq, c("mitotic", "projection neurons", "interneurons"), name = 'gRNA',
  width_ratio = 0.3, min_size = 10, min_degree = 1,
  annotations = list(
    'gRNA' = list(
      aes = aes(x = intersection, fill = gRNA),
      geom = list(
        geom_bar(stat = 'count', position = 'fill', na.rm = TRUE),
        geom_text(
          aes(label = !!aes_percentage(relative_to = 'group'), group = gRNA),
          stat = 'count',
          position = position_fill(vjust = 0.5)
        ),
        scale_y_continuous(labels = scales::percent_format()),
        scale_color_manual(values = c('show' = 'black', 'hide' = 'transparent'), guide = 'none'),
        scale_fill_manual(values = c('glacZ' = '#377EB8', 'gMeis2' = '#FF7F00'))
      )
    )
  )
)

# Visualize individual clones on UMAP
clonesviz <- c(114, 2910, 4132, 7744)
for (i in clonesviz) {
  DimPlot(
    Eminence.combined.sct, reduction = "umap",
    cells.highlight = which(Eminence.combined.sct@meta.data$cloneID == i),
    sizes.highlight = 4, pt.size = 0.5,
    cols.highlight = "royalblue4"
  ) + NoLegend()
  ggsave(file = paste("E16_clones/Clone_", i, ".png", sep = ""), width = 2, height = 2, scale = 4, dpi = 600)
}




########################## Z-score analysis ##################################
# Author: Florian Neuhaus 

# load dataset -  the important metadata are: gRNA, cloneID and neuron_class 
E16_lin_subset <- readRDS("E16_lin_subset.rds")

## subset data into gMeis2 and gLacz:
gLacZ_e16_line_sub_seurat <- subset(E16_lin_subset, subset = gRNA == "glacZ")
gMeis2_e16_line_sub_seurat <- subset(E16_lin_subset, subset = gRNA == "gMeis2")

## save cell- clone- and cluster-ID in table:
gLacZ_e16_line_sub_df <- data.frame(
  "cellID" = colnames(gLacZ_e16_line_sub_seurat),
  "cloneID" = gLacZ_e16_line_sub_seurat$cloneID,
  "ident" = gLacZ_e16_line_sub_seurat$neuron_class
)
rownames(gLacZ_e16_line_sub_df) <- 1:nrow(gLacZ_e16_line_sub_df)
gLacZ_e16_line_sub_df <- gLacZ_e16_line_sub_df[!is.na(gLacZ_e16_line_sub_df$cloneID), ]
clone_dist <- table(gLacZ_e16_line_sub_df$cloneID)
gLacZ_e16_line_sub_df <- gLacZ_e16_line_sub_df[clone_dist[gLacZ_e16_line_sub_df$cloneID] > 1, ]
write.table(gLacZ_e16_line_sub_df, "results/data_for_weinreb/gLacZ_e16_line_sub_table_neuron_class.csv", sep = ",", row.names = F)
###

gMeis2_e16_line_sub_df <- data.frame(
  "cellID" = colnames(gMeis2_e16_line_sub_seurat),
  "cloneID" = gMeis2_e16_line_sub_seurat$cloneID,
  "ident" = gMeis2_e16_line_sub_seurat$neuron_class
)
rownames(gMeis2_e16_line_sub_df) <- 1:nrow(gMeis2_e16_line_sub_df)
gMeis2_e16_line_sub_df <- gMeis2_e16_line_sub_df[!is.na(gMeis2_e16_line_sub_df$cloneID), ]
clone_dist <- table(gMeis2_e16_line_sub_df$cloneID)
gMeis2_e16_line_sub_df <- gMeis2_e16_line_sub_df[clone_dist[gMeis2_e16_line_sub_df$cloneID] > 1, ]
write.table(gMeis2_e16_line_sub_df, "results/data_for_weinreb/gMeis2_e16_line_sub_table_neuron_class.csv", sep = ",", row.names = F)

## run python script ##

system("conda activate lineage_coupling")
## n_iteration=2000
system("python lineage_coupling_analysis_weinreb_way_5.py -N 2000 -u 0.0 -v 12.0 -f '/datastore_share/Users/neuhaus/Meis2_z_score/results/data_for_weinreb/gLacZ_e16_line_sub_table_neuron_class.csv' -L '/datastore_share/Users/neuhaus/Meis2_z_score/results/weinreb/gLacZ_e16_lin_sub/lineage_coupling_scores_matrix_Method_Weinreb_et_al_N_2000_gLacZ_e16_lin_sub_table_neuron_class.csv' -M '/datastore_share/Users/neuhaus/Meis2_z_score/results/weinreb/gLacZ_e16_lin_sub/metric_values_per_state_pair_matrix_N_2000_gLacZ_e16_lin_sub_table_neuron_class.csv'")
system("python lineage_coupling_analysis_weinreb_way_5.py -N 2000 -u 0.0 -v 12.0 -f '/datastore_share/Users/neuhaus/Meis2_z_score/results/data_for_weinreb/gMeis2_e16_line_sub_table_neuron_class.csv' -L '/datastore_share/Users/neuhaus/Meis2_z_score/results/weinreb/gMeis2_e16_lin_sub/lineage_coupling_scores_matrix_Method_Weinreb_et_al_N_2000_gMeis2_e16_lin_sub_table_neuron_class.csv' -M '/datastore_share/Users/neuhaus/Meis2_z_score/results/weinreb/gMeis2_e16_lin_sub/metric_values_per_state_pair_matrix_N_2000_gMeis2_e16_lin_sub_table_neuron_class.csv'")

## visualize results:
gLacZ_WE_mtx2 <- read.table("/datastore_share/Users/neuhaus/Meis2_z_score/results/weinreb/gLacZ_e16_lin_sub/lineage_coupling_scores_matrix_Method_Weinreb_et_al_N_2000_gLacZ_e16_lin_sub_table_neuron_class.csv", h=T, sep=",",row.names = 1)
gLacZ_WE_mtx2 <- gLacZ_WE_mtx2[row_order, col_order]
p11 <- pheatmap(gLacZ_WE_mtx2, cluster_rows = F, cluster_cols = F, display_numbers = T, color = hcl.colors(50, "Peach", rev = T), breaks = seq(0.78,1.60, length.out = 50))[[4]]

gMeis2_WE_mtx2 <- read.table("/datastore_share/Users/neuhaus/Meis2_z_score/results/weinreb/gMeis2_e16_lin_sub/lineage_coupling_scores_matrix_Method_Weinreb_et_al_N_2000_gMeis2_e16_lin_sub_table_neuron_class.csv", h=T, sep=",",row.names = 1)
gMeis2_WE_mtx2 <- gMeis2_WE_mtx2[row_order, col_order]
p21 <- pheatmap(gMeis2_WE_mtx2, cluster_rows = F, cluster_cols = F, display_numbers = T, color = hcl.colors(50, "Peach", rev = T), breaks = seq(0.78,1.60, length.out = 50))[[4]]

pdf("/datastore_share/Users/neuhaus/Meis2_z_score/results/weinreb/gLacZ_gMeis2_lin_coupling_heat_N2000.pdf", width = 10, height = 5)
grid.arrange(p11,p21, ncol = 2, left = "gLacZ", right = "gMeis2")
dev.off()


