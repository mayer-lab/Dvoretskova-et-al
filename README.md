# Dvoretskova-et-al

This repository includes the necessary code to reproduce the analysis presented in the paper titled 'Spatial enhancer activation determines inhibitory neuron identity.'

## Data and Code Associated with the Study

**Title:** Spatial Enhancer Activation Determines Inhibitory Neuron Identity

**DOI:** [https://doi.org/10.1101/2023.01.30.525356](https://doi.org/10.1101/2023.01.30.525356)

## Authors

Elena Dvoretskova, May C. Ho, Volker Kittke, Ilaria Vitali, Daniel D. Lam, Irene Delgado, Chao Feng, Miguel Torres, Juliane Winkelmann, Christian Mayer

**Correspondence:** [christian.mayer@bi.mpg.de](mailto:christian.mayer@bi.mpg.de)

## Data

The sequencing datasets generated for this study are available in the Gene Expression Omnibus (GEO). Preprocessed data to reproduce the analysis shown here can be downloaded from: [https://keeper.mpdl.mpg.de/d/147c56f5bcf74f03a8c0/](https://keeper.mpdl.mpg.de/d/147c56f5bcf74f03a8c0/)

## Scripts

There are individual scripts for different parts of the analysis.

### tCropSeq

This directory contains code used to process and analyze tCROP-seq and TrackerSeq data. Below is a list of the files and a brief explanation of what each script does:

- `ChipSeq_tCropSeq_intersect.R`: This R script performs the intersection analysis between ChIP-seq data and tCrop-seq data. It reads data from an Excel file, filters the data based on adjusted p-values and fold change values, and extracts the upregulated and downregulated genes. The script also reads ChIP-seq data from a text file and creates a Venn diagram to visualize the intersections between the upregulated genes, downregulated genes, and ChIP-seq genes. The script calculates the lengths of the four intersections (upregulated and downregulated genes, upregulated genes and ChIP-seq genes, downregulated genes and ChIP-seq genes, and upregulated genes, downregulated genes, and ChIP-seq genes) and prints the numbers of the intersections.

- `DE_analysis.R`: This R script performs differential expression (DE) analysis using Seurat and Libra packages. It reads Seurat objects and assigns cell types based on the assigned labels. The script then runs DE analysis using the pseudobulk method and the edgeR package. The DE analysis is performed for different cell types separately. The script calculates the total number of differentially expressed genes (DEGs) for each inhibitory cluster and plots a lollipop chart to visualize the DEGs. The script also includes functions to create EnhancedVolcano plots for DEGs in different cell classes and saves the plots as PDF files.

- `GO_Term.Rmd`: This R Markdown file performs gene ontology (GO) enrichment analysis using the enrichR package. It reads differential expression data for upregulated genes in E16 projection neurons, filters the data based on adjusted p-values, and performs GO enrichment analysis for biological processes. The script identifies significantly enriched biological processes based on an adjusted p-value threshold and stores the results in the `sig_bp` variable.

- `ineage_coupling_analysis_weinreb.py`: This Python script calculates the lineage coupling for different cell types to quantify the strength of lineage coupling between different cell types. This method is adapted from Weinreb et al., 2020 (Science), 10.1126/science.aaw3381; and taken from Bandler et al., 2022 (Nature), 10.1038/s41586-021-04237-0; 

- `module.ipynb`: This Jupyter Notebook file contains the code for module analysis. It uses the Hotspot package to identify modules and explore co-expression patterns within the data. The script follows the steps outlined in the documentation of the Hotspot package (https://github.com/neurorestore/Libra) to perform module analysis and generate visualizations.

- `P7_module_analysis.Rmd`: This R Markdown file performs module analysis on P7 data, calculating the effect size of perturbations on module scores, and generating heatmaps of top genes within the modules.

- `Proportion_change.R`: This script calculates the proportion change of perturbed cell clusters compared to the lacZ control. It uses regression analysis as described in [Github - ivPerturbSeq](https://github.com/klarman-cell-observatory/ivPerturbSeq) to model the relationship between cell cluster proportions and perturbations (gRNAs) while accounting for the batch effect. The script then generates a plot showing the effect size and statistical significance of the perturbations on the cell composition of different clusters.

- `Remove_background_decontX.R`: This script uses the celda package to remove background noise from single-cell RNA sequencing data and applies decontX II to improve data quality. It also performs DotPlot visualization for selected gene markers before and after the decontX transformation.

- `TrackerSeq_lineage_analysis.R`: This script performs lineage analysis using the TrackerSeq approach on single-cell RNA sequencing data. It calculates clone sizes, filters out single-cell clones, visualizes the distribution of clone sizes between "glacZ" and "gMeis2" groups, and generates UpSet plots to visualize clonal intersection between different cell classes. The script also performs lineage coupling analysis on the two groups separately.

Please refer to the individual files for more detailed information on the code implementation and usage.

### ChipSeq_Meis2

This directory contains code used to process and analyze Meis2 ChipSeq data. Below is a list of the files and a brief explanation of what each script does:

- `preprocessing_peak_calling_MEIS2_mouse_GE.sh`: MEIS1/2 ChIP-seq analysis starting from .fastq files to peak calling.

- `peak overlaps.R`: Overlap of ChIP-seq peak regions for MEIS1/2, DLX5, and LHX6 datasets with promoters, enhancers, and Vista enhancers.

- `peak-TSS_distance.R`: Calculation of the distance from MEIS1/2 peaks to the nearest TSS; plot for Fig. 3A.

- `plot_enhancer-promoter_vista.R`: Upset plots for overlap of MEIS1/2 and DLX5 peaks with promoters from EPD, enhancers from Gorkin et al., and Vista enhancers.

- `HOMER_motif_analysis_preprocessing.R`: Extraction of peak sequences (coordinates) of MEIS1/2 binding sites.

- `HOMER_motif_anaylsis1.sh`: Motif discovery and enrichment in MEIS1/2 peaks.

- `HOMER_motif_anaylsis2.R`: Plots for motif density and motif enrichment.

- `SpaMo_input_sequences.R`: Extraction of peak sequences for motif spacing analysis of MEIS1/2-DLX5 overlapping peaks.

- `SpaMo_motif_analysis.txt`: Motif spacing analysis parameters.

Please note that additional descriptions and usage details can be found within each script.
