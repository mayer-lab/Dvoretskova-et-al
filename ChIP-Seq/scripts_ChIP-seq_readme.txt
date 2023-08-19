# preprocessing_peak_calling_MEIS2_mouse_GE.sh
MEIS1/2 ChIP-seq analysis starting from .fastq files to peak calling

# peak overlaps.R
overlap of ChIP-seq peak regions for MEIS1/2, DLX5 and LHX6 datasets.
overlap of peaks with promoters from EPD, enhancers from Gorkin et al. and with Vista enhancers
plot Fig S3e

# peak-TSS_distance.R
calculation of distance from MEIS1/2 peaks to the nearest TSS; plot for Fig. 3A

# plot_enhancer-promoter_vista.R 
Upset plots for overlap of MEIS1/2 and DLX5 peaks with promoters from EPD, enhancers from Gorkin et al. and with Vista enhancers
Fig. 3e, S3d

# HOMER_motif_analysis_preprocessing.R
extraction of peak sequences (coordinates) of MEIS1/2 binding sites
 
# HOMER_motif_anaylsis1.sh
motif discovery and enrichment in MEIS1/2 peaks

# HOMER_motif_anaylsis2.R
plots for motif density and motif enrichment, Figs. 3c&d, S3b

# SpaMo_input_sequences.R
extraction of peak sequences for motif spacing analysis of MEIS1/2-DLX5 overlapping peaks

# SpaMo_motif_analysis.txt
motif spacing analysis parameters
Fig S3c