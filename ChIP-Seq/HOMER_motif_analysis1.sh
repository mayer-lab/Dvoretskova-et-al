#!/bin/bash


#HOMER motif analysis
# data generation for Fig3 c & d, S3 b & c



# all MEIS2 peaks, centered on peak summit (determined by MACS2)
# MEIS2, 200bp window, 6,8,10,12 bp motifs
nohup findMotifsGenome.pl ./meis2_e14_mm10_500bp_summit.bed mm10 ./HOMER/MEIS2_200bp_new -size 200 -len 6,8,10,12 > output.txt &

# MEIS2_promoter-overlapping peaks, 200bp window(default), 6,8,10,12 bp motifs
nohup findMotifsGenome.pl ./meis2_e14_mm10_promoter_peaks_500bp.bed mm10 ./MEIS2_200bp_promoter_peaks_new -size 200 -len 6,8,10,12 & > output.txt

# MEIS2_enhancer-overlapping peaks, 200bp window(default), 6,8,10,12 bp motifs
nohup findMotifsGenome.pl ./meis2_e14_mm10_enhancer_peaks_500bp.bed mm10 ./MEIS2_200bp_enhancer_peaks_new -size 200 -len 6,8,10,12 & > output.txt

# density plot: using hexa and deca motif (HOMER findMotifsGenome output) as input; 
annotatePeaks.pl ./meis2_e14_mm10.bed mm10 -size 500 -hist 5 -m ./MEIS2_200bp_new/homerResults/hexa_deca_motifs.motif -o "./MEIS2_200bp_new/density/motif_density_meis2peaks.txt"



