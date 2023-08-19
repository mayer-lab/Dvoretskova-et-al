#########
## motif analysis MEIS/DLX

library(plyranges)
library(tidyverse)
library(data.table)
library("BSgenome.Mmusculus.UCSC.mm10")
library(seqinr)
library(ChIPpeakAnno)
genome.mm10 <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10")


merged_peaks.df5 <- read_rds("./Table_S5_ChIP-seq.RDS")

meis2_peaks.df <- merged_peaks.df5 %>% 
  filter(meis2_peak == "y") %>%
  mutate(strand = ".") %>%
  dplyr::select(12:14, 16, 17, strand, 18:21) %>%
  data.table::setnames(c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" )) 


############################################
##### extract peak sequences for downstream analysis
# of +- 250bp around peak summit (for motif enrichment within peaks)

# MEIS2 all peaks, 500 bp

meis2_e14_mm10_summit_500bp <- meis2_peaks.df %>%
  mutate(start2 = start + rel_summit_pos - 250, end2 = start + rel_summit_pos + 250) %>%
  mutate(width = end2 - start2) %>%
  dplyr::select(1, start=start2, end= end2, 4:6)

write.table(meis2_e14_mm10_summit_500bp,"./meis2_e14_mm10_500bp_summit.bed", quote=F, col.names=F, row.names = F, sep= "\t")


# MEIS2 promoter overlapping peaks, 500 bp

meis2_peaks_prom.df <- merged_peaks.df5 %>% 
  filter(meis2_peak == "y" & !is.na(TSS_overlap)) %>%
  mutate(strand = ".") %>%
  dplyr::select(12:14, 16, 17, strand, 18:21)  %>%
  data.table::setnames(c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" ))  %>%
  mutate(start2 = start + rel_summit_pos - 250, end2 = start + rel_summit_pos + 250) %>%
  mutate(width = end2 - start2) %>%
  dplyr::select(1, start=start2, end= end2, 4:6)

write.table(meis2_peaks_prom.df,"./meis2_e14_mm10_promoter_peaks_500bp.bed", quote=F, col.names=F, row.names = F, sep= "\t")


# MEIS2 enhancer overlapping peaks, 500 bp

meis2_peaks_enh.df <- merged_peaks.df5 %>% 
  filter(meis2_peak == "y" & !is.na(enhancer_interacting_genes)) %>%
  mutate(strand = ".") %>%
  dplyr::select(12:14, 16, 17, strand, 18:21)  %>%
  data.table::setnames(c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" ))  %>%
  mutate(start2 = start + rel_summit_pos - 250, end2 = start + rel_summit_pos + 250) %>%
  mutate(width = end2 - start2) %>%
  dplyr::select(1, start=start2, end= end2, 4:6)

write.table(meis2_peaks_enh.df,"./meis2_e14_mm10_enhancer_peaks_500bp.bed", 
            quote=F, col.names=F, row.names = F, sep= "\t")


# overlap between enhancer & promoter peaks:

intersect(meis2_peaks_prom.df$name,meis2_peaks_enh.df$name )
#35 peaks



# MEIS2 & DLX5 common peaks, overlapping enhancers, 500 bp

meis2_dlx5_peaks_enh.df <- merged_peaks.df8 %>% 
  filter(meis2_peak == "y" & enhancer_overlap == "y" & dlx5_peak == "y") %>%
  mutate(strand = ".") %>% 
  dplyr::select(22:24, 26, 27, strand, 28:31) %>%
  data.table::setnames(c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" ))  %>%
  mutate(start2 = start + rel_summit_pos - 250, end2 = start + rel_summit_pos + 250) %>%
  mutate(width = end2 - start2) %>%
  dplyr::select(1, start=start2, end= end2, 4:6)

write.table(meis2_dlx5_peaks_enh.df,"/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peaksets/processed_peaksets/meis2_dlx5_enhancer_peaks_500bp.bed", 
            quote=F, col.names=F, row.names = F, sep= "\t")


# MEIS2 & DLX5 common peaks, overlapping enhancers & target genes are DE in PN, 500 bp
PN_de_genes <- read.table("/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/target_genes/E16_proj_DE_filter.csv", header=T, sep=",")

meis2_dlx5_peaks_enh_DE.df <- merged_peaks.df8 %>% 
  filter(meis2_peak == "y" & enhancer_overlap == "y" & dlx5_peak == "y") %>%
  filter(grepl(paste(PN_de_genes$gene, collapse="|"),enhancer_interacting_genes )) 
write.table(meis2_dlx5_peaks_enh_DE.df, "/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peaksets/processed_peaksets/gene_overlap/meis2_dlx5_peaks_enhancer_DE.df.txt", quote=F, col.names=T,row.names=F, sep="\t")

# 19 peaks only
meis2_dlx5_peaks_enh_DE.df2 <- meis2_dlx5_peaks_enh_DE.df %>%
  mutate(strand = ".") %>% 
  dplyr::select(22:24, 26, 27, strand, 28:31) %>%
  data.table::setnames(c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" ))  %>%
  mutate(start2 = start + rel_summit_pos - 250, end2 = start + rel_summit_pos + 250) %>%
  mutate(width = end2 - start2) %>%
  dplyr::select(1, start=start2, end= end2, 4:6)

write.table(meis2_dlx5_peaks_enh_DE.df2,"/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peaksets/processed_peaksets/meis2_dlx5_enhancer_peaks_DE_500bp.bed", 
            quote=F, col.names=F, row.names = F, sep= "\t")

## LHX6 500bp 

lhx6_peaks_all.df <- merged_peaks.df8 %>% 
  filter(lhx6_peak == "y") %>%
  mutate(strand = ".") %>%
  dplyr::select(42:44,46,47,strand,48:51) %>%
  data.table::setnames(c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" ))  %>%
  mutate(start2 = start + rel_summit_pos - 250, end2 = start + rel_summit_pos + 250) %>%
  mutate(width = end2 - start2) %>%
  dplyr::select(1, start=start2, end= end2, 4:6)

write.table(lhx6_peaks_all.df,"/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peaksets/processed_peaksets/lhx6_all_peaks_summit_500bp.bed", quote=F, col.names=F, row.names = F, sep= "\t")


