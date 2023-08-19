#################
# extraction of MEIS1/2-DLX5 overlapping peak sequences for input in SpaMo
# using MEIS1/2 peak coordinates, summit position +/- 250 bp

genome.mm10 <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Mmusculus.UCSC.mm10")
library(seqinr)


merged_peaks.df5 <- read_rds("./Table_S5_ChIP-seq.RDS")

meis2_dlx5_peaks_500bp_summit.df <- merged_peaks.df5 %>% 
  filter(meis2_peak == "y", dlx5_peak =="y") %>%
  mutate(strand = ".") %>%
  dplyr::select(12:14, 16, 17, strand, 18:21) %>%
  data.table::setnames(c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" )) %>%
  mutate(start2 = start + rel_summit_pos - 250, end2 = start + rel_summit_pos + 250) %>%
  mutate(width = end2 - start2) %>%
  dplyr::select(1, start=start2, end= end2, 4:6) %>%
  dplyr::mutate(strand = "*") %>%
  dplyr::filter(seqnames != "chrM") 

meis2_dlx5_peaks_500bp_summit.GR <- as_granges(meis2_dlx5_peaks_500bp_summit.df)

seq_meis2_dlx5_peaks_500bp_summit <- as.data.frame(getAllPeakSequence(myPeakList = meis2_dlx5_peaks_500bp_summit.GR, upstream=0, downstream=0,genome=genome.mm10)) 
write.fasta(sequences=as.list(dplyr::pull(seq_meis2_dlx5_peaks_500bp_summit,sequence)), names=pull(seq_meis2_dlx5_peaks_500bp_summit,name), "./seq_meis2_dlx5_overlap_500bp_summit_seqs.fa", open = "w", nbchar = 60, as.string = FALSE)

