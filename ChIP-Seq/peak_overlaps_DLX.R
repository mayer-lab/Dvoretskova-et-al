#DLX and MEIS peaks in mouse

## overlap meis2 peaks with dlx-merged table form from Lindtner et al., 2019
####### keep nearest gene info just to compare for now
## do overlaps with Vista enhancers etc. afterwards, based on new merged peak coordinates 
## in cases where several close peaks in dlx dataset overlap the same peak in meis2 dataset, keep only the one with highest fold change in dlx2 (has the most peaks)

library(tidyverse) #v1.3.1
library(plyranges, lib.loc= "/home/volker.kittke/R/x86_64-pc-linux-gnu-library/4.0") #v1.8.0
library(liftOver)  #v1.12.0
library(readxl)    #v1.4.0
library(UpSetR)

colnames <- c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" )

## read in the ChIP-seq data: from Lindtner et al.(Dlx5), ING/Torres lab (Meis2)
setwd("/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/")

########################################
##### import peaks from Lindtner et al., supplementary table S3
#
dlx_suppT3 <- read_xls("./scripts/final_scripts_for_Christian/input_data/Lindnter2019_Table_S3.xls", skip=2)
dlx_suppT3.GR <- as_granges(dlx_suppT3)
#names(dlx_suppT3)
# liftOver merged peak coordinates from mm9 to mm10
mm9tom10.chain <- import.chain("./mm9ToMm10.over.chain")
dlx_all_lifted_mm10.df <- liftOver(dlx_suppT3.GR, mm9tom10.chain) %>% 
  as.data.frame() %>% 
  distinct(PeakID,.keep_all=T) %>% # remove overlapping peaks 
  dplyr::select(3:24, 29:32,41:44, 53:56 )

names(dlx_all_lifted_mm10.df)
###########################
#########Meis2
MEIS_GE_e14_mm10 <- read.table("./scripts/final_scripts_for_Christian/input_data/GE_meis2_IP_q0.01_peaks.narrowPeak", col.names=colnames)  %>% mutate(strand="*")
MEIS_GE_e14_mm10.GR <- as_granges(MEIS_GE_e14_mm10) 

# blacklist regions from ENCODE
mm10_blacklist <- read.table("./scripts/final_scripts_for_Christian/input_data/mm10_blacklist_ENCFF547MET.bed", header=F, col.names=c("seqnames","start","end"))

# remove blacklisted regions
mm10_blacklist.GR <- mm10_blacklist %>% as_granges()
meis2_e14_mm10_filtered.GR <-  filter_by_non_overlaps(MEIS_GE_e14_mm10.GR,mm10_blacklist.GR)

## save as bed files

meis2_e14_mm10.df <- as.data.frame(meis2_e14_mm10_filtered.GR)  %>% dplyr::select(1,2,3,6,7,5,8:11)
#write.table(mutate(meis2_e14_mm10.df, strand = ".")[,1:6], "./output/peaksets/processed_peaksets/meis2_e14_mm10.bed", col.names = F,row.names = F,quote = F, sep="\t")
#write_rds(meis2_e14_mm10.df, "./output/peaksets/processed_peaksets/meis2_e14_mm10.rds")


##################### peak overlaps
###################
#################
#### create merged peaklist by merging combining all peaks in one dataset, merging overlapping peaks and assigning "merged peak name"
names(meis2_e14_mm10.df)
names(dlx1_e13_lifted_mm10.df)


all_peaks <- bind_rows(meis2 = meis2_e14_mm10.df[,1:4],dlx= dlx_all_lifted_mm10.df[,c(1:3,7)], .id="ID") %>% mutate(row = 1:nrow(.))
all_peaks.GR <- as_granges(all_peaks)
#View(as.data.frame(all_peaks.GR))
merged_peaks.GR <- GenomicRanges::reduce(all_peaks.GR,with.revmap=T) # revmap = T: indicates which rows of the original table were merged to create the new ranges
length(merged_peaks.GR) #18564

####
merged_peaks.GR$merged_peak_ID <- paste("merged_peak_",1:18564, sep="") 
#View(as.data.frame(merged_peaks.GR))



# overlap idividual peakset ranges with merged peak list to assign merged peak IDs to the narrowPeak tables
meis2_e14_mm10.GR <- as_granges(meis2_e14_mm10.df)
dlx_all_lifted_mm10.GR <-as_granges(dlx_all_lifted_mm10.df) 

meis2_peaks_IDs <- plyranges::find_overlaps(merged_peaks.GR, meis2_e14_mm10.GR, minoverlap = 1) %>% as.data.frame() %>% 
  left_join(dplyr::select(as.data.frame(meis2_e14_mm10.GR), seqnames_meis2 = seqnames, start_meis2 = start, end_meis2 = end, width_meis2=width, name),by="name") %>% dplyr::select(-c(1:6))

length(pull(meis2_peaks_IDs, merged_peak_ID))
length(unique(pull(meis2_peaks_IDs, merged_peak_ID))) # 37 peaks are "duplicates", i.e. 2 or more MEIS2 peaks overlap with the same "merged peak" -> keep the MEIS2 peak with higher fold-change by sorting according to FC beforehand

meis2_peaks_IDs_2 <- meis2_peaks_IDs %>%
  dplyr::arrange(desc(fold_change)) %>%
  distinct(merged_peak_ID, .keep_all=T)


dlx_peaks_IDs <- plyranges::find_overlaps(merged_peaks.GR, dlx_all_lifted_mm10.GR, minoverlap = 1)  %>% as.data.frame() %>%
  left_join(dplyr::select(as.data.frame(dlx_all_lifted_mm10.GR), seqnames_dlx = seqnames, start_dlx = start, end_dlx = end, width_dlx=width, PeakID),by="PeakID")  %>% dplyr::select(-c(1:6))
length(pull(dlx_peaks_IDs, merged_peak_ID)) 
length(unique(dlx_peaks_IDs$merged_peak_ID)) #3 duplicates

dlx_peaks_IDs_2 <- dlx_peaks_IDs %>%
  dplyr::arrange(desc(Dlx1.BG.e13.5.peakHeight)) %>%
  distinct(merged_peak_ID, .keep_all=T)

# convert back to df and join tables
merged_peaks.df <- as.data.frame(merged_peaks.GR) %>%
  left_join(meis2_peaks_IDs_2, by = "merged_peak_ID") %>%
  left_join(dlx_peaks_IDs_2,  by = "merged_peak_ID")  %>%
  dplyr::select(- c(5,6)) %>%
  dplyr::rename("seqnames_merged"= "seqnames", "start_merged" = "start", "end_merged" = "end", "width_merged"="width", 
                "peakname_meis2" = "name", "score_meis2" = "score", "fold_change_meis2" = "fold_change", 
                "X.log10pvalue_meis2" = "X.log10pvalue", "X.log10qvalue_meis2" = "X.log10qvalue", 
                "rel_summit_pos_meis2" = "rel_summit_pos")
  

names(merged_peaks.df)

#######################  pause
####################### 



# add ID column for each ChIP y/n
#########
merged_peaks.df2 <- merged_peaks.df %>%
  mutate(meis2_peak = ifelse(! is.na(peakname_meis2), "y", "n"), 
         dlx1_peak = ifelse(! is.na(Dlx1.BG.e13.5.PeakID), "y", "n"),
         dlx2_peak = ifelse(! is.na(Dlx2.BG.e13.5.PeakID), "y", "n"),
         dlx5_peak = ifelse(! is.na(Dlx5.BG.e13.5.PeakID), "y", "n"))


#### Upset plot
# create names list 


peaks_overlap_list <- list(MEIS12 = pull(filter(merged_peaks.df2, meis2_peak == "y"), merged_peak_ID), 
                           DLX1 = pull(filter(merged_peaks.df2, dlx1_peak == "y"), merged_peak_ID),
                           DLX2 = pull(filter(merged_peaks.df2, dlx2_peak == "y"), merged_peak_ID),
                           DLX5 = pull(filter(merged_peaks.df2, dlx5_peak == "y"), merged_peak_ID))
str(peaks_overlap_list)

peaks_overlap_list$MEIS2

upset(fromList(peaks_overlap_list),set_size.angles = 90) 

# text_scale= c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)

upset_peaks_overlap1  <- upset(fromList(peaks_overlap_list), order.by="freq", set_size.show=F, set_size.scale_max=15000, set_size.angles = 90,
                              sets.bar.color=c("black","black","black", "black"), sets.x.label = "peak count",
                              text.scale=1, shade.alpha =  1,shade.color="darkgrey",nintersects = 10, mb.ratio = c(0.7,0.3),
                              sets = c("MEIS12","DLX5", "DLX2", "DLX1" ), keep.order = T)

upset_peaks_overlap1

png(file="/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peak_overlaps/peak_overlap_meis_dlx_1.png", res=400,  width=80, height=80 , units = "mm") # or other device
upset_peaks_overlap1
dev.off()

pdf(file="/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peak_overlaps/peak_overlap_meis_dlx_1.pdf", onefile=FALSE,width=3.149, height=3.149 ) # or other device
upset_peaks_overlap1
dev.off()


####### keep all intersects

upset_peaks_overlap2  <- upset(fromList(peaks_overlap_list), order.by="freq", set_size.show=F, set_size.scale_max=15000, set_size.angles = 90,
                               sets.bar.color=c("black","black","black", "black"), sets.x.label = "peak count",
                               text.scale=1, shade.alpha =  1,shade.color="darkgrey",nintersects = 20, mb.ratio = c(0.7,0.3),
                               sets = c("MEIS12","DLX5", "DLX2", "DLX1" ), keep.order = T)

upset_peaks_overlap2

png(file="/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peak_overlaps/peak_overlap_meis_dlx_2.png", res=400,  width=80, height=80 , units = "mm") # or other device
upset_peaks_overlap2
dev.off()

pdf(file="/home/volker.kittke/DLX_and_MEIS_binding_mouse_embryo/output/peak_overlaps/peak_overlap_meis_dlx_2.pdf",   width=80, height=80 ) # or other device
upset_peaks_overlap2
dev.off()

text.scale=c(2,1.5,2,1.5,2,1.5)
upset()
#######################
###########################
  
## save

write.table(merged_peaks.df2, "./output/peaksets/processed_peaksets/meis2_dlx125_overlap.txt", col.names = T,row.names = F,quote = F, sep="\t")
write_rds(merged_peaks.df2, "./output/peaksets/processed_peaksets/meis2_dlx125_overlap.RDS")

