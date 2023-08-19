#DLX and MEIS peaks in mouse

## data from Lindtner et al., 2019

library(tidyverse) #v1.3.1
library(plyranges) #v1.8.0
library(liftOver)  #v1.12.0
library(readxl)    #v1.4.0


colnames <- c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" )

## read in the ChIP-seq data: from Lindtner et al.(Dlx5), Sandberg et al.(Lhx6), ING/Torres lab (Meis2)

####################
#Dlx5
dlx5_e13_mm9 <- read.table("./GSM3559643_Dlx5_e13.5_BG_WT-vs.Input.no_model_peaks.narrowPeak",col.names = colnames)%>% mutate(strand="*")
dlx5_e13_mm9.GR <- as_granges(dlx5_e13_mm9)
#dlx5_e13_mm9.GR <- as_granges(dlx5_e13_mm9)

# liftOver Dlx5 to mm10
mm9tom10.chain <- import.chain("./mm9ToMm10.over.chain")
dlx5_e13_lifted_mm10.GR <- liftOver(dlx5_e13_mm9.GR, mm9tom10.chain)
#View(as.data.frame(dlx5_e13_lifted_mm10.GR))
nrow(as.data.frame(dlx5_e13_lifted_mm10.GR))
#3110 (= unfiltered peaklist) -> paper lists 2429 peaks

##### import peaks from Lindtner et al., supplementary table S3
# filter rows containing peak names of dlx5_e13.5; remove blacklisted, satellite & repeat regions
dlx_suppT3 <- read_xls("./Lindnter2019_Table_S3.xls", skip=2)
#names(dlx_suppT3)
dlx5_e13_suppT3 <- dlx_suppT3 %>%
  dplyr::filter(! is.na(Dlx5.BG.e13.5.PeakID)) %>%
  dplyr::select(1:21, 50:57) %>%
  dplyr::filter(Blacklist == "NA" & gaps == "NA" & Repeats == "NA")
# returns 2428 peaks -> paper states 2429  -> good enough 

# filter lifted peaks for peak names from supplement 
dlx5_e13_peaknames <- pull(dlx5_e13_suppT3,Dlx5.BG.e13.5.PeakID)

dlx5_e13_lifted_mm10.df <- dlx5_e13_lifted_mm10.GR %>%
  as.data.frame() %>%
  dplyr::filter(name %in% dlx5_e13_peaknames) %>%
  dplyr::select(3:5, 8, 9,7, 10:13)
#2422 difference due to lifover from mm9 to mm10

#write.table(mutate(dlx5_e13_lifted_mm10.df, strand = ".")[,1:6], "./dlx5_e13.5_lifted_mm10.bed", col.names = F,row.names = F,quote = F)
#write_rds(dlx5_e13_lifted_mm10.df, "./dlx5_e13.5_lifted_mm10.rds")


########################################
## #Lhx6
lhx6_e13_mm9 <- read.table("./GSM2281993_Lhx6_ChIPvsInput_p0.00001_peaks_mm9.txt", skip=28, header=T) %>%
  mutate(strand="*") %>%
  mutate(rel_summit_pos = abs_summit - start) %>%
  dplyr::select(1,2,3,10,6,11,8,7,9,12)
  
names(lhx6_e13_mm9) <- c("seqnames", "start", "end", "name", "score" ,"strand" ,"fold_change" ,"-log10pvalue" ,"-log10qvalue" ,"rel_summit_pos" )


lhx6_e13_supp_tbl1 <- read_xlsx("./Sandberg2016_Data_S1.xlsx", sheet=7)
lhx6_e13_peaknames <- lhx6_e13_supp_tbl1 %>% pull(name)

lhx6_e13_mm9.GR <- as_granges(lhx6_e13_mm9)
lhx6_e13_lifted_mm10.GR <- liftOver(lhx6_e13_mm9.GR, mm9tom10.chain)

lhx6_e13_lifted_mm10.df <- lhx6_e13_lifted_mm10.GR %>%
  as.data.frame() %>%
  filter(name %in%  lhx6_e13_peaknames) %>%
  dplyr::select(3:5, 8, 9,7,10:13)

#write bed and RDS files
#write.table(mutate(lhx6_e13_lifted_mm10.df, strand = ".")[,1:6], "./lhx6_e13.5_lifted_mm10.bed", col.names = F,row.names = F,quote = F)
#write_rds(lhx6_e13_lifted_mm10.df, "./lhx6_e13.5_lifted_mm10.rds")

###########################
#########Meis2
MEIS_GE_e14_mm10 <- read.table("./GE_meis2_IP_q0.01_peaks.narrowPeak", col.names=colnames)  %>% mutate(strand="*")
MEIS_GE_e14_mm10.GR <- as_granges(MEIS_GE_e14_mm10)

# blacklist regions from ENCODE
mm10_blacklist <- read.table("./mm10_blacklist_ENCFF547MET.bed", header=F, col.names=c("seqnames","start","end"))

# remove blacklisted regions
mm10_blacklist.GR <- mm10_blacklist %>% as_granges()
meis2_e14_mm10_filtered.GR <-  filter_by_non_overlaps(MEIS_GE_e14_mm10.GR,mm10_blacklist.GR)

## save as bed files

meis2_e14_mm10.df <- as.data.frame(meis2_e14_mm10_filtered.GR)  %>% dplyr::select(1,2,3,6,7,5,8:11)
#write.table(mutate(meis2_e14_mm10.df, strand = ".")[,1:6], "./meis2_e14_mm10.bed", col.names = F,row.names = F,quote = F, sep="\t")
#write_rds(meis2_e14_mm10.df, "./meis2_e14_mm10.rds")


##################### peak overlaps
###################
#################
#### create merged peaklist by merging combining all peaks in one dataset, merging overlapping peaks and assigning "merged peak name"
names(meis2_e14_mm10.df)
names(dlx5_e13_lifted_mm10.df)
names(lhx6_e13_lifted_mm10.df)


all_peaks <- bind_rows(meis2 = meis2_e14_mm10.df[,1:4],dlx5= dlx5_e13_lifted_mm10.df[,1:4],lhx6=lhx6_e13_lifted_mm10.df[,1:4],.id="ID") %>% mutate(row = 1:nrow(.))
all_peaks.GR <- as_granges(all_peaks)
#View(as.data.frame(all_peaks.GR))
merged_peaks.GR <- GenomicRanges::reduce(all_peaks.GR,with.revmap=T) # revmap = T: indicates which rows of the original table were merged to create the new ranges
length(merged_peaks.GR) #6100
merged_peaks.GR$merged_peak_ID <- paste("merged_peak_",1:6100, sep="") 
#View(as.data.frame(merged_peaks.GR))
# write_rds(merged_peaks.GR, "./merged_peaks.GR.rds")


# overlap idividual peakset ranges with merged peak list
meis2_e14_mm10.GR <- as_granges(meis2_e14_mm10.df)
dlx5_e13_lifted_mm10.GR <- as_granges(dlx5_e13_lifted_mm10.df)
lhx6_e13_lifted_mm10.GR <- as_granges(lhx6_e13_lifted_mm10.df)

#View(as.data.frame(meis2_e14_mm10.GR))

meis2_peaks_IDs <- plyranges::find_overlaps(merged_peaks.GR, meis2_e14_mm10.GR, minoverlap = 1) %>% as.data.frame() %>% 
  left_join(dplyr::select(as.data.frame(meis2_e14_mm10.GR), seqnames_meis2 = seqnames, start_meis2 = start, end_meis2 = end, width_meis2=width, name),by="name") %>% dplyr::select(-c(1:5))
dlx5_peaks_IDs <- plyranges::find_overlaps(merged_peaks.GR, dlx5_e13_lifted_mm10.GR, minoverlap = 1)  %>% as.data.frame() %>% 
  left_join(dplyr::select(as.data.frame(dlx5_e13_lifted_mm10.GR), seqnames_dlx5 = seqnames, start_dlx5 = start, end_dlx5 = end, width_dlx5=width, name),by="name")  %>% dplyr::select(-c(1:5))
lhx6_peaks_IDs <- plyranges::find_overlaps(merged_peaks.GR, lhx6_e13_lifted_mm10.GR, minoverlap = 1)  %>% as.data.frame() %>%
  left_join(dplyr::select(as.data.frame(lhx6_e13_lifted_mm10.GR), seqnames_lhx6 = seqnames, start_lhx6 = start, end_lhx6 = end, width_lhx6=width, name),by="name") %>% dplyr::select(-c(1:5))


length(pull(meis2_peaks_IDs, merged_peak_ID))
length(unique(pull(meis2_peaks_IDs, merged_peak_ID))) #5 peaks are duplicates
a<-meis2_peaks_IDs$merged_peak_ID[duplicated((pull(meis2_peaks_IDs, merged_peak_ID)))]
filter(meis2_peaks_IDs, merged_peak_ID %in% a)

length(pull(dlx5_peaks_IDs, merged_peak_ID))
length(unique(pull(dlx5_peaks_IDs, merged_peak_ID))) #1 peak is duplicates
b<-dlx5_peaks_IDs$merged_peak_ID[duplicated((pull(dlx5_peaks_IDs, merged_peak_ID)))]
filter(dlx5_peaks_IDs, merged_peak_ID %in% b)

length(pull(lhx6_peaks_IDs, merged_peak_ID))
length(unique(pull(lhx6_peaks_IDs, merged_peak_ID))) #1 peak is duplicates
c<-lhx6_peaks_IDs$merged_peak_ID[duplicated((pull(lhx6_peaks_IDs, merged_peak_ID)))]
filter(lhx6_peaks_IDs, merged_peak_ID %in% c)

# for dual peaks, remove the row containing the peak with lower enrichment
dual_peaks <- c("GE_meis_IP_q0.01_peak_856", "GE_meis_IP_q0.01_peak_869", "GE_meis_IP_q0.01_peak_1832", "GE_meis_IP_q0.01_peak_2667","GE_meis_IP_q0.01_peak_2686", 
                "Dlx5_e13.5_BG_WT_mrgdf_trim-vs.Input.no_model_peak_394", "Lhx6_ChIPvsInput_p0.00001_peak_1312" )

# convert back to df and join tables
merged_peaks.df <- as.data.frame(merged_peaks.GR) %>%
  left_join(meis2_peaks_IDs, by = "merged_peak_ID") %>%
  left_join(dlx5_peaks_IDs,  by = "merged_peak_ID") %>%
  left_join(lhx6_peaks_IDs,  by = "merged_peak_ID") 


names(merged_peaks.df)
names(merged_peaks.df) <- c("seqnames_merged", "start_merged", "end_merged", "width_merged", "strand_merged", "revmap.merged", "peakname_merged",
                            "revmap_meis2", "peakname_meis2", "score_meis2", "fold_change_meis2", "neg_log10pvalue_meis2", "neg_log10_qvalue_meis2", "rel_summit_pos_meis2",
                            "seqnames_meis2", "start_meis2", "end_meis2", "width_meis2",
                            "revmap_dlx5", "peakname_dlx5", "score_dlx5", "fold_change_dlx5", "neg_log10pvalue_dlx5", "neg_log10_qvalue_dlx5", "rel_summit_pos_dlx5",
                            "seqnames_dlx5", "start_dlx5", "end_dlx5", "width_dlx5", 
                            "revmap_lhx6", "peakname_lhx6", "score_lhx6", "fold_change_lhx6", "neg_log10pvalue_lhx6", "neg_log10_qvalue_lhx6", "rel_summit_pos_lhx6",
                            "seqnames_lhx6", "start_lhx6", "end_lhx6", "width_lhx6")
#View(as.data.frame(merged_peaks.GR))

merged_peaks.df2 <- merged_peaks.df %>%
  filter(! peakname_meis2 %in% dual_peaks , ! peakname_dlx5 %in% dual_peaks , ! peakname_lhx6 %in% dual_peaks) %>%
  mutate(binding_type = ifelse(! is.na(peakname_meis2) & ! is.na(peakname_dlx5) & ! is.na(peakname_lhx6), "meis2_dlx5_lhx6", 
                               ifelse(! is.na(peakname_meis2) & ! is.na(peakname_dlx5) &  is.na(peakname_lhx6), "meis2_dlx5",
                                      ifelse(! is.na(peakname_meis2) &  is.na(peakname_dlx5) & !  is.na(peakname_lhx6), "meis2_lhx6",
                                             ifelse(! is.na(peakname_meis2) &  is.na(peakname_dlx5) &  is.na(peakname_lhx6), "meis2",
                                                    ifelse( is.na(peakname_meis2) & ! is.na(peakname_dlx5) &  ! is.na(peakname_lhx6), "dlx5_lhx6",
                                                            ifelse( is.na(peakname_meis2) & ! is.na(peakname_dlx5) &  is.na(peakname_lhx6), "dlx5",
                                                                    ifelse( is.na(peakname_meis2) &  is.na(peakname_dlx5) &  ! is.na(peakname_lhx6), "lhx6","error")))))))) %>%
  mutate(meis2_peak = ifelse(! is.na(peakname_meis2), "y", "n"), 
         dlx5_peak = ifelse(! is.na(peakname_dlx5), "y", "n"),
         lhx6_peak = ifelse(! is.na(peakname_lhx6), "y", "n")) %>%
  dplyr::select(1:4,7, 41:44,15:18, 9:14, 26:29, 20:25, 37:40, 31:36) 

########### overlap with Vista Enhancers
# read in bed file of enhancers (hg19_lifted_to_mm10)
# run plyranges::find_overlaps() with merged peak list to append merged peak ID
# left_join the resulting table by merged peak ID
# check if everything worked properly

######### list extracted from Vista homepage (from fasta output)
# human (hg19) amd mouse (mm9) enhancers were downloaded separately and lifted to mm10 using the UCSC liftOver webtool 

vista_full_mm10 <- read_rds("./vista_enhancers_all_mm10.rds")
table(vista_full_mm10$species)
# human mouse 
# 1931    1371  
# on vista website (14.7.2022): human - 1942, mouse - 1372  

vista_full_mm10.GR <- as_granges(vista_full_mm10, seqnames = chr)
#View(as.data.frame(vista_full_mm10.GR))

vista_full_mm10.GR2 <- find_overlaps(vista_full_mm10.GR, merged_peaks.GR)
#View(as.data.frame(vista_full_mm10.GR2))
vista_full_mm10_2 <- as.data.frame(vista_full_mm10.GR2) %>%
  dplyr::select(6,20)
x <- pull(filter(vista_full_mm10_2,duplicated(vista_full_mm10_2$merged_peak_ID)), merged_peak_ID)
filter(vista_full_mm10_2, merged_peak_ID %in% x) %>% dplyr::arrange(merged_peak_ID)
# some peaks overlap multiple vista enhancers -> these enhancers are already overlapping
# add to column by hand


vista_full_mm10_3 <- vista_full_mm10_2 %>%
  mutate(vista_enhancer = ifelse(merged_peak_ID == "merged_peak_601", "hs1354;hs1540;hs998", 
                                 ifelse(merged_peak_ID == "merged_peak_951", "mm876;mm999",
                                        ifelse(merged_peak_ID == "merged_peak_1220", "mm1676;hs1466", 
                                               ifelse(merged_peak_ID == "merged_peak_1436", "mm1605;hs2234",
                                                      ifelse(merged_peak_ID == "merged_peak_1517", "mm14;hs1323 ",
                                                          ifelse(merged_peak_ID == "merged_peak_1599", "mm1660;hs268",
                                                                 ifelse(merged_peak_ID == "merged_peak_2123", "mm17;hs1403",
                                                                        ifelse(merged_peak_ID == "merged_peak_2321", "mm817;hs1188",
                                                                               ifelse(merged_peak_ID == "merged_peak_2630", "mm15;hs1399" ,
                                                                                      ifelse(merged_peak_ID == "merged_peak_3540","mm1541;hs1681",
                                                                                             ifelse(merged_peak_ID == "merged_peak_3869","hs174;hs322", 
                                                                                                    ifelse(merged_peak_ID == "merged_peak_4044", "hs1597;hs1717",
                                                                                                           ifelse(merged_peak_ID == "merged_peak_4244", "mm1684;hs1388",
                                                                                                                  ifelse(merged_peak_ID == "merged_peak_4330", "hs2094;mm1532",
                                                                                                                         ifelse(merged_peak_ID == "merged_peak_4922", "mm123;hs1652", 
                                                                                                                                ifelse(merged_peak_ID == "merged_peak_5814", "mm1608;hs872",
                                                                                                                                       ifelse(merged_peak_ID == "merged_peak_5815", "mm1608;hs872",
                                                                                                                                              ifelse(merged_peak_ID == "merged_peak_5842", "mm845;mm870", vista_ID))))))))))))))))))) %>%
                                 
                                 
dplyr::distinct()


## ad information on enhancer activity

vista_full_activity_mm10 <- vista_full_mm10 %>%
  unite(tissue, tissue1:tissue10, sep=";") %>%
  mutate(enhancer_activity = ifelse(grepl("forebrain*", tissue) & activity == "positive","forebrain_positive", 
                                    ifelse(! grepl("forebrain*", tissue) & activity == "positive","other_positive", "negative"))) %>%
  dplyr::select(4, 8)

vista_full_mm10_4 <- vista_full_mm10_3 %>%
  tidyr::separate(col=vista_enhancer, into = c("vista1", "vista2", "vista3"), sep=";", remove=F) %>%
  left_join(vista_full_activity_mm10, by = c("vista1" = "vista_ID")) %>%
  left_join(vista_full_activity_mm10, by = c("vista2" = "vista_ID")) %>%
  left_join(vista_full_activity_mm10, by = c("vista3" = "vista_ID")) %>%
  tidyr::unite(vista_activity, enhancer_activity.x, enhancer_activity.y, enhancer_activity, sep=";", na.rm=T) %>%
  dplyr::select(-vista1,-vista2, -vista3, -vista_ID) %>%
  distinct()
  

merged_peaks.df3 <- merged_peaks.df2 %>%
  left_join(vista_full_mm10_4, by=c("peakname_merged" = "merged_peak_ID")) %>%
  dplyr::select(names(.)[c(1:9,40,41,10:39)]) 

#write_rds(merged_peaks.df3, "./merged_peaks_meis2_dlx5_lhx6_df3.RDS")
#merged_peaks.df3 <- read_rds("./merged_peaks_meis2_dlx5_lhx6_df3.RDS")



########################### overlap with promoters & enhancers
#
##### EPDnew v003 (downloaded from EPD server "Promoters from EPDnew mouse version 003")
# extend promoter region to +- 5kb from TSS

EPDnewV3 <- read.table("./EPD_new_v3_mm10.bed", header=F) 
names(EPDnewV3) <- c("seqnames","start","end", "name", "score", "strand", "thickStart", "thickEnd")

EPDnewV3_2 <- EPDnewV3 %>% 
  separate(col = name, into = c("name", "promnr"),  sep="_") %>%
  mutate(tssStart = ifelse(strand == "+", start + 49 -5000, start + 10 -5000)) %>%
  mutate(tssEnd = ifelse(strand == "+", start + 49 +5000, start + 10 +5000)) %>%
  mutate(tssWidth = tssEnd - tssStart)

EPDnewV3.GR <- as_granges(EPDnewV3_2, start = tssStart, end=tssEnd )

## overlap

EPDnewV3.GR2 <- find_overlaps(EPDnewV3.GR, merged_peaks.GR)
EPDnewV3.df <- as.data.frame(EPDnewV3.GR2)


length(EPDnewV3.df$merged_peak_ID) # 1426
length(unique(EPDnewV3.df$merged_peak_ID)) #1004

##### clean up list, to get only one peak name per row
# some peaks bind multiple promoters of the same gene -> remove first
# some peaks bind promoters of multiple genes -> report genes in multiple columns
# separate duplicate overlaps into separate tables, then left_join in 3 steps to merged_peaks.df3; keep only the gene name
# combine gene names in one column, separated by semicolon



EPDnewV3.df_b <- EPDnewV3.df %>% # remove duplicate promoters per peak
  distinct(merged_peak_ID, name)

EPDnewV3.df_1 <- filter(EPDnewV3.df_b,  ! duplicated(merged_peak_ID)) %>% dplyr::select(name, merged_peak_ID)
EPDnewV3.df_2 <- filter(EPDnewV3.df_b, duplicated(merged_peak_ID)) %>% dplyr::select(name, merged_peak_ID)
EPDnewV3.df_3 <- filter(EPDnewV3.df_2, duplicated(merged_peak_ID))



merged_peaks.df4 <- merged_peaks.df3 %>% 
  left_join(EPDnewV3.df_1, by = c("peakname_merged"= "merged_peak_ID")) %>%
  left_join(EPDnewV3.df_2, by = c("peakname_merged"= "merged_peak_ID")) %>%
  left_join(EPDnewV3.df_3, by = c("peakname_merged"= "merged_peak_ID")) %>%
  tidyr::unite(EPDnewV3_promoter, name.x, name.y, name, remove=T, na.rm =T, sep="; ") %>%
  mutate_at(vars(EPDnewV3_promoter), na_if,"") %>%
  mutate(promoter_overlap = ifelse(is.na(EPDnewV3_promoter), "n", "y")) %>%
  distinct(peakname_merged, .keep_all=T)  %>%
  dplyr::select(c(1:11, 42:43, 12:41))

#write_rds(merged_peaks.df4, "./merged_peaks_meis2_dlx5_lhx6_df4.RDS")



############### promoter-interacting enhancers
##### from Gorkin et al. 2020, Nature
# using all promoter-enhancer interactions

prom_enh_interaction_rep1 <- read_excel("./Gorkin_2020_Nature_Supp_table_5-11_41586_2020_2093_MOESM6_ESM.xlsx",
                                        skip=1, sheet = 4 )
prom_enh_interaction_rep2 <- read_excel("./Gorkin_2020_Nature_Supp_table_5-11_41586_2020_2093_MOESM6_ESM.xlsx",
                                        skip=1, sheet = 5 )

# combine replicates
prom_enh_interaction_all <- bind_rows(prom_enh_interaction_rep1, prom_enh_interaction_rep2, .id = "replicate") %>%
  dplyr::arrange(chrom, start) %>%
  mutate(enhancer_name = paste("e", 1:nrow(.), sep="")) 

# remove duplicates in enhancer coordinates + gene 
prom_enh_interaction_all_II <- prom_enh_interaction_all %>% 
  distinct(chrom, start, end, ensembl, geneName = symbol,  .keep_all=T) %>% 
  mutate(interaction_name = paste(symbol,enhancer_name, sep="/"))
#43555 interactions left

prom_enh_interaction_all.bed <- prom_enh_interaction_all_II %>%
  dplyr::select( chrom, start, end, interaction_name) %>%
  mutate(end = format(end, scientific=F), start= format(start, scientific = F)) %>%
  mutate(start= as.character(start), end= as.character(end))

# convert to ranged object
prom_enh_interaction_all.GR <- as_granges(prom_enh_interaction_all_II, seqnames = chrom)
#View(as.data.frame(prom_enh_interaction_all.GR))

### overlap
enhancersAll.GR2 <- find_overlaps(prom_enh_interaction_all.GR, merged_peaks.GR)
#View(as.data.frame(enhancersEnc3.GR2))
enhancersAll.df <- as.data.frame(enhancersAll.GR2)

length(enhancersAll.df$merged_peak_ID) # 1311
length(unique(enhancersAll.df$merged_peak_ID)) #1193


# 118 are listed twice -> multiple enhancers per gene reported
# separate duplicate overlaps into separate tables, then left_join in 2 steps to merged_peaks.df3; keep only the promoter ID
enhancersAll.df_1 <- filter(enhancersAll.df,  ! duplicated(merged_peak_ID)) %>% dplyr::select(geneName, merged_peak_ID)
enhancersAll.df_2 <- filter(enhancersAll.df,   duplicated(merged_peak_ID) ) %>% dplyr::select(geneName, merged_peak_ID)


merged_peaks.df5 <- merged_peaks.df4 %>% 
  left_join(enhancersAll.df_1, by = c("peakname_merged"= "merged_peak_ID")) %>%
  left_join(enhancersAll.df_2, by = c("peakname_merged"= "merged_peak_ID")) %>%
  tidyr::unite(enhancer_interacting_genes, geneName.x, geneName.y, remove=T, na.rm =T, sep="; ") %>%
  mutate_at(vars(enhancer_interacting_genes), na_if,"") %>%
  mutate(enhancer_overlap = ifelse(is.na(enhancer_interacting_genes), "n", "y")) %>%
  dplyr::select(c(1:13, 44:45, 14:43)) %>%
  dplyr::select(-c(4,6, 13,15)) %>%
  rename("chr_merged"="seqnames_merged","TSS_overlap" = "EPDnewV3_promoter", "vista_overlap" = "vista_enhancer", "chr_meis2" = "seqnames_meis2", "chr_dlx5" = "seqnames_dlx5", "chr_lhx6" = "seqnames_lhx6")


write.table(merged_peaks.df5, "./Table_S5_ChIP-seq.txt",col.names = T,row.names = F,quote = F, sep="\t")
write_rds(merged_peaks.df5, "./Table_S5_ChIP-seq.RDS")

