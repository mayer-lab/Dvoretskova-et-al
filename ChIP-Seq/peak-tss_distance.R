## Fig 3a


library(tidyverse) #v1.3.1
library(plyranges) #v1.18.0

# import peak table
merged_peaks_meis2.df <- read.table("./Table_S5_ChIP-seq.txt",header= T, sep="\t") %>%
  filter(meis2_peak == "y") %>%
  dplyr::select(12:17) %>%
  mutate(strand = "+") %>%
  dplyr::rename("seqnames" = "chr_meis2", start = "start_meis2", end = "end_meis2" )

merged_peaks_meis2.GR <- merged_peaks_meis2.df %>%
  as_granges()

# import promoter coordinates from Ekaryotic Promoter Database
EPDnewV3 <- read.table("./EPD_new_v3_mm10.bed", header=F, sep="\t")
names(EPDnewV3) <- c("seqnames","start","end", "name", "score", "strand", "thickStart", "thickEnd")
# extract TSS coordinates
EPDnewV3_TSS <- EPDnewV3 %>%
  mutate(tssStart = ifelse(strand == "+", start + 49 , start + 10 )) %>%
  mutate(tssEnd = ifelse(strand == "+", start + 49 , start + 10 )) %>%
  separate(col = name, into = c("geneName", "promnr"),  sep="_", remove=F) %>%
  rename("promoter_name" = "name") %>%
  dplyr::select(seqnames,start = tssStart, end =tssEnd, geneName, promoter_name, score, strand) 

EPDnew.GR <- as_granges(EPDnewV3_TSS)

meis2_peaks_tss_dist_nearest <- join_nearest(merged_peaks_meis2.GR, EPDnew.GR,distance=T ) %>% as.data.frame()


# join_nearest does not consider up- or downstream
# use join_nearest_left/right to assign all up- or downstream genes, not just the ones on the same strand

meis2_peaks_tss_dist_right.GR <- join_nearest_right(merged_peaks_meis2.GR, EPDnew.GR,distance=T )
meis2_peaks_tss_dist_right <- as.data.frame(meis2_peaks_tss_dist_right.GR) %>% dplyr::rename("distance_downstream" = "distance")

meis2_peaks_tss_dist_left.GR <- join_nearest_left(merged_peaks_meis2.GR, EPDnew.GR,distance=T )
meis2_peaks_tss_dist_left <- as.data.frame(meis2_peaks_tss_dist_left.GR) %>% dplyr::rename("distance_upstream" = "distance") %>% mutate(distance_upstream = distance_upstream*-1)

meis2_peaks_tss_dist_left2 <- meis2_peaks_tss_dist_left %>% select(c(peakname_meis2,geneName, promoter_name,distance_upstream))
meis2_peaks_tss_dist_right2 <- meis2_peaks_tss_dist_right %>% select(c(peakname_meis2,geneName, promoter_name,distance_downstream))


# add both upstream and downstream genes (left & right); 
# add correct sign depending on gene strand and direction

meis2_peaks_tss_dist_LR <- meis2_peaks_tss_dist_nearest %>%
  select(c(peakname_meis2,geneName, promoter_name,distance)) %>%
  left_join(meis2_peaks_tss_dist_left2,  by="peakname_meis2") %>%
  left_join(meis2_peaks_tss_dist_right2,  by="peakname_meis2") %>%
  replace_na(list(distance_upstream =0, distance_downstream=0)) %>%
  mutate(nearest_promoter_distance = ifelse(abs(distance_upstream) > distance_downstream, distance_upstream, distance_downstream)) %>%
  left_join(select(EPDnew, promoter_name, strand), by=c( "promoter_name.x"="promoter_name"))  %>%
  mutate(nearest_promoter_distance2 = ifelse(nearest_promoter_distance > 0 & strand =="+", distance*-1,
                                             ifelse(nearest_promoter_distance > 0 & strand =="-", distance,
                                                    ifelse(nearest_promoter_distance < 0 & strand =="+", distance,
                                                           ifelse(nearest_promoter_distance < 0 & strand =="-", distance*-1, 0)))))


# assign to distance categories for plotting

meis2_peaks_tss_dist_LR_2 <- meis2_peaks_tss_dist_LR %>%
  mutate(dist_cat = ifelse(nearest_promoter_distance2 < -500000, "< -500",
                           ifelse(nearest_promoter_distance2 >= -500000 & nearest_promoter_distance2 < -50000, "-500 to -50",
                                  ifelse(nearest_promoter_distance2 >= -50000 & nearest_promoter_distance2 < -5000, "-50 to -5",
                                         ifelse(nearest_promoter_distance2 >= -5000 & nearest_promoter_distance2 < 0, "-5 to 0",
                                                ifelse(nearest_promoter_distance2 == 0, "TSS overlap",
                                                       ifelse(nearest_promoter_distance2 > 0 & nearest_promoter_distance2 <= 5000, "0 to 5",
                                                             ifelse(nearest_promoter_distance2 > 5000 & nearest_promoter_distance2 <= 50000, "5 to 50",
                                                                   ifelse(nearest_promoter_distance2 > 50000 & nearest_promoter_distance2 <= 500000, "50 to 500", "> 500"))))))))) %>%
  mutate(dist_cat2 = as.factor(dist_cat)) %>%
  mutate(dist_cat3 = factor(dist_cat2, levels = c("< -500","-500 to -50", "-50 to -5", "-5 to 0" , "TSS overlap" ,"0 to 5" ,"5 to 50" , "50 to 500","> 500" ))) %>%
  select(c(1,2,3,12,13,16))%>%
  rename("geneName" = "geneName.x", "promoteName" = "promoter_name.x", "gene_strand" = "strand", "distance_cat" = "dist_cat3", "peak_promoter_distance" ="nearest_promoter_distance2")

meis2_peaks_tss_dist_LR_2  

#write.table(meis2_peaks_tss_dist_LR_2, "/meis2_peak_gene_distance.txt", col.names = T, row.names = F, quote=F, sep="\t")
#write_rds(meis2_peaks_tss_dist_LR_2, "/meis2_peak_gene_distance.rds")


# barplot

peak_gene_distance_plot <- ggplot(meis2_peaks_tss_dist_LR_2 , aes(distance_cat))+
  geom_bar(fill="#3A53A4", color="black",width=0.7) +
  geom_text(aes(label=after_stat(count)), stat="count", vjust=-0.4, size = 7*0.36)+
  scale_y_continuous(expand=c(0,0), limits=c(0,1300))+
  theme_bw()+
  theme(panel.grid= element_blank(), axis.text.x = element_text(size=7, angle=45, hjust=1),
        axis.text.y=element_text(size=7),
        axis.title.x = element_text(size=9, margin=margin(t=5)),
        axis.title.y = element_text(size=9, margin=margin(r=5))) +
  theme(title=element_text(family = "Arial"))+
  labs(y = "Peak-gene associations", x = "distance to nearest TSS [kb]", title="MEIS1/2 binding sites in GE at E14.5")
peak_gene_distance_plot

ggsave(plot=peak_gene_distance_plot, path="./", filename="peak_tss_distance.eps", width= 65, height=60, units = "mm", device = "eps")





