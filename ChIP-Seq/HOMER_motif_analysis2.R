#########################
# Visualization of HOMER output
# Fig 3 c,d, S3 c



### density plot
library(tidyverse)

# import output from HOMER
meis2_all_density <- read.table("./input_data/motif_density_meis2peaks.txt", skip=1,header=F, sep="\t")

columns <- c("distance" , "hexa_total_sites"	,"hexa_pos_sites", "hexa_neg_sites",	"deca_total_sites",	"deca_pos_sites", "deca_sites",
             	"A_frequency", "C_frequency", "G_frequency",	"T_frequency")

names(meis2_all_density) <- columns


meis2_motif_enrichment2 <- ggplot(meis2_all_density, aes(x=distance)) +
  geom_smooth(aes(y=hexa_total_sites), se=F, method="gam", color="#F47E70")+
  geom_smooth(aes(y=deca_total_sites), se=F, method="gam", color="#8FD1C4") +
  theme_classic()+
  labs(x="distance to peak summit [bp]", y="relative_motif_enrichment") +
  theme(axis.text = element_text(size=7), axis.title=element_text(size=8))
meis2_motif_enrichment2




### motif enrichment plot
library(tidyverse)
library(RColorBrewer)

#### results from HOMER
## include all DLX motifs

############ all MEIS2 peaks

homer_enrichment_meis2_all <- read.table("./output/motif_analysis/HOMER/MEIS2_200bp_new/knownResults.txt", 
                                   header=F,row.names=NULL,sep="\t", skip=1) 

homer_enrichment_meis2_all[1,4]
names(homer_enrichment_meis2_all) <- c("motif_name", "consensus", "p_value", "log_pvalue", "q_value_Benjamini", "no_of_targets_seqs_with_motif", "percent_of_target_seqs_with_motif", "no_of_background_seqs_with_motif", "percent_of_background_seqs_with_motif")

#convert ln value to log10
homer_enrichment_meis2_all <- homer_enrichment_meis2_all %>%
  mutate(n_log10_pval = -round((log_pvalue / 2.303 ),1))
  
# "log p-value" column is natural logarithm!

homer_enrichment_meis2_all2 <- homer_enrichment_meis2_all %>%
    mutate(percent_of_target_seqs_with_motif = as.numeric(gsub("%","",percent_of_target_seqs_with_motif)), percent_of_background_seqs_with_motif = as.numeric(gsub("%","",percent_of_background_seqs_with_motif))) %>%
    mutate(fold_enrichment = percent_of_target_seqs_with_motif/percent_of_background_seqs_with_motif) %>%
    tidyr::separate(motif_name,c( "motif_short","motif_info1", "motif_info2"), sep="\\(", remove=F) %>%
    dplyr::select(-c("motif_info1","motif_info2")) %>%
    mutate(rownumber = seq(1:nrow(.)))
  
  
  
  
# manually selected sequences from HOMER knownResults, to avoid redundancy
rowselectionB <- c(1,4, 7,8,12, 16,20,13,69,15 )
homer_enrichment_all_selectedB <- homer_enrichment_meis2_all2 %>%
  filter(rownumber %in% rowselectionB) %>%
  mutate() %>%
  arrange(factor(rownumber, levels = rev(rowselectionB))) %>%
  mutate(motif_short = as.factor(motif_short)) %>%
  mutate(motif_short = factor(motif_short, 
                              levels=c("Nanog" , "Sp5" ,   "Tbx5"  , "Oct6" ,  "Isl1"  , "Nkx6.1" ,"Lhx2"  , "Dlx3" ,  "Pbx3"  , "Meis1" )))

# homer_enrichment_plot 
homer_enrichment_plot_selectionB <- ggplot(homer_enrichment_all_selectedB, aes(x= motif_short))+
  theme_bw()+
  geom_col(aes(y=percent_of_background_seqs_with_motif),fill="orange",alpha=0.5, color="red")+
  geom_col(aes(y=percent_of_target_seqs_with_motif),fill="lightgrey", alpha=0.5, color="black")+
  theme(axis.line = element_line(colour = "black"),panel.grid = element_blank(), panel.border =  element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(face="bold", size=9), axis.text.x= element_text(size=9), axis.title.x=element_text(size=10))+
  geom_text(aes(y=percent_of_target_seqs_with_motif, label=paste("(",round(fold_enrichment,2),")", sep=""),hjust=-0.1),size=9*0.36)+
  ylim(0,100)+
  labs(x=NULL, y="percent of peaks with motif")+
  coord_flip()
homer_enrichment_plot_selectionB

ggsave(plot=homer_enrichment_plot_selectionB, 
path= ".",
filename ="enrichment_plot_selectionB.png",
device="png", width = 60, height=60, units="mm")

ggsave(plot=homer_enrichment_plot_selectionB, 
path= ".",
filename ="enrichment_plot_selectionB.pdf",
device="png", width = 60, height=60, units="mm", dpi=300)

###########################################
############ divide by enhancers and promoters


############# MEIS2 peaks in promoters 


homer_enrichment_promoter <- read.table("/MEIS2_200bp_promoter_peaks_new/knownResults.txt", 
                                        header=F,row.names=NULL,sep="\t", skip=1)
names(homer_enrichment_promoter) <- c("motif_name", "consensus", "p_value", "log_pvalue", "q_value_Benjamini", "no_of_targets_seqs_with_motif", "percent_of_target_seqs_with_motif", "no_of_background_seqs_with_motif", "percent_of_background_seqs_with_motif")

homer_enrichment_promoter2 <- homer_enrichment_promoter %>%
  mutate(percent_of_target_seqs_with_motif = as.numeric(gsub("%","",percent_of_target_seqs_with_motif)), percent_of_background_seqs_with_motif = as.numeric(gsub("%","",percent_of_background_seqs_with_motif))) %>%
  mutate(fold_enrichment = percent_of_target_seqs_with_motif/percent_of_background_seqs_with_motif) %>%
  tidyr::separate(motif_name,c( "motif_short","motif_info1", "motif_info2"), sep="\\(", remove=F) %>%
  dplyr::select(-c("motif_info1","motif_info2")) %>%
  mutate(rownumber = seq(1:nrow(.)))

# log p-value is natural logarithm!

rowselection_promoter <- c(1,3,7,8,9,13,15,17,18,27)
homer_enrichment_promoter_selected <- homer_enrichment_promoter2 %>%
  filter(rownumber %in% rowselection_promoter) %>%
  mutate()



##############  MEIS2 peaks in ennhancers 

homer_enrichment_enhancer <- read.table("/HOMER/MEIS2_200bp_enhancer_peaks_new/knownResults.txt", 
                                        header=F,row.names=NULL,sep="\t", skip=1)
names(homer_enrichment_enhancer) <- c("motif_name", "consensus", "p_value", "log_pvalue", "q_value_Benjamini", "no_of_targets_seqs_with_motif", "percent_of_target_seqs_with_motif", "no_of_background_seqs_with_motif", "percent_of_background_seqs_with_motif")

# log p-value is natural logarithm!

homer_enrichment_enhancer2 <- homer_enrichment_enhancer %>%
  mutate(percent_of_target_seqs_with_motif = as.numeric(gsub("%","",percent_of_target_seqs_with_motif)), percent_of_background_seqs_with_motif = as.numeric(gsub("%","",percent_of_background_seqs_with_motif))) %>%
  mutate(fold_enrichment = percent_of_target_seqs_with_motif/percent_of_background_seqs_with_motif) %>%
  tidyr::separate(motif_name,c( "motif_short","motif_info1", "motif_info2"), sep="\\(", remove=F) %>%
  dplyr::select(-c("motif_info1","motif_info2")) %>%
  mutate(rownumber = seq(1:nrow(.)))


# enhancers-manually selected sequences
rowselection_enhancer <- c(2,4,6,7,10,13,15,16,17,19)
homer_enrichment_enhancer_selected <- homer_enrichment_enhancer2 %>%
  filter(rownumber %in% rowselection_enhancer) %>%
  mutate()


##################################
##### grouped barplot for easier comparison between enhancers and promoters
# 1) rename columns
# 2) bind rows to generate one table, add ID


library(RColorBrewer)
display.brewer.pal(12, "Paired")
brewer.pal(12, "Paired")
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
#"#33A02C", "#FF7F00"
homer_enrichment_promoter_enhancer_combined <- homer_enrichment_promoter2 %>%
  bind_rows(homer_enrichment_enhancer2,.id = "ID") %>%
  mutate(ID = ifelse(ID == "1", "promoter", "enhancer"))


# select motifs by name

motifs <- c("Meis1", "Pbx3", "Dlx3", "Lhx2", "Nkx6.1", "Isl1", "Oct6", "Tbx5", "Sp5", "Nanog")


homer_enrichment_promoter_enhancer_selected <- homer_enrichment_promoter_enhancer_combined %>%
  filter(motif_short %in% motifs)

# convert to wide
homer_enrichment_promoter_enhancer_wide <- homer_enrichment_promoter_enhancer_selected %>%
  dplyr::select(c(1:4,12)) %>%
  pivot_wider(names_from = ID, values_from = fold_enrichment) %>%
  mutate(promoterVSenhancer = enhancer/promoter)



homer_enrichment_plot_promoter_enhancer_selected <- ggplot(homer_enrichment_promoter_enhancer_selected, 
                                                           aes(x= factor(motif_short, levels=rev(motifs)), fill= ID))+
  theme_bw()+
  geom_col(aes(y=percent_of_background_seqs_with_motif ),width=0.8,alpha=1, position="dodge", color="red")+
  geom_col(aes(y=percent_of_target_seqs_with_motif),width = 0.8, alpha=0.5, position="dodge", color="black")+
  geom_text(aes(y=percent_of_target_seqs_with_motif, label=paste("(",round(fold_enrichment,2),")", sep="")),size=8*0.36, hjust=-0.1, 
            position=position_dodge(width=0.8))+
  ylim(0,100)+
  labs(x=NULL, y="% of peaks with motif", fill="peak category")+
  scale_fill_manual(values = c("#FDBF6F", "#CAB2D6"))+
  theme(axis.line = element_line(colour = "black"),panel.grid = element_blank(), panel.border =  element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(face="bold", size=9), axis.text.x= element_text(size=9), axis.title.x=element_text(size=9),
        legend.position = "right", legend.direction = "vertical")+
  coord_flip()
homer_enrichment_plot_promoter_enhancer_selected


ggsave(plot=homer_enrichment_plot_promoter_enhancer_selected, 
       path= "",
       filename ="enrichment_plot_enhancer_and_promoter_comparison_new5_12.png",
       device="png", width = 100, height=100, units="mm")

ggsave(plot=homer_enrichment_plot_promoter_enhancer_selected, 
       path= "",
       filename ="enrichment_plot_enhancer_and_promoter_comparison_new5_12.pdf",
       device="pdf", width = 100, height=100, units="mm", dpi=300)