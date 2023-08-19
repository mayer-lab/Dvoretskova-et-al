# Fig 3e, S3d
####### plot enhancer & promoter overlaps

library(tidyverse)#v1.3.1
library(plyranges)#v1.8.0
library(UpSetR)#v1.4.0
library(ComplexUpset)# v1.3.3

merged_peaks.df5 <- read.table( "./Table_S5_ChIP-seq.txt", header=T, sep="\t")

####################################################################################################################
################################
########  Gorkin et al. Enhancers/EPD promoters
# convert data frame
# meis2_peak, dlx5_peak and lhx6_peak to logical instead of character

merged_peaks_plot_enh_prom <- merged_peaks.df5 %>%
  filter(meis2_peak == "y" | dlx5_peak == "y") %>%
  mutate(dplyr::across(.cols = c(meis2_peak, dlx5_peak), .fns=~if_else(. == "y", TRUE, FALSE))) %>%
  mutate(reg_element = ifelse(! is.na(TSS_overlap) & ! is.na(enhancer_interacting_genes), "enhancer & promoter", 
                              ifelse(! is.na(TSS_overlap) & is.na(enhancer_interacting_genes), "promoter",
                                     ifelse(is.na(TSS_overlap) & !is.na(enhancer_interacting_genes), "enhancer", "none")))) %>%
  dplyr::select(-c(32:41)) %>%
  mutate(reg_element = as.factor(reg_element)) %>%
  mutate(reg_element = factor(reg_element,levels=c("enhancer", "promoter","enhancer & promoter" ,"none")))  %>%
  rename("MEIS2" = "meis2_peak", "DLX5" = "dlx5_peak")

enh_prom_cols <- c("MEIS2", "DLX5")  

complexUpset_enhancer_promoter_meis2_dlx5_Helvetica <- ComplexUpset::upset(merged_peaks_plot_enh_prom,  enh_prom_cols,
                                                                         name="binding site",
                                                                         base_annotations = list("Intersection size" = intersection_size(text=list(size=8*0.36))),
                                                                         themes=upset_modify_themes(list("overall_sizes"=theme(axis.text=element_text(size=8, angle = 90, family = "Helvetica"), axis.title=element_text(size=9, family = "Helvetica"), panel.grid.minor = element_blank()), # set sizes, left panel
                                                                                                         "Intersection size"=theme(axis.text=element_text(size=8, family = "Helvetica"),axis.title=element_text(size=9, family = "Helvetica")),
                                                                                                         "intersections_matrix"=theme(axis.text=element_text(size=8, family = "Helvetica"),axis.title=element_text(size=9, family = "Helvetica"),panel.grid = element_blank()))), # barplot of intersections
                                                                         annotations=list('peak category'=(
                                                                           ggplot(mapping=aes(x=intersection, fill= reg_element))+
                                                                             geom_bar(stat="count", position = "fill")+
                                                                             scale_y_continuous(labels = scales::percent_format())+
                                                                             scale_fill_manual(values= c("#FDBF6F", "#CAB2D6","#98DD99", "black"),breaks=c("enhancer","promoter","enhancer & promoter","none"))+
                                                                             labs(y= "reg. Element fraction", fill = "regulatory element\nin binding site")+
                                                                             theme(axis.title.y = element_text(size=9, family = "Helvetica"), axis.text = element_text(size=8, family = "Helvetica"),legend.text = element_text(size=8, family = "Helvetica"), legend.title=element_text(size=9, family = "Helvetica"), legend.margin = margin(l=1) ))),
                                                                         set_sizes=upset_set_size(position="right"), guides="over", width_ratio = 0.5)
complexUpset_enhancer_promoter_meis2_dlx5_Helvetica  

####### statistical comparisons

#1) are types of reg. elements independent from class of MEIS1/2-DLX5 binding?

meis2_dlx5_binding <- merged_peaks_plot_enh_prom %>%
  mutate(binding_type = ifelse(MEIS2 == TRUE & DLX5 == FALSE, "meis2", 
                               ifelse(MEIS2 == TRUE & DLX5 == TRUE, "meis2_dlx5",
                                      ifelse(MEIS2==FALSE & DLX5 == TRUE, "dlx5", "error")))) %>%
  dplyr::select(peakname_merged, binding_type, reg_element)  %>%
  mutate(reg_element = ifelse(reg_element == "enhancer", "enhancer", 
                              ifelse(reg_element == "enhancer & promoter", "enhancer", "no enhancer"))) 

meis2_dlx5_binding_summary <- meis2_dlx5_binding %>%
  group_by(binding_type, reg_element) %>%
  summarise(count = dplyr::n()) %>%
  pivot_wider(names_from = reg_element, values_from = count )

meis2_dlx5_binding_summary_matrix <- meis2_dlx5_binding_summary %>%
  column_to_rownames(., var = "binding_type") %>%
  as.matrix()
chisq.test(meis2_dlx5_binding_summary_matrix)

#Pearson's Chi-squared test
#data:  meis2_dlx5_binding_summary_matrix
#X-squared = 37.084, df = 2, p-value = 8.856e-09




##################################################
#### vista enhancers
###############################################


str(merged_peaks.df5)

merged_peaks.df5.vista <-  merged_peaks.df5 %>% 
  dplyr::filter(meis2_peak == "y" | dlx5_peak == "y") %>%
  mutate(vista_activity = ifelse(is.na(vista_activity), "no_vista_overlap", vista_activity)) %>%
  mutate(vista_activity = ifelse(grepl("forebrain_positive", vista_activity), "forebrain_positive", 
                                 ifelse(! grepl("forebrain_positive", vista_activity) &  grepl("other_positive", vista_activity), "other_positive",
                                        ifelse(! grepl("forebrain_positive", vista_activity) &  !grepl("other_positive", vista_activity) & grepl("negative", vista_activity), "negative", "no_vista_overlap")))) %>%
  mutate(vista_activity = factor(vista_activity, levels = c("forebrain_positive", "other_positive", "negative", "no_vista_overlap"))) %>%
  mutate(meis2_peak = ifelse(meis2_peak == "y", "TRUE", "FALSE"), dlx5_peak = ifelse(dlx5_peak == "y", "TRUE", "FALSE"), lhx6_peak = ifelse(lhx6_peak == "y", "TRUE", "FALSE")) %>%
  rename("MEIS2" = "meis2_peak", "DLX5" = "dlx5_peak") %>%
  mutate(MEIS2 = as.logical(MEIS2), DLX5 = as.logical(DLX5))

# top part of the plot has to be created separate from the upset plot to only display different Vista categories, else most peaks overlap no Vista enhancer

merged_peaks.df5.vista.top <- merged_peaks.df5.vista %>% filter(vista_activity != "no_vista_overlap")

enhancer_cols2 <- c("MEIS2", "DLX5")  

ComplexUpset_vista_count_top <- ComplexUpset::upset(
  merged_peaks.df5.vista.top,
  enhancer_cols2,
  name="binding site",
  base_annotations = list("Intersection size" = intersection_size(text=list(size=8*0.36))),
  themes=upset_modify_themes(list("overall_sizes"=theme(axis.text=element_text(size=8, angle = 90, family = "Helvetica"), axis.title=element_text(size=9, family = "Helvetica"), panel.grid.minor = element_blank()), # set sizes, left panel
                                  "Intersection size"=theme(axis.text=element_text(size=8, family = "Helvetica"),axis.title=element_text(size=9, family = "Helvetica")),
                                  "intersections_matrix"=theme(axis.text=element_text(size=8, family = "Helvetica"),axis.title=element_text(size=9, family = "Helvetica"),panel.grid = element_blank()))),
  annotations = list("vista"=(
    ggplot(mapping=aes(x=intersection,fill=vista_activity))+
      geom_bar()+
      scale_fill_manual(values = c("#9CB3D3", "#E9A860", "black"), labels = c( "forebrain positive","other positive","negative"))+
      theme(axis.title.y = element_text(size=9), axis.text = element_text(size=8),legend.position = "right",legend.text = element_text(size=8), legend.title=element_text(size=9), )+
      labs(fill="activity", y="Vista enhancers [count]")+
      theme(axis.title.y = element_text(size=9, family = "Helvetica"), axis.text = element_text(size=8, family = "Helvetica"),legend.text = element_text(size=8, family = "Helvetica"), legend.title=element_text(size=9, family = "Helvetica"), legend.margin = margin(l=1) ))),
  set_sizes=upset_set_size(position="right"), guides="over", width_ratio = 0.5)

ComplexUpset_vista_count_top

#### lower part of the plot to contain all peaks 


ComplexUpset_vista_count_bottom <- ComplexUpset::upset(
  merged_peaks.df5.vista,
  enhancer_cols2,
  name="binding site",
  base_annotations = list("Intersection size" = intersection_size(text=list(size=8*0.36))),
  themes=upset_modify_themes(list("overall_sizes"=theme(axis.text=element_text(size=8, angle = 90, family = "Helvetica"), axis.title=element_text(size=9, family = "Helvetica"), panel.grid.minor = element_blank()), # set sizes, left panel
                                  "Intersection size"=theme(axis.text=element_text(size=8, family = "Helvetica"),axis.title=element_text(size=9, family = "Helvetica")),
                                  "intersections_matrix"=theme(axis.text=element_text(size=8, family = "Helvetica"),axis.title=element_text(size=9, family = "Helvetica"),panel.grid = element_blank()))),
  annotations = list("vista"=(
    ggplot(mapping=aes(x=intersection,fill=vista_activity))+
      geom_bar()+
      scale_fill_manual(values = c("#9CB3D3", "#E9A860", "black", "grey"), labels = c( "forebrain positive","other positive","negative", "no Vista overlap"))+
      theme(axis.title.y = element_text(size=9), axis.text = element_text(size=8),legend.position = "right",legend.text = element_text(size=8), legend.title=element_text(size=9), )+
      labs(fill="activity", y="Vista enhancers [count]")+
      theme(axis.title.y = element_text(size=9, family = "Helvetica"), axis.text = element_text(size=8, family = "Helvetica"),legend.text = element_text(size=8, family = "Helvetica"), legend.title=element_text(size=9, family = "Helvetica"), legend.margin = margin(l=1) ))),
  set_sizes=upset_set_size(position="right"), guides="over", width_ratio = 0.5)

ComplexUpset_vista_count_bottom



