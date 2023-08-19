---
title: "Proportion"
---

## Proportion Change 

# Load necessary libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

# Set working directory
setwd("/data/mayerlab/cmayer/meis2_paper/results")

# Read Seurat objects
E16 <- readRDS("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/E16.rds")
P7GE_sub <- readRDS("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/P7GE_sub.rds")
E16_CCA_inhib <- readRDS("/data/mayerlab/Github/Dvoretskova_et_al_2023-DATASETS/E16_CCA_inhib.rds")

seur_obj <- E16

# Extract unique gRNAs and inhibitory cell types
gRNAs <- unique(seur_obj@meta.data$gRNA)
inhibitory <- unique(seur_obj@meta.data$CCA_Assignment)

# Calculate proportion change
prop <- data.frame(class=character(), dataset=character(), value=numeric(), gRNA=character(), stringsAsFactors=FALSE)

for (i in gRNAs) {
  guide_object <- subset(seur_obj, subset = gRNA == i)
  data <- guide_object@meta.data[c("CCA_Assignment", "dataset", "gRNA")]
  prop_table <- prop.table(table(data$CCA_Assignment, data$dataset), margin = 2) %>% reshape2::melt()
  colnames(prop_table) <- c("class", "dataset", "value")
  prop_table$gRNA <- rep(i, length(prop_table$dataset))
  prop <- rbind(prop, prop_table)
}

lacZ_prop <- prop[prop$gRNA == "glacZ", ]
lacZ_prop <- data.table(lacZ_prop)
lacZ_prop[, Mean := mean(value), by = class]

lacZ_AvgProp <- data.frame(
  clusters = unique(lacZ_prop$class),
  mean_prop = unique(lacZ_prop$Mean),
  stringsAsFactors = FALSE
)

clusters <- unique(prop$class)
prop_change_total <- data.frame(
  class = character(),
  dataset = character(),
  value = numeric(),
  gRNA = character(),
  prop_change = numeric(),
  stringsAsFactors = FALSE
)

for (i in clusters) {
  single_change <- prop[prop$class == i, ]
  single_change <- single_change %>% mutate(prop_change = single_change$value / lacZ_AvgProp[lacZ_AvgProp$clusters == i, ]$mean_prop)
  prop_change_total <- rbind(prop_change_total, single_change)
}
prop_change_total <- filter(prop_change_total, value > 0)

# Categorize proportion change
prop_change_total$category <- ifelse(prop_change_total$prop_change > 1, "high", "low")
prop_change_total <- prop_change_total %>%
  dplyr::group_by(class, gRNA) %>%
  dplyr::mutate(mean_prop_change = mean(prop_change, na.rm = TRUE))
prop_change_total$category <- ifelse(prop_change_total$mean_prop_change < 1, "decrease", "increase")

# Filter unwanted classes
unwanted_class <- c("Ependymal", "Unknown")
prop_filter <- prop_change_total[!prop_change_total$class %in% unwanted_class, ]
prop_filter <- prop_filter[!prop_filter$gRNA %in% "glacZ", ] # remove lacZ

# Plot the proportion change
prop_filter <- prop_filter[, c("class", "dataset", "prop_change")]
prop_filter$class <- factor(prop_filter$class, levels = c("IN:Cck/Reln", "IN:Nfib/Tcf4", "IN:Tiam2/Zfp704", "IN:Tshz1/Pbx1", "IN:Calb2/Nxph1", "Mitotic", "IN:Lhx6/Npy", "PN:Foxp1/Six3", "PN:Foxp1/Isl1", "PN:Isl1/Bcl11b", "PN:Meis2/Bcl11b", "PN:Ebf1/Zfp503", "PN:Isl1/Meis2", "IN:Nr2f2/Nnat"))

# Set the desired dimensions (2:1 aspect ratio) for the plot
plot_width <- 6  # Inches
plot_height <- 4  # Inches

# Create the plot
plot <- ggplot(prop_filter, aes(x = class, y = prop_change)) +
  geom_boxplot(color = "black", fill = "lightgray", outlier.shape = NA) +
  labs(x = "Cell cluster", y = "log10(Proportion change)") +
  ggtitle("Proportion change of perturbed clusters compared to lacZ control") +
  scale_y_log10() +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")
plot




## Regression Analysis

# Load necessary libraries for regression analysis (already loaded earlier)
library(datasets)
library(foreign)
library(MASS)
library(broom)
library(lmtest)

# Function to calculate gRNA effect on cell composition. Taken from Di https://github.com/klarman-cell-observatory/ivPerturbSeq

CellComp_Poisson<-function(seur, seur2, celltype="CellType",perturbations="perturbation",batch="batch",cutoff=10)
{
print("Clean Data")
meta=seur@meta.data[,c(celltype,perturbations,batch)]
meta2 =seur2@meta.data[c(perturbations, batch)]
colnames(meta)=c("CellType","Pert","Batch")
colnames(meta2)= c("Pert", "Batch")
meta <- meta %>% group_by(CellType,Pert,Batch) %>% summarise(Num=length(CellType)) %>% as.data.frame()
meta2<-meta %>% group_by(Pert,Batch) %>% summarise(Tot=sum(Num)) %>% as.data.frame()
meta=left_join(meta,meta2)
##so meta is dataframe of 5 columns: celltype, perturbation, batch, Number of total cells of that celltype/pert/batch, and number of total cells of that pert/batch
meta=meta[meta[,"Tot"]>cutoff,]
# meta = meta[, "Num" > 10]
meta["Pert"]=relevel(factor(meta[,"Pert"]),ref="glacZ")
lst=list()
for(i in unique(meta[,"CellType"])){lst[[i]]=meta[meta[,"CellType"]==i,]}
print("Fit model!")
out<-lapply(lst,function(cur){
celltype=cur[1,"CellType"]
print(celltype)
cur["logTot"]=log(cur[,"Tot"])
fit <-glm(Num~offset(logTot)+Batch+Pert,data=cur,family=poisson())
tab=summary(fit)
tab=tab$coefficients
tab=data.frame(tab)
tab=tab[grep("Pert",rownames(tab)),]
tab["Gene"]=sub("Pert","",rownames(tab))
tab["CellType"]=celltype
tab=tab[,c(5,6,4,1,2,3)]
colnames(tab)[3]="pval"
return(tab)
})
tab=do.call(rbind,out)
rownames(tab)=NULL
tab=tab[order(tab[,"pval"]),]
tab["padj"]=p.adjust(tab[,"pval"],"fdr")
print("Done!")
return(tab)
}

#Calculate P7 proportion composition effect.
P7_Poi <- CellComp_Poisson(seur =  P7GE_sub, seur2 = P7GE_sub, celltype = "level4", perturbations = "gRNA", batch = "dataset", cutoff = 10)

#Calculate E16 proportion composition effect. 
E16_Poi <- CellComp_Poisson(seur =  E16_CCA_inhib, seur2 = E16, celltype = "CCA_Assignment", perturbations = "gRNA", batch = "dataset", cutoff = 0)
E16_Poi <- CellComp_Poisson(seur =  E16_CCA_inhib, seur2 = E16_CCA_inhib, celltype = "CCA_Assignment", perturbations = "gRNA", batch = "dataset", cutoff = 0)


E16_Poi$CellType <- factor(E16_Poi$CellType, levels = c("Mitotic", "IN:Calb2/Nxph1", "IN:Nr2f2/Nnat", "IN:Tiam2/Zfp704", "IN:Nfib/Tcf4", "IN:Lhx6/Npy", "IN:Cck/Reln", "IN:Tshz1/Pbx1", "PN:Isl1/Meis2", "PN:Isl1/Bcl11b", "PN:Ebf1/Zfp503", "PN:Meis2/Bcl11b", "PN:Foxp1/Six3", "PN:Foxp1/Isl1"))



poisson_df <- E16_Poi # or P7_Poit
poisson_df <- P7_Poi

mod_plot <- ggplot(poisson_df) +
  geom_point(mapping = aes(x = CellType, y = Gene, colour = Estimate, size = -log(padj))) +
 # scale_colour_viridis_b() +  
  scale_colour_gradient2(low = muted("blue"), high = muted("red")) +
  scale_size_continuous("-log(Pvalue)", range = c(2,9), breaks = c(1,2,3,4), labels = c("1", "2", "3", ">4")) +
  geom_point(mapping = aes(x = CellType, y = Gene, size = -log(padj)), shape = 21, colour = "black", stroke =  1, data = poisson_df[which(-log(poisson_df$padj)> 3),]) + 
  labs(size = "-log(Pvalue)", colour = "Effect size") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(panel.grid.major = element_blank()) +
  ggtitle("Perturbation effect on proportion compared to lacZ control") 
mod_plot

seurat_obj <- SetIdent(P7GE_sub, value = "level4")

