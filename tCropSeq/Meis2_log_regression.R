## logistic regression: interneuron vs PN ##

library(Seurat)
library(dplyr)
library(glmnet)
library(ggrepel)
library(MASS)
library(bmrm)
library(ggsignif)

set.seed(4)

## helper functions ##
train_model <- function(expr_mtx, target_vec, regr_model = "binomial", type.multinomial = "no") {
  ## account for balanced design:
  min_num <- min(table(target_vec))
  ids_to_keep <- sapply(unique(target_vec), function(el) {
    sample(names(target_vec)[target_vec == el], min_num)
  })
  ids_to_keep <- as.vector(ids_to_keep)
  
  ## subset for sampled cellIDs:
  target_vec <- target_vec[ids_to_keep]
  expr_mtx <- expr_mtx[, ids_to_keep]
  
  ## randomly sample 2/3 of the data to be training data:
  train_ids <- sample(colnames(expr_mtx), size = round(2/3*ncol(expr_mtx)))
  test_ids <- colnames(expr_mtx)[!colnames(expr_mtx) %in% train_ids]
  
  train_mtx <- expr_mtx[, train_ids]
  train_vec <- target_vec[train_ids]
  test_mtx <- expr_mtx[, test_ids]
  test_vec <- target_vec[test_ids]
  
  ## run regression:
  if(regr_model == "binomial") {
    cv_fit <- cv.glmnet(x = t(train_mtx), y = train_vec, family = regr_model, type.measure = "class", alpha = 1)
  } else if(regr_model == "multinomial") {
    if(type.multinomial == "grouped") {
      cv_fit <- cv.glmnet(x = t(train_mtx), y = train_vec, family = regr_model, type.measure = "class", type.multinomial = "grouped", alpha = 1)
    } else {
      cv_fit <- cv.glmnet(x = t(train_mtx), y = train_vec, family = regr_model, type.measure = "class", alpha = 1)
    }
  } 
  else {
    print("Please input proper model family, see  argument in cv.glmnet")
  }
  ## predict test sample:
  test_prediction <- predict(cv_fit, newx = t(test_mtx), s = "lambda.min", type = "class")
  correct_test <- sum(test_prediction == test_vec) / length(test_prediction)
  print(paste0("Percentage of correct predictions in test set: ", correct_test))
  
  ## return trained model:
  return(cv_fit)
}


get_gene_coefs <- function(cv_fit, seurat_obj, filter_genes = TRUE, variable_features, num_active_th = 1000) {
  coef_mtx <- as.matrix(coef(cv_fit))[2:nrow(coef(cv_fit)),]
  coef_mtx <- coef_mtx[coef_mtx != 0]
  coef_mtx_prob <- exp(coef_mtx)/(1 + exp(coef_mtx))
  coef_mtx_df <- data.frame("coef" = coef_mtx_prob, "gene" = names(coef_mtx), "idx" = seq(length(coef_mtx_prob)))
  coef_mtx_tb <- as_tibble(coef_mtx_df)
  coef_mtx_tb %>% 
    arrange(coef) -> coef_mtx_tb
  coef_mtx_tb$coef_order <- seq(1:nrow(coef_mtx_tb))
  ## filter genes for abundancy:
  if(filter_genes == TRUE) {
    num_active_cells <- apply(seurat_obj@assays$RNA@counts[variable_features, ], 1, function(vec) {
      sum(vec > 0)
    })
    filtered_genes <- coef_mtx_tb$gene[sapply(coef_mtx_tb$gene, function(gene){
      num_active_cells[gene] > num_active_th
    })]
    coef_mtx_tb <- coef_mtx_tb[coef_mtx_tb$gene %in% filtered_genes, ]
    coef_mtx_tb$coef_order <- seq(1,nrow(coef_mtx_tb))
  }
  return(coef_mtx_tb)
}

## load data from bandler et al
## you can download the pre-processed seuratobject (seurobj_Bandler_et_al.RDS) from here: https://keeper.mpdl.mpg.de/f/6eb4e412d05f4a18b80d/?dl=1

bandler_seurat <- readRDS("seurobj_Bandler_et_al.RDS")

## visualize:
DimPlot(bandler_seurat, group.by = "LR")
table(bandler_seurat$LR)

## train model on post-mitotic cells:
bandler_pm_sub <- subset(bandler_seurat, subset = LR %in% c("IN","PN"))

GEX_PM_mtx <- bandler_pm_sub@assays$integrated@scale.data
dim(GEX_PM_mtx)

fate_vec <- as.character(bandler_pm_sub$LR)
names(fate_vec) <- colnames(bandler_pm_sub)

log_fate_fit <- train_model(expr_mtx = GEX_PM_mtx, target_vec = fate_vec, regr_model = "binomial")
# [1] "Percentage of correct predictions in test set: 0.991500934897161"

plot(log_fate_fit)
log_fate_fit$lambda.min
#log_fate_fit$nzero

coef_log_tb  <- get_gene_coefs(cv_fit = log_fate_fit, seurat_obj = bandler_pm_sub, filter_genes = TRUE,
                               variable_features = VariableFeatures(bandler_seurat), num_active_th = 1000)

ggplot(coef_log_tb, aes(x = coef_order,y = coef)) + 
  geom_point(size = 0.8) +
  geom_text_repel(mapping = aes(x = coef_order, y = coef, label = gene), max.overlaps = 200)  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Index") +
  geom_hline(yintercept = 0.50, color = "red")
ggsave("/datastore_share/Users/neuhaus/Meis2_z_score/results/regression/regression_binom_fate_gene_coefs.pdf", width = 7, height = 4)

## interneuron features:
FeaturePlot(bandler_seurat, features = head(coef_log_tb$gene, n = 10))
ggsave("/datastore_share/Users/neuhaus/Meis2_z_score/results/regression/interneuron_top10_features.pdf", width = 15, height = 15)

## PN features:
FeaturePlot(bandler_seurat, features = tail(coef_log_tb$gene, n = 10))
ggsave("/datastore_share/Users/neuhaus/Meis2_z_score/results/regression/PN_top10_features.pdf", width = 15, height = 15)



## predict for mitotic cells:
bandler_m_sub <- subset(bandler_seurat, subset = LR == "mitotic")

GEX_M_mtx <- bandler_m_sub@assays$integrated@scale.data
mitotic_prediction <- predict(log_fate_fit, newx = t(GEX_M_mtx), s = "lambda.min")
head(mitotic_prediction)

mitotic_prediction_prob <- sapply(as.numeric(mitotic_prediction), function(el) {exp(el)/(1+exp(el))})
names(mitotic_prediction_prob) <- rownames(mitotic_prediction)

pred_vec <- rep(NA, ncol(bandler_seurat))
names(pred_vec) <- colnames(bandler_seurat)
pred_vec[names(mitotic_prediction_prob)] <- mitotic_prediction_prob

hist(mitotic_prediction_prob, breaks = 20)
## most cells predicted to become interneurons..

bandler_seurat$predicted_fate <- pred_vec
FeaturePlot(bandler_seurat, features = "predicted_fate", cols = c("chartreuse4", "green", "yellow", "orange","red")) +
  ggtitle("Predicted Fate of Mitotic Cells")
ggsave("/datastore_share/Users/neuhaus/Meis2_z_score/results/regression/mitotic_fate_prediction_umap.pdf", width = 9, height = 6)

save.image("/datastore_share/Users/neuhaus/Meis2_z_score/regression_workspace_05_08.RData")
