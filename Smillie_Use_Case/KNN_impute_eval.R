#!/usr/bin/env Rscript 

library('Seurat')
library('dplyr')
library(tidyverse)
library(gridExtra)
library('ggplot2')
library(neighbr)
library(caret)

smillie_qc_file <-readRDS(file='/mnt/rstudio/smohammed/CNA/smillie/Smillie_CNA_fibroblasts_cells.RDS')
temp_qc_file <- readRDS(file='/mnt/rstudio/smohammed/Impute/Smillie_IBD/smillie_temp_qc_file_imputed.rds')


seurat_knn_ranks<-function(seurat_obj, k, g, c, indexes) {
  
   prc_n <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@data[g, c])))

  prc_n$id <- c
  prc_n$celltype <-factor(seurat_obj@meta.data[c,'Cluster'])

  ## Split data for train-test 
  prc_train = prc_n[indexes, ]
  prc_test = prc_n[-indexes, -c(dim(prc_n)[2],dim(prc_n)[2]-1)]

  fit <- neighbr::knn(train_set=prc_train, test_set=prc_test,
           k=k,
           categorical_target="celltype",
           continuous_target= NULL,
           comparison_measure="euclidean",
           return_ranked_neighbors=k,
           id="id")

  return(list(fit$test_set_scores,prc_n))
}

set.seed(11)
n=100
g=10
# create partition of randomly sampled cells for training and testing
cell_of_interest <- sample(ncol(smillie_qc_file), n, replace = FALSE)
cell.indexes = createDataPartition(cell_of_interest, p = .70, list = FALSE)

# sampling of gene indexes without replacement
gene.indexes <- sample(nrow(smillie_qc_file), g, replace = FALSE)

# Inflammatory fibroblasts cells 
#cell_of_interest<-c(intersect(smillie_qc_file@meta.data %>% subset(seurat_clusters==2) %>% rownames(), temp_qc_file@meta.data %>% subset(seurat_clusters==5) %>% rownames()))
#cell.indexes = createDataPartition(1:length(cell_of_interest), p = .70, list = FALSE)

                    
acc<-c()
raw_knn <- list()
impute_knn <- list()
k<-c(10, 20, 40, 50)
for (i in 1:length(k)) {
  cat("Running k -", k[i], "(raw data)","\n") 
  raw_knn[[i]]<-seurat_knn_ranks(seurat_obj=smillie_qc_file,k=k[i], g=gene.indexes, c=cell_of_interest, indexes=cell.indexes)
  
    
  cat("Running k -", k[i], "(imputed data)","\n") 
  impute_knn[[i]]<-seurat_knn_ranks(seurat_obj=temp_qc_file,k=k[i], g=gene.indexes, c=cell_of_interest, indexes=cell.indexes)
  
  acc[i] <- mean(raw_knn[[i]][[1]]==impute_knn[[i]][[1]])
}

png(filename="~/Impute/Smillie_IBD/knn_acc.png",width=1000,height=1000,res=150)
plot(k,acc,type='b')
dev.off()

# https://towardsdatascience.com/calculate-similarity-the-most-relevant-metrics-in-a-nutshell-9a43564f533e
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jacc_score_list<-list()
jacc_score_vector<-c()

for (i in 1:length(k)) {
  cat("Running k -", k[i],"\n")
  
  for (row in 1:nrow(raw_knn[[i]][[1]])) {
   #cat("Running row -", row,"\n") 
    jacc_score_vector[row] <- jaccard( (raw_knn[[i]][[1]][row,-1] %>% remove_rownames(.) %>% unlist()), 
                                       (impute_knn[[i]][[1]][row,-1] %>% remove_rownames(.)%>% unlist()) )
  }
 jacc_score_list[[i]] <- jacc_score_vector
}


library(gridExtra)
library(ggplot2)

stat_box_data <- function(y, upper_limit = max(df$Jaccard_index) *1.05) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n'))
  )
}

p <- list()
for (i in 1:length(k)){
  
  ## Test data indexes are used to extract cell types and associated neighbors jaccard index. 
  df=data.frame(Celltype=impute_knn[[i]][[2]][-cell.indexes,'celltype'], Jaccard_index=jacc_score_list[[i]])
  
  p[[i]] <- df %>% group_by(Celltype) %>%
          ggplot(aes(x=Celltype,y=Jaccard_index)) + 
          geom_boxplot() +
          stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9) +
          scale_x_discrete(guide = guide_axis(angle = 90)) +
          ylim(0,1) + labs(title = paste0("knn=",k[i])) +
         geom_hline(yintercept=0.5, linetype='dashed', col = 'red')

}

png(filename="~/Impute/Smillie_IBD/knn_jaccard_index.png",width=1000,height=1000,res=150)
do.call(grid.arrange,p)
dev.off()
