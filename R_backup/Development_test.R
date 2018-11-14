#CORE clustering 
day5 <- sample2
cellnames<-colnames(day5$dat5_counts)
cluster <-day5$dat5_clusters
cellnames <- data.frame("cluster" = cluster, "cellBarcodes" = cellnames)
#day5$dat5_counts needs to be in a matrix format
mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
                       CellMetadata = day5$dat5_clusters)


day2 <- sample1
cellnames<-colnames(day2$dat2_counts)
cluster <-day2$dat2_clusters
cellnames <- data.frame("cluster" = cluster, "cellBarcodes" = cellnames)
#day5$dat5_counts needs to be in a matrix format
mixedpop1 <-NewscGPS_SME(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
                       CellMetadata = day2$dat2_clusters)

#For sample 2 bagging
test2 <- CORE_scGPS_bagging(mixedpop2, remove_outlier = c(0), PCA=FALSE, nPCs=50, ngenes=1500, 
	bagging_run = 20, subsample_proportion = 0.6)

#For sample 1 bagging
test <- CORE_scGPS_bagging(mixedpop1, remove_outlier = c(1), PCA=FALSE, bagging_run = 20, subsample_proportion = .8)

plot_CORE(test$tree, list_clusters = test$Cluster)
#see all counts clusters
table(test$optimalClust)
#see the optimal clusters
test$optimalMax

test2 <- CORE_scGPS_bagging(mixedpop2, remove_outlier = c(0), PCA=FALSE, nPCs=50, ngenes=1500,
bagging_run = 10, subsample_proportion = 0.8)
table(test2$optimalClust)
#see the optimal clusters
test2$optimalMax
plot_CORE(test2$tree, list_clusters = test2$Cluster)
#see all counts clusters
table(test2$optimalClust)
#see the optimal clusters
test2$optimalMax

#For without bagging
test2 <- CORE_scGPS(mixedpop2, remove_outlier = c(0), PCA=FALSE, nPCs=50, ngenes=1500)
plot_CORE(test2$tree, list_clusters = test2$Cluster)


test2$optimalClust
test2$optimalMax

mean(test1$optimalClust)
median(test1$optimalClust)

##################################scGPS prediction model development############################################
library(SingleCellExperiment)

c_selectID<-1
out_idx<-1
day2 <- sample1
mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
                      CellMetadata = day2$dat2_clusters)
day5 <- sample2
mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
                      CellMetadata = day5$dat5_clusters)
genes <-GeneList
genes <-genes$Merged_unique
listData  <- training_scGPS(genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID, listData =list(), out_idx=out_idx)
names(listData)
listData$Accuracy

#CORE clustering
day5 <- dat_exprs
cellnames<-colnames(dat_exprs)
cluster <-day5$dat5_clusters
cellnames <- data.frame("cluster" = cluster, "cellBarcodes" = cellnames)
#day5$dat5_counts needs to be in a matrix format
mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
CellMetadata = day5$dat5_clusters)

#New dataset
path_dat="/Users/quan.nguyen/Documents/Powell_group_MacQuan/scGPSmanuscript/BenchMarking/Datasets/scMapDataset/"
dat <-readRDS(paste0(path_dat,"yan.rds"))
dat_exprs <- assay(dat)



##################################sankey#########################################

c_selectID <- 1
genes = DEgenes$DE_Subpop1vsRemaining$id[1:200] #top 200 gene markers distinguishing cluster 1
genes <- gsub("_.*", "", genes)

LSOLDA_dat1 <- bootstrap_scGPS(nboots = 1,mixedpop1 = mixedpop2, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())

#>

c_selectID <- 2
genes = DEgenes$DE_Subpop2vsRemaining$id[1:200]
genes <- gsub("_.*", "", genes)
LSOLDA_dat2 <- bootstrap_scGPS(nboots = 1,mixedpop1 = mixedpop2, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#>


c_selectID <- 3
genes = DEgenes$DE_Subpop3vsRemaining$id[1:200]
genes <- gsub("_.*", "", genes)
LSOLDA_dat3 <- bootstrap_scGPS(nboots = 1,mixedpop1 = mixedpop2, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#>

row_cluster <-length(unique(colData(mixedpop2)[,1]))
#>

reformat_LASSO <-function(c_selectID = NULL, s_selectID = NULL, LSOLDA_dat = NULL,
nPredSubpop = row_cluster, Nodes_group = "#7570b3"){
    LASSO_out <- summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop = nPredSubpop)
    LASSO_out <-as.data.frame(LASSO_out)
    temp_name <- gsub("LASSO for subpop", "C", LASSO_out$names)
    temp_name <- gsub(" in target mixedpop", "S", temp_name)
    LASSO_out$names <-temp_name
    source <-rep(paste0("C",c_selectID,"S",s_selectID), length(temp_name))
    LASSO_out$Source <- source
    LASSO_out$Node <- source
    LASSO_out$Nodes_group <- rep(Nodes_group, length(temp_name))
    colnames(LASSO_out) <-c("Value", "Target", "Source", "Node", "NodeGroup")
    LASSO_out$Value <- as.numeric(as.vector(LASSO_out$Value))
    return(LASSO_out)
}

LASSO_C1S2  <- reformat_LASSO(c_selectID=1, s_selectID =2, LSOLDA_dat=LSOLDA_dat1,
nPredSubpop = row_cluster, Nodes_group = "#7570b3")

LASSO_C2S2  <- reformat_LASSO(c_selectID=2, s_selectID =2, LSOLDA_dat=LSOLDA_dat2,
nPredSubpop = row_cluster, Nodes_group = "#1b9e77")

LASSO_C3S2  <- reformat_LASSO(c_selectID=3, s_selectID =2, LSOLDA_dat=LSOLDA_dat3,
nPredSubpop = row_cluster, Nodes_group = "#e7298a")


combined <- rbind(LASSO_C1S2,LASSO_C2S2,LASSO_C3S2 )
combined <- combined[is.na(combined$Value) != TRUE,]
combined_D3obj <-list(Nodes=combined[,4:5], Links=combined[,c(3,2,1)])

library(networkD3)

Node_source <- as.vector(sort(unique(combined_D3obj$Links$Source)))
Node_target <- as.vector(sort(unique(combined_D3obj$Links$Target)))
Node_all <-unique(c(Node_source, Node_target))

#assign IDs for Source (start from 0)
Source <-combined_D3obj$Links$Source
Target <- combined_D3obj$Links$Target

for(i in 1:length(Node_all)){
    Source[Source==Node_all[i]] <-i-1
    Target[Target==Node_all[i]] <-i-1
}

combined_D3obj$Links$Source <- as.numeric(Source)
combined_D3obj$Links$Target <- as.numeric(Target)
combined_D3obj$Links$LinkColor <- combined$NodeGroup

#prepare node info
node_df <-data.frame(Node=Node_all)
node_df$id <-as.numeric(c(0, 1:(length(Node_all)-1)))

suppressMessages(library(dplyr))
Color <- combined %>% count(Node, color=NodeGroup) %>% select(2)
node_df$color <- Color$color

suppressMessages(library(networkD3))
p1<-sankeyNetwork(Links =combined_D3obj$Links, Nodes = node_df,  Value = "Value", NodeGroup ="color", LinkGroup = "LinkColor", NodeID="Node", Source="Source", Target="Target",
fontSize = 22 )
p1

LASSO_out <- summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat1, nPredSubpop = 4)












