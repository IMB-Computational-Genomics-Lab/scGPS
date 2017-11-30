#' Main training function for a subpopulation
#'
#' @description  Training a haft of all cells to find optimal LASSO and LDA
#' models to predict a subpopulation
#' @param mixedpop1 is a \linkS4class{SingleCellExperiment} object from the train mixed population
#' @param mixedpop2 is a \linkS4class{SingleCellExperiment} object from the target mixed population
#' @param genes a vector of gene names (for LASSO shrinkage); gene symbols must be
#' in the same format with gene names in subpop2. Note that genes are listed by the order
#' of importance, e.g. differentially expressed genes that are most significan, so that
#' if the gene list contains too many genes, only the top 500 genes are used.
#' @param c_selectID a selected number to specify which subpopulation to be used for training
#' @param out_idx a number to specify index to write results into the list output.
#' This is needed for running bootstrap.
#' @return a \code{list} with prediction results written in to the index \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#'
#' c_selectID<-1
#' out_idx<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                      CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                      CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' listData  <- training_scGPS(genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID, listData =list(), out_idx=out_idx)
#'

training_scGPS <-function(genes, mixedpop1 = NULL,
                         mixedpop2 = NULL, c_selectID, listData =list(), out_idx) {
  # subsammpling-----------------------------------------------------------------
  # taking a subsampling size of 50% of the cluster_select out for training
  subpop1cluster_indx <- which(colData(mixedpop1)[,1] == c_selectID)
  subremaining1_indx <- which(colData(mixedpop1)[,1] != c_selectID)

  subsampling <- round(length(subpop1cluster_indx)/2)
  subpop1_train_indx <-   sample(subpop1cluster_indx, subsampling , replace = F)
  # check if there is a very big cluster present in the dataset
  # C/2 = SubSampling >length(total - C, which is cluster_compare)
  if (ncol(mixedpop1) > 2*subsampling) {
    subremaining1_train_indx <- sample(subremaining1_indx, subsampling, replace = F)
  } else {
    subremaining1_train_indx <- subremaining1_indx
  }
  # done subsampling------------------------------------------------------------

  # making sure the gene list contains fewer than 500 genes---------------------
  # select top 500, the DE_result is an already sorted list
  if(length(genes) >500){genes <- genes[1:500]}
  #select genes
  names1 <- elementMetadata(mixedpop1)[,1]
  subpop1_selected_genes <- names1[which(names1 %in% genes)]
  names2 <-elementMetadata(mixedpop2)[,1]
  genes_in_both <-names2[which(names2 %in% subpop1_selected_genes)]
  genes_in_both_idx1 <-which(names1 %in%  genes_in_both)
  # done selecting genes--------------------------------------------------------

  #prepare LASSO training matrices----------------------------------------------

  #prepare predictor matrix containing both clutering classes
  predictor_S1 <- mixedpop1[genes_in_both_idx1,c(subpop1_train_indx,subremaining1_train_indx)]
  predictor_S1 <-assay(predictor_S1)

  #generate categorical response

  #assign class labels############
  c_compareID <- 1:length(unique(colData(mixedpop1)[,1]))
  c_compareID <- paste0(c_compareID[-which(c_compareID==c_selectID)], collapse = "")
  # change cluster number to character
  c_selectID <- as.character(c_selectID)
  #done assigning class labels####

  #set all y values to cluster select (character type) length = #cells
  y_cat = rep(c_selectID,ncol(predictor_S1))

  #replace values for cluster compare
  #get indexes for cells in predictor_S1 belong to "remaining class"
  remainingClass_Indx_in_y <-which(colnames(predictor_S1) %in%
    colnames(mixedpop1[,subremaining1_train_indx]))

  #change value of the cluster id
  y_cat[remainingClass_Indx_in_y] <-rep(c_compareID, length(remainingClass_Indx_in_y))
  #Done prepare LASSO training matrices-----------------------------------------

  #prepare LDA training matrices------------------------------------------------
  #fitting with lda, also with cross validation
  dataset <- t(predictor_S1) #(note predictor_S1 =t(gene_S1))
  dataset <- as.data.frame(dataset)

  #remove genes with no variation across cells
  Zero_col<-which(colSums(dataset) == 0)
  duplicated_col <-which(duplicated(colnames(dataset)) == TRUE)
  if(length(c(Zero_col, duplicated_col)) !=0){dataset <-dataset[,-c(Zero_col, duplicated_col)]}
  dataset$Cluster_class <- as.character(y_cat)
  #Done LDA training matrices---------------------------------------------------

  #Fit the LASSO and LDA models-------------------------------------------------

  #fitting with cross validation to find the best LASSO model
  cvfit = cv.glmnet(t(predictor_S1), y_cat, family = "binomial", type.measure = "class")

  #fit LDA
  trainControl <- trainControl(method="repeatedcv", number=10, repeats=3)
  metric <- "Accuracy"
  fit.lda <- train(Cluster_class ~., data=dataset, method="lda", metric=metric,
    trControl=trainControl, na.action=na.omit)

  #Done fitting the LASSO and LDA models----------------------------------------

  #Extract coefficient Beta for a gene for an optimized lambda value------------

  cvfit_out <-as.matrix(coef(cvfit, s = cvfit$lambda.min))
  cvfit_out <-as.data.frame(cvfit_out)
  #find number of genes with coefficient different to 0
  cvfit_out$name <-row.names(cvfit_out)
  sub_cvfit_out <-cvfit_out[cvfit_out$`1` != 0,]
  #extract deviance explained
  t_DE <- as.matrix(print(cvfit$glmnet.fit))
  dat_DE <-as.data.frame(t_DE)
  colnames(dat_DE) <-c('Dfd', 'Deviance', 'lambda')
  #to get the coordinate for lambda that produces minimum error
  dat_DE_Lambda_idx <- which(round(dat_DE$lambda,digit=3) == round(cvfit$lambda.min,digits=3))
  dat_DE <-dat_DE[1:dat_DE_Lambda_idx[1],]

  dat_DE %>% group_by(Dfd) %>% summarise(Deviance = max(Deviance)) -> dat_DE_fm_DE
  dat_DE_fm_DE <-as.data.frame(dat_DE_fm_DE)
  dat_DE_fm_DE$DEgenes <-paste0('genes_cluster',c_selectID)
  remaining <-c('remaining', 1, 'DEgenes')
  dat_DE_fm_DE <-rbind(dat_DE_fm_DE, remaining)
  #Done extracting the coefficients---------------------------------------------

  #Cross validation test to estimate accuracy-----------------------------------

  #keep all cells except for those used in the training set
  #recall that: subpop1_train_indx = subpop1cluster_indx[sample(1:ncol(mixedpop1),subsampling , replace = F)]

  cluster_select_indx_Round2 <- subpop1cluster_indx[- which(subpop1cluster_indx %in% subpop1_train_indx)]

  subremaining1_train_order <-  which(subremaining1_indx %in% subremaining1_train_indx)

  #check the SubSampling for the target population
  if (length(subremaining1_indx) - length(subremaining1_train_indx) > subsampling) {
    cluster_compare_indx_Round2 <- sample(subremaining1_indx[-subremaining1_train_order], subsampling, replace=F )
  } else {
    cluster_compare_indx_Round2 <- subremaining1_indx[-subremaining1_train_order]
  }

  temp_S2 <-mixedpop1[genes_in_both_idx1, c(cluster_select_indx_Round2  , cluster_compare_indx_Round2 )]
  temp_S2 <-assay(temp_S2)
  predictor_S2 <-t(temp_S2)

  #start prediction for estimating accuracy
  predict_clusters<-predict(cvfit, newx = predictor_S2,  type = "class", s = cvfit$lambda.min)
  #Done cross validation test to estimate accuracy------------------------------

  #estimate accuracy------------------------------------------------------------
  # reference classification for comparing predicted results
  #this cellNames_cluster is to calculate the accuracy
  cellNames_cluster <- cbind(colnames(mixedpop1), colData(mixedpop1)[,1] )
  #from here compare to original clusering classes to check for accuracy (just for the training set)
  #it uses the initial cellNames_cluster
  predictor_S2_name <-row.names(predict_clusters)
  predict_label <-predict_clusters[,1]

  predict_index <-which(cellNames_cluster[,1] %in%   predictor_S2_name)

  original_cluster <- cellNames_cluster[predict_index,]
  original_cluster <-original_cluster[order(original_cluster[,2], decreasing = T),]
  original_cluster <-as.data.frame(original_cluster)
  predict_cluster_dat <-as.data.frame(predict_label)
  predict_cluster_dat$cellnames <-row.names(predict_cluster_dat)

  compare <-merge(original_cluster, predict_cluster_dat, by.x='V1', by.y='cellnames')
  #this is the accurate prediction
  cluster_select_predict <- subset(compare,(as.numeric(compare$V2) ==
     c_selectID & compare$predict_label== c_selectID) | ( as.numeric(compare$V2) !=
      c_selectID & compare$predict_label != c_selectID))

  accurate <-dim(cluster_select_predict)[1]
  inaccurate <-dim(compare)[1] - dim(cluster_select_predict)[1]

  list_acc_inacc<- list(accurate, inaccurate)

  #done estimate accuracy-------------------------------------------------------

  #to write the 5 lists into the object

  listData$Accuracy[[out_idx]] <- list(list_acc_inacc)
  listData$LassoGenes[[out_idx]] <- list(sub_cvfit_out)
  listData$Deviance[[out_idx]] <- list(dat_DE_fm_DE)
  listData$LassoFit[[out_idx]] <- list(cvfit)
  listData$LDAFit[[out_idx]] <- list(fit.lda)
  listData$predictor_S1[[out_idx]] <- list(predictor_S1)

  return(listData)
}


#' Main prediction function applying the optimal LASSO and LDA models
#'
#' @description  Predict a new mixed population after training the model for a
#' subpopulation in the first mixed population.
#' All subpopulations in the new target mixed population will be predicted, where
#' each targeted subpopulation will have a transition score from the orginal
#' subpopulation to the new subpopulation.
#' @param listData a \code{list} object containing trained results for the
#' selected subpopulation in the first mixed population
#' @param mixedpop2 a \linkS4class{SingleCellExperiment} object from the target
#' mixed population of importance, e.g. differentially expressed genes that are most significant
#' @param out_idx a number to specify index to write results into the list output.
#' This is needed for running bootstrap.
#' @return a \code{list} with prediction results written in to the index \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' c_selectID<-1
#' out_idx<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                      CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                      CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' listData  <- training_scGPS(genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID, listData =list(), out_idx=out_idx)
#' listData  <- predicting_scGPS(listData =listData,  mixedpop2 = mixedpop2, out_idx=out_idx)

predicting_scGPS <-function(listData = NULL,  mixedpop2 = NULL, out_idx=NULL){
  #predictor_S1 is the dataset used for the training phase (already transposed)
  predictor_S1 <- listData$predictor_S1[[out_idx]][[1]] #1 for extract matrix
  cvfit_best <- listData$LassoFit[[out_idx]][[1]]
  fit.lda<- listData$LDAFit[[out_idx]][[1]]

  my.clusters <- colData(mixedpop2)[,1]
  ori_dat_2 <- assay(mixedpop2)
  names <-elementMetadata(mixedpop2)$GeneSymbol

  list_predict_clusters_LASSO <- list()
  list_predict_clusters_LDA <- list()

  for (clust in unique(my.clusters)){
    c_selectID_2 <- clust
    cluster_select <- which(my.clusters == as.numeric(c_selectID_2)) #select cells
    dataset <- predictor_S1

        #for getting genes
    gene_cvfit <- cvfit_best$glmnet.fit$beta@Dimnames[[1]]

    #reformat the names (to do: may write functions to check format)
    gene_cvfit <- gsub("_.*", '', gene_cvfit )
    cvfitGenes_idx <-which(names %in% gene_cvfit)
    #check how many genes in gene_cvfit  but not in mixedpop2
    #(should be 0, as this is checked in the training step)
    to_add <-length(gene_cvfit) - length(cvfitGenes_idx)
    to_add_idx <- c()

    if(to_add >0) {
      to_add_idx <- sample(cvfitGenes_idx, to_add, replace=F)} else if (to_add < 0) {
        cvfitGenes_idx <-sample(cvfitGenes_idx, length(gene_cvfit), replace=F)}

    predictor_S2_temp <-ori_dat_2[c(to_add_idx, cvfitGenes_idx),cluster_select]

    #predict LASSO:
    predict_clusters_LASSO <- predict(cvfit_best, newx = t(predictor_S2_temp),  type = "class", s = cvfit_best$lambda.min)
    LASSO_result <- as.data.frame(table(predict_clusters_LASSO )) #convert table() to 2x2 dataframe, it will always have 2 variable names: $name, Freq
    LASSO_cluster_idx <- which(LASSO_result[,1]== c_selectID)
    predict_clusters_LASSO <-as.numeric(LASSO_result[LASSO_cluster_idx,2])/sum(as.numeric(LASSO_result[,2]))*100
    #print out the prediction output
    report_LASSO <- paste0('LASSO for subpop', c_selectID_2, ' in target mixedpop2')
    predict_clusters_LASSO <- list( report_LASSO, predict_clusters_LASSO)
    list_predict_clusters_LASSO <-c(list_predict_clusters_LASSO, predict_clusters_LASSO )

    #predict LDA:
    newdataset <-t(predictor_S2_temp)
    newdataset <-as.data.frame(newdataset)
    predict_clusters_LDA <- predict(fit.lda, newdataset)
    LDA_result <- as.data.frame(table(predict_clusters_LDA )) #convert table() to 2x2 dataframe, it will always have 2 variable names: $name, Freq
    LDA_cluster_idx <- which(LDA_result[,1] == c_selectID)
    predict_clusters_LDA <-as.numeric(LDA_result[LDA_cluster_idx,2])/sum(as.numeric(LDA_result[,2]))*100

    #write prediction output
    report_LDA <- paste0('LDA for subpop ', c_selectID_2, ' in target mixedpop2')
    predict_clusters_LDA <- list(report_LDA, predict_clusters_LDA)
    list_predict_clusters_LDA <-c(list_predict_clusters_LDA, predict_clusters_LDA)
  }

  #to write prediction result
  listData$LassoPredict[out_idx] <-  list(list_predict_clusters_LASSO)
  listData$LDAPredict[out_idx] <-  list(list_predict_clusters_LDA)

  return(listData)
}


#' BootStrap runs for both scGPS training and prediction
#'
#' @description  LASSO and LDA prediction for each of all the subpopulations in
#' the new mixed population after training the model for a subpopulation in the
#' first mixed population. The number of bootstraps to be run can be specified.
#' @seealso \code{\link{bootstrap_scGPS_parallel}} for parallel options
#' @param listData  a \code{list} object, which contains trained results for the first mixed population
#' @param mixedpop1 a \linkS4class{SingleCellExperiment} object from a mixed population for training
#' @param mixedpop2 a \linkS4class{SingleCellExperiment} object from a target mixed population for prediction
#' @param genes a gene list to build the model
#' @param nboots a number specifying how many bootstraps to be run
#' @return a \code{list} with prediction results written in to the index \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' c_selectID<-1
#' out_idx<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                      CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                      CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' test <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#' names(test)
#' test$LassoPredict
#' test$LDAPredict

bootstrap_scGPS <- function(nboots = 1, genes = genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID, listData =list()){

  for (out_idx in 1:nboots){
    listData  <- training_scGPS(genes =genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID, listData =listData, out_idx=out_idx)
    print(paste0("done ", out_idx))
    listData  <- predicting_scGPS(listData =listData,  mixedpop2 = mixedpop2, out_idx=out_idx)
  }
  return(listData)
}


#' BootStrap runs for both scGPS training and prediction with parallel option
#'
#' @description  same as bootstrap_scGPS, but with an multicore option
#' @param listData  a \code{list} object, which contains trained results for the first mixed population
#' @param mixedpop1 a \linkS4class{SingleCellExperiment} object from a mixed population for training
#' @param mixedpop2 a \linkS4class{SingleCellExperiment} object from a target mixed population for prediction
#' @param genes a gene list to build the model
#' @param nboots a number specifying how many bootstraps to be run
#' @param ncores a number specifying how many cpus to be used for running
#' @return a \code{list} with prediction results written in to the index \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#'
#'


bootstrap_scGPS_parallel <- function(ncores=4, nboots = 1, genes = genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID, listData =list()) {

  cl <- makeCluster(ncores)

  #clusterExport(cl, c("listData1", "mixedpop1", "mixpop2", "out_idx"))
  clusterExport(cl, c( "mixedpop1","mixedpop2", "c_selectID", "listData1", "genes",
                       "training_scGPS", "predicting_scGPS"))

  clusterEvalQ(cl, library("glmnet", "caret", "scGPS",
                           "SingleCellExperiment"))
  parLapply(cl, 1:nboots, function(out_idx){

    listData  <- training_scGPS(genes =genes, mixedpop1 = mixedpop1,
    mixedpop2 = mixedpop2, c_selectID, listData =listData1, out_idx=out_idx)

    print(paste0("done ", out_idx))

    listData  <- predicting_scGPS(listData =listData,  mixedpop2 = mixedpop2, out_idx=out_idx)

  }
)
  stopCluster(cl)
  return(listData)
}





