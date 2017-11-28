#' Prediction
#'
#' @description  this function takes cvfit as the input and predict cells in another subpopulation
#' @param cvfit the LASSO fit object
#' @param predictor_S1 a predictor matrix
#' @param dayID a character specifying day ID
#' @return a list containing:	list_predict_clusters_LASSO and list_predict_clusters_LDA
#' @keywords Input
#' @author QN
#' @Example
#' prepare_2_subpop(mat1, mat2, geneNames)
#' \dontrun{
#'
#' predict_subpop(
#'                  dayID = NULL,
#'                  predictor_S1 = NULL,
#'                  cvfit = NULL)
#' }
#'

predict_subpop <-function(dayID, predictor_S1, cvfit ){
  dayIDtarget <- dayID
  my.clusters <-readRDS(paste0("my.clusters_0.45_day",dayIDtarget, ".RDS"))
  ori_dat_2 <- readRDS(paste0("../Exprs_DCVLnorm_unlog_minus1_pos_Day", dayIDtarget, ".RDS"))
  names <-rownames(ori_dat_2)
  names <-gsub("_.*", '', names)
  rownames(ori_dat_2) <-names

  list_predict_clusters_LASSO <- list()
  list_predict_clusters_LDA <- list()

  for (clust in unique(my.clusters)){
    c_selectID_2 <-clust
    cluster_select <- which(my.clusters == as.numeric(c_selectID_2)) #select cells
    dataset <- predictor_S1 #predictor_S1 is the dataset used for the training phase (already transposed)
    gene_cvfit <- colnames(dataset)
    cvfitGenes_idx <-which(rownames(ori_dat_2) %in% gene_cvfit) #gene idx

    #predictor_S2_temp <-ori_dat_2[c(cvfitGenes_idx),cluster_select]

    #for getting genes
    cvfit_best <- cvfit
    gene_cvfit <- cvfit_best$glmnet.fit$beta@Dimnames[[1]]
    cvfitGenes_idx <-which(rownames(ori_dat_2) %in% gene_cvfit)

    to_add <-length(gene_cvfit) - length(cvfitGenes_idx)
    to_add_idx <- c()
    if(to_add >0) {
      to_add_idx <- sample(cvfitGenes_idx, to_add, replace=F)} else if (to_add <0) {
        cvfitGenes_idx <-sample(cvfitGenes_idx, length(gene_cvfit),replace=F)}

    predictor_S2_temp <-ori_dat_2[c(to_add_idx, cvfitGenes_idx),cluster_select]

    #predict LASSO:
    predict_clusters_LASSO <- predict(cvfit_best, newx = t(predictor_S2_temp),  type = "class", s = cvfit_best$lambda.min)
    LASSO_result <- as.data.frame(table(predict_clusters_LASSO )) #convert table() to 2x2 dataframe, it will always have 2 variable names: $name, Freq
    LASSO_cluster_idx <- which(LASSO_result[,1]== c_selectID)
    predict_clusters_LASSO <-as.numeric(LASSO_result[LASSO_cluster_idx,2])/sum(as.numeric(LASSO_result[,2]))*100
    #print out the prediction output
    report_LASSO <- paste0('LASSO From cluster', c_selectID, ' Day ', dayID, ' to predict cluster ',c_selectID_2, 'Day ', dayIDtarget)
    predict_clusters_LASSO <- list(predict_clusters_LASSO, report_LASSO)
    list_predict_clusters_LASSO <-c(list_predict_clusters_LASSO, predict_clusters_LASSO )
    #cluster_select_predict <- subset(compare, (as.numeric(compare$my.clusters) == c_selectID &  compare$`1`==c_selectID) | (compare$my.clusters!=c_selectID &  compare$`1`!=c_selectID))

    #predict LDA:
    newdataset <-t(predictor_S2_temp)
    newdataset <-as.data.frame(newdataset)

    predict_clusters_LDA <- predict(fit.lda, newdataset)
    LDA_result <- as.data.frame(table(predict_clusters_LDA )) #convert table() to 2x2 dataframe, it will always have 2 variable names: $name, Freq
    LDA_cluster_idx <- which(LDA_result[,1] == c_selectID)
    predict_clusters_LDA <-as.numeric(LDA_result[LDA_cluster_idx,2])/sum(as.numeric(LDA_result[,2]))*100

    #write prediction output
    report_LDA <- paste0('LDA From cluster', c_selectID, ' Day ', dayID, ' to predict cluster ',c_selectID_2, 'Day ', dayIDtarget)
    predict_clusters_LDA <- list(predict_clusters_LDA, report_LDA)
    list_predict_clusters_LDA <-c(list_predict_clusters_LDA, predict_clusters_LDA)
  }
  return(list(predict_percent_LASSO =	list_predict_clusters_LASSO, predict_clusters_LDA = list_predict_clusters_LDA))
}
