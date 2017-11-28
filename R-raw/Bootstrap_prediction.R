
#' Bootstrap run
#'
#' @description  this function takes cvfit as the input and predict cells in another subpopulation
#' @param cvfit the LASSO fit object
#' @param predictor_S1 a predictor matrix
#' @param dayID a character specifying day ID
#' @return a list containing:list_acc_inacc,  list_SigGenes, list_Deviance, list_cvFit,  list_predict_sameDay, list_predict_DifferentDay
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

bootstrap_run <- function(bootstrap_number=100){
  #Lists to be saved. Note that only list_predict_sameDay and list_predict_DifferentDay are day-specific; all others are shared for the same models LASSO or LDA
  list_SigGenes = list() #for extracting genes post LASSO
  list_Deviance = list() #for plotting deviance explained
  list_cvFit = list() #for saving lasso model
  list_lda = list() #for saving lda model

  list_acc_inacc = list() #for plotting accuracy of a model
  list_predict_sameDay = list() #for percent predicted by LASSO
  list_predict_DifferentDay = list() #for percent predicted by LDA

  for (i in 1:bootstrap_number){

    cluster_select_indx_S1 <-sample(cluster_select, subsampling, replace=F)
    predict_marker<- Lit_New_Lasso(cluster_select_indx_S1 = cluster_select_indx_S1)

    predictor_S2_name <-predict_marker[[2]]
    predict_label <-predict_marker[[1]]

    #from here compare to original clusering classes to check for accuracy (just for the training set)
    #it uses the initial cellNames_cluster
    predict_index <-which(cellNames_cluster[,1] %in% row.names(predictor_S2_name))
    original_cluster <- cellNames_cluster[predict_index,]
    original_cluster <-original_cluster[order(original_cluster[,2], decreasing = T),]
    original_cluster <-as.data.frame(original_cluster)
    predict_clusters <-as.data.frame(predict_label)
    predict_clusters$cellnames <-row.names(predict_clusters)

    compare <-merge(original_cluster, predict_clusters, by.x='V1', by.y='cellnames')
    cluster_select_predict <- subset(compare, (as.numeric(compare$my.clusters) == c_selectID &  compare$`1`==c_selectID) | (compare$my.clusters!=c_selectID &  compare$`1`!=c_selectID)) #this is the accurate prediction
    accurate <-dim(cluster_select_predict)[1]
    inaccurate <-dim(compare)[1] - dim(cluster_select_predict)[1]

    list_acc_inacc[[i]] <- list(accurate, inaccurate)
    list_SigGenes[[i]] <-predict_marker[[3]]
    list_Deviance[[i]] <-predict_marker[[4]]
    list_cvFit[[i]]  <-predict_marker[[5]]
    list_lda[[i]] = predict_marker[[7]]

    #More new prediction
    cvfit=predict_marker[[5]]
    fit.lda =predict_marker[[7]]
    #predict lasso
    list_predict_sameDay[[i]] <- predict_subpop(dayID = dayID, predictor_S1 = predict_marker[[6]], cvfit= predict_marker[[5]])
    #predict LDA
    list_predict_DifferentDay[[i]] <- predict_subpop(dayID = dayID2, predictor_S1 = predict_marker[[6]], cvfit= predict_marker[[5]])

  }
  list_all <-list(list_acc_inacc =list_acc_inacc,  list_SigGenes = list_SigGenes, list_Deviance = list_Deviance,  list_cvFit = list_cvFit,  list_predict_sameDay = list_predict_sameDay, list_predict_DifferentDay = list_predict_DifferentDay)
  saveRDS(list_all, paste0('LASSO_LDA_DEgenes_cluster', c_selectID,'fromDay',dayID,'_toDay_', dayID2, '.RDS'))
}





