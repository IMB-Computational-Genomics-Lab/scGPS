#' @name summary_accuracy
#' @title get percent accuracy for Lasso model, from \code{n} bootstraps
#' @description The training results from \code{training_scGPS} were written to
#' @param LSOLDA_dat is a list containing the training results from \code{training_scGPS}
#' the object \code{LSOLDA_dat}, the \code{summary_accuracy} summarise \code{n} bootstraps
#' @return a vector of percent accuracy for the selected subpopulation
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                     CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#' summary_accuracy(LSOLDA_dat)
#' summary_deviance(LSOLDA_dat)
#' @author Quan Nguyen, 2017-11-25


summary_accuracy <-function(object=LSOLDA_dat){
  acc_inacc <- object$Accuracy
  pcAcc <-as.vector(unlist(lapply(acc_inacc, function(x){x[[1]][[1]]/(x[[1]][[1]] +x[[1]][[2]])*100})))
  return(pcAcc)
}

#' @name summary_deviance
#' @title get percent deviance explained for Lasso model, from \code{n} bootstraps
#' @description the training results from \code{training_scGPS} were written to
#' the object \code{LSOLDA_dat}, the \code{summary_devidance} summarises deviance explained
#' for \code{n} bootstrap runs and also returns the best deviance matrix for plotting, as
#' well as the best matrix with Lasso genes and coefficients
#' @param LSOLDA_dat is a list containing the training results from \code{training_scGPS}
#' @return a \code{list} containing three elements, with a vector of percent maximum
#' deviance explained, a dataframe containg information for the full deviance, and a
#' dataframe containing gene names and coefficients
#' of the best model
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                     CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#' summary_deviance(LSOLDA_dat)
#' @author Quan Nguyen, 2017-11-25

summary_deviance <-function(object=LSOLDA_dat){
  deviDat <- object$Deviance
  deviVec <-as.vector(unlist(lapply(deviDat, function(x){
    temp <- x[[1]]$Deviance
    temp_max <-temp[length(temp)-1]}
    )))
  deviVec_max <- which(deviVec== max(deviVec,na.rm = T))
  #get gene info
  genesSig <- LSOLDA_dat$LassoGenes

  GeneNames_max <- genesSig[[deviVec_max]][[1]]
  return(list("allDeviance" = deviVec, "DeviMax" = deviDat[[deviVec_max]][[1]],
              "LassoGenesMax" = GeneNames_max))
}


#' @name summary_prediction_lasso
#' @title get percent deviance explained for Lasso model, from \code{n} bootstraps
#' @description the training results from \code{training_scGPS} were written to
#' the object \code{LSOLDA_dat}, the \code{summary_prediction} summarises prediction
#' for \code{n} bootstrap runs
#' @param LSOLDA_dat is a list containing the training results from \code{training_scGPS}
#' @param nPredSubpop is the number of subpopulations in the target mixed population
#' @return a dataframe containg information for the Lasso prediction results, each column
#' contains prediction results for all subpopulations from each bootstrap run
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                     CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#' summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop=4)


summary_prediction_lasso <-function(LSOLDA_dat=NULL, nPredSubpop=NULL){
  pred_lasso <- LSOLDA_dat$LassoPredict

  pred_lasso_tranformed <- as.vector(unlist(pred_lasso))

  toremove <- grep("target", pred_lasso_tranformed)

  pred_lasso_percentOnly <- pred_lasso_tranformed[-toremove]

  pred_lasso_mtrx <-matrix(pred_lasso_percentOnly, nrow=nPredSubpop, byrow=F)

  row_names <- pred_lasso_tranformed[toremove[c(1:nPredSubpop)]]

  pred_lasso_mtrx <-as.data.frame(pred_lasso_mtrx)

  pred_lasso_mtrx$names <-row_names

  return(pred_lasso_mtrx)
}

#' @name summary_prediction_lda
#' @title get percent deviance explained for LDA model, from \code{n} bootstraps
#' @description the training results from \code{training_scGPS} were written to
#' the object \code{LSOLDA_dat}, the \code{summary_prediction} summarises prediction explained
#' for \code{n} bootstrap runs and also returns the best deviance matrix for plotting, as
#' well as the best matrix with Lasso genes and coefficients
#' @param LSOLDA_dat is a list containing the training results from \code{training_scGPS}
#' @param nPredSubpop is the number of subpopulations in the target mixed population
#' @return a dataframe containg information for the LDA prediction results, each column
#' contains prediction results for all subpopulations from each bootstrap run
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                     CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#' summary_prediction_lda(LSOLDA_dat=LSOLDA_dat, nPredSubpop=4)
#'

summary_prediction_lda <-function(LSOLDA_dat=NULL, nPredSubpop=NULL){
  pred_lda <- LSOLDA_dat$LDAPredict

  pred_lda_tranformed <- as.vector(unlist(pred_lda))

  toremove <- grep("target", pred_lda_tranformed)

  pred_lda_percentOnly <- pred_lda_tranformed[-toremove]

  pred_lda_mtrx <-matrix(pred_lda_percentOnly, nrow=nPredSubpop, byrow=F)

  row_names <- pred_lda_tranformed[toremove[c(1:nPredSubpop)]]

  pred_lda_mtrx <-as.data.frame(pred_lda_mtrx)

  pred_lda_mtrx$names <-row_names

  return(pred_lda_mtrx)
}


