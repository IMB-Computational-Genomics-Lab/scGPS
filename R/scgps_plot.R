#' @name summary_accuracy
#' @title get percent accuracy for Lasso model, from \code{n} bootstraps
#' @description The training results from \code{training_scGPS} were written to
#' the object \code{listData}, the \code{summary_accuracy} summarise \code{n} bootstraps
#' @return a vector of percent accuracy for the selected subpopulation
#' @docType data
#' @usage sample2
#' @author Quan Nguyen, 2017-11-25


summary_accuracy <-function(object=LSOLDA_dat){
  acc_inacc <- object$Accuracy
  pcAcc <-as.vector(unlist(lapply(acc_inacc, function(x){x[[1]][[1]]/(x[[1]][[1]] +x[[1]][[2]])*100})))
  return(pcAcc)
}


