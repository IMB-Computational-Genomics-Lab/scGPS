#' @name summary_accuracy
#' @title get percent accuracy for Lasso model, from \code{n} bootstraps
#' @description The training results from \code{training_scGPS} were written to
#' @param object is a list containing the training results from 
# \code{training_scGPS} the object \code{LSOLDA_dat},
#' the \code{summary_accuracy} summarise \code{n} bootstraps
#' @return a vector of percent accuracy for the selected subpopulation
#' @export
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, 
#'     mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
#'     cluster_mixedpop1 = colData(mixedpop1)[,1],
#'     cluster_mixedpop2=colData(mixedpop2)[,1])
#' summary_accuracy(LSOLDA_dat)
#' summary_deviance(LSOLDA_dat)
#' @author Quan Nguyen, 2017-11-25


summary_accuracy <- function(object = NULL) {
    acc_inacc <- object$Accuracy
    pcAcc <- as.vector(unlist(lapply(acc_inacc, function(x) {
        x[[1]][[1]]/(x[[1]][[1]] + x[[1]][[2]]) * 100
    })))
    return(pcAcc)
}

#' @name summary_deviance
#' @title get percent deviance explained for Lasso model, 
#' from \code{n} bootstraps
#' @description the training results from \code{training_scGPS} were written to
#' the object \code{LSOLDA_dat}, the \code{summary_devidance} summarises 
#' deviance explained for \code{n} bootstrap runs and also returns the best
#' deviance matrix for plotting, as well as the best matrix with Lasso genes 
#' and coefficients
#' @param object is a list containing the training results from 
#' \code{training_scGPS}
#' @return a \code{list} containing three elements, with a vector of percent
#' maximum deviance explained, a dataframe containg information for the full 
#' deviance, and a dataframe containing gene names and coefficients of the best 
#' model
#' @export
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo,
#'                     CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, 
#'     mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
#'     cluster_mixedpop1 = colData(mixedpop1)[,1],
#'     cluster_mixedpop2=colData(mixedpop2)[,1])
#' summary_deviance(LSOLDA_dat)
#' @author Quan Nguyen, 2017-11-25

summary_deviance <- function(object = NULL) {
    deviDat <- object$Deviance
    deviVec <- as.vector(unlist(lapply(deviDat, function(x) {
        temp <- x[[1]]$Deviance
        temp_max <- temp[length(temp) - 1]
    })))
    deviVec_max <- which(deviVec == max(deviVec, na.rm = TRUE))
    # get gene info
    genesSig <- object$LassoGenes
    
    GeneNames_max <- genesSig[[deviVec_max]][[1]]
    return(list(allDeviance = deviVec, DeviMax = deviDat[[deviVec_max]][[1]], LassoGenesMax = GeneNames_max))
}


#' @name summary_prediction_lasso
#' @title get percent deviance explained for Lasso model, from \code{n} 
#' bootstraps
#' @description the training results from \code{training_scGPS} were written to
#' the object \code{LSOLDA_dat}, the \code{summary_prediction} summarises 
#' prediction for \code{n} bootstrap runs
#' @param LSOLDA_dat is a list containing the training results from 
#' \code{training_scGPS}
#' @param nPredSubpop is the number of subpopulations in the target mixed 
#' population
#' @return a dataframe containg information for the Lasso prediction 
#' results, each column
#' contains prediction results for all subpopulations from each bootstrap run
#' @export
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, 
#'     mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
#'     cluster_mixedpop1 = colData(mixedpop1)[,1],
#'     cluster_mixedpop2=colData(mixedpop2)[,1])
#' summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop=4)
#'


summary_prediction_lasso <- function(LSOLDA_dat = NULL, nPredSubpop = NULL) {
    pred_lasso <- LSOLDA_dat$ElasticNetPredict
    for (i in 1:length(pred_lasso)) {
        for (j in 1:length(pred_lasso[[i]])) {
            if (identical(pred_lasso[[i]][[j]], numeric(0))) {
                pred_lasso[[i]][[j]] <- "NA"
            }
        }
    }
    
    pred_lasso_tranformed <- as.vector(unlist(pred_lasso))
    
    
    toremove <- grep("target", pred_lasso_tranformed)
    
    pred_lasso_percentOnly <- pred_lasso_tranformed[-toremove]
    
    pred_lasso_mtrx <- matrix(pred_lasso_percentOnly, nrow = nPredSubpop, byrow = FALSE)
    
    row_names <- pred_lasso_tranformed[toremove[c(1:nPredSubpop)]]
    
    pred_lasso_mtrx <- as.data.frame(pred_lasso_mtrx)
    
    pred_lasso_mtrx$names <- row_names
    
    return(pred_lasso_mtrx)
}

#' @name summary_prediction_lda
#' @title get percent deviance explained for LDA model, from \code{n} bootstraps
#' @description the training results from \code{training_scGPS} were written to
#' the object \code{LSOLDA_dat}, the \code{summary_prediction} summarises 
#' prediction explained for \code{n} bootstrap runs and also returns the best
#' deviance matrix for plotting, as well as the best matrix with Lasso genes
#' and coefficients
#' @param LSOLDA_dat is a list containing the training results from 
#' \code{training_scGPS}
#' @param nPredSubpop is the number of subpopulations in the target mixed 
#' population
#' @return a dataframe containg information for the LDA prediction 
#' results, each column contains prediction results for all subpopulations from
#' each bootstrap run
#' @export
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#' GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, 
#' mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
#'     cluster_mixedpop1 = colData(mixedpop1)[,1],
#'     cluster_mixedpop2=colData(mixedpop2)[,1])
#' summary_prediction_lda(LSOLDA_dat=LSOLDA_dat, nPredSubpop=4)
#'

summary_prediction_lda <- function(LSOLDA_dat = NULL, nPredSubpop = NULL) {
    pred_lda <- LSOLDA_dat$LDAPredict
    for (i in 1:length(pred_lda)) {
        for (j in 1:length(pred_lda[[i]])) {
            if (identical(pred_lda[[i]][[j]], numeric(0))) {
                pred_lda[[i]][[j]] <- "NA"
            }
        }
    }
    
    pred_lda_tranformed <- as.vector(unlist(pred_lda))
    
    toremove <- grep("target", pred_lda_tranformed)
    
    pred_lda_percentOnly <- pred_lda_tranformed[-toremove]
    
    pred_lda_mtrx <- matrix(pred_lda_percentOnly, nrow = nPredSubpop, byrow = FALSE)
    
    row_names <- pred_lda_tranformed[toremove[c(1:nPredSubpop)]]
    
    pred_lda_mtrx <- as.data.frame(pred_lda_mtrx)
    
    pred_lda_mtrx$names <- row_names
    
    return(pred_lda_mtrx)
}


#' @name reformat_LASSO
#' @title summarise bootstrap runs for Lasso model, from \code{n} bootstraps
#' @description the training and prediction results from \code{bootstrap_scGPS}
#' were written to the object \code{LSOLDA_dat}, the \code{reformat_LASSO}
#' summarises prediction for \code{n} bootstrap runs
#' @param c_selectID  is the original cluster to be projected
#' @param mp_selectID  is the target mixedpop to project to
#' @param LSOLDA_dat is the results from the bootstrap_scGPS
#' @param nPredSubpop is the number of clusters in the target mixedpop 
#' \code{row_cluster <-length(unique(target_cluster))}
#' @param Nodes_group string representation of hexidecimal color code for node
#' @param nboots is an integer for how many bootstraps are run
#' @return a dataframe containg information for the Lasso prediction results, 
#' each column
#' contains prediction results for all subpopulations from each bootstrap run
#' @export
#' @examples
#' c_selectID<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, 
#'     mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
#'     cluster_mixedpop1 = colData(mixedpop1)[,1],
#'     cluster_mixedpop2=colData(mixedpop2)[,1])
#' reformat_LASSO(LSOLDA_dat=LSOLDA_dat, 
#'     nPredSubpop=length(unique(colData(mixedpop2)[,1])), c_selectID = 1, 
#'     mp_selectID =2)
#'
reformat_LASSO <- function(c_selectID = NULL, mp_selectID = NULL, LSOLDA_dat = NULL, 
    nPredSubpop = NULL, Nodes_group = "#7570b3", nboots = 2) {
    LASSO_out <- summary_prediction_lasso(LSOLDA_dat = LSOLDA_dat, nPredSubpop = nPredSubpop)
    LASSO_out <- as.data.frame(LASSO_out)
    temp_name <- gsub("ElasticNet for subpop", "C", LASSO_out$names)
    temp_name <- gsub(" in target mixedpop", "_MP", temp_name)
    LASSO_out$names <- temp_name
    source <- rep(paste0("C", c_selectID, "_MP", mp_selectID), length(temp_name))
    LASSO_out$Source <- source
    LASSO_out$Node <- source
    LASSO_out$Nodes_group <- rep(Nodes_group, length(temp_name))
    colnames(LASSO_out) <- c(paste0("Boostrap", 1:nboots), "Target", "Source", "Node", 
        "NodeGroup")
    matrx_mean <- apply(LASSO_out[, c(1:nboots)], 1, function(x) {
        mean(as.numeric(x))
    })
    
    LASSO_out$Value <- matrx_mean
    return(LASSO_out)
}



