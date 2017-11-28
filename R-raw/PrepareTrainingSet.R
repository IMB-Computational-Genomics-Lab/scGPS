#' prepares two subpopulations
#'
#' @description  given two matrix inputs, prepare a training subpopulation to build a prediction model and
#' evaluate the model - subpop1. The resulting model is applied to predict the second subpopulation - subpop2
#' Two matrixes are used instead of one object because scGPS allows flexible prediction between any two expression matrix.
#' @param MixedPop1 the training subpop, a data frame, rows contain gene names,
#' columns contain cell IDs, can be an  \code{ExpressionSet} object
#' @param MixedPop2 the prediction subpop rows contain gene names,
#' columns contain cell IDs, can be an \code{ExpressionSet} object
#' @param genes a vector of gene symbols used for prediction, gene symbols must be
#' in the same format with gene names in subpop2
#' genes are listed by the order of importance
#' e.g. differentially expressed genes that are most significant
#' @return Two dataframes and a vector of genes in the right format ready for training 1 vs 2
#' @author QN
#' @Example
#' prepare_2_subpop(mat1, mat2, geneNames)
#' \dontrun{
#'
#' prepare_2_subpop(subpop1 = NULL,
#'                  subpop2 = NULL,
#'                  genes = "")if(!require(installr)) { install.packages("installr"); require(installr)} #load / install+load installr

#' }
#'

prepare_2_subpop <-function(subpop1,cluster1, subpop2, cluster2, genes){
  # making sure the gene list contains fewer than 500 genes
  # select top 500, the DE_result is an already sorted list
  if(length(genes) >500){genes <- genes[1:500]}

  names1 <-rownames(subpop1)
  subpop1_selected_genes <- names1[which(names1 %in% genes)]
  names2 <-rownames(subpop2)
  selected_genes_in_both <-names2[which(names2 %in% subpop1_selected_genes)]

  subpop1_train <- subpop1[which(names1 %in% subpop1_selected_genes)]
  subpop2_target <- subpop2[which(names2 %in% subpop1_selected_genes)]

  # reference classification for comparing predicted results
  cellNames_subpop <- cbind(colnames(subpop1), colnames(subpop2)) #note: may edit this to maje 2 genetic classes 1 and 2
}

