#' PCA
#'
#' @description Select top variable genes and perform prcomp
#' @param expression.matrix An expression matrix, with genes in rows
#' @param ngenes number of genes used for clustering calculations.
#' @param npcs an integer specifying the number of principal components to use.
#' @param scaling a logical of whether we want to scale the matrix
#' @export
#' @return a list containing PCA results and variance explained
#' @examples
#' day2 <- sample1
#' mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' t <-PCA(expression.matrix=assay(mixedpop1))
#'


PCA <- function(expression.matrix = NULL, ngenes = 1500, scaling = TRUE, 
    npcs = 50) {
    
    message(paste0("Preparing PCA inputs using the top ", as.integer(ngenes),
        " genes ..."))
    subset.matrix <- top_var(expression.matrix = expression.matrix,
        ngenes = ngenes)
    
    # transpose to perform pca for cells in rows
    pca.input.matrix <- t(subset.matrix)
    # pca.input.matrix <- transposed.matrix[, 
    #     CalcColVariance(transposed.matrix) > 0]
    
    # Compute PCA
    message("Computing PCA values...")
    pca.result <- stats::prcomp(pca.input.matrix, scale = scaling)
    pca.percent.var <- pca.result$sdev^2/sum(pca.result$sdev^2)

    return(list(reduced_dat = pca.result$x[, 1:npcs], 
        variance = pca.percent.var))
}


#' tSNE
#'
#' @description calculate tSNE from top variable genes
#' @param expression.mat An expression matrix, with genes in rows
#' @param topgenes number of genes used for clustering calculations.
#' @param perp numeric; Perplexity parameter
#' (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for 
#' interpretation)
#' @param thet numeric; Speed/accuracy trade-off (increase for less accuracy)
#' @param scale a logical of whether we want to scale the matrix
#' @export
#' @return a tSNE reduced matrix containing three tSNE dimensions
#' @examples
#' day2 <- sample1
#' mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' t <-tSNE(expression.mat = assay(mixedpop1))
#'
tSNE <- function(expression.mat = NULL, topgenes = 1500, scale = TRUE, 
    thet = 0.5, perp = 30) {
    reducedDat <- PCA(expression.matrix = expression.mat, 
        ngenes = topgenes, scaling = scale)
    reducedDat <- reducedDat$reduced_dat
    message("Running tSNE ...")
    tsne <- Rtsne::Rtsne(reducedDat, pca = TRUE, perplexity = perp, 
        theta = thet, ignore_duplicates = TRUE)
    return(tsne$Y)
}

