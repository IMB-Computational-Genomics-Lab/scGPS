#' PCA
#'
#' @description Select top variable genes and perform prcomp
#' @param expression_matrix An expression matrix, with genes in rows
#' @export
#' @return a list containing PCA results and variance explained
#' @examples
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' t <-PCA_scGPS(expression.matrix=assay(mixedpop1))


PCA_scGPS <- function(expression.matrix = NULL, ngenes = 1500, scaling = TRUE, npcs = 50) {
    print(paste0("Preparing PCA inputs using the top ", as.integer(ngenes), " genes ..."))
    subset.matrix <- topvar_scGPS(expression.matrix = expression.matrix, ngenes = ngenes)

    # transpose to perform pca for cells in rows
    pca.input.matrix <- t(subset.matrix)
    # pca.input.matrix <- transposed.matrix[, CalcColVariance(transposed.matrix) > 0]

    # Compute PCA
    print("Computing PCA values...")
    if (scaling == TRUE) {
        pca.result <- stats::prcomp(pca.input.matrix, scale = TRUE)
    } else {
        pca.result <- stats::prcomp(pca.input.matrix, scale = FALSE)
    }

    pca.percent.var <- pca.result$sdev^2/sum(pca.result$sdev^2)
    return(list(reduced_dat = pca.result$x[, 1:npcs], variance = pca.percent.var))
}


#' CIDR
#'
#' @description calculate CIDR using top variable genes
#' @param expression_matrix An expression matrix, with genes in rows
#' @return a CIDR reduced matrix for the top 20 components
#' @export
#' @examples
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' t <-CIDR_scGPS(expression.matrix=assay(mixedpop1))
#'
CIDR_scGPS <- function(expression.matrix = NULL, ngenes = 1500, scaling = TRUE, theta = 0.5,
    perplexity = 30) {
    subset.matrix <- topvar_scGPS(expression.matrix = expression.matrix, ngenes = 1500)
    print("building cidr object...")
    sData <- scDataConstructor(subset.matrix)
    print("determine dropout candidates...")
    sData <- determineDropoutCandidates(sData)
    print("determine the imputation weighting threshold...")
    sData <- wThreshold(sData)
    print("computes the _CIDR_ dissimilarity matrix...")
    sData <- scDissim(sData)
    print("PCA plot with proportion of variance explained...")
    sData <- scPCA(sData)  #plotPC = FALSE)
    print("find the number of PC...")
    sData <- nPC(sData)
    print("perform clustering...")
    sData <- scCluster(sData)
    cidr_PCA_20 <- sData@PC[, 1:20]
    return(cidr_PCA_20)
}

#' tSNE
#'
#' @description calculate tSNE from top variable genes
#' @param expression_matrix An expression matrix, with genes in rows
#' @export
#' @return a tSNE reduced matrix containing three tSNE dimensions
#' @examples
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' t <-tSNE_scGPS(expression.matrix=assay(mixedpop1))
#'
tSNE_scGPS <- function(expression.matrix = NULL, ngenes = 1500, scaling = TRUE, theta = 0.5,
    perplexity = 30) {
    reducedDat <- PCA_scGPS(expression.matrix = expression.matrix, ngenes = ngenes,
        scaling = scaling)
    reducedDat <- reducedDat$reduced_dat
    print("Running tSNE ...")
    tsne <- Rtsne::Rtsne(reducedDat, pca = TRUE, perplexity = perplexity, theta = theta,
        ignore_duplicates = TRUE)
    return(tsne$Y)
}

