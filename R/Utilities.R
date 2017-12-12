
#' select top variable genes
#' @description subset a matrix by top variable genes
#' @param expression.matrix is a matrix with genes in rows and cells in columns
#'
topvar_scGPS <- function(expression.matrix=NULL,ngenes = 1500){
  CalcRowVariance <- function(x) {
  row.variance <- rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
  return(row.variance)
}
  gene.variance <- CalcRowVariance(expression.matrix)
  names(gene.variance) <- rownames(expression.matrix)
  sorted.gene.variance <- gene.variance[order(gene.variance, decreasing = TRUE)]
  top.genes <- sorted.gene.variance[1:ngenes]
  subset.matrix <- expression.matrix[names(top.genes), ]
  return(subset.matrix)
}

#' plot reduced data
#' @description plot PCA, tSNE, and CIDR reduced datasets
#' @param reduced_dat is a matrix with genes in rows and cells in columns
#' @examples
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
#'                     CellMetadata = day2$dat2_clusters)
#' CIDR_dim <-CIDR_scGPS(expression.matrix=assay(mixedpop1))
#' p <-plotReduced_scGPS(CIDR_dim)
#' plot(p)
#' tSNE_dim <-tSNE_scGPS(expression.matrix=assay(mixedpop1))
#' p2 <-plotReduced_scGPS(tSNE_dim)
#' plot(p2)
#'
#'
plotReduced_scGPS <- function(reduced_dat, dims = c(1,2),
  dimNames = c("Dim 1", "Dim 2")){
  reduced_dat_toPlot <-as.data.frame(reduced_dat[,dims])
  colnames(reduced_dat_toPlot) <- dimNames
  p <- qplot(reduced_dat[,dims[1]], reduced_dat[,dims[2]], geom = "point")
  p <- p + ylab(dimNames[1]) + xlab(dimNames[2])
  return(p)
}


#' build disance matrice V2.0
#'
#' @description  Calculate distance matrix and perform hclust
#' @param mixedpop1 is a \linkS4class{SingleCellExperiment} object from the train mixed population
#' @return a \code{matrix} with Eucleadean distance used for clustreting
#' @export
#' @author Quan Nguyen, 2017-11-25
#'

distance_scGPS <-function(object = NULL){
  print("Calculating distance matrix")
  exprs_mat <- assay(object)
  exprs_mat_topVar <- topvar_scGPS(exprs_mat, ngenes = 1500)
  #Take the top variable genes
  #make a transpose
  exprs_mat_t <-t(exprs_mat_topVar)
  dist_mat <- dist(exprs_mat_t)
  print("Performing hierarchical clustering")
  return(list("tree"=original.tree, "ref_clust" = original.clusters))
}
