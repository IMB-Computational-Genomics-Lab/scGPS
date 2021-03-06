% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimReduction.R
\name{PCA}
\alias{PCA}
\title{PCA}
\usage{
PCA(expression.matrix = NULL, ngenes = 1500, scaling = TRUE, npcs = 50)
}
\arguments{
\item{expression.matrix}{An expression matrix, with genes in rows}

\item{ngenes}{number of genes used for clustering calculations.}

\item{scaling}{a logical of whether we want to scale the matrix}

\item{npcs}{an integer specifying the number of principal components to use.}
}
\value{
a list containing PCA results and variance explained
}
\description{
Select top variable genes and perform prcomp
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
t <-PCA(expression.matrix=assay(mixedpop1))

}
