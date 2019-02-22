
# Constructor function for the scgps object of the SingleCellExperiment class
#'
#' NewscGPS
#'
#' \code{\link{NewscGPS}} generates a scGPS object in the 
#' \linkS4class{SingleCellExperiment} class for use with the scGPS package. This
#' object contains an expression matrix, associated metadata (cells, genes,
#' clusters). The data are expected to be normalised counts.
#' @param ExpressionMatrix An expression matrix in data.frame or matrix format.
#' Rows should represent a transcript and its normalised counts,
#' while columns should represent individual cells.
#' @param GeneMetadata A data frame or vector containing gene identifiers used 
#' in the expression matrix. The first column should hold the gene identifiers
#' you are using in the expression matrix. Other columns contain information 
#' about the genes, such as their corresponding ENSEMBL transcript identifiers.
#' @param CellMetadata A data frame containing cell identifiers 
#' (usually barcodes) and an integer representing which batch they belong to.
#' The column containing clustering information needs to be the first column in 
#' the CellMetadata dataframe If clustering information is not available, users
#' can run CORE function and add the information to the scGPS before running
#' scGPS prediction
#' @return This function generates an scGPS object belonging to the 
#' \linkS4class{SingleCellExperiment}.
#' @seealso \linkS4class{SingleCellExperiment}
#' @export
#' @examples
#' day2 <- sample1
#' t <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' colData(t); show(t); colnames(t)
#' @author Quan Nguyen, 2018-04-06
#'
#'
NewscGPS <- function(ExpressionMatrix = NULL, GeneMetadata = NULL, CellMetadata = NULL) {
    # Check that we have the essential arguments - an expression matrix
    arg.check <- list(ExpressionMatrix = missing(ExpressionMatrix), GeneMetadata = missing(GeneMetadata), 
        CellMetadata = missing(CellMetadata))
    if (any(arg.check == TRUE)) {
        missing.args <- names(which(arg.check == TRUE))
        msg <- sprintf("Please supply the following arguments: %s\n", as.character(unlist(missing.args)))
        stop(msg)
    }
    
    if (is.data.frame(ExpressionMatrix) == FALSE & is.matrix(ExpressionMatrix) == 
        FALSE) {
        # stop('Please supply an expression matrix in one of the following formats:
        # data.frame or matrix')
        stop(paste0("Please supply an expression matrix in one of the ", "following formats: data.frame or matrix"))
    }
    
    # Create a new scGPS object.
    scGPSset <- SingleCellExperiment(assays = list(counts = ExpressionMatrix), rowData = GeneMetadata, 
        colData = CellMetadata)
    
    # All clear, return the object
    return(scGPSset)
}



# Constructor function for the scgps object of the SummarizedExperiment class
#'
#' NewscGPS_SME
#'
#' \code{\link{NewscGPS}} generates a scGPS object in the 
#' \linkS4class{SingleCellExperiment} class for use with the scGPS package. This
#' object contains an expression matrix, associated metadata (cells, genes,
#' clusters). The data are expected to be normalised counts.
#' @param ExpressionMatrix An expression dataset in matrix format.
#' Rows should represent a transcript and its normalised counts,
#' while columns should represent individual cells.
#' @param GeneMetadata A data frame or vector containing gene identifiers used 
#' in the expression matrix. The first column should hold the cell identifiers
#' you are using in the expression matrix. Other columns contain information 
#' about the genes, such as their corresponding ENSEMBL transcript identifiers.
#' @param CellMetadata A data frame containing cell identifiers 
#' (usually barcodes) and clustering information (the first column of the data
#' frame contains clustering information). The column containing clustering
#' information needs to be named as 'Cluster'. If clustering information is not
#' available, users can run CORE function and add the information to the scGPS
#' before running scGPS prediction
#' @return This function generates an scGPS object belonging to the
#' \linkS4class{SingleCellExperiment}.
#' @seealso \linkS4class{SingleCellExperiment}
#' @export
#' @examples
#' day2 <- sample1
#' t <-NewscGPS_SME(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' colData(t); show(t); colnames(t)
#' @author Quan Nguyen, 2017-11-25
#'
#'
#'
NewscGPS_SME <- function(ExpressionMatrix = NULL, GeneMetadata = NULL, CellMetadata = NULL) {
    # Check that we have the essential arguments - an expression matrix
    arg.check <- list(ExpressionMatrix = missing(ExpressionMatrix), GeneMetadata = missing(GeneMetadata), 
        CellMetadata = missing(CellMetadata))
    if (any(arg.check == TRUE)) {
        missing.args <- names(which(arg.check == TRUE))
        msg <- sprintf("Please supply the following arguments: %s\n", as.character(unlist(missing.args)))
        stop(msg)
    }
    
    # Check data formats (to do: can add a series of checking here)
    if (is.matrix(ExpressionMatrix) == FALSE) {
        stop("Please supply an expression data in the matrix format")
    }
    
    # Create a new scGPS object.
    scGPSset <- SummarizedExperiment(assays = list(ExpressionMatrix), rowData = GeneMetadata, 
        colData = CellMetadata)
    
    # All clear, return the object
    return(scGPSset)
}

