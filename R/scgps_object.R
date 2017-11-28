
# Constructor function for scgps obje SingleCellExperiment class to contain scgps inputs and outputs
#'
#' NewscGPS
#'
#' \code{\link{NewscGPS}} generates a \linkS4class{SingleCellExperiment} object for use with the scGPS package. This object contains an expression matrix, associated metadata (cells, genes, clusters), and bootstrap results. elementMetadata refer to genes/genomic features.
#' @param ExpressionMatrix An expression matrix in data.frame or matrix format. Rows should represent a transcript and its normalised counts, while columns should represent individual cells.
#' @param GeneMetadata A data frame or vector containing gene identifiers used in the expression matrix. The first column should hold the cell identifiers you are using in the expression matrix. Other columns contain information about the genes, such as their corresponding ENSEMBL transcript identifiers, whether or not they are a control and any additional information supplied by the user. This is an optional field.
#' @param CellMetadata A data frame containing cell identifiers (usually barcodes) and an integer representing which batch they belong to. This is an optional field, and it best used for experiments that contain data from multiple samples. This data frame can also hold additional information supplied by the user.
#' @return This function generates an object belonging to the \linkS4class{SingleCellExperiment}.
#' @seealso \linkS4class{SingleCellExperiment}
#' @export
#' @examples
#' day2 <- sample1
#' t <-NewscGPS(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' colData(t); show(t); colnames(t)
#' @author Quan Nguyen, 2017-11-25
#'
NewscGPS <- function(ExpressionMatrix = NULL, GeneMetadata= NULL, CellMetadata = NULL) {
  # Check that we have the essential arguments - an expression matrix
  arg.check <- list(ExpressionMatrix = missing(ExpressionMatrix),
                    GeneMetadata =  missing(GeneMetadata),
                    CellMetadata = missing(CellMetadata))
  if (any(arg.check == TRUE)) {
    missing.args <- names(which(arg.check == TRUE))
    msg <- sprintf("Please supply the following arguments: %s\n", as.character(unlist(missing.args)))
    stop(msg)
  }

  # Convert the data.frame input into a sparse matrix
  if (is.data.frame(ExpressionMatrix) == FALSE & is.matrix(ExpressionMatrix) == FALSE) {
    stop("Please supply an expression matrix in one of the following formats: data.frame or matrix")
  }

  # Create a new scGPS object.
  scGPSset <- SingleCellExperiment(assays = list(ExpressionMatrix), rowData = GeneMetadata, colData = CellMetadata)

  # All clear, return the object
  return(scGPSset)
  }
