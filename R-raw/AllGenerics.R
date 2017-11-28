#' @name QC_metrics
#' @aliases QC_metrics
#' @export
#' @docType methods
#' @return a dataframe of quality control matrics
#' @rdname QC_metrics

setGeneric("QC_metrics", function(object) standardGeneric("QC_metrics"))

#' @name QC_metrics
#' @aliases QC_metrics<-
#' @export
#' @docType methods
#' @rdname QC_metrics

setGeneric("QC_metrics<-", function(object, value) standardGeneric("QC_metrics<-"))


#' SubsetBatch
#'
#' Subset a specific batch from a \linkS4class{EMSet} object. This data is already normalised, but if you wish to recluster the data, you will need to use the RunPCA function again.
#'
#' @param object A \linkS4class{EMSet} object
#' @param batches Name or number of the batch(es) you would like to subset.
#'
#' @return An \linkS4class{EMSet} containing only this batch.
#' @include ascend_objects.R
#' @export
setGeneric(name = "SubsetBatch", def = function(object, batches) {
  standardGeneric("SubsetBatch")
})

setMethod("SubsetBatch", signature("EMSet"), function(object, batches = c()) {
  # Check if selected batches are present in the batches column
  if (!any(batches %in% object@CellInformation[, 2])) {
    stop("None of your selected batches are present in the dataset.")
  }

  # Create a new object to output, ensures original object does not get overwritten.
  subset.obj <- object

  # Retrieve data from relevant slots
  expression.matrix <- GetExpressionMatrix(subset.obj, format = "data.frame")
  cell.info <- subset.obj@CellInformation

  # Get barcodes for selected clusters
  subset.cell.info <- cell.info[which(cell.info[, 2] %in% batches), ]

  # Subset the expression matrix.
  subset.matrix <- expression.matrix[, unlist(subset.cell.info[, 1])]
  subset.obj <- ReplaceExpressionMatrix(subset.obj, subset.matrix)
  subset.obj <- SyncSlots(subset.obj)

  # Clean up the object
  subset.obj@PCA <- list()
  subset.obj@Clusters <- list()
  subset.obj@Log <- c(subset.obj@Log, list(SubsetByBatches = TRUE, SubsettedBatches = batches))
  return(subset.obj)
})

