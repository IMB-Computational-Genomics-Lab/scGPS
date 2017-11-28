
#' Get or set \code{gene_id_type} from a SingleCellExperiment object
#' @rdname gene_id_type
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return gene id type string
#' @author Luyi Tian
#'
#' @export
#'
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' sce = SingleCellExperiment(assays = list(counts =as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#'
#' gene_id_type(sce)
#'
gene_id_type.sce <- function(object) {
  return(object@metadata$Biomart$gene_id_type)
}


#' @rdname gene_id_type
#' @aliases gene_id_type
#' @export
setMethod("gene_id_type", signature(object = "SingleCellExperiment"),
          gene_id_type.sce)


#' @aliases gene_id_type
#' @rdname gene_id_type
#' @export
setReplaceMethod("gene_id_type",signature="SingleCellExperiment",
                 function(object, value) {
                   if(is.null(value)){
                     object@metadata$Biomart$gene_id_type = NA
                   }else if(value == "NA"){
                     object@metadata$Biomart$gene_id_type = NA
                   }else{
                     object@metadata$Biomart$gene_id_type = value
                   }
                   return(object)
                 })


#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the t-SNE Default is \code{500}, but any \code{ntop} argument is
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values),
#' or any other named element of the \code{assayData} slot of the \code{SingleCellExperiment}
#' object that can be accessed with the \code{assay} function.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the t-SNE calculation. If character, entries must all be
#' in \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param use_dimred character(1), use named reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object instead of recomputing (e.g. "PCA").
#'  Default is \code{NULL}, no reduced dimension values are provided to \code{Rtsne}.
#' @param n_dimred integer(1), number of components of the reduced dimension slot
#' to use. Default is \code{NULL}, in which case (if \code{use_dimred} is not \code{NULL})
#' all components of the reduced dimension slot are used.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param rand_seed (optional) numeric scalar that can be passed to
#' \code{set.seed} to make plots reproducible.
#' @param perplexity numeric scalar value defining the "perplexity parameter"
#' for the t-SNE plot. Passed to \code{\link[Rtsne]{Rtsne}} - see documentation
#' for that package for more details.
#'
#' @rdname plotTSNE
#' @export
runTSNE <- function(object, ntop = 500, ncomponents = 2, exprs_values = "logcounts",
                    feature_set = NULL, use_dimred = NULL, n_dimred = NULL, scale_features = TRUE,
                    rand_seed = NULL, perplexity = floor(ncol(object) / 5), ...) {

  if (!is.null(use_dimred)) {
    ## Use existing dimensionality reduction results (turning off PCA)
    dr <- reducedDim(object, use_dimred)
    if (!is.null(n_dimred)) {
      dr <- dr[,seq_len(n_dimred),drop = FALSE]
    }
    vals <- dr
    do_pca <- FALSE
    pca_dims <- ncol(vals)

  } else {
    ## Define an expression matrix depending on which values we're
    ## using
    exprs_mat <- assay(object, i = exprs_values)

    ## Define features to use: either ntop, or if a set of features is
    ## defined, then those
    if ( is.null(feature_set) ) {
      rv <- .general_rowVars(exprs_mat)
      ntop <- min(ntop, length(rv))
      feature_set <- order(rv, decreasing = TRUE)[seq_len(ntop)]
    }

    ## Drop any features with zero variance
    vals <- exprs_mat[feature_set,,drop = FALSE]
    keep_feature <- .general_rowVars(vals) > 0.001
    keep_feature[is.na(keep_feature)] <- FALSE
    vals <- vals[keep_feature,,drop = FALSE]

    ## Standardise expression if stand_exprs(object) is null
    vals <- t(vals)
    if (scale_features) {
      vals <- scale(vals, scale = TRUE)
    }
    do_pca <- TRUE
    pca_dims <- max(50, ncol(object))
  }

  # Actually running the Rtsne step.
  if ( !is.null(rand_seed) )
    set.seed(rand_seed)
  tsne_out <- Rtsne::Rtsne(vals, initial_dims = pca_dims, pca = do_pca,
                           perplexity = perplexity, dims = ncomponents,...)
  reducedDim(object, "TSNE") <- tsne_out$Y
  return(object)
}


#' Plot t-SNE for an SingleCellExperiment object
#'
#' Produce a t-distributed stochastic neighbour embedding (t-SNE) plot of two
#' components for an \code{SingleCellExperiment} dataset.
#'
#' @param object an \code{SingleCellExperiment} object
#' @param ncomponents numeric scalar indicating the number of t-SNE
#' components to plot, starting from the first t-SNE component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of component 1 vs component
#' 2 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#' top components is produced. NB: computing more than two components for t-SNE
#' can become very time consuming.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively,
#' a data frame with one column containing values to map to colours for all cells.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' Alternatively, a data frame with one column containing values to map to shapes.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' Alternatively, a data frame with one column containing values to map to sizes.
#' @param return_SCE logical, should the function return an \code{SingleCellExperiment}
#' object with principal component values for cells in the
#' \code{reducedDims} slot. Default is \code{FALSE}, in which case a
#' \code{ggplot} object is returned.
#' @param rerun logical, should PCA be recomputed even if \code{object} contains a
#' "PCA" element in the \code{reducedDims} slot?
#' @param draw_plot logical, should the plot be drawn on the current graphics
#' device? Only used if \code{return_SCE} is \code{TRUE}, otherwise the plot
#' is always produced.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param ... further arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @details The function \code{\link[Rtsne]{Rtsne}} is used internally to
#' compute the t-SNE. Note that the algorithm is not deterministic, so different
#' runs of the function will produce differing plots (see \code{\link{set.seed}}
#' to set a random seed for replicable results). The value of the
#' \code{perplexity} parameter can have a large effect on the resulting plot, so
#' it can often be worthwhile to try multiple values to find the most appealing
#' visualisation.
#'
#' @return If \code{return_SCE} is \code{TRUE}, then the function returns a
#' \code{SingleCellExperiment} object, otherwise it returns a \code{ggplot} object.
#' @name plotTSNE
#' @rdname plotTSNE
#' @aliases plotTSNE plotTSNE,SingleCellExperiment-method
#'
#' @export
#' @seealso
#' \code{\link[Rtsne]{Rtsne}}
#' @references
#' L.J.P. van der Maaten. Barnes-Hut-SNE. In Proceedings of the International
#' Conference on Learning Representations, 2013.
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' ## Examples plotting t-SNE
#' plotTSNE(example_sce, perplexity = 10)
#' plotTSNE(example_sce, colour_by = "Cell_Cycle", perplexity = 10)
#' plotTSNE(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
#' size_by = "Mutation_Status", perplexity = 10)
#' plotTSNE(example_sce, shape_by = "Treatment", size_by = "Mutation_Status",
#' perplexity = 5)
#' plotTSNE(example_sce, feature_set = 1:100, colour_by = "Treatment",
#' shape_by = "Mutation_Status", perplexity = 5)
#'
#' plotTSNE(example_sce, shape_by = "Treatment", return_SCE = TRUE,
#' perplexity = 10)
#'
#'
plotTSNE <- function(object, colour_by = NULL, shape_by = NULL, size_by = NULL,
                     return_SCE = FALSE, draw_plot = TRUE,
                     theme_size = 10, legend = "auto",
                     rerun = FALSE, ncomponents = 2, ...) {

  if ( !("TSNE" %in% names(reducedDims(object))) || rerun) {
    object <- runTSNE(object, ncomponents = ncomponents, ...)
  }

  plot_out <- plotReducedDim(object, ncomponents = ncomponents, use_dimred = "TSNE",
                             colour_by = colour_by, shape_by = shape_by, size_by = size_by,
                             theme_size = theme_size, legend = legend)

  if (return_SCE) {
    if ( draw_plot )
      print(plot_out)
    return(object)
  } else {
    return(plot_out)
  }
}
