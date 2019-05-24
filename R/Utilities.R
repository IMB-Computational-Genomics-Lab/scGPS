
#' select top variable genes
#' @description subset a matrix by top variable genes
#' @param expression.matrix is a matrix with genes in rows and cells in columns
#' @return a subsetted expression matrix with the top n most variable genes
#' @param ngenes number of genes used for clustering calculations.
#' @export
#' @examples
#' day2 <- day_2_cardio_cell_sample
#' mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' SortedExprsMat <-top_var(expression.matrix=assay(mixedpop1))

top_var <- function(expression.matrix = NULL, ngenes = 1500) {
    calc_row_variance <- function(x) {
        row.variance <- rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
        return(row.variance)
    }
    gene.variance <- calc_row_variance(expression.matrix)
    names(gene.variance) <- rownames(expression.matrix)
    sorted.gene.variance <- gene.variance[order(gene.variance, 
        decreasing = TRUE)]
    top.genes <- sorted.gene.variance[1:ngenes]
    subset.matrix <- expression.matrix[names(top.genes), ]
    return(subset.matrix)
}

#' plot reduced data
#' @description plot PCA, tSNE, and CIDR reduced datasets
#' @param reduced_dat is a matrix with genes in rows and cells in columns
#' @param color_fac is a vector of colors corresponding to clusters to determine
#' colors of scattered plots 
#' @param palletes can be a customised color pallete that determine colors for 
#' density plots, if NULL it will use RColorBrewer 
#' colorRampPalette(RColorBrewer::brewer.pal(sample_num, 'Set1'))(sample_num)
#' @param dims an integer of the number of dimestions
#' @param dimNames a vector of the names of the dimensions
#' @param legend_title title of the plot's legend
#' @export
#' @return a matrix with the top 20 CIDR dimensions
#' @examples
#' day2 <- day_2_cardio_cell_sample
#' mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' #CIDR_dim <-CIDR(expression.matrix=assay(mixedpop1))
#' #p <- plot_reduced(CIDR_dim, color_fac = factor(colData(mixedpop1)[,1]),
#' #     palletes =1:length(unique(colData(mixedpop1)[,1])))
#' #plot(p)
#' tSNE_dim <-tSNE(expression.mat=assay(mixedpop1))
#' p2 <- plot_reduced(tSNE_dim, color_fac = factor(colData(mixedpop1)[,1]),
#'     palletes =1:length(unique(colData(mixedpop1)[,1])))
#' plot(p2)
#'

plot_reduced <- function(reduced_dat, color_fac = NULL, dims = c(1, 2), 
    dimNames = c("Dim1", "Dim2"), palletes = NULL, legend_title = "Cluster") {
    reduced_dat_toPlot <- as.data.frame(reduced_dat[, dims])
    sample_num <- length(unique(color_fac))
    if (is.null(palletes)) {
        palletes <- colorRampPalette(RColorBrewer::brewer.pal(sample_num, 
            "Set1"))(sample_num)
    }
    reduced_dat_toPlot <- as.data.frame(reduced_dat[, dims])
    sample_num <- length(unique(color_fac))
    colnames(reduced_dat_toPlot) <- dimNames
    reduced_dat_toPlot$color_fac <- color_fac
    p <- qplot(x = reduced_dat[, dims[1]], y = reduced_dat[, dims[2]], 
        alpha = I(0.7), geom = "point", color = color_fac) + theme_bw()
    p <- p + ylab(dimNames[2]) + xlab(dimNames[1]) + 
        scale_color_manual(name = legend_title, values = palletes[1:sample_num],
        limits = sort(as.character(as.vector(unique(color_fac)))))
    p <- p + theme(panel.border = element_rect(colour = "black", fill = NA, 
        size = 1.5)) + theme(legend.position = "bottom") + 
        theme(text = element_text(size = 20))
    
    yaxis <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE) + 
        geom_density(data = reduced_dat_toPlot, aes(reduced_dat_toPlot$Dim2, 
            ..count.., fill = color_fac), size = 0.2, alpha = 0.7) + 
        coord_flip() + scale_fill_manual(name = "Samples", 
            values = palletes[1:sample_num], 
            limits = sort(as.character(as.vector(unique(color_fac)))))
    
    xaxis <- cowplot::axis_canvas(p, axis = "x") + 
        geom_density(data = reduced_dat_toPlot, aes(reduced_dat_toPlot$Dim1,
            ..count.., fill = color_fac), size = 0.4, alpha = 0.7) + 
        scale_fill_manual(name = "Samples", values = palletes[1:sample_num], 
            limits = sort(as.character(as.vector(unique(color_fac)))))
    
    p1_x <- cowplot::insert_xaxis_grob(p, xaxis, grid::unit(0.2, "null"), 
        position = "top")
    p1_x_y <- cowplot::insert_yaxis_grob(p1_x, yaxis, grid::unit(0.2, "null"), 
        position = "right")
    p2 <- cowplot::ggdraw(p1_x_y)
    return(p2)
}


#' find marker genes
#'
#' @description  Find DE genes from comparing one clust vs remaining
#' @param expression_matrix is  a normalised expression matrix.
#' @param cluster corresponding cluster information in the expression_matrix
#' by running CORE clustering or using other methods.
#' @param selected_cluster a vector of unique cluster ids to calculate
#' @param fitType string specifying 'local' or 'parametric' for DEseq dispersion
#' estimation
#' @param dispersion_method one of the options c( 'pooled', 'pooled-CR', 
#' per-condition', 'blind' )
#' @param sharing_Mode one of the options c("maximum", "fit-only", 
#' "gene-est-only")
#' @return a \code{list} containing sorted DESeq analysis results
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day2 <- day_2_cardio_cell_sample
#' mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' # depending on the data, the DESeq::estimateDispersions function requires
#' # suitable fitType
#'# and dispersion_method options
#'DEgenes <- find_markers(expression_matrix=assay(mixedpop1),
#'                         cluster = colData(mixedpop1)[,1],
#'                         selected_cluster=c(1), #can also run for more
#'                         #than one clusters, e.g.selected_cluster = c(1,2)
#'                         fitType = "parametric", 
#'                         dispersion_method = "blind",
#'                         sharing_Mode="fit-only"
#'                         )
#'names(DEgenes)



find_markers <- function(expression_matrix = NULL, cluster = NULL, 
    selected_cluster = NULL, fitType = "local", 
    dispersion_method = "per-condition",
    sharing_Mode = "maximum") {
    DE_exprsMat <- round(expression_matrix + 1)
    
    DE_results <- vector(mode = "list", 
        length = length(unique(selected_cluster)))
    for (cl_id in unique(selected_cluster)) {
        # arrange clusters and exprs matrix
        cl_index <- which(as.character(cluster) == as.character(cl_id))
        mainCl_idx <- which(as.character(cluster) != as.character(cl_id))
        diff_mat <- DE_exprsMat[, c(mainCl_idx, cl_index)]
        # start DE
        
        condition_cluster = as.vector(cluster)
        condition_cluster[1:length(mainCl_idx)] <- rep("Others", 
            length(mainCl_idx))
        condition_cluster[(length(mainCl_idx) + 1):ncol(diff_mat)] <- 
            rep(as.character(cl_id), length(cl_index))
        
        message(paste0("Start estimate dispersions for cluster ", 
            as.character(cl_id), "..."))
        # Note: the local fit option requires the library
        cds = DESeq::newCountDataSet(diff_mat, condition_cluster)
        cds = DESeq::estimateSizeFactors(cds)
        #library(locfit)
        cds = DESeq::estimateDispersions(cds, method = dispersion_method, 
            fitType = fitType, sharingMode = sharing_Mode)
        message(paste0("Done estimate dispersions. ",
            "Start nbinom test for cluster ", as.character(cl_id), "..."))
        res1 = DESeq::nbinomTest(cds, "Others", as.character(cl_id))
        message(paste0("Done nbinom test for cluster ", as.character(cl_id),
            " ..."))
        # adjust folchange
        message(paste0("Adjust foldchange by subtracting basemean to 1..."))
        res1 <- mutate(res1, AdjustedFC = (res1$baseMeanB - 1)/
            (res1$baseMeanA - 1))
        res1 <- mutate(res1, AdjustedLogFC = log2((res1$baseMeanB - 1)/
            (res1$baseMeanA - 1)))
        # order
        res1_order <- arrange(res1, res1$pval, desc(abs(res1$AdjustedLogFC)))
        # write to list
        DE_results <- c(DE_results, list(res1_order))
        name_list = paste0("DE_Subpop", cl_id, "vsRemaining")
        names(DE_results)[length(DE_results)] <- name_list
    }
    
    return(DE_results)
}

#' annotate functionally annotates the identified clusters
#'
#' @description often we need to label clusters with unique biological 
#' characters. One of the common approach to annotate a cluster is to perform 
#' functional enrichment analysis. The annotate implements ReactomePA and
#' clusterProfiler for this analysis type in R. The function require 
#' installation of several databases as described below.
#' @param DEgeneList is a vector of gene symbols, convertable to ENTREZID
#' @param pvalueCutoff is a numeric of the cutoff p value
#' @param gene_symbol logical of whether the geneList is a gene symbol
#' @param species is the selection of 'human' or 'mouse', default to 'human' 
#' genes
#' @export
#' @return write enrichment test output to a file and an enrichment test object 
#' for plotting
#' @examples
#' genes <-training_gene_sample
#' genes <-genes$Merged_unique[1:50]
#' enrichment_test <- annotate(genes, pvalueCutoff=0.05, 
#'     gene_symbol=TRUE, species = 'human')
#' clusterProfiler::dotplot(enrichment_test, showCategory=15)
#'


# Installation needed for reactome pathway analysis reactome in R---------------
# source('https://bioconductor.org/biocLite.R') biocLite('ReactomePA') Package
# Genome wide annotation for Human
#http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
# biocLite('org.Hs.eg.db'); biocLite('org.Mm.eg.db');
# biocLite('clusterProfiler'); install.packages('xlsx') Note: users may need to
# download and install clusterProfiler from source clusterProfiler_3.6.0.tgz'
# use: manual installing install.packages(path_to_file, repos = NULL,
# type='source') Done installation needed for reactome pathway analysis reactome
# in R----------------------------

annotate <- function(DEgeneList, pvalueCutoff = 0.05, gene_symbol = TRUE,
    species = "human") {
    # assumming the geneList is gene symbol (common for 10X data)
    if (species == "human") {
        if (gene_symbol == TRUE) {
            convert_to_gene_ID = clusterProfiler::bitr(DEgeneList, 
                fromType = "SYMBOL", toType = "ENTREZID", 
                OrgDb = "org.Hs.eg.db")
            message("Original gene number in geneList")
            message(length(DEgeneList))
            message("Number of genes successfully converted")
            message(nrow(convert_to_gene_ID))
        } else {
            stop("The list must contain human gene symbols")
        }
    } else if (species == "mouse") {
        if (gene_symbol) {
            convert_to_gene_ID = clusterProfiler::bitr(DEgeneList, 
                fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = "org.Mm.eg.db")
            message("Original gene number in geneList")
            message(length(DEgeneList))
            message("Number of genes successfully converted")
            message(nrow(convert_to_gene_ID))
        } else {
            stop("The list must contain mouse gene symbols")
        }
    }
    
    Reactome_pathway_test <- ReactomePA::enrichPathway(
        gene = convert_to_gene_ID$ENTREZID, pvalueCutoff = 0.05,
        readable = TRUE)
    
    # plot some results: note Reactome_pathway_test is a reactomePA object write
    # Reactome_pathway_test results, need to convert to data.frame
    output_df <- as.data.frame(Reactome_pathway_test)
    # xlsx::write.xlsx(output_df, paste0(output_path, output_filename))
    return(Reactome_pathway_test)
    # note can conveniently plot the outputs by running the followings
    # dotplot(Reactome_pathway_test, showCategory=15) 
    # barplot(Reactome_pathway_test, showCategory=15)
}


#' add_import
#'
#' @description temp function to import packages to namespace using devtools/roxygen
#' @name add_import
#' @useDynLib scGPS
#' @importFrom Rcpp evalCpp
#' @import glmnet
#' @import caret
#' @import dplyr
#' @import dynamicTreeCut
#' @import ggplot2
#' @import RcppArmadillo
#' @import RcppParallel
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import DESeq
#' @import locfit
#' @importFrom graphics barplot lines rect strheight strwidth text
#' @importFrom stats as.dist coef na.omit prcomp predict sd as.dendrogram
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline layout par plot
#' @importFrom utils head capture.output

NULL
