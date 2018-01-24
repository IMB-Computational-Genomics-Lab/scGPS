
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

plotReduced_scGPS <- function(reduced_dat, color_fac =factor(Sample_id),
    dims = c(1,2),dimNames = c("Dim 1", "Dim 2"), palletes =NULL, legend_title ="Samples" ){
  library(cowplot)
  reduced_dat_toPlot <-as.data.frame(reduced_dat[,dims])
  sample_num <-length(unique(color_fac))
  colnames(reduced_dat_toPlot) <- dimNames
  reduced_dat_toPlot$color_fac <-color_fac
  p <- qplot(reduced_dat[,dims[1]], reduced_dat[,dims[2]], alpha=I(0.7), geom = "point",color=color_fac) +theme_bw()
  p <- p + ylab(dimNames[1]) + xlab(dimNames[2]) +
    scale_color_manual(name= legend_title, values=palletes[1:sample_num], limits=as.character(as.vector(unique(color_fac))))
  p<- p + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +  theme(legend.position="bottom") +
    theme(text=element_text(size=20))

  yaxis <- axis_canvas(p, axis="y", coord_flip = TRUE) + geom_density(
    data=reduced_dat_toPlot, aes(`Dim 2`, ..count.., fill=color_fac), size=.2,
    alpha=0.7) + coord_flip() + scale_fill_manual(name="Samples",
     values=palletes[1:sample_num], limits=as.character(as.vector(unique(color_fac))))

  xaxis <- axis_canvas(p,axis="x") + geom_density(data=reduced_dat_toPlot, aes(`Dim 1`, ..count..,  fill=color_fac), size=.4, alpha=0.7) +
    scale_fill_manual(name="Samples", values=palletes[1:sample_num], limits=as.character(as.vector(unique(color_fac))))

  p1_x <- insert_xaxis_grob(p, xaxis, grid::unit(.2, "null"), position = "top")
  p1_x_y <- insert_yaxis_grob(p1_x, yaxis, grid::unit(.2, "null"), position = "right")
  p2 <- ggdraw(p1_x_y)
  return(p2)
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

#' find DE genes
#'
#' @description  Find DE genes from comparing one clust vs remaining
#' @param expression_matrix is  a normalised expression matrix.
#' @param cluster corresponding cluster information in the expression_matrix
#' by going through CORE or from other method).
#' @return a \code{list} containing DE analysis results
#' @export
#' @author Quan Nguyen, 2017-11-25
#'


findMarkers_scGPS <- function(expression_matrix=NULL, cluster = NULL) {
  library(DESeq)

  DE_exprsMat <- round(expression_matrix+1)

  DE_results <- list()
  for (cl_id in unique(cluster)) {
    #arrange clusters and exprs matrix
    cl_index <- which(as.character(cluster) == as.character(cl_id))
    mainCl_idx <- which(as.character(cluster) != as.character(cl_id))
    condition_cluster = cluster
    condition_cluster[mainCl_idx] <- rep("Others", length(mainCl_idx))
    condition_cluster[cl_index] <- rep(as.character(cl_id), length(cl_index))
    diff_mat <- DE_exprsMat[,c(mainCl_idx, cl_index)]
    #start DE
    print(paste0("Start estimate dispersions for cluster ", as.character(cl_id) , "..."))
    cds = newCountDataSet(diff_mat, condition_cluster)
    cds = estimateSizeFactors( cds )
    cds = estimateDispersions( cds, method="per-condition" , fitType = "local" )
    print(paste0("Done estimate dispersions. Start nbinom test for cluster ",as.character(cl_id) , "..."))
    res1 = nbinomTest( cds, "Others", as.character(cl_id))
    print(paste0("Done nbinom test for cluster ",as.character(cl_id) , " ..."))
    #adjust folchange
    print(paste0("Adjust foldchange by subtracting basemean to 1..."))
    res1 <- mutate(res1,  AdjustedFC = (baseMeanB-1)/(baseMeanA-1))
    res1 <- mutate(res1,  AdjustedLogFC = log2((baseMeanB-1)/(baseMeanA-1)))
    #order
    res1_order <- arrange(res1, pval, desc(abs(AdjustedLogFC)))
    #write to list
    DE_results <- c(DE_results, list(res1))
    name_list =paste0("DE_Subpop", cl_id, "vsRemaining")
    names(DE_results)[length(DE_results)] <- name_list
  }

  return(DE_results)
}

#' annotate_scGPS functionally annotates the identified clusters
#'
#' @description often we need to label clusters with unique biological characters.
#' One of the common approach to annotate a cluster is to perform functional enrichment
#' analysis. The annotate_scGPS implements ReactomePA and clusterProfiler for this analysis
#' type in R. The function require installation of several databases as described below.
#' @param DEgeneList is a vector of gene symbols, convertable to ENTREZID



#Installation needed for reactome pathway analysis reactome in R----------------
#source("https://bioconductor.org/biocLite.R")
#biocLite("ReactomePA")
#Package Genome wide annotation for Human http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
#biocLite("org.Hs.eg.db")
#biocLite("clusterProfiler")
#install.packages("xlsx")
#Note: users may need to download and install clusterProfiler from source clusterProfiler_3.6.0.tgz"
#use: manual installing install.packages(path_to_file, repos = NULL, type="source")
#Done installation needed for reactome pathway analysis reactome in R-----------



annotate_scGPS <-function(DEgenelist, pvalueCutoff=0.05, gene_symbol=TRUE,
    output_filename = "PathwayEnrichment.xlsx", output_path = NULL ){
    library(ReactomePA)
    library(org.Hs.eg.db)
    library(org.Hs.eg.db)
    library(xlsx)
    #assumming the geneList is gene symbol (common for 10X data)
    if(gene_symbol==TRUE){
      convert_to_gene_ID = bitr(geneList, fromType="SYMBOL",
                              toType="ENTREZID", OrgDb="org.Hs.eg.db")
      print("Orignial gene number in geneList")
      print(length(geneList))
      print("Number of genes successfully converted")
      print(nrow(convert_to_gene_ID))
      } else {
        stop("The list must contain gene symbols")
        }

    Reactome_pathway_test <- enrichPathway(gene=convert_to_gene_ID$ENTREZID,
                                         pvalueCutoff=0.05, readable=T)

    #plot some results: note Reactome_pathway_test is a reactomePA object
    #write Reactome_pathway_test results, need to convert to data.frame

    dotplot(Reactome_pathway_test, showCategory=15)
    barplot(Reactome_pathway_test, showCategory=15)
    output_df <- as.data.frame(Reactome_pathway_test)
    write.xlsx(output_df, paste0(output_path, output_filename))
    return(Reactome_pathway_test)
}

