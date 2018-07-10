#' Main clustering CORE V2.0 with bagging analysis
#'
#' @description  CORE is an algorithm to generate reproduciable clustering,
#' CORE is first implemented in ascend R package. Here, CORE V2.0 uses bagging analysis
#' to find a stable clustering result and detect rare clusters mixed population
#' @param windows a numeric specifying the number of windows to test
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons.
#' @param PCA logical specifying if PCA is used before calculating distance matrix
#' @return a \code{list} with clustering results of all iterations, and a selected
#' optimal resolution
#' @examples
#' day5 <- sample2
#' cellnames<-colnames(day5$dat5_counts)
#' cluster <-day5$dat5_clusters
#' cellnames <- data.frame("cluster" = cluster, "cellBarcodes" = cellnames)
#' #day5$dat5_counts needs to be in a matrix format
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                          CellMetadata = day5$dat5_clusters)
#' test <- CORE_scGPS_bagging(mixedpop2, remove_outlier = c(1), PCA=FALSE, bagging_run = 20, subsample_proportion = .8)
#' plot_CORE(test$tree, list_clusters = test$clustering_param)
#'
#' @export
#' @author Quan Nguyen, 2018-05-11

CORE_scGPS_bagging <- function(mixedpop = NULL, bagging_run = 10, subsample_proportion = 0.8, windows = seq(0.025:1, by = 0.025), remove_outlier = c(0),
                       nRounds = 1, PCA=FALSE, nPCs=20, ngenes=1500) {
  cluster_all <- clustering_scGPS_bagging(object = mixedpop, windows = windows,
                                          bagging_run = bagging_run , subsample_proportion = subsample_proportion,
                                          remove_outlier = remove_outlier,
                                  nRounds = nRounds,PCA=PCA, nPCs=nPCs)

  optimal_stab <- list()
  for(i in 1:bagging_run){
  	stab_df <- FindStability(list_clusters = cluster_all$bootstrap_clusters[[i]][[2]],
  	                         cluster_ref = unname(unlist(cluster_all$bootstrap_clusters[[i]][[2]][[1]])))
    #stab_df <- FindStability(list_clusters = cluster_all$bootstrap_list[[i]], cluster_ref = cluster_all$bootstrap_list[[i]][[1]])


    #optimal_stab[[i]] <- FindOptimalStabilityBagging(list_clusters = cluster_all$bootstrap_list[[i]], stab_df)
    optimal_stab[[i]] <- FindOptimalStability(list_clusters = cluster_all$bootstrap_clusters[[i]][[2]], stab_df, bagging = TRUE)
  }

   OptimalCluster_bagging <-vector()
   #to find stable cluster
   for(i in 1:bagging_run){
     OptimalCluster_bagging[i] <- optimal_stab[[i]]$OptimalClust
   }
  # #to find rare cluster
   RareCluster_bagging <-vector()
   for(i in 1:bagging_run){
     RareCluster_bagging[i] <- optimal_stab[[i]]$HighestRes
   }


  #  #OptimalCluster_bagging_count<- max(table(OptimalCluster_bagging))[1]
   OptimalCluster_bagging_count<- which.max(tabulate(OptimalCluster_bagging))[1]


   #to find rare cluster
   #NEED TO ADD HERE TO HANDLE THE HIGH RESOLUTION CLUSTERS
   #HighestRes = list_clusters[[1]]

  return(list(Cluster = cluster_all$list_clusters, tree = cluster_all$tree, optimalClust = OptimalCluster_bagging,
              cellsRemoved = cluster_all$cellsRemoved, cellsForClustering = cluster_all$cellsForClustering,
              optimalMax = OptimalCluster_bagging_count, baggingClusters = cluster_all$bootstrap_clusters,
              highResCluster=RareCluster_bagging,
              clustering_param = cluster_all$clustering_param))
}

#' HC clustering for a number of resolutions
#'
#' @description  performs 40 clustering runs or more depending on windows
#' @param mixedpop1 is a \linkS4class{SingleCellExperiment} object from the
#' train mixed population
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons
#' @param ngenes number of top variable genes to be used
#' @param nPCs number of principal components from PCA dimensional reduction to be used
#' @param nRounds number of iterations to remove a selected clusters
#' @return clustering results
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                          CellMetadata = day5$dat5_clusters)
#' test <-clustering_scGPS(mixedpop2, remove_outlier = c(1))

clustering_scGPS_bagging <- function(object = NULL, ngenes = 1500, bagging_run = 10, subsample_proportion = 0.1, windows = seq(0.025:1,
                            by = 0.025), remove_outlier = c(0), nRounds = 1, PCA=FALSE, nPCs=20) {


  # function for the highest resolution clustering (i.e. no window applied)
  firstRoundClustering <- function(object = NULL) {
    exprs_mat <- assay(object)
    # take the top variable genes
    print("Identifying top variable genes")
    exprs_mat_topVar <- topvar_scGPS(exprs_mat, ngenes = ngenes)
    # tranpose so that cells are in rows
    exprs_mat_t <- t(exprs_mat_topVar)
    #-------------------------------------Work in progress--------#
    if(PCA==TRUE){
      # perform PCA dimensionality reduction
      print("Performing PCA analysis (Note: the variance for each cell needs to be >0)")
      exprs_mat_topVar_PCA <-prcomp(t(exprs_mat_topVar))
      exprs_mat_t <- as.data.frame(exprs_mat_topVar_PCA$x[,1:nPCs])

    }


    # calculate distance matrix for the rows
    print("Calculating distance matrix")
    dist_mat <- rcpp_parallel_distance(as.matrix(exprs_mat_t))
    #-------------------------------------Work in progress--------#

    print("Performing hierarchical clustering")
    original.tree <- fastcluster::hclust(as.dist(dist_mat), method = "ward.D2")
    # the original clusters to be used as the reference
    print("Finding clustering information")
    original.clusters <- unname(cutreeDynamic(original.tree, distM = as.matrix(dist_mat),
                                              verbose = 0))
    original.tree$labels <- original.clusters
    return(list(tree = original.tree, cluster_ref = original.clusters,
                  dist_mat = dist_mat, exprs_mat_t = exprs_mat_t))
  }

  removeOutlierCluster <- function(object = object, remove_outlier = remove_outlier,
                                   nRounds = nRounds) {


    # Initial Message to the user
    if (nRounds == 1) {
      print(paste0("Performing ", nRounds, " round of filtering"))
    } else {
      print(paste0("Performing ", nRounds, " rounds of filtering"))
    }

    # loop for the number of filtering rounds
    i = 1
    objectTemp <- object
    cells_to_remove <- c()

    while (i <= nRounds) {
            filter_out <- firstRoundClustering(objectTemp)
            cluster_toRemove <- which(filter_out$cluster_ref %in% remove_outlier)
            if (length(cluster_toRemove) > 0) {
                print(paste0("Found ", length(cluster_toRemove), " cells as outliers at round ",  i, " ..."))
                cells_to_remove <- c(cells_to_remove, cluster_toRemove)
                objectTemp <- object[, -cells_to_remove]
                i <- i + 1
            } else {
                print(paste0("No more outliers detected in filtering round ", i))
                i <- nRounds + 1
            }
        }

        filter_out <- firstRoundClustering(objectTemp)
        cluster_toRemove <- which(filter_out$cluster_ref %in% remove_outlier)
        if (length(cluster_toRemove) > 0) {
                print(paste0("Found ", length(cluster_toRemove), " cells as outliers at round ",  i, " ..."))
                print(paste0("Select ", i , " removal rounds if you want to remove these cells"))
        }

        exprs_mat_t <- filter_out$exprs_mat_t
        if(length(cells_to_remove) > 0){
          output <- list(firstRound_out = filter_out, cellsRemoved = colnames(object[,cells_to_remove]),
                       cellsForClustering = colnames(object[, -cells_to_remove]), exprs_mat_t = exprs_mat_t)
        }else{
          output <- list(firstRound_out = filter_out, cellsRemoved = c("No outliers found"),
                              cellsForClustering = "All cells are kept for clustering", exprs_mat_t = exprs_mat_t)
        }
        print(paste0(dim(exprs_mat_t)[1], " cells left after filtering"))
        return(output)
  }


  firstRoundPostRemoval <- removeOutlierCluster(object = object, remove_outlier = remove_outlier,
                                                nRounds = nRounds)
  firstRound_out <- firstRoundPostRemoval$firstRound_out

  exprs_mat_t = firstRoundPostRemoval$exprs_mat_t

  # return variables for the next step
  original.tree <- firstRound_out$tree
  original.clusters <- firstRound_out$cluster_ref
  dist_mat <- firstRound_out$dist_mat

  clustering_windows <- function(tree=original.tree, dist_mat=dist_mat){

    clustering_param <- list()

    for (i in 1:length(windows)) {

      namelist = paste0("window", windows[i])
      toadd <- as.vector(cutreeDynamic(tree, distM = as.matrix(dist_mat),
                                       minSplitHeight = windows[i], verbose = 0))

      #print(paste0("writing clustering result for run ", i))
      clustering_param[[i]] <- list(toadd)
      names(clustering_param[[i]]) <- namelist
    }

    return(clustering_param)
  }

  clustering_param_all <- clustering_windows(tree=original.tree, dist_mat=dist_mat)
  print(paste0("Running ",bagging_run, " bagging runs, with ", subsample_proportion, " subsampling..."))
  bootstrap_list <- list()
  for(i in 1:bagging_run){

    #dist_mat <- rcpp_parallel_distance(as.matrix(exprs_mat_t))
    #dist_mat_bootstrap_idx <- sample(1:ncol(dist_mat), subsample_proportion * ncol(dist_mat), replace=TRUE)

   # dist_mat_bootstrap_column_idx <- sample(1:ncol(exprs_mat_t),subsample_proportion * ncol(exprs_mat_t), replace=FALSE)
    dist_mat_bootstrap_row_idx <- sample(1:nrow(exprs_mat_t),subsample_proportion * nrow(exprs_mat_t), replace=TRUE)

    #dist_mat_bootstrap <- dist_mat[dist_mat_bootstrap_idx, dist_mat_bootstrap_idx]

    #dist_mat_bootstrap <- rcpp_parallel_distance(exprs_mat_t[dist_mat_bootstrap_row_idx,])
    dist_mat_bootstrap <- dist_mat[dist_mat_bootstrap_row_idx,dist_mat_bootstrap_row_idx]
    iter_tree <- fastcluster::hclust(as.dist(dist_mat_bootstrap), method = "ward.D2")
    iter_temp <- clustering_windows(tree=iter_tree, dist_mat=dist_mat_bootstrap)
    iter_write <-list(cellnames,iter_temp)
    bootstrap_list[[i]] <- iter_write
  }

  print("Done clustering, moving to stability calculation...")

  return(list(tree = original.tree, cluster_ref = original.clusters,
             bootstrap_clusters = bootstrap_list,
             clustering_param = clustering_param_all))

}

#'Find the optimal cluster
#'
#' @description from calculated stability based on Rand indexes for consecutive
#' clustering run, find the resolution (window), where the stability is the highest
#' @param run_RandIdx is a \code{data frame} object from iterative clustering runs
#' @param list_clusters is a \code{list} object containing 40 clustering results
#' @return a \code{list} with optimal stability, cluster count and summary stats
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                         CellMetadata = day5$dat5_clusters)
#' cluster_all <-clustering_scGPS(object=mixedpop2)
#' stab_df <- FindStability(list_clusters=cluster_all$list_clusters, cluster_ref = cluster_all$cluster_ref)
#' optimal_stab <- FindOptimalStability(list_clusters = cluster_all$list_clusters, stab_df)


FindOptimalStabilityBagging <- function(list_clusters, run_RandIdx) {
  library(reshape2)
  print("Start finding optimal clustering...")
  # get the number of cluster
  t <- lapply(list_clusters, function(x) {
    length(unique(unlist(x)))
  })
  run_RandIdx$cluster_count <- as.vector(as.numeric(t))

  #-----------------------------------------------------------------------------
  # Diagnostic plot for comparing clustering results
  #-----------------------------------------------------------------------------
  run_RandIdx$stability_count <- run_RandIdx$stability_count/40

  KeyStats <- as.data.frame(cbind(as.numeric(run_RandIdx$order) * 0.025, run_RandIdx$stability_count,
                                  run_RandIdx$cluster_index_ref, run_RandIdx$cluster_index_consec))

  colnames(KeyStats) <- c("Height", "Stability", "RandIndex", "ConsecutiveRI")

  KeyStats$Height <- as.character(KeyStats$Height)
  day_melt <- melt(KeyStats, id = "Height")
  day_melt$Height <- as.numeric(day_melt$Height)
  p <- ggplot(day_melt)
  p <- p + geom_line(aes(x = Height, y = value, colour = variable), size = 2) +
    theme_bw() + theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24)) + theme(legend.title = element_blank()) +
    xlab("Parameter from 0.025 to 1") + ylab("Scores") + theme(panel.border = element_rect(colour = "black",
         fill = NA, size = 1.5)) + guides(colour = FALSE)

  #-----------------------------------------------------------------------------
  # End diagnostic plot for comparing clustering results
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # Find the optimal parameter for clustering
  #-----------------------------------------------------------------------------

  KeyStats$Cluster_count <- run_RandIdx$cluster_count

  ###INITIALISING THE HEIGHT###
  stabilities <- KeyStats$Stability
  clusters <- KeyStats$Cluster_count
  max_cluster <- clusters[1]
  quality <-vector()

  for (i in 1:40){
    quality[i] <- stabilities[i]
  }
  optimal_param <- which(quality == max(quality))[1]
  optimal_cluster <- KeyStats$Cluster_count[optimal_param]

  # Final result
  return(list(HighestRes = max_cluster, Qualities = quality,
              OptimalClust = optimal_cluster,  KeyStats = KeyStats))
  #-----------------------------------------------------------------------------
  # Done finding the optimal parameter for clustering
  #-----------------------------------------------------------------------------
}
