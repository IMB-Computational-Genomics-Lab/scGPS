#' Main clustering SCORE (CORE V2.0) Stable Clustering
#' at Optimal REsolution with bagging and bootstrapping
#'
#' @description  CORE is an algorithm to generate reproduciable clustering,
#' CORE is first implemented in ascend R package. Here, CORE V2.0 uses bagging 
#' analysis to find a stable clustering result and detect rare clusters mixed
#' population.
#' @param mixedpop is a \linkS4class{SingleCellExperiment} object from the train
#' mixed population.
#' @param bagging_run an integer specifying the number of bagging runs to be 
#' computed.
#' @param subsample_proportion a numeric specifying the proportion 
#' of the tree to be chosen in subsampling.
#' @param windows a numeric vector specifying the ranges of each window.
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons.
#' @param nRounds an integer specifying the number rounds to attempt to remove 
#' outliers.
#' @param PCA logical specifying if PCA is used before calculating distance 
#' matrix.
#' @param nPCs an integer specifying the number of principal components to use.
#' @param ngenes number of genes used for clustering calculations.
#' @return a \code{list} with clustering results of all iterations, and a 
#' selected
#' optimal resolution
#' @examples
#' day5 <- day_5_cardio_cell_sample
#' cellnames<-colnames(day5$dat5_counts)
#' cluster <-day5$dat5_clusters
#' cellnames <- data.frame('cluster' = cluster, 'cellBarcodes' = cellnames)
#' #day5$dat5_counts needs to be in a matrix format
#' mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' test <- CORE_bagging(mixedpop2, remove_outlier = c(0), PCA=FALSE,
#'     bagging_run = 2, subsample_proportion = .7)

#' @export
#' @author Quan Nguyen, 2018-05-11

CORE_bagging <- function(mixedpop = NULL, bagging_run = 20, 
    subsample_proportion = 0.8, windows = seq(from = 0.025, to = 1, by = 0.025),
    remove_outlier = c(0), nRounds = 1, PCA = FALSE, nPCs = 20, ngenes = 1500) {
    
    # perform the clustering runs
    cluster_all <- clustering_bagging(object = mixedpop, 
        windows = windows, bagging_run = bagging_run, 
        subsample_proportion = subsample_proportion, 
        remove_outlier = remove_outlier, ngenes = ngenes, nRounds = nRounds, 
        PCA = PCA, nPCs = nPCs)
    
    # find the optimal stability for each of the bagging runs
    optimal_stab <- vector(mode = "list", length = bagging_run)
    for (i in seq_len(bagging_run)) {
        stab_df <- find_stability(list_clusters = 
            cluster_all$bootstrap_clusters[[i]][[1]], 
            cluster_ref = unname(unlist(
            cluster_all$bootstrap_clusters[[i]][[1]][[1]])))
        
        optimal_stab[[i]] <- find_optimal_stability(list_clusters = 
            cluster_all$bootstrap_clusters[[i]][[1]], 
            stab_df, bagging = TRUE, windows = windows)
    }
    
    # record the optimal and highest resolutions to find stable cluster
    OptimalCluster_bagging <- vector(mode = "double", length = bagging_run)
    for (i in seq_len(bagging_run)) {
        OptimalCluster_bagging[i] <- optimal_stab[[i]]$OptimalClust
    }
    
    # to find rare cluster
    RareCluster_bagging <- vector(mode = "double", length  = bagging_run)
    for (i in seq_len(bagging_run)) {
        RareCluster_bagging[i] <- optimal_stab[[i]]$HighestRes
    }
    
    # record the most frequently occurring result, favour higher resolution
    OptimalCluster_bagging_count <- which.max(tabulate(
        OptimalCluster_bagging))[1]
    
    NumberClusters <- vector(mode = "double", length  = length(windows))
    for (i in seq_len(length(windows))) {
        NumberClusters[i] <- max(unlist(cluster_all$clustering_param[[i]]))
    }
    
    #-----------REWRITE THIS CODE SNIPPET---------------------#
    # Check to see if the optimal is valid in the original tree
    if (!(OptimalCluster_bagging_count %in% NumberClusters)) {
        if (OptimalCluster_bagging_count > max(NumberClusters)) {
            OptimalCluster_bagging_count <- max(NumberClusters)
        } else {
            above <- OptimalCluster_bagging_count + 1
            below <- OptimalCluster_bagging_count - 1
            above_count <- length(which(NumberClusters == above))
            below_count <- length(which(NumberClusters == below))
            if (above_count >= below_count) {
                OptimalCluster_bagging_count <- above
            } else {
                OptimalCluster_bagging_count <- below
            }
        }
    }
    
    optimal_index <- which(NumberClusters == OptimalCluster_bagging_count)[1]
    
    
    return(list(Cluster = cluster_all$clustering_param, tree = cluster_all$tree,
        optimalClust = OptimalCluster_bagging, 
        cellsRemoved = cluster_all$cellsRemoved, 
        cellsForClustering = cluster_all$cellsForClustering, 
        optimalMax = OptimalCluster_bagging_count, 
        highResCluster = RareCluster_bagging, optimal_index = optimal_index))
}



#' HC clustering for a number of resolutions
#'
#' @description  subsamples cells for each bagging run and performs 40 
#' clustering runs or more depending on windows.
#' @param object is a \linkS4class{SingleCellExperiment} object from the train
#' mixed population.
#' @param bagging_run an integer specifying the number of bagging runs to be 
#' computed.
#' @param subsample_proportion a numeric specifying the proportion of the tree 
#' to be chosen in subsampling.
#' @param windows a numeric vector specifying the rages of each window.
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons.
#' @param nRounds a integer specifying the number rounds to attempt to remove 
#' outliers.
#' @param PCA logical specifying if PCA is used before calculating distance 
#' matrix.
#' @param nPCs an integer specifying the number of principal components to use.
#' @param ngenes number of genes used for clustering calculations.
#' @return a list of clustering results containing each bagging run
#' as well as the clustering of the original tree and the tree itself.
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day5 <- day_5_cardio_cell_sample
#' mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' test <-clustering_bagging(mixedpop2, remove_outlier = c(0),
#'     bagging_run = 2, subsample_proportion = .7)

clustering_bagging <- function(object = NULL, ngenes = 1500, 
    bagging_run = 20, subsample_proportion = 0.8, 
    windows = seq(from = 0.025, to = 1, by = 0.025), remove_outlier = c(0),
    nRounds = 1, PCA = FALSE, nPCs = 20) {
    
    
    # function for the highest resolution clustering (i.e. no window applied)
    first_round_clustering <- function(object = NULL) {
        exprs_mat <- assay(object)
        # take the top variable genes
        message("Identifying top variable genes")
        exprs_mat_topVar <- top_var(exprs_mat, ngenes = ngenes)
        # tranpose so that cells are in rows
        exprs_mat_t <- t(exprs_mat_topVar)
        #-------------------------------------Work in progress--------#
        if (PCA) {
            # perform PCA dimensionality reduction
            message(paste0("Performing PCA analysis ", 
                "(Note: the variance for each cell needs to be > 0)"))
            exprs_mat_topVar_PCA <- prcomp(t(exprs_mat_topVar))
            exprs_mat_t <- as.data.frame(exprs_mat_topVar_PCA$x[,seq_len(nPCs)])
            
        }
        
        
        # calculate distance matrix for the rows
        message("Calculating distance matrix")
        dist_mat <- rcpp_parallel_distance(as.matrix(exprs_mat_t))
        #-------------------------------------Work in progress--------#
        
        message("Performing hierarchical clustering")
        original.tree <- fastcluster::hclust(as.dist(dist_mat), 
            method = "ward.D2")
        # the original clusters to be used as the reference
        message("Finding clustering information")
        original.clusters <- unname(cutreeDynamic(original.tree, 
            distM = as.matrix(dist_mat), verbose = 0))
        original.tree$labels <- original.clusters
        return(list(tree = original.tree, cluster_ref = original.clusters,
            dist_mat = dist_mat, exprs_mat_t = exprs_mat_t))
    }
    
    # function to remove outlier clusters
    remove_outlier_cluster <- function(object = object, 
        remove_outlier = remove_outlier, nRounds = nRounds) {
        
        
        # Initial Message to the user
        if (nRounds == 1) {
            message(paste0("Performing ", nRounds, " round of filtering"))
        } else {
            message(paste0("Performing ", nRounds, " rounds of filtering"))
        }
        
        # loop for the number of filtering rounds
        i = 1
        objectTemp <- object
        cells_to_remove <- c()
        
        while (i <= nRounds) {
            filter_out <- first_round_clustering(objectTemp)
            cluster_toRemove <- which(filter_out$cluster_ref %in% 
                remove_outlier)
            if (length(cluster_toRemove) > 0) {
                message(paste0("Found ", length(cluster_toRemove), 
                    " cells as outliers at round ", i, " ..."))
                cells_to_remove <- c(cells_to_remove, cluster_toRemove)
                objectTemp <- object[, -cells_to_remove]
                i <- i + 1
            } else {
                message(paste0("No more outliers detected in filtering round ",
                    i))
                i <- nRounds + 1
            }
        }
        
        filter_out <- first_round_clustering(objectTemp)
        cluster_toRemove <- which(filter_out$cluster_ref %in% remove_outlier)
        if (length(cluster_toRemove) > 0) {
            message(paste0("Found ", length(cluster_toRemove), 
                " cells as outliers at round ", i, " ..."))
            message(paste0("Select ", i, 
                " removal rounds if you want to remove these cells"))
        }
        
        exprs_mat_t <- filter_out$exprs_mat_t
        if (length(cells_to_remove) > 0) {
            output <- list(firstRound_out = filter_out, 
                cellsRemoved = colnames(object[, cells_to_remove]), 
                cellsForClustering = colnames(object[, -cells_to_remove]), 
                exprs_mat_t = exprs_mat_t)
        } else {
            output <- list(firstRound_out = filter_out, 
                cellsRemoved = c("No outliers found"), 
                cellsForClustering = "All cells are kept for clustering", 
                exprs_mat_t = exprs_mat_t)
        }
        message(paste0(dim(exprs_mat_t)[1], " cells left after filtering"))
        return(output)
    }
    
    
    firstRoundPostRemoval <- remove_outlier_cluster(object = object, 
        remove_outlier = remove_outlier, nRounds = nRounds)
    firstRound_out <- firstRoundPostRemoval$firstRound_out
    
    exprs_mat_t = firstRoundPostRemoval$exprs_mat_t
    
    # return variables for the next step
    original.tree <- firstRound_out$tree
    original.clusters <- firstRound_out$cluster_ref
    dist_mat <- firstRound_out$dist_mat
    
    # function that performs clustering for each window
    clustering_windows <- function(tree = original.tree, dist_mat = dist_mat) {
        
        clustering_param <- vector(mode = "list", length = length(windows))
        
        for (i in seq_len(length(windows))) {
            
            namelist = paste0("window", windows[i])
            toadd <- as.vector(cutreeDynamic(tree, distM = as.matrix(dist_mat), 
                minSplitHeight = windows[i], verbose = 0))
            clustering_param[[i]] <- list(toadd)
            names(clustering_param[[i]]) <- namelist
        }
        
        return(clustering_param)
    }
    
    # save the whole tree
    clustering_param_all <- clustering_windows(tree = original.tree, 
        dist_mat = dist_mat)
    
    # cluster a subsample for each of the bagging runs
    message(paste0("Running ", bagging_run, " bagging runs, with ", 
        subsample_proportion, " subsampling..."))
    bootstrap_list <- vector(mode = "list", length = bagging_run)
    for (i in seq_len(bagging_run)) {
        # subsample the distance matrix for every bagging run
        dist_mat_bootstrap_row_idx <- sample(seq_len(nrow(exprs_mat_t)), 
            subsample_proportion * nrow(exprs_mat_t), replace = TRUE)
        dist_mat_bootstrap <- dist_mat[dist_mat_bootstrap_row_idx, 
            dist_mat_bootstrap_row_idx]
        
        # tempoarily store clustering results for each run
        iter_tree <- fastcluster::hclust(as.dist(dist_mat_bootstrap), 
            method = "ward.D2")
        iter_temp <- clustering_windows(tree = iter_tree, 
            dist_mat = dist_mat_bootstrap)
        iter_write <- list(iter_temp)  #NEED TO PARSE CELLNAMES TO HERE
        bootstrap_list[[i]] <- iter_write
    }
    
    message("Done clustering, moving to stability calculation...")
    
    return(list(tree = original.tree, cluster_ref = original.clusters, 
        bootstrap_clusters = bootstrap_list, 
        clustering_param = clustering_param_all, 
        cellsRemoved = firstRoundPostRemoval$cellsRemoved, 
        cellsForClustering = firstRoundPostRemoval$cellsForClustering))
    
}

