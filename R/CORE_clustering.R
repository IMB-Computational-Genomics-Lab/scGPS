#' Main clustering CORE V2.0 updated
#'
#' @description  CORE is an algorithm to generate reproduciable clustering,
#' CORE is first implemented in ascend R package. Here, CORE V2.0 introduces 
#' several new functionalities, including three key features:
#' fast (and more memory efficient) implementation with C++ and paralellisation
#' options allowing clustering of hundreds of thousands of cells 
#' (ongoing development), outlier revomal important if singletons exist (done),
#' a number of dimensionality reduction methods including the imputation
#' implementation (CIDR) for confirming clustering results (done), and an option
#' to select the number of optimisation tree height windows for increasing
#' resolution
#' @param mixedpop is a \linkS4class{SingleCellExperiment} object from the train
#' mixed population
#' @param windows a numeric specifying the number of windows to test
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons.
#' @param PCA logical specifying if PCA is used before calculating distance
#' matrix
#' @param nRounds an integer specifying the number rounds to attempt to remove 
#' outliers.
#' @param nPCs an integer specifying the number of principal components to use.
#' @param ngenes number of genes used for clustering calculations.
#' @return a \code{list} with clustering results of all iterations, and a 
#' selected optimal resolution
#' @examples
#' day5 <- sample2
#' #day5$dat5_counts needs to be in a matrix format
#' cellnames <- colnames(day5$dat5_counts)
#' cluster <-day5$dat5_clusters
#' cellnames <-data.frame('Cluster'=cluster, 'cellBarcodes' = cellnames)
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' test <- CORE_scGPS(mixedpop2, remove_outlier = c(1), PCA=FALSE, nPCs=20,
#'     ngenes=1500)
#' @export
#' @author Quan Nguyen, 2017-11-25

CORE_scGPS <- function(mixedpop = NULL, windows = seq(0.025:1, by = 0.025), 
    remove_outlier = c(0), nRounds = 1, PCA = FALSE, nPCs = 20, ngenes = 1500) {
    cluster_all <- clustering_scGPS(object = mixedpop, windows = windows, 
        remove_outlier = remove_outlier, 
        nRounds = nRounds, PCA = PCA, nPCs = nPCs)
    
    stab_df <- FindStability(list_clusters = cluster_all$list_clusters,
        cluster_ref = cluster_all$cluster_ref)
    optimal_stab <- FindOptimalStability(
        list_clusters = cluster_all$list_clusters, stab_df, windows = windows)
    
    return(list(Cluster = cluster_all$list_clusters, tree = cluster_all$tree,
        optimalClust = optimal_stab, cellsRemoved = cluster_all$cellsRemoved, 
        cellsForClustering = cluster_all$cellsForClustering))
}

#' Subclustering (optional) after running CORE 'test'
#'
#' @description  CORE_Subcluster_scGPS allows re-cluster the CORE clustering 
#' result
#' @param mixedpop is a \linkS4class{SingleCellExperiment} object from the train
#' mixed population
#' @param windows a numeric specifying the number of windows to test
#' @param  select_cell_index a vector containing indexes for cells in selected 
#' clusters to be reclustered
#' @param ngenes number of genes used for clustering calculations.
#' @return a \code{list} with clustering results of all iterations, and a 
#' selected optimal resolution
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts,
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' test <- CORE_scGPS(mixedpop2,remove_outlier= c(0))
#' @export
#' @author Quan Nguyen, 2017-11-25

CORE_Subcluster_scGPS <- function(mixedpop = NULL, windows = seq(0.025:1, 
    by = 0.025), select_cell_index = NULL, ngenes = 1500) {
    
    cluster_all <- SubClustering_scGPS(object = mixedpop, windows = windows, 
        select_cell_index = select_cell_index, ngenes = ngenes)
    
    stab_df <- FindStability(list_clusters = cluster_all$list_clusters, 
        cluster_ref = cluster_all$cluster_ref)
    
    optimal_stab <- FindOptimalStability(
        list_clusters = cluster_all$list_clusters, stab_df)
    
    return(list(Cluster = cluster_all$list_clusters, tree = cluster_all$tree, 
        optimalClust = optimal_stab, cellsRemoved = cluster_all$cellsRemoved, 
        cellsForClustering = cluster_all$cellsForClustering))
}


#' HC clustering for a number of resolutions
#'
#' @description  performs 40 clustering runs or more depending on windows
#' @param object is a \linkS4class{SingleCellExperiment} object from the
#' train mixed population
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons
#' @param ngenes number of top variable genes to be used
#' @param PCA logical specifying if PCA is used before calculating distance 
#' matrix
#' @param nPCs number of principal components from PCA dimensional reduction to
#' be used
#' @param nRounds number of iterations to remove a selected clusters
#' @param windows a numeric specifying the number of windows to test
#' @return clustering results
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' test <-clustering_scGPS(mixedpop2, remove_outlier = c(1))

clustering_scGPS <- function(object = NULL, ngenes = 1500, 
    windows = seq(0.025:1, by = 0.025), remove_outlier = c(0), nRounds = 1, 
    PCA = FALSE, nPCs = 20) {
    
    # function for the highest resolution clustering (i.e. no window applied, 
    # no cell removal)
    firstRoundClustering <- function(object = NULL) {
        exprs_mat <- assay(object)
        # take the top variable genes
        print("Identifying top variable genes")
        exprs_mat_topVar <- topvar_scGPS(exprs_mat, ngenes = ngenes)
        # exprs_mat_t <- t(exprs_mat_topVar)
        if (PCA == TRUE) {
            # perform PCA dimensionality reduction
            print(paste0("Performing PCA analysis (Note: the variance for ",
                "each cell needs to be >0)"))
            # print('Performing PCA analysis (Note: the variance for each cell 
            # needs to be >0)')
            exprs_mat_topVar_PCA <- prcomp(t(exprs_mat_topVar))
            exprs_mat_t <- as.data.frame(exprs_mat_topVar_PCA$x[, 1:nPCs])
            
        } else {
            exprs_mat_t <- t(exprs_mat_topVar)
        }  # tranpose so that cells are in rows
        
        # calculate distance matrix for the rows
        print("Calculating distance matrix")
        dist_mat <- rcpp_parallel_distance(as.matrix(exprs_mat_t))
        print("Performing hierarchical clustering")
        original.tree <- fastcluster::hclust(as.dist(dist_mat), 
            method = "ward.D2")
        # the original clusters to be used as the reference
        print("Finding clustering information")
        original.clusters <- unname(cutreeDynamic(original.tree, 
            distM = as.matrix(dist_mat), verbose = 0, 
            minClusterSize = round(ncol(dist_mat)/100)))
        original.tree$labels <- original.clusters
        return(list(tree = original.tree, cluster_ref = original.clusters, 
            dist_mat = dist_mat))
    }
    
    # function to remove outlier clusters
    removeOutlierCluster <- function(object = object, 
        remove_outlier = remove_outlier, nRounds = nRounds) {
        
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
            cluster_toRemove <- which(
                filter_out$cluster_ref %in% remove_outlier)
            if (length(cluster_toRemove) > 0) {
                print(paste0("Found ", length(cluster_toRemove), 
                    " cells as outliers at round ", i, " ..."))
                cells_to_remove <- c(cells_to_remove, cluster_toRemove)
                objectTemp <- object[, -cells_to_remove]
                i <- i + 1
            } else {
                print(paste0("No more outliers detected in filtering round ", 
                    i))
                i <- nRounds + 1
            }
        }
        
        filter_out <- firstRoundClustering(objectTemp)
        cluster_toRemove <- which(filter_out$cluster_ref %in% remove_outlier)
        if (length(cluster_toRemove) > 0) {
            print(paste0("Found ", length(cluster_toRemove), 
                " cells as outliers at round ", i, " ..."))
            print(paste0("Select ", i, 
                " removal rounds if you want to remove these cells"))
        }
        
        if (length(cells_to_remove) > 0) {
            output <- list(firstRound_out = filter_out, 
                cellsRemoved = colnames(object[, cells_to_remove]), 
                cellsForClustering = colnames(object[, -cells_to_remove]))
        } else {
            output <- list(firstRound_out = filter_out, 
                cellsRemoved = "No cells removed", 
                cellsForClustering = colnames(object))
        }
        return(output)
    }
    
    
    firstRoundPostRemoval <- removeOutlierCluster(object = object, 
        remove_outlier = remove_outlier, nRounds = nRounds)
    firstRound_out <- firstRoundPostRemoval$firstRound_out
    # return variables for the next step
    original.tree <- firstRound_out$tree
    original.clusters <- firstRound_out$cluster_ref
    dist_mat <- firstRound_out$dist_mat
    
    clustering_param <- list()
    for (i in 1:length(windows)) {
        
        namelist = paste0("window", windows[i])
        toadd <- as.vector(cutreeDynamic(original.tree, 
            distM = as.matrix(dist_mat), minSplitHeight = windows[i], 
            verbose = 0))
        
        print(paste0("writing clustering result for run ", i))
        clustering_param[[i]] <- list(toadd)
        names(clustering_param[[i]]) <- namelist
    }
    
    names(clustering_param[[i]]) <- "cluster_ref"
    print("Done clustering, moving to stability calculation...")
    return(list(list_clusters = clustering_param, tree = original.tree,
        cluster_ref = original.clusters, 
        cellsRemoved = firstRoundPostRemoval$cellsRemoved, 
        cellsForClustering = firstRoundPostRemoval$cellsForClustering))
    
}

#' Subclustering for selected cells
#'
#' @description  performs 40 clustering runs or more depending on windows
#' @param object is a \linkS4class{SingleCellExperiment} object from the
#' train mixed population
#' @param select_cell_index a vector containing indexes for cells in selected 
#' clusters to be reclustered
#' @param ngenes number of genes used for clustering calculations.
#' @param windows a numeric vector specifying the ranges of each window.
#' @return clustering results
#' @export
#' @author Quan Nguyen, 2018-01-31
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' test_SubClustering <-SubClustering_scGPS(mixedpop2,
#'     select_cell_index = c(1:100))


SubClustering_scGPS <- function(object = NULL, ngenes = 1500, 
    windows = seq(0.025:1, by = 0.025), select_cell_index = NULL) {
    
    print("Calculating distance matrix")
    # function for the highest resolution clustering (i.e. no window applied)
    Clustering <- function(object = NULL) {
        exprs_mat <- assay(object)
        exprs_mat_topVar <- topvar_scGPS(exprs_mat, ngenes = ngenes)
        # take the top variable genes
        exprs_mat_t <- t(exprs_mat_topVar)
        dist_mat <- rcpp_parallel_distance(exprs_mat_t)
        print("Performing hierarchical clustering")
        original.tree <- fastcluster::hclust(as.dist(dist_mat), 
            method = "ward.D2")
        # the original clusters to be used as the reference
        print("Finding clustering information")
        original.clusters <- unname(cutreeDynamic(original.tree, 
            distM = as.matrix(dist_mat), verbose = 0))
        original.tree$labels <- original.clusters
        return(list(tree = original.tree, cluster_ref = original.clusters, 
            dist_mat = dist_mat))
    }
    
    Select_out <- function(object = NULL, select_cell_index = NULL) {
        
        Clustering_out <- Clustering(object[, select_cell_index])
        
        return(list(SubClustering_out = Clustering_out, 
            cellsRemoved = colnames(object[, -select_cell_index]), 
            cellsForClustering = colnames(object[, select_cell_index])))
    }
    
    
    SelectCluster_out <- Select_out(object = object, 
        select_cell_index = select_cell_index)
    
    # return variables for the next step
    original.tree <- SelectCluster_out$SubClustering_out$tree
    original.clusters <- SelectCluster_out$SubClustering_out$cluster_ref
    dist_mat <- SelectCluster_out$SubClustering_out$dist_mat
    
    clustering_param <- list()
    for (i in 1:length(windows)) {
        
        namelist = paste0("window", windows[i])
        toadd <- as.vector(cutreeDynamic(original.tree, 
            distM = as.matrix(dist_mat), minSplitHeight = windows[i],
            verbose = 0))
        
        print(paste0("writing clustering result for run ", i))
        clustering_param[[i]] <- list(toadd)
        names(clustering_param[[i]]) <- namelist
    }
    
    names(clustering_param[[i]]) <- "cluster_ref"
    print("Done clustering, moving to stability calculation...")
    return(list(list_clusters = clustering_param, tree = original.tree, 
        cluster_ref = original.clusters, 
        cellsRemoved = SelectCluster_out$cellsRemoved, 
        cellsForClustering = SelectCluster_out$cellsForClustering))
    
}


#' Calculate rand index
#'
#' @description  Comparing clustering results Function for calculating randindex
#' (adapted from the function by Steve Horvath and Luohua Jiang, UCLA, 2003)
#' @param tab a table containing different clustering results in rows
#' @param adjust a logical of whether to use the adjusted rand index
#' @return a randIndex value
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts,
#' GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' cluster_all <-clustering_scGPS(object=mixedpop2)
#'
#' randIndex(table(unlist(cluster_all$list_clusters[[1]]), 
#' cluster_all$cluster_ref))
#'
#' @export
#' @author Quan Nguyen and Michael Thompson, 2018-05-11
#'


randIndex <- function(tab, adjust = TRUE) {
    choosenew <- function(n, k) {
        n <- c(n)
        out1 <- rep(0, length(n))
        for (i in c(1:length(n))) {
            if (n[i] < k) {
                out1[i] <- 0
            } else {
                out1[i] <- choose(n[i], k)
            }
        }
        out1
    }
    a <- 0
    b <- 0
    c <- 0
    d <- 0
    nn <- 0
    m <- nrow(tab)
    n <- ncol(tab)
    for (i in 1:m) {
        c <- 0
        for (j in 1:n) {
            a <- a + choosenew(tab[i, j], 2)
            nj <- sum(tab[, j])
            c <- c + choosenew(nj, 2)
        }
        ni <- sum(tab[i, ])
        b <- b + choosenew(ni, 2)
        nn <- nn + ni
    }
    if (adjust) {
        d <- choosenew(nn, 2)
        adrand <- (a - (b * c)/d)/(0.5 * (b + c) - (b * c)/d)
        adrand
    } else {
        b <- b - a
        c <- c - a
        d <- choosenew(nn, 2) - a - b - c
        rand <- (a + d)/(a + b + c + d)
        rand
    }
}

#' Calculate stability index
#'
#' @description  from clustering results, compare similarity between clusters by
#' adjusted Randindex
#' @param list_clusters is a object from the iterative clustering runs
#' @param cluster_ref is a object from the reference cluster
#' @return a \code{data frame} with stability scores and randIndex results
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' cluster_all <-clustering_scGPS(object=mixedpop2)
#' stab_df <- FindStability(list_clusters=cluster_all$list_clusters,
#'                          cluster_ref = cluster_all$cluster_ref)
#' @export
#' @author Quan Nguyen, 2017-11-25
#'

FindStability <- function(list_clusters = NULL, cluster_ref = NULL) {
    
    #---------------------------------------------------------------------------
    # End function for calculating randindex
    #---------------------------------------------------------------------------
    cluster_index_consec <- list()
    cluster_index_ref <- list()
    
    cluster_index_consec[[1]] <- 1
    cluster_index_ref[[1]] <- randIndex(table(unlist(list_clusters[[1]]), 
        cluster_ref))
    
    for (i in 2:(length(list_clusters))) {
        cluster_index_consec[[i]] <- randIndex(table(unlist(list_clusters[[i]]),
            unlist(list_clusters[[i - 1]])))
        cluster_index_ref[[i]] <- randIndex(table(unlist(list_clusters[[i]]), 
            cluster_ref))
    }
    
    cluster_index_consec <- unlist(cluster_index_consec)
    cluster_index_ref <- unlist(cluster_index_ref)
    run_RandIdx <- as.data.frame(cbind(cluster_index_consec, cluster_index_ref))
    run_RandIdx$order <- row.names(run_RandIdx)
    stability <- run_RandIdx$cluster_index_consec
    
    #---------------------------------------------------------------------------
    # find stability score of the cluster results
    #---------------------------------------------------------------------------
    
    # first get the general counter
    counter = rep(0, length(stability))
    
    counter[1] = 1
    for (i in 2:length(stability) - 1) {
        if (stability[i] == stability[i + 1]) {
            counter[i + 1] <- counter[i] + 1
        } else {
            counter[i + 1] <- 1
        }
    }
    
    # second get the counter location where there is no change
    counter_0 = counter
    
    for (i in 1:(length(counter) - 1)) {
        if (counter_0[i] == 1 & counter_0[i + 1] == 1) {
            counter_0[i] = 0
        }
    }
    
    index_0 <- which(counter_0 == 0)
    # third reset the counter to the count values
    counter_adjusted <- counter
    
    # set counter_0 on the last right or left or middle
    if (length(index_0) > 0) {
        
        # setup counter 0 on the last right
        if (index_0[1] == 1) {
            counter_adjusted[1] = 1
        } else {
            counter_adjusted[1:index_0[1] - 1] <- counter[index_0[1] - 1]
        }
        
        # setup counter 0 on the last left
        if (index_0[length(index_0)] == length(counter)) {
            length_id0 <- length(index_0)
            counter_adjusted[index_0[length_id0]] = 1
        } else {
            length_id0 <- length(index_0)
            counter_adjusted[(index_0[length_id0] + 1):length(counter)] <- 
                counter[length(counter)]
            counter_adjusted[index_0[length_id0]] = 1
        }
        
        # setup counter 0 in the middle (index_0 >=3)
        if (length(index_0) > 2) {
            for (i in 1:(length(index_0) - 1)) {
                if (index_0[i + 1] - index_0[i] > 1) {
                    # assign value to the last count
                    counter_adjusted[(index_0[i] + 1):(index_0[i + 1] - 1)] = 
                        counter[index_0[i + 1] - 1]
                    counter_adjusted[index_0[i]] = 1  #reset
                } else {
                    counter_adjusted[(index_0[i] + 1)] = 1
                }
            }
        }
        
    }
    
    
    run_RandIdx$stability_count <- counter_adjusted
    print("Done calculating stability...")
    return(run_RandIdx)
    #---------------------------------------------------------------------------
    # Done stability score of the cluster results
    #---------------------------------------------------------------------------
}

#'Find the optimal cluster
#'
#' @description from calculated stability based on Rand indexes for consecutive
#' clustering run, find the resolution (window), where the stability is the 
#' highest
#' @param run_RandIdx is a \code{data frame} object from iterative clustering 
#' runs
#' @param list_clusters is a \code{list} object containing 40 clustering results
#' @param bagging is a logical that is true if bagging is to be performed, 
#' changes return
#' @param windows a numeric vector specifying the ranges of each window.
#' @return bagging == FALSE => a \code{list} with optimal stability, cluster 
#' count and summary stats bagging == TRUE => a \code{list} with high res 
#' cluster count, optimal cluster count and keystats
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' cluster_all <-clustering_scGPS(object=mixedpop2)
#' stab_df <- FindStability(list_clusters=cluster_all$list_clusters,
#'                          cluster_ref = cluster_all$cluster_ref)
#' optimal_stab <- FindOptimalStability(list_clusters = 
#'     cluster_all$list_clusters, stab_df, bagging = FALSE)
#'


FindOptimalStability <- function(list_clusters, run_RandIdx, bagging = FALSE, 
    windows = seq(0.025:1, by = 0.025)) {
    print("Start finding optimal clustering...")
    
    window_param <- length(windows)
    
    # get the number of cluster
    t <- lapply(list_clusters, function(x) {
        length(unique(unlist(x)))
    })
    run_RandIdx$cluster_count <- as.vector(as.numeric(t))
    
    #---------------------------------------------------------------------------
    # Diagnostic plot for comparing clustering results
    #---------------------------------------------------------------------------
    run_RandIdx$stability_count <- run_RandIdx$stability_count/window_param
    
    
    
    # CHANGED HERE FOR .025 -->> based on the seq
    KeyStats <- as.data.frame(cbind(as.numeric(run_RandIdx$order) * windows[1], 
        run_RandIdx$stability_count, run_RandIdx$cluster_index_ref, 
        run_RandIdx$cluster_index_consec))
    
    colnames(KeyStats) <- c("Height", "Stability", "RandIndex", "ConsecutiveRI")
    
    KeyStats$Height <- as.character(KeyStats$Height)
    day_melt <- reshape2::melt(KeyStats, id = "Height")
    day_melt$Height <- as.numeric(day_melt$Height)
    p <- ggplot(day_melt)
    p <- p + geom_line(aes(x = day_melt$Height, y = day_melt$value, 
        colour = day_melt$variable), size = 2) + theme_bw() + 
        theme(axis.text = element_text(size = 24), 
            axis.title = element_text(size = 24)) + 
        theme(legend.text = element_text(size = 24)) + 
        theme(legend.title = element_blank()) + 
        xlab("Parameter from 0.025 to 1") + ylab("Scores") + 
        theme(panel.border = element_rect(colour = "black", 
            fill = NA, size = 1.5)) + 
        guides(colour = FALSE)
    
    #---------------------------------------------------------------------------
    # End diagnostic plot for comparing clustering results
    #---------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------
    # Find the optimal parameter for clustering
    #---------------------------------------------------------------------------
    
    #------NEW OPTIMAL METHOD ALLOWING FOR MAX RESOLUTION-----------------------
    KeyStats$Cluster_count <- run_RandIdx$cluster_count
    Cluster_count <- KeyStats$Cluster_count
    optimal_param = 1
    St <- KeyStats$Stability
    St_max <- St[which.max(St)]
    
    if ((St[window_param] > 0.5) && (St[window_param] == St_max)) {
        optimal_param = window_param
    } else {
        concat_St <- vector()
        for (i in 1:window_param) {
            if (Cluster_count[i] != min(Cluster_count)) {
                concat_St <- c(concat_St, St[i])
            }
        }
        optimal_param = which.max(concat_St)[1]
        if (optimal_param != 1) {
            optimal_param = optimal_param - 1
        }
    }
    
    if (bagging == TRUE) {
        output <- list(HighestRes = KeyStats$Cluster_count[1], 
            OptimalClust = KeyStats$Cluster_count[optimal_param], 
            KeyStats = KeyStats)
    } else {
        print("Done finding optimal clustering...")
        # Final result
        output <- list(StabilityPlot = p, KeyStats = KeyStats, 
            OptimalRes = KeyStats$Height[optimal_param], 
            OptimalClust = KeyStats$Cluster_count[optimal_param])
    }
    
    return(output)
    #---------------------------------------------------------------------------
    # Done finding the optimal parameter for clustering
    #---------------------------------------------------------------------------
    
}


#' Plot dendrogram tree for CORE result
#'
#' @description This function plots CORE and all clustering results underneath
#' @param original.tree the original dendrogram before clustering
#' @param list_clusters a list containing clustering results for each of the 
# resolution run
#' @param color_branch is a vector containing user-specified colors (the number
#' of unique colors should be equal or larger than the number of clusters). This
#' parameter allows better selection of colors for the display.
#' @return a plot with clustering bars underneath the tree
#' @examples
#' day5 <- sample2
#' cellnames <- colnames(day5$dat5_counts)
#' cluster <-day5$dat5_clusters
#' cellnames <-data.frame('Cluster'=cluster, 'cellBarcodes' = cellnames)
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts,
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = cellnames)
#' CORE_cluster <- CORE_scGPS(mixedpop2, remove_outlier = c(0))
#' plot_CORE(CORE_cluster$tree, CORE_cluster$Cluster)

plot_CORE <- function(original.tree, list_clusters = NULL, 
    color_branch = NULL) {
    
    n <- length(unique(unlist(list_clusters[[1]])))
    qual_col_pals = 
        RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 
        "qual", ] 
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, 
        qual_col_pals$maxcolors, rownames(qual_col_pals)))
    if (is.null(color_branch)) {
        color_branch <- col_vector
    }
    #---------------------------------------------------------------------------
    # Function to plot dendrogram and color The plot_CORE function implements 
    # the code by Steve Horvarth, Peter Langelder, and Tal Galili (used in the 
    # WGCNA package, version 1.61) The code were customised to plot scGPS 
    # object and clustering colors
    #---------------------------------------------------------------------------
    col_all <- matrix(unlist(list_clusters), ncol = length(list_clusters))
    col_all <- as.data.frame(col_all)
    # remove branch labels
    original.tree$labels <- rep("", length(original.tree$labels))
    
    plotDendroAndColors <- function(dendro, colors, groupLabels = NULL, 
        rowText = NULL, rowTextAlignment = c("left", "center", "right"), 
        rowTextIgnore = NULL, textPositions = NULL, setLayout = TRUE, 
        autoColorHeight = TRUE, colorHeight = 0.2, rowWidths = NULL, 
        dendroLabels = NULL, addGuide = FALSE, guideAll = FALSE, 
        guideCount = 50, guideHang = 0.2, addTextGuide = FALSE, 
        cex.colorLabels = 0.8, cex.dendroLabels = 0.9, cex.rowText = 0.8, 
        marAll = c(1, 5, 3, 1), saveMar = TRUE, abHeight = NULL, 
        abCol = "red", ...) {
        oldMar = par("mar")
        if (!is.null(dim(colors))) {
            nRows = dim(colors)[2]
        } else nRows = 1
        if (!is.null(rowText)) 
            nRows = nRows + if (is.null(textPositions)) 
                nRows else length(textPositions)
        if (autoColorHeight) 
            colorHeight = 0.2 + 0.3 * (1 - exp(-(nRows - 1)/6))
        if (setLayout) 
            layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight, 
                colorHeight))
        par(mar = c(0, marAll[2], marAll[3], marAll[4]))
        plot(dendro, labels = dendroLabels, cex = cex.dendroLabels, ...)
        if (addGuide) 
            WGCNA::addGuideLines(dendro, count = if (guideAll) 
                length(dendro$height) + 1 else guideCount, hang = guideHang)
        if (!is.null(abHeight)) 
            abline(h = abHeight, col = abCol)
        par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
        plotColorUnderTree(dendro, colors, groupLabels, 
            cex.rowLabels = cex.colorLabels, rowText = rowText, 
            rowTextAlignment = rowTextAlignment, rowTextIgnore = rowTextIgnore, 
            textPositions = textPositions, cex.rowText = cex.rowText, 
            rowWidths = rowWidths, addTextGuide = addTextGuide)
        if (saveMar) 
            par(mar = oldMar)
    }
    #---------------------------------------------------------------------------
    # plotColorUnderTree
    #---------------------------------------------------------------------------
    plotColorUnderTree <- function(dendro, colors, rowLabels = NULL, 
        rowWidths = NULL, rowText = NULL, 
        rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL, 
        textPositions = NULL, addTextGuide = TRUE, cex.rowLabels = 1, 
        cex.rowText = 0.8, ...) {
        
        plotOrderedColors(dendro$order, colors = colors, rowLabels = rowLabels, 
            rowWidths = rowWidths, rowText = rowText, 
            rowTextAlignment = rowTextAlignment, rowTextIgnore = rowTextIgnore, 
            textPositions = textPositions, addTextGuide = addTextGuide, 
            cex.rowLabels = cex.rowLabels, cex.rowText = cex.rowText, 
            startAt = 0, ...)
    }
    #---------------------------------------------------------------------------
    # plotOrderedColors
    #---------------------------------------------------------------------------
    plotOrderedColors <- function(order, colors, rowLabels = NULL, 
        rowWidths = NULL, rowText = NULL, 
        rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL, 
        textPositions = NULL, addTextGuide = TRUE, cex.rowLabels = 1, 
        cex.rowText = 0.8, startAt = 0, ...) {
        colors = as.matrix(colors)
        dimC = dim(colors)
        if (is.null(rowLabels) & (length(dimnames(colors)[[2]]) == dimC[2])) 
            rowLabels = colnames(colors)
        sAF = options("stringsAsFactors")
        options(stringsAsFactors = FALSE)
        on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)
        nColorRows = dimC[2]
        if (length(order) != dimC[1]) 
            stop(paste0("ERROR: length of colors vector not compatible with ",
                "number of objects in 'order'."))
        C = colors[order, , drop = FALSE]
        step = 1/(dimC[1] - 1 + 2 * startAt)
        barplot(height = 1, col = "white", border = FALSE, space = 0, 
            axes = FALSE)
        charWidth = strwidth("W")/2
        if (!is.null(rowText)) {
            if (is.null(textPositions)) 
                textPositions = c(1:nColorRows)
            if (is.logical(textPositions)) 
                textPositions = c(1:nColorRows)[textPositions]
            nTextRows = length(textPositions)
        } else nTextRows = 0
        nRows = nColorRows + nTextRows
        ystep = 1/nRows
        if (is.null(rowWidths)) {
            rowWidths = rep(ystep, nColorRows + nTextRows)
        } else {
            if (length(rowWidths) != nRows) 
                stop(paste0("plotOrderedColors: Length of 'rowWidths' must ",
                    "equal the total number of rows."))
            rowWidths = rowWidths/sum(rowWidths)
        }
        hasText = rep(0, nColorRows)
        hasText[textPositions] = 1
        csPosition = cumsum(c(0, hasText[-nColorRows]))
        colorRows = c(1:nColorRows) + csPosition
        rowType = rep(2, nRows)
        rowType[colorRows] = 1
        physicalTextRow = c(1:nRows)[rowType == 2]
        yBottom = c(0, cumsum(rowWidths[nRows:1]))
        yTop = cumsum(rowWidths[nRows:1])
        if (!is.null(rowText)) {
            rowTextAlignment = match.arg(rowTextAlignment)
            rowText = as.matrix(rowText)
            textPos = list()
            textPosY = list()
            textLevs = list()
            for (tr in 1:nTextRows) {
                charHeight = max(strheight(rowText[, tr], cex = cex.rowText))
                width1 = rowWidths[physicalTextRow[tr]]
                nCharFit = floor(width1/charHeight/1.7/par("lheight"))
                if (nCharFit < 1) 
                    stop(paste0("Rows are too narrow to fit text. Consider ",
                        "decreasing cex.rowText."))
                set = textPositions[tr]
                textLevs[[tr]] = sort(unique(rowText[, tr]))
                textLevs[[tr]] = textLevs[[tr]][!textLevs[[tr]] %in% 
                    rowTextIgnore]
                nLevs = length(textLevs[[tr]])
                textPos[[tr]] = rep(0, nLevs)
                orderedText = rowText[order, tr]
                for (cl in 1:nLevs) {
                    ind = orderedText == textLevs[[tr]][cl]
                    sind = ind[-1]
                    ind1 = ind[-length(ind)]
                    starts = c(if (ind[1]) 1 else NULL, which(!ind1 & sind) + 1)
                    ends = which(c(ind1 & !sind, ind[length(ind)]))
                    if (length(starts) == 0) 
                        starts = 1
                    if (length(ends) == 0) 
                        ends = length(ind)
                    if (ends[1] < starts[1]) 
                        starts = c(1, starts)
                    if (ends[length(ends)] < starts[length(starts)]) 
                        ends = c(ends, length(ind))
                    lengths = ends - starts
                    long = which.max(lengths)
                    textPos[[tr]][cl] = switch(rowTextAlignment, 
                        left = starts[long], center = (starts[long] + 
                        ends[long])/2 + 0.5, right = ends[long] + 1)
                }
                if (rowTextAlignment == "left") {
                    yPos = seq(from = 1, to = nCharFit, by = 1)/(nCharFit + 1)
                } else {
                    yPos = seq(from = nCharFit, to = 1, by = -1)/(nCharFit + 1)
                }
                textPosY[[tr]] = rep(yPos, ceiling(nLevs/nCharFit) + 5)[
                    1:nLevs][rank(textPos[[tr]])]
            }
        }
        jIndex = nRows
        if (is.null(rowLabels)) 
            rowLabels = c(1:nColorRows)
        C[is.na(C)] = "grey"
        for (j in 1:nColorRows) {
            jj = jIndex
            ind = (1:dimC[1])
            xl = (ind - 1.5 + startAt) * step
            xr = (ind - 0.5 + startAt) * step
            yb = rep(yBottom[jj], dimC[1])
            yt = rep(yTop[jj], dimC[1])
            if (is.null(dim(C))) {
                rect(xl, yb, xr, yt, col = as.character(C), 
                    border = as.character(C))
            } else {
                rect(xl, yb, xr, yt, col = as.character(C[, j]), 
                    border = as.character(C[, j]))
            }
            text(rowLabels[j], pos = 2, x = -charWidth/2 + xl[1], 
                y = (yBottom[jj] + yTop[jj])/2, cex = cex.rowLabels, xpd = TRUE)
            textRow = match(j, textPositions)
            if (is.finite(textRow)) {
                jIndex = jIndex - 1
                xt = (textPos[[textRow]] - 1.5) * step
                xt[xt < par("usr")[1]] = par("usr")[1]
                xt[xt > par("usr")[2]] = par("usr")[2]
                yt = yBottom[jIndex] + (yTop[jIndex] - yBottom[jIndex]) * 
                    (textPosY[[textRow]] + 1/(2 * nCharFit + 2))
                nt = length(textLevs[[textRow]])
                if (addTextGuide) 
                    for (l in 1:nt) lines(c(xt[l], xt[l]), 
                        c(yt[l], yTop[jIndex]), col = "darkgrey", lty = 3)
                textAdj = c(0, 0.5, 1)[match(rowTextAlignment, c("left", 
                    "center", "right"))]
                text(textLevs[[textRow]], x = xt, y = yt, adj = c(textAdj, 1), 
                    xpd = TRUE, cex = cex.rowText)
            }
            jIndex = jIndex - 1
        }
        for (j in 0:(nColorRows + nTextRows)) lines(x = c(0, 1), 
            y = c(yBottom[j + 1], yBottom[j + 1]))
    }
    
    #---------------------------------------------------------------------------
    # Start selecting colors
    #---------------------------------------------------------------------------
    
    # this color range assumes a maximum 15 clusters color_range <-
    # c(,'#f4215a','#00c9ce','#824f00','#714973','#006837','#ffa5a1','#e4c27f')
    color_range <- color_branch
    number_colors <- length(unique(unlist(col_all)))
    col_all2 <- col_all
    for (i in 1:number_colors) {
        col_all2[col_all2 == i] <- as.character(color_range[i])
    }
    colnames(col_all2) <- gsub("V", "", colnames(col_all2))
    
    plotDendroAndColors(original.tree, col_all2)
    
}


#' plot one single tree with the optimal clustering result
#'
#' @description after an optimal cluster has been identified, users may use this
#' function to plot the resulting dendrogram with the branch colors represent 
#' clutering results
#' @param original_tree a dendrogram object
#' @param optimal_cluster a vector of cluster IDs for cells in the dendrogram
#' @param shift a numer specifying the gap between the dendrogram and the 
#' colored
#' @param values a vector containing color values of the branches and the
#' colored bar underneath the tree bar underneath the dendrogram. This
#' parameter allows better selection of colors for the display.
#' @export
#' @return a plot with colored braches and bars for the optimal clustering 
#' result
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' CORE_cluster <- CORE_scGPS(mixedpop2, remove_outlier = c(0))
#' key_height <- CORE_cluster$optimalClust$KeyStats$Height
#' optimal_res <- CORE_cluster$optimalClust$OptimalRes
#' optimal_index = which(key_height == optimal_res)
#' plot_optimal_CORE(original_tree= CORE_cluster$tree, 
#'     optimal_cluster = unlist(CORE_cluster$Cluster[optimal_index]), shift = -2000)
#'

plot_optimal_CORE <- function(original_tree, optimal_cluster = NULL, 
    shift = -100, values = NULL) {
    if (is.null(values)) {
        n <- length(unique(optimal_cluster))
        qual_col_pals = RColorBrewer::brewer.pal.info[
            RColorBrewer::brewer.pal.info$category == "qual", ]
        col_vector = unlist(mapply(RColorBrewer::brewer.pal, 
            qual_col_pals$maxcolors, rownames(qual_col_pals)))
        values <- col_vector
    }
    
    print("Ordering and assigning labels...")
    dendro.obj <- as.dendrogram(original_tree)
    # Sort clusters by order in dendrogram
    ordered.clusters <- optimal_cluster[stats::order.dendrogram(dendro.obj)]
    # Count table
    cluster.df <- as.data.frame(table(optimal_cluster))
    
    # Sort cluster sizes in same order.
    dendro.labels <- cluster.df$Freq[unique(ordered.clusters)]
    
    index_labels <- unique(ordered.clusters)
    # for the first index
    index_to_overwriteNA <- c(rep(NA, length(dendro.labels)))
    index_to_overwriteNA[1] <- c(round(dendro.labels[1]/2))
    # the loop only uses information from dendro.labels
    for (i in 2:length(dendro.labels)) {
        index_to_overwriteNA[i] <- sum(dendro.labels[1:i - 1]) + 
            round(dendro.labels[i]/2)
        print(i)
        print(index_to_overwriteNA)
    }
    
    branch_names <- rep(NA, length(optimal_cluster))
    
    # index_labels_cluster_names <- paste0('C', index_labels,) not
    # dendro_labels_names <- paste( dendro.labels,index_labels_cluster_names)
    # branch_names[index_to_overwriteNA] <- index_labels_cluster_names
    
    # Apply labels directly to dendrogram
    coloured.dendro <- dendextend::branches_attr_by_clusters(dendro.obj, 
        clusters = ordered.clusters, attr = "col", values = values)
    # need to specify the name space dendextend::
    coloured.dendro <- dendextend::set(coloured.dendro, "labels", branch_names)
    
    # make branch lines bigger
    coloured.dendro <- dendextend::set(coloured.dendro, "branches_lwd", 2)
    print("Plotting the colored dendrogram now....")
    plot(coloured.dendro)
    print("Plotting the bar underneath now....")
    dendro.colours <- unique(dendextend::get_leaves_branches_col(
        coloured.dendro))
    coloured.order <- stats::order.dendrogram(coloured.dendro)
    sorted.levels <- dendextend::sort_levels_values(as.vector(optimal_cluster)[
        coloured.order])
    
    sorted.levels <- sorted.levels[match(seq_along(coloured.order),
        coloured.order)]
    dendextend::colored_bars(dendro.colours[sorted.levels], coloured.dendro,
        rowLabels = "Cluster", y_shift = shift)
}
