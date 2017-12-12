#' Main clustering CORE V2.0
#'
#' @description  CORE is an algorithm to generate reproduciable clustering,
#' CORE is first implemented in ascend R package. Here, CORE V2.0 introduces several new
#' functionalities, including three key feature:
#' fast implementation with C++ and paralellisation options allowing clustering
#' of hundreds of thousands of cells (ongoing development), outlier revomal important if singletons
#' exist (done), a number of dimensionality reduction methods including the imputation
#' implementation (CIDR) for confirming clustering results (done), and an option
#' to select the number of optimisation tree height windows for increasing resolution
#' @param mixedpop1 is a \linkS4class{SingleCellExperiment} object from the train
#' mixed population
#' @param windows a numeric specifying the number of windows to test
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons
#' @return a \code{list} with clustering results of all iterations, and a selected
#' optimal resolution
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                          CellMetadata = day5$dat5_clusters)
#' test <- CORE_scGPS(mixedpop2, remove_outlier = c(1))
#' @export
#' @author Quan Nguyen, 2017-11-25

CORE_scGPS <-function(mixedpop = NULL, windows = seq(0.025:1, by=0.025),
                      remove_outlier = c(0)){
  cluster_all <-clustering_scGPS(object=mixedpop, windows = windows, remove_outlier = remove_outlier)

  stab_df <- FindStability(list_clusters=cluster_all$list_clusters,
    cluster_ref = cluster_all$cluster_ref)

  optimal_stab <- FindOptimalStability(list_clusters = list_clusters, stab_df)

  return(list("Cluster" = cluster_all$list_clusters,
              "tree" = cluster_all$tree, "optimalClust" = optimal_stab))
}

#' Iterative HC clustering
#'
#' @description  performs 40 clustering runs or more depending on windows
#' @param mixedpop1 is a \linkS4class{SingleCellExperiment} object from the
#' train mixed population
#' @param remove_outlier a vector containing IDs for clusters to be removed
#' the default vector contains 0, as 0 is the cluster with singletons
#' @return a \code{matrix} with Eucleadean distance used for clustreting
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                          CellMetadata = day5$dat5_clusters)
#' test <-clustering_scGPS(mixedpop2, remove_outlier = c(1))
#' @export
#' @author Quan Nguyen, 2017-11-25
#'

clustering_scGPS <- function(object = NULL, windows = seq(0.025:1, by=0.025),
  remove_outlier = c(0)){

  print("Calculating distance matrix")
  #function for the highest resolution clustering (i.e. no window applied)
  firstRoundClustering <- function(object = NULL){
    exprs_mat <- assay(object)
    exprs_mat_topVar <- topvar_scGPS(exprs_mat, ngenes = 1500)
    #Take the top variable genes
    #make a transpose
    exprs_mat_t <-t(exprs_mat_topVar)
    dist_mat <- dist(exprs_mat_t)
    print("Performing hierarchical clustering")
    original.tree <- hclust(dist_mat, method="ward.D2")
    #The original clusters to be used as the reference
    print("Finding clustering information")
    original.clusters <- unname(cutreeDynamic(original.tree, distM=as.matrix(dist_mat), verbose=0))
    original.tree$labels <- original.clusters
    return(list("tree" = original.tree, "cluster_ref" = original.clusters, "dist_mat" = dist_mat))
    }

  removeOutlierCluster <-function(object = NULL,object_rmOutlier = object_rmOutlier){
    #check for singletons
    firstRound_out <- firstRoundClustering(object)
    firstRound_cluster <- as.data.frame(table(firstRound_out$cluster_ref))
    cluster_toRemove <- which(firstRound_out$cluster_ref %in% remove_outlier)

    if(length(cluster_toRemove) > 0){
      print("Removing outlier clusters...")
      object_rmOutlier <- object[,-cluster_toRemove]
      firstRound_out <- firstRoundClustering(object_rmOutlier)
    }
    return(firstRound_out)
  }

  firstRound_out <- removeOutlierCluster(object, object_rmOutlier)
  #return variables for the next step
  original.tree <- firstRound_out$tree
  original.clusters <-firstRound_out$cluster_ref
  dist_mat <-firstRound_out$dist_mat

  clustering_param <-list()
  for (i in 1:length(windows)){

  namelist =paste0("window",windows[i])
  toadd <-as.vector(cutreeDynamic(original.tree, distM=as.matrix(dist_mat),
                        minSplitHeight=windows[i], verbose=0))

  print(paste0("writing clustering result for run ", i))
  clustering_param[[i]] <-  list(toadd)
  names(clustering_param[[i]]) <- namelist
  }

  names(clustering_param[[i]]) <- "cluster_ref"
  print("Done clustering, moving to stability calculation...")
  return(list("list_clusters" = clustering_param, "tree" = original.tree,
              "cluster_ref" = original.clusters))

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
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                         CellMetadata = day5$dat5_clusters)
#' cluster_all <-clustering_scGPS(object=mixedpop2)
#' stab_df <- FindStability(list_clusters=cluster_all$list_clusters, cluster_ref = cluster_all$cluster_ref)
#' @export
#' @author Quan Nguyen, 2017-11-25
#'
FindStability <- function(list_clusters=NULL, cluster_ref =NULL){
  choosenew <- function(n,k) {
    n <- c(n); out1 <- rep(0,length(n));
    for (i in c(1:length(n)) ){
      if ( n[i]<k ) {out1[i] <- 0}
      else {out1[i] <- choose(n[i],k) }
    }
    out1
  }
  #-----------------------------------------------------------------------------
  #Comparing clustering results
  #Function for calculating randindex (adapted from the function by Steve Horvath
  #and Luohua Jiang, UCLA, 2003)
  #-----------------------------------------------------------------------------

  randIndex <- function (tab, adjust = TRUE)
  {
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
    }
    else {
      b <- b - a
      c <- c - a
      d <- choosenew(nn, 2) - a - b - c
      rand <- (a + d)/(a + b + c + d)
      rand
    }
  }
  #-----------------------------------------------------------------------------
  #End function for calculating randindex
  #-----------------------------------------------------------------------------
  cluster_index_consec <-list()
  cluster_index_ref <-list()

  cluster_index_consec[[1]] <-1
  cluster_index_ref[[1]] <-randIndex(table(unlist(list_clusters[[1]]),cluster_ref))

  for (i in 2:(length(list_clusters))){
    cluster_index_consec[[i]] <-randIndex(table(unlist(list_clusters[[i]]),
      unlist(list_clusters[[i-1]])))
    cluster_index_ref[[i]] <-randIndex(table(unlist(list_clusters[[i]]),cluster_ref))
  }

  cluster_index_consec <-unlist(cluster_index_consec)
  cluster_index_ref <-unlist(cluster_index_ref)
  run_RandIdx <-as.data.frame(cbind(cluster_index_consec,cluster_index_ref))
  run_RandIdx$order <-row.names(run_RandIdx)
  stability <-run_RandIdx$cluster_index_consec

  #-----------------------------------------------------------------------------
  #find stability score of the cluster results
  #-----------------------------------------------------------------------------

  #first get the general counter
  counter=rep(0, length(stability))

  counter[1]=1
  for (i in 2:length(stability)-1){
    if(stability[i]==stability[i+1]){counter[i+1]<-counter[i]+1} else {
      counter[i+1] <- 1}
  }

  #second get the counter location where there is no change
  counter_0 = counter
  index_0 <- which(counter_0 == 0)

  for (i in 1:length(counter)){
    if(counter_0[i] == 1 & counter_0[i+1] == 1){counter_0[i] = 0 }
  }

  #third reset the counter to the count values
  counter_adjusted <- counter
  #setup counter 0 on the last right
  if (length(index_0) > 0){
    if(index_0[1] == 1){counter_adjusted[1] = 1}else{
    counter_adjusted[1:index_0[1] - 1] <- counter[index_0[1] - 1]
    }

  #setup counter 0 on the last left
  if (index_0[length(index_0)]==length(counter)){
    length_id0 <- length(index_0)
    counter_adjusted[index_0[length_id0]]=1} else {
    counter_adjusted[(index_0[length_id0]+1):length(counter)] <- counter[length(counter)]
    counter_adjusted[index_0[length_id0]]=1
    } #to be considered
  #setup counter 0 in the middle
  for (i in 2:length(index_0)-1){
    counter_adjusted[(index_0[i]+1):(index_0[i+1]-1)]= counter[index_0[i+1]-1]
    counter_adjusted[index_0[i]]=1
  }
  }

  run_RandIdx$stability_count <-counter_adjusted
  print("Done calculating stability...")
  return(run_RandIdx)
  #-----------------------------------------------------------------------------
  #Done stability score of the cluster results
  #-----------------------------------------------------------------------------
}

#'Find the optimal cluster
#'
#' @description from calculated stability based on Rand indexes for consecutive
#' clustering run, find the resolution (window), where the stability is the highest
#' @param run_RandIdx is a \code{data frame} object from iterative clustering runs
#' @param list_clusters is a \code{list} object containing 40 clustering results
#' @return a \code{list} with optimal stability, cluster count and summary stats
#' @examples
#' day5 <- sample2
#' mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
#'                         CellMetadata = day5$dat5_clusters)
#' cluster_all <-clustering_scGPS(object=mixedpop2)
#' stab_df <- FindStability(list_clusters=cluster_all$list_clusters, cluster_ref = cluster_all$cluster_ref)
#' optimal_stab <- FindOptimalStability(list_clusters = list_clusters, stab_df)
#' @export
#' @author Quan Nguyen, 2017-11-25
#'

FindOptimalStability <- function(list_clusters, run_RandIdx){

  print("Start finding optimal clustering...")
  #get the number of cluster
  t <-lapply(list_clusters, function(x){length(unique(unlist(x)))})
  run_RandIdx$cluster_count <-as.vector(as.numeric(t))

  #-----------------------------------------------------------------------------
  #Diagnostic plot for comparing clustering results
  #-----------------------------------------------------------------------------
  run_RandIdx$stability_count <-run_RandIdx$stability_count/40

  KeyStats <-as.data.frame(cbind(as.numeric(run_RandIdx$order)*0.025,
                                 run_RandIdx$stability_count, run_RandIdx$cluster_index_ref,
                                 run_RandIdx$cluster_index_consec))

  colnames(KeyStats) <-c('Height', 'Stability', 'RandIndex', 'ConsecutiveRI')
  library(ggplot2)
  library(reshape2)

  KeyStats$Height <-as.character(KeyStats$Height)
  day_melt <-melt(KeyStats, id='Height')
  day_melt$Height <-as.numeric(day_melt$Height)
  p<-ggplot(day_melt)
  p <- p+geom_line(aes(x=Height, y=value,  colour=variable), size=2) + theme_bw() +
    theme(axis.text=element_text(size=24), axis.title=element_text(size=24)) +
    theme(legend.text=element_text(size=24))+theme(legend.title=element_blank()) +
    xlab('Parameter from 0.025 to 1') +ylab('Scores') +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    guides(colour=FALSE)

  #-----------------------------------------------------------------------------
  #End diagnostic plot for comparing clustering results
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  #Find the optimal parameter for clustering
  #-----------------------------------------------------------------------------

  KeyStats$Cluster_count <-run_RandIdx$cluster_count
  optimal_param = 0
  St <- KeyStats$Stability
  St_unique <-unique(St)
  St_max <-St[which.max(St)]

  St_max_middle <- St_unique[-which(St_unique == St[1])]
  St_max_middle <- St_max_middle[-which(St_max_middle == St[40])]
  St_max_middle <-St_max_middle[which.max(St_max_middle)]

  if(St[40] >= St[1]){St_Minus_max =  St - St[40]} else {
    St_Minus_max =  St_max - St[1]
    }

  if(St[40] > St[1]){ if(St[1] >0.5){optimal_param = 1; break} else {
    for (i in 2:39){
      if((St_Minus_max[i] == 0) && ((St_Minus_max[i+1] <0) )) {optimal_param = i; break}
    }}}
  if(St[1] > St[40]){ if(St[1] >0.5){optimal_param = 1; break} else {
    for (i in 2:39){
      if((St_Minus_max[i] == 0) && ((St_Minus_max[i-1] <0) )) {optimal_param = i; break}
    }}}

  if (optimal_param == 0){ for (i in 2:39){ if (St[i] == St_max_middle) {optimal_param = i; break}}}

  print("Done finding optimal clustering...")
  #Final result
   return(list("StabilityPlot" = p, "KeyStats" = KeyStats,
     "OptimalRes" = KeyStats$Height[optimal_param],
     "OptimalClust" = KeyStats$Cluster_count[optimal_param])
     )
  #-----------------------------------------------------------------------------
  #Done finding the optimal parameter for clustering
  #-----------------------------------------------------------------------------

}
