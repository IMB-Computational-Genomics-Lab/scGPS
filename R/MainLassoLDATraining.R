#' Main model training function for finding the best model that characterises a
#' subpopulation
#'
#' @description  Training a haft of all cells to find optimal ElasticNet and LDA
#' models to predict a subpopulation
#' @param mixedpop1 is a \linkS4class{SingleCellExperiment} object from the 
#' train mixed population 
#' @param mixedpop2 is a \linkS4class{SingleCellExperiment} object from the 
#' target mixed population
#' @param genes a vector of gene names (for ElasticNet shrinkage); gene symbols 
#' must be in the same format with gene names in subpop2. Note that genes are 
#' listed by the order of importance, e.g. differentially expressed genes that 
#' are most significan, so that if the gene list contains too many genes, only 
#' the top 500 genes are used.
#' @param cluster_mixedpop1 a vector of cluster assignment in mixedpop1
#' @param c_selectID a selected number to specify which subpopulation to be used
#' for training
#' @param out_idx a number to specify index to write results into the list 
#' output. This is needed for running bootstrap.
#' @param standardize a logical value specifying whether or not to standardize 
#' the train matrix
#' @param listData list to store output in
#' @param trainset_ratio a number specifying the proportion of cells to be part
#' of the training subpopulation
#' @param LDA_run logical, if the LDA run is added to compare to ElasticNet 
#' @return a \code{list} with prediction results written in to the indexed 
#' \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#'
#' c_selectID<-1
#' out_idx<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts,
#' GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' listData  <- training_scGPS(genes, cluster_mixedpop1 = colData(mixedpop1)[, 1],
#'                             mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID,
#'                             listData =list(), out_idx=out_idx, trainset_ratio = 0.5)
#' names(listData)
#' listData$Accuracy

training_scGPS <- function(genes = NULL, cluster_mixedpop1 = NULL, mixedpop1 = NULL, 
    mixedpop2 = NULL, c_selectID = NULL, listData = list(), out_idx = 1, standardize = TRUE, 
    trainset_ratio = 0.5, LDA_run = FALSE) {
    # subsammpling-----------------------------------------------------------------
    # taking a subsampling size of trainset_ratio of the cluster_select out for
    # training
    subpop1cluster_indx <- which(cluster_mixedpop1 == c_selectID)  #class 1 
    print(paste0("Total ", length(subpop1cluster_indx), " cells as source subpop"))
    subremaining1_indx <- which(cluster_mixedpop1 != c_selectID)  #remaining class 
    print(paste0("Total ", length(subremaining1_indx), " cells in remaining subpops"))
    subsampling <- round(length(subpop1cluster_indx) * trainset_ratio)
    print(paste0("subsampling ", subsampling, " cells for training source subpop"))
    subpop1_train_indx <- sample(subpop1cluster_indx, subsampling, replace = FALSE)
    # check if there is a very big cluster present in the dataset C/2 = SubSampling
    # >length(total - C, which is cluster_compare)
    if (length(subremaining1_indx) > subsampling) {
        print(paste0("subsampling ", subsampling, " cells in remaining subpops for training"))
        subremaining1_train_indx <- sample(subremaining1_indx, subsampling, replace = FALSE)
    } else {
        print(paste0("no subsampling, using ", length(subremaining1_indx), " cells in remaining subpops for training"))
        subremaining1_train_indx <- subremaining1_indx
    }
    # done subsampling----------------------------------------------------------
    
    # select top 500 genes for the gene list used in the model------------------ e.g.
    # The DE_result is a sorted list by p-values
    if (length(genes) > 500) {
        genes <- genes[1:500]
    }
    # select genes in both mixedpop1 and mixedpop2
    names1 <- elementMetadata(mixedpop1)[, 1]
    subpop1_selected_genes <- names1[which(names1 %in% genes)]
    names2 <- elementMetadata(mixedpop2)[, 1]
    genes_in_both <- names2[which(names2 %in% subpop1_selected_genes)]
    genes_in_both_idx1 <- which(names1 %in% genes_in_both)
    print(paste0("use ", length(genes_in_both_idx1), " genes for training model"))
    # done selecting genes--------------------------------------------------------
    
    # prepare ElasticNet training
    # matrices--------------------------------------------- prepare predictor matrix
    # containing both clustering classes
    predictor_S1 <- mixedpop1[genes_in_both_idx1, c(subpop1_train_indx, subremaining1_train_indx)]
    predictor_S1 <- assay(predictor_S1)
    print(paste0("use ", dim(predictor_S1)[1], " genes ", dim(predictor_S1)[2], " cells for testing model"))
    # generate categorical response
    
    # assign class labels############ rename the group of remaining clusters
    c_compareID <- unique(cluster_mixedpop1)
    c_compareID <- paste0(c_compareID[-which(c_compareID == c_selectID)], collapse = "_")
    print(paste0("rename remaining subpops to ", c_compareID))
    # change cluster number to character
    c_selectID <- as.character(c_selectID)
    # done assigning class labels####
    
    # set all y values to cluster select (character type) length = #cells
    y_cat = rep(c_selectID, ncol(predictor_S1))
    
    # replace values for cluster compare get indexes for cells in predictor_S1 belong
    # to 'remaining class'
    remainingClass_Indx_in_y <- which(colnames(predictor_S1) %in% colnames(mixedpop1[, 
        subremaining1_train_indx]))
    
    # change value of the cluster id for the remaining class
    y_cat[remainingClass_Indx_in_y] <- rep(c_compareID, length(remainingClass_Indx_in_y))
    print(paste0("there are ", length(remainingClass_Indx_in_y), " cells in class ", 
        c_compareID, " and ", length(y_cat) - length(remainingClass_Indx_in_y), " cells in class ", 
        c_selectID))
    # Done prepare ElasticNet training matrices------------------------------------
    
    # prepare training matrices------------------------------------------------
    
    # fitting with lda, also with cross validation
    dataset <- t(predictor_S1)  #(note predictor_S1 =t(gene_S1))
    # remove genes with no variation across cells and duplicated genes
    Zero_col <- which(colSums(dataset) == 0)
    duplicated_col <- which(duplicated(colnames(dataset)) == TRUE)
    if (length(c(Zero_col, duplicated_col)) != 0) {
        print(paste0("removing ", length(c(Zero_col, duplicated_col)), " genes with no variance"))
        dataset <- dataset[, -c(Zero_col, duplicated_col)]
    }
    
    # scaled and centered data before training standardizing data (centered and
    # scaled) - this step is not necessary in the function glmnet
    standardizing <- function(X) {
        X <- X - mean(X)
        X <- X/sd(X)
        return(X)
    }
    
    if (standardize == TRUE) {
        print("standardizing prediction/target dataset")
        dataset <- t(apply(dataset, 1, standardizing))  #dataset is transposed after standardised and scaled
    }
    dataset <- as.data.frame(dataset)
    dataset$Cluster_class <- as.character(as.vector(y_cat))
    # remove NA before prediction (if a gene has zero variance, it will have NA
    # values after standardization)
    dataset <- na.omit(dataset)
    # fitting with cross validation to find the best ElasticNet model
    print("performning elasticnet model training...")
    cvfit = cv.glmnet(as.matrix(dataset[, -which(colnames(dataset) == "Cluster_class")]), 
        as.vector(dataset$Cluster_class), family = "binomial", type.measure = "class")
    
    # fit LDA
    if (LDA_run == TRUE) {
        print("performning LDA model training...")
        trainControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
        fit.lda <- train(Cluster_class ~ ., preProcess = NULL, data = dataset, method = "lda", 
            metric = "Accuracy", trControl = trainControl, na.action = na.omit, standardize = FALSE)
    }
    
    # Done fitting the ElasticNet and LDA models-----------------------------------
    
    # Extract coefficient Beta for a gene for an optimized lambda value------------
    print("extracting deviance and best gene features...")
    cvfit_out <- as.matrix(coef(cvfit, s = cvfit$lambda.min))
    cvfit_out <- as.data.frame(cvfit_out)
    # Find number of genes with coefficient different to 0
    cvfit_out$name <- row.names(cvfit_out)
    sub_cvfit_out <- cvfit_out[cvfit_out$`1` != 0, ]
    # Extract deviance explained
    t_DE <- as.matrix(print(cvfit$glmnet.fit))
    dat_DE <- as.data.frame(t_DE)
    colnames(dat_DE) <- c("Dfd", "Deviance", "lambda")
    # Get the coordinate for lambda that produces minimum error
    dat_DE_Lambda_idx <- c()
    for (i in (1:nrow(dat_DE))) {
        if (dat_DE$lambda[i] == round(cvfit$lambda.min, nchar(dat_DE$lambda[i]) - 
            2)) {
            dat_DE_Lambda_idx <- c(dat_DE_Lambda_idx, i)
        }
    }
    if (length(dat_DE_Lambda_idx) > 0) {
        dat_DE <- dat_DE[1:dat_DE_Lambda_idx[1], ]
        print(paste0("lambda min is at location ", dat_DE_Lambda_idx[1]))
    } else {
        print("no lambda min found, please check output ...")
    }
    
    dat_DE_fm_DE <- dat_DE %>% group_by(dat_DE$Dfd) %>% summarise(Deviance = max(dat_DE$Deviance))
    dat_DE_fm_DE <- as.data.frame(dat_DE_fm_DE)
    dat_DE_fm_DE$DEgenes <- paste0("genes_cluster", c_selectID)
    remaining <- c("remaining DEgenes")
    dat_DE_fm_DE <- rbind(dat_DE_fm_DE, remaining)
    # Done extracting the coefficients---------------------------------------------
    
    # fit the model on the leave-out data to estimate accuracy---------------------
    
    # Keep all cells except for those used in the training set note that:
    # subpop1cluster_indx is for the total cells in the source subpop
    # subpop1_train_indx is for the subsampled cells as in
    # subpop1cluster_indx[sample(1:ncol(mixedpop1),subsampling , replace = F)]
    
    cluster_select_indx_Round2 <- subpop1cluster_indx[-which(subpop1cluster_indx %in% 
        subpop1_train_indx)]
    print(paste0("the leave-out cells in the source subpop is ", length(cluster_select_indx_Round2)))
    
    # subremaining1_train_indx is the index for remaining subpops used for training
    # subremaining1_indx is the total cells in the remaining subpops
    
    subremaining1_train_order <- which(subremaining1_indx %in% subremaining1_train_indx)  #to remove remaininng class already used for training 
    
    # Select cluster_compare_indx_Round2 as the remaining cells in class 2 (opposite
    # to class 1)
    if ((length(subremaining1_indx) - length(subremaining1_train_indx)) > subsampling) {
        cluster_compare_indx_Round2 <- sample(subremaining1_indx[-subremaining1_train_order], 
            subsampling, replace = FALSE)
        print(paste0("use ", subsampling, " target subpops cells for leave-out test set"))
    } else {
        cluster_compare_indx_Round2 <- subremaining1_indx[-subremaining1_train_order]  #use all of the class 2 cells after taken those for training 
        print(paste0("use ", length(cluster_compare_indx_Round2), " target subpop cells for leave-out test set"))
    }
    
    # the cells present in select (class 1) vs compare (remaining class, i.e.  class
    # 2) c(cluster_select_indx_Round2, cluster_compare_indx_Round2) use the leave-out
    # dataset in mixed pop1 to evaluate the model
    temp_S2 <- mixedpop1[, c(cluster_select_indx_Round2, cluster_compare_indx_Round2)]  #genes_in_trainset_idx_ordered
    temp_S2 <- assay(temp_S2)
    predictor_S2 <- t(temp_S2)
    
    # before prediction, check genes in cvfit and predictor_S2 are compatible Get
    # gene names
    gene_cvfit <- cvfit$glmnet.fit$beta@Dimnames[[1]]
    names <- colnames(predictor_S2)
    cvfitGenes_idx <- match(gene_cvfit, as.character(names))  #names are genes in target subpop
    if (length(cvfitGenes_idx) == 0) {
        print("No genes in the model present in the data...check your input")
    }
    
    # check how many genes in gene_cvfit but not in mixedpop2
    to_add <- which(is.na(cvfitGenes_idx) == TRUE)
    if (length(to_add) > 0) {
        print(paste0("add ", length(to_add), " random indexes for genes in model but not in target subpop, later to be replaced by 0"))
        to_add_idx <- sample(cvfitGenes_idx[-to_add], length(to_add), replace = FALSE)
        cvfitGenes_idx[to_add] <- to_add_idx  #replace NA indexes by random           
    }
    
    if (length(cvfitGenes_idx) > 0) {
        predictor_S2 <- predictor_S2[, cvfitGenes_idx]
    }
    # replace gene values for those in model but not in prediction dataset by NA
    if (length(to_add) > 0) {
        print("Replacing missing genes by NA...")
        predictor_S2[, to_add] <- NA
    }
    
    # standardizing the leave-out target and source subpops
    if (standardize == TRUE) {
        print("standardizing the leave-out target and source subpops...")
        predictor_S2 <- t(apply(predictor_S2, 1, standardizing))
    }
    # remove NA before prediction (if a gene has zero variance, it will have NA
    # values after standardization)
    predictor_S2 <- na.omit(predictor_S2)
    # Evaluation: predict ElasticNet Start prediction for estimating accuracy
    print("start ElasticNet prediction for estimating accuracy...")
    predict_clusters <- predict(cvfit, newx = predictor_S2, type = "class", s = cvfit$lambda.min)
    
    # Evaluation: predict LDA
    if (LDA_run == TRUE) {
        print(paste0("start LDA prediction for estimating accuracy for ", nrow(predictor_S2), 
            " cells and ", ncol(predictor_S2), " genes..."))
        lda_predict <- predict(fit.lda, predictor_S2)
    }
    # Done validation test to estimate accuracy------------------------------------
    
    # Estimate accuracy------------------------------------------------------------
    # Reference classification for comparing predicted results this cellNames_cluster
    # is to calculate the accuracy
    cellNames_cluster <- cbind(colnames(mixedpop1), cluster_mixedpop1)
    # Compare to original clusering classes to check for accuracy (just for the
    # source subpop), using the initial cellNames_cluster
    predictor_S2_name <- row.names(predict_clusters)
    predict_label <- predict_clusters[, 1]
    
    predict_index <- which(cellNames_cluster[, 1] %in% predictor_S2_name)
    
    original_cluster <- cellNames_cluster[predict_index, ]
    original_cluster <- original_cluster[order(original_cluster[, 2], decreasing = TRUE), 
        ]
    original_cluster <- as.data.frame(original_cluster)
    # for ElasticNet
    predict_cluster_dat <- as.data.frame(predict_label)
    predict_cluster_dat$cellnames <- row.names(predict_cluster_dat)
    # merge original and predicted datasets
    compare <- merge(original_cluster, predict_cluster_dat, by.x = "V1", by.y = "cellnames")
    # fing how many predictions are accurate
    cluster_select_predict <- subset(compare, (as.numeric(compare[, 2]) == c_selectID & 
        compare$predict_label == c_selectID) | (as.numeric(compare[, 2]) != c_selectID & 
        compare$predict_label != c_selectID))
    accurate <- dim(cluster_select_predict)[1]
    inaccurate <- dim(compare)[1] - dim(cluster_select_predict)[1]
    list_acc_inacc <- list(accurate, inaccurate)
    print(paste0("evaluation accuracy ElasticNet ", accurate/(accurate + inaccurate)))
    # for LDA
    if (LDA_run == TRUE) {
        predict_cluster_dat_LDA <- as.data.frame(lda_predict)
        predict_cluster_dat_LDA$cellnames <- row.names(predictor_S2)
        # merge original and predicted datasets
        compare <- merge(original_cluster, predict_cluster_dat_LDA, by.x = "V1", 
            by.y = "cellnames")
        # fing how many predictions are accurate
        cluster_select_predict_LDA <- subset(compare, (as.numeric(compare[, 2]) == 
            c_selectID & compare$lda_predict == c_selectID) | (as.numeric(compare[, 
            2]) != c_selectID & compare$lda_predict != c_selectID))
        accurate_LDA <- dim(cluster_select_predict_LDA)[1]
        inaccurate_LDA <- dim(compare)[1] - dim(cluster_select_predict_LDA)[1]
        list_acc_inacc <- list(accurate_LDA, inaccurate_LDA)
        print(paste0("evaluation accuracy LDA ", accurate_LDA/(accurate_LDA + inaccurate_LDA)))
    }
    
    # done estimate accuracy----------------------------------------------------
    
    # write the 5 lists into the output object
    
    listData$Accuracy[[out_idx]] <- list(list_acc_inacc)
    listData$ElasticNetGenes[[out_idx]] <- list(sub_cvfit_out)
    listData$Deviance[[out_idx]] <- list(dat_DE_fm_DE)
    listData$ElasticNetFit[[out_idx]] <- list(cvfit)
    if (LDA_run == TRUE) {
        listData$LDAFit[[out_idx]] <- list(fit.lda)
    } else {
        listData$LDAFit[[out_idx]] <- list(NA)
    }
    listData$predictor_S1[[out_idx]] <- list(t(dataset))
    
    return(listData)
}


#' Main prediction function applying the optimal ElasticNet and LDA models
#'
#' @description  Predict a new mixed population after training the model for a
#' subpopulation in the first mixed population.
#' All subpopulations in the new target mixed population will be predicted, 
#' where each targeted subpopulation will have a transition score from the 
#' orginal subpopulation to the new subpopulation.
#' @param listData a \code{list} object containing trained results for the
#' selected subpopulation in the first mixed population
#' @param mixedpop2 a \linkS4class{SingleCellExperiment} object from the target
#' mixed population of importance, e.g. differentially expressed genes that are
#' most significant
#' @param out_idx a number to specify index to write results into the list 
#' output. This is needed for running bootstrap.
#' @param cluster_mixedpop2 a vector of cluster assignment for mixedpop2
#' @param standardize a logical of whether to standardize the data 
#' @param LDA_run logical, if the LDA prediction is added to compare to ElasticNet, 
#' the LDA model needs to be trained from the training_scGPS before inputting to this 
#' prediction step
#' @param c_selectId a number to specify the trained cluster used for prediction 
#' @return a \code{list} with prediction results written in to the index
#' \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' c_selectID<-1
#' out_idx<-1
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' listData  <- training_scGPS(genes, 
#'     cluster_mixedpop1 = colData(mixedpop1)[, 1], mixedpop1 = mixedpop1, 
#'     mixedpop2 = mixedpop2, c_selectID, listData =list(), out_idx=out_idx)
#' listData  <- predicting_scGPS(listData =listData,  mixedpop2 = mixedpop2, 
#'     out_idx=out_idx, cluster_mixedpop2 = colData(mixedpop2)[, 1], c_selectID = c_selectID)
#'


predicting_scGPS <- function(listData = NULL, cluster_mixedpop2 = NULL, mixedpop2 = NULL, 
    out_idx = NULL, standardize = TRUE, LDA_run = FALSE, c_selectID = NULL) {
    # predictor_S1 is the dataset used for the training phase (already transposed)
    predictor_S1 <- listData$predictor_S1[[out_idx]][[1]]  #1 for extract matrix
    cvfit_best <- listData$ElasticNetFit[[out_idx]][[1]]
    fit.lda <- listData$LDAFit[[out_idx]][[1]]
    
    my.clusters <- cluster_mixedpop2
    ori_dat_2 <- assay(mixedpop2)
    
    # standardizing data (centered and scaled)
    standardizing <- function(X) {
        X <- X - mean(X)
        X <- X/sd(X)
        return(X)
    }
    
    if (standardize == TRUE) {
        print("standardizing target subpops before prediction...")
        ori_dat_2 <- t(apply(ori_dat_2, 1, standardizing))
    }
    
    list_predict_clusters_ElasticNet <- list()
    list_predict_clusters_LDA <- list()
    
    for (clust in unique(my.clusters)) {
        print(paste0("predicting from source to target subpop ", clust, "..."))
        c_selectID_2 <- clust
        cluster_select <- which(my.clusters == c_selectID_2)  #select cells
        print(paste0("number of cells in the target subpop ", clust, " is ", length(cluster_select)))
        # Get gene names from the trained model
        target_genes <- elementMetadata(mixedpop2)[, 1]
        
        # A function to match target dataset
        matching_genes <- function(target_dataset, model_genes, cluster_idx = cluster_select) {
            target_genes <- row.names(target_dataset)
            targetGenes_remove <- which(!(as.character(target_genes) %in% model_genes))  #genes in target but not in modelGenes
            print(paste0("Number of genes in the target data, but not in model genes is ", 
                length(targetGenes_remove)))
            target_dataset <- target_dataset[-targetGenes_remove, ]
            # update target genes
            target_genes <- row.names(target_dataset)
            # assumming the same gename formats between target and model
            modelGenes_idx <- match(model_genes, as.character(target_genes))  #which indexes of the target genes that match model genes
            print(paste0("Number of genes in the model present in the target data is ", 
                length(na.omit(modelGenes_idx))))
            if (length(na.omit(modelGenes_idx)) == 0) {
                print("No genes in the model present in the data...check your input")
            }
            # check how many genes in model_genes but not in mixedpop2 and add those genes
            # into the target
            to_add <- which(is.na(modelGenes_idx) == TRUE)
            print(paste0("There are ", length(to_add), " genes that are in the model, but not in target subpopulations"))
            # update target_dataset
            predictor_S2_temp <- target_dataset[na.omit(modelGenes_idx), cluster_select]
            
            if (length(to_add) > 0) {
                print(paste0("for dim conformation, add ", length(to_add), " random indexes for genes in model but not in target subpop, later to be replaced by 0"))
                # if the genes not in the target dataset, set the expression values to 0
                to_add_df <- as.data.frame(matrix(0, nrow = length(to_add), ncol = length(cluster_select)))  #columns = number of cells
                row.names(to_add_df) <- gsub("`", "", model_genes[to_add])  #format special names with '`MCM3AP-AS1`'
                print(paste0("genes to be added to target subpop are ", as.vector(row.names(to_add_df))))
                colnames(to_add_df) <- colnames(predictor_S2_temp)
                predictor_S2_temp <- rbind(predictor_S2_temp, to_add_df)
            }
            
            print(paste0("the prediction (target) subop has ", dim(predictor_S2_temp)[1], 
                " genes and ", dim(predictor_S2_temp)[2], " cells. The trained model has ", 
                length(model_genes), " genes"))
            print("first 10 genes in model ")
            print(head(model_genes, n = 10))
            print("first 10 genes in target ")
            print(head(row.names(predictor_S2_temp), n = 10))
            return(predictor_S2_temp)
        }
        
        # update predictor dataset
        predictor_S2_temp <- matching_genes(target_dataset = ori_dat_2, model_genes = cvfit_best$glmnet.fit$beta@Dimnames[[1]], 
            cluster_idx = cluster_select)
        # predict ElasticNet:
        print("running elasticNet classification...")
        predict_clusters_ElasticNet <- predict(cvfit_best, newx = t(predictor_S2_temp), 
            type = "class", s = cvfit_best$lambda.min)
        
        ElasticNet_result <- as.data.frame(table(predict_clusters_ElasticNet))  #convert table() to 2x2 dataframe, it will always have 2 variable names: $name, Freq
        report_ElasticNet <- paste0("ElasticNet for subpop", c_selectID_2, " in target mixedpop2")
        # report class probability
        if (ncol(ElasticNet_result) == 2) {
            ElasticNet_cluster_idx <- which(ElasticNet_result[, 1] == c_selectID)
            if (length(ElasticNet_cluster_idx) > 0) {
                predict_clusters_ElasticNet_proportion <- as.numeric(ElasticNet_result[ElasticNet_cluster_idx, 
                  2])/sum(as.numeric(ElasticNet_result[, 2])) * 100
                print(paste0("class probability prediction ElasticNet for target subpop ", 
                  clust, " is ", predict_clusters_ElasticNet_proportion))
                
            } else {
                predict_clusters_ElasticNet_proportion = 0
                print(paste0("class probability prediction ElasticNet for target subpop ", 
                  clust, " is 0"))
            }
            predict_clusters_ElasticNet <- list(report_ElasticNet, predict_clusters_ElasticNet_proportion)
        } else {
            print(paste0("ElasticNet for target subpop ", clust, " has no solution"))
            predict_clusters_ElasticNet <- list(report_ElasticNet, "NA")
        }
        
        # write to the next list level in the list if bootstrap is used
        list_predict_clusters_ElasticNet <- c(list_predict_clusters_ElasticNet, predict_clusters_ElasticNet)
        
        # done checking---------------------------------------------------------- predict
        # LDA:
        if (LDA_run == TRUE) 
            {
                # newdataset <- as.data.frame(t(ori_dat_2[,cluster_select])) for better LDA
                # conversion, the target data should not be standardised
                ori_dat_2 <- assay(mixedpop2)
                
                newdataset <- matching_genes(target_dataset = ori_dat_2, model_genes = fit.lda$finalModel$xNames, 
                  cluster_idx = cluster_select)
                
                print("running LDA classification...")
                predict_clusters_LDA <- predict(fit.lda, t(newdataset), na.action = na.omit)
                
                LDA_result <- as.data.frame(table(predict_clusters_LDA))  #convert
                # table() to 2x2 dataframe, it will always have 2 variable names: $name,
                report_LDA <- paste0("LDA for subpop ", c_selectID_2, " in target mixedpop2")
                if (ncol(LDA_result) == 2) {
                  LDA_cluster_idx <- which(LDA_result[, 1] == c_selectID)
                  predict_clusters_LDA <- as.numeric(LDA_result[LDA_cluster_idx, 
                    2])/sum(as.numeric(LDA_result[, 2])) * 100
                  print(paste0("class probability prediction LDA for target subpop ", 
                    clust, " is ", predict_clusters_LDA))
                  # write prediction output
                  predict_clusters_LDA <- list(report_LDA, predict_clusters_LDA)
                  
                } else {
                  
                  print(paste0("prediction LDA has no solution for target subpop ", 
                    clust))
                  # write prediction output
                  predict_clusters_LDA <- list(report_LDA, "NA")
                }
                
                list_predict_clusters_LDA <- c(list_predict_clusters_LDA, predict_clusters_LDA)
            }  # end LDA run  
    }  # end the loop through all subpops
    
    # write to the next list level in the list if bootstrap is used to write
    # prediction result
    listData$ElasticNetPredict[out_idx] <- list(list_predict_clusters_ElasticNet)
    if (LDA_run == TRUE) {
        listData$LDAPredict[out_idx] <- list(list_predict_clusters_LDA)
    } else {
        listData$LDAPredict[out_idx] <- list(NA)
    }
    
    return(listData)
}


#' BootStrap runs for both scGPS training and prediction
#'
#' @description  ElasticNet and LDA prediction for each of all the 
#' subpopulations in the new mixed population after training the model for a 
#' subpopulation in the first mixed population. The number of bootstraps to be 
#' run can be specified.
#' @seealso \code{\link{bootstrap_scGPS_parallel}} for parallel options
#' @param listData  a \code{list} object, which contains trained results for the
#' first mixed population
#' @param mixedpop1 a \linkS4class{SingleCellExperiment} object from a mixed 
#' population for training
#' @param mixedpop2 a \linkS4class{SingleCellExperiment} object from a target 
#' mixed population for prediction
#' @param cluster_mixedpop1 a vector of cluster assignment for mixedpop1 
#' @param cluster_mixedpop2 a vector of cluster assignment for mixedpop2 
#' @param c_selectID the root cluster in mixedpop1 to becompared to clusters in
#' mixedpop2 
#' @param genes a gene list to build the model
#' @param nboots a number specifying how many bootstraps to be run
#' @param trainset_ratio a number specifying the proportion of cells to be part 
#' of the training subpopulation
#' @param LDA_run logical, if the LDA prediction is added to compare to ElasticNet
#' @return a \code{list} with prediction results written in to the index 
#' \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' cluster_mixedpop1 <- colData(mixedpop1)[,1]
#' cluster_mixedpop2 <- colData(mixedpop2)[,1]
#' c_selectID <- 2
#' test <- bootstrap_scGPS(nboots = 2, mixedpop1 = mixedpop1, 
#'     mixedpop2 = mixedpop2, genes=genes, listData =list(), 
#'     cluster_mixedpop1 = cluster_mixedpop1, 
#'     cluster_mixedpop2 = cluster_mixedpop2, c_selectID = c_selectID)
#' names(test)
#' test$ElasticNetPredict
#' test$LDAPredict

bootstrap_scGPS <- function(nboots = 1, genes = genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, 
    c_selectID = NULL, listData = list(), cluster_mixedpop1 = NULL, cluster_mixedpop2 = NULL, 
    trainset_ratio = 0.5, LDA_run = TRUE) {
    
    for (out_idx in 1:nboots) {
        listData <- training_scGPS(genes = genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, 
            trainset_ratio = trainset_ratio, c_selectID, listData = listData, out_idx = out_idx, 
            cluster_mixedpop1 = cluster_mixedpop1, standardize = TRUE, LDA_run = LDA_run)
        print(paste0("done training for bootstrap ", out_idx, ", moving to prediction..."))
        listData <- predicting_scGPS(listData = listData, mixedpop2 = mixedpop2, 
            out_idx = out_idx, standardize = TRUE, cluster_mixedpop2 = cluster_mixedpop2, 
            LDA_run = LDA_run, c_selectID = c_selectID)
    }
    return(listData)
}


#' BootStrap runs for both scGPS training and prediction with parallel option
#'
#' @description  same as bootstrap_scGPS, but with an multicore option
#' @param listData  a \code{list} object, which contains trained results for the
#' first mixed population
#' @param mixedpop1 a \linkS4class{SingleCellExperiment} object from a mixed 
#' population for training
#' @param mixedpop2 a \linkS4class{SingleCellExperiment} object from a target 
#' mixed population for prediction
#' @param cluster_mixedpop1 a vector of cluster assignment for mixedpop1 
#' @param cluster_mixedpop2 a vector of cluster assignment for mixedpop2 
#' @param genes a gene list to build the model
#' @param nboots a number specifying how many bootstraps to be run
#' @param ncores a number specifying how many cpus to be used for running
#' @param c_selectID the root cluster in mixedpop1 to becompared to clusters in
#' mixedpop2 
#' @return a \code{list} with prediction results written in to the index 
#' \code{out_idx}
#' @export
#' @author Quan Nguyen, 2017-11-25
#' @examples
#' day2 <- sample1
#' mixedpop1 <-NewscGPS(ExpressionMatrix = day2$dat2_counts, 
#'     GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#' day5 <- sample2
#' mixedpop2 <-NewscGPS(ExpressionMatrix = day5$dat5_counts, 
#'     GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
#' genes <-GeneList
#' genes <-genes$Merged_unique
#' #prl_boots <- bootstrap_scGPS_parallel(ncores = 4, nboots = 2, genes=genes,
#' #    mixedpop1 = mixedpop2, mixedpop2 = mixedpop2,  c_selectID=1,
#' #    listData =list())
#' #prl_boots[[1]]$ElasticNetPredict
#' #prl_boots[[1]]$LDAPredict
#'

bootstrap_scGPS_parallel <- function(ncores = 4, nboots = 1, genes = genes, mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, c_selectID, listData = list(), cluster_mixedpop1 = NULL, 
    cluster_mixedpop2 = NULL) {
    
    bootstrap_single <- function(genes = genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, 
        c_selectID = c_selectID, out_idx = 1, listData = list(), cluster_mixedpop1 = NULL, 
        cluster_mixedpop2 = NULL) {
        listData <- training_scGPS(genes = genes, mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, 
            c_selectID, listData = listData, out_idx = 1, cluster_mixedpop1 = cluster_mixedpop1)
        listData <- predicting_scGPS(listData = listData, mixedpop2 = mixedpop2, 
            out_idx = 1, standardize = TRUE, cluster_mixedpop2 = cluster_mixedpop2)
        return(listData)
    }
    
    BiocParallel::register(BiocParallel::MulticoreParam(workers = ncores, progressbar = TRUE))
    
    
    listData <- BiocParallel::bplapply(1:nboots, bootstrap_single, genes = genes, 
        mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID = c_selectID, listData = list())
    return(listData)
}





