#' Main training
#'
#' @description  Training LASSO and LDA
#' @param mixedpop1 is a \linkS4class{SingleCellExperiment} object from the train mixed population
#' @param mixedpop2 is a \linkS4class{SingleCellExperiment} object from the target mixed population
#' @param genes a vector of gene names (for LASSO shrinkage); gene symbols must be
#' in the same format with gene names in subpop2; #' genes are listed by the order
#' of importance, e.g. differentially expressed genes that are most significant
#' @param c_selectID a selected number to specify which subpop to be used for training
#' (original subpop)
#' @return a list containing: predict_clusters, predictor_S2, sub_cvfit_out, dat_DE_fm_DE,
#' cvfit, predictor_S1,fit.lda
#' @export
#' @author QN
#' @Example
#' \dontrun { #to be updated
#' prepare_2_SubPop(mat1, mat2, geneNames)
#'
#' prepare_2_SubPop(SubPop1 = NULL,
#'                  SubPop2 = NULL,
#'                  genes = "")
#' }
#'

training_scGPS <-function(genes, mixedpop1 = NULL,
                         mixedpop2 = NULL, c_selectID) {
  # subsammpling-----------------------------------------------------------------
  # taking a subsampling size of 50% of the cluster_select out for training
  subpop1cluster_indx <- which(colData(mixedpop1)[,1] == c_selectID)
  subremaining1_indx <- which(colData(mixedpop1)[,1] != c_selectID)

  subsampling <- round(length(subpop1_cluster_indx)/2)
  subpop1_train_indx <-   subpop1cluster_indx[sample(1:ncol(mixedpop1),
                                                     subsampling , replace = F)]
  # check if there is a very big cluster present in the dataset
  # C/2 = SubSampling >length(total - C, which is cluster_compare)
  if (ncol(mixedpop1) > 2*subsampling) {
    subremaining1_train_indx <- subremaining1_indx[sample(1:length(subremaining1_indx),
                                                          subsampling, replace = F)]
  } else {
    subremaining_train_indx <- subremaining1_indx
  }
  # done subsampling------------------------------------------------------------

  # making sure the gene list contains fewer than 500 genes---------------------
  # select top 500, the DE_result is an already sorted list
  if(length(genes) >500){genes <- genes[1:500]}
  #select genes
  names1 <-names(mixedpop1)
  subpop1_selected_genes <- names1[which(names1 %in% genes)]
  names2 <-names(mixedpop2)
  genes_in_both <-names2[which(names2 %in% subpop1_selected_genes)]
  # done selecting genes--------------------------------------------------------

  #prepare LASSO training matrices----------------------------------------------

  #prepare predictor matrix containing both clutering classes
  predictor_S1 <-cbind(mixedpop1[genes_in_both,subpop1_train_indx],
                       mixedpop1[genes_in_both,subremaining_train_indx])

  #generate categorical response

  #assign class labels############
  c_compareID <- 1:length(unique(colData(t)[,1]))
  c_compareID <- paste0(c_compareID[-which(c_compareID==c_selectID)], collapse = "")
  # change cluster number to character
  c_selectID <- as.character(c_selectID)
  #done assigning class labels####

  #set all y values to cluster select (character type) length = #cells
  y_cat = rep(c_selectID,ncol(predictor_S1))

  #replace values for cluster compare
  #get indexes for cells in predictor_S1 belong to "remaining class"
  remainingClass_Indx_in_y <-which(colnames(predictor_S1) %in%
                                     colnames(mixedpop1[,subremaining_train_indx]))
  #change value of the cluster id
  y_cat[remainingClass_Indx_in_y] <-rep(c_compareID, length(remainingClass_Indx_in_y))
  #Done prepare LASSO training matrices-----------------------------------------

  #prepare LDA training matrices------------------------------------------------
  #fitting with lda, also with cross validation
  dataset <- t(predictor_S1) #(note predictor_S1 =t(gene_S1))
  dataset <- as.data.frame(dataset)

  #remove genes with no variation across cells
  Zero_col<-which(colSums(dataset) == 0)
  duplicated_col <-which(duplicated(colnames(dataset)) == TRUE)
  if(length(c(Zero_col, duplicated_col)) !=0){dataset <-dataset[,-c(Zero_col, duplicated_col)]}
  dataset$Cluster_class <- as.character(y_cat)
  #Done LDA training matrices---------------------------------------------------

  #Fit the LASSO and LDA models-------------------------------------------------

  #fitting with cross validation to find the best LASSO model
  cvfit = cv.glmnet(t(predictor_S1), y_cat, family = "binomial", type.measure = "class")

  #fit LDA
  trainControl <- trainControl(method="repeatedcv", number=10, repeats=3)
  metric <- "Accuracy"
  fit.lda <- train(Cluster_class ~., data=dataset, method="lda", metric=metric,
                   trControl=trainControl, na.action=na.omit)

  #Done fitting the LASSO and LDA models----------------------------------------

  #Extract coefficient Beta for a gene for an optimized lambda value------------

  cvfit_out <-as.matrix(coef(cvfit, s = cvfit$lambda.min))
  cvfit_out <-as.data.frame(cvfit_out)
  #find number of genes with coefficient different to 0
  cvfit_out$name <-row.names(cvfit_out)
  sub_cvfit_out <-cvfit_out[cvfit_out$`1` != 0,]
  #extract deviance explained
  t_DE <- as.matrix(print(cvfit$glmnet.fit))
  dat_DE <-as.data.frame(t_DE)
  colnames(dat_DE) <-c('Dfd', 'Deviance', 'lambda')
  #to get the coordinate for lambda that produces minimum error
  dat_DE_Lambda_idx <- which(round(dat_DE$lambda,digit=3) == round(cvfit$lambda.min,digits=3))
  dat_DE <-dat_DE[1:dat_DE_Lambda_idx[1],]
  require(dplyr)
  dat_DE %>% group_by(Dfd) %>% summarise(Deviance = max(Deviance)) -> dat_DE_fm_DE
  dat_DE_fm_DE <-as.data.frame(dat_DE_fm_DE)
  dat_DE_fm_DE$DEgenes <-paste0('DEgenes_C',c_selectID,'_day_', dayID)
  remaining <-c('remaining', 1, 'DEgenes')
  dat_DE_fm_DE <-rbind(dat_DE_fm_DE, remaining)
  #Done extracting the coefficients---------------------------------------------

  #Cross validation test to estimate accuracy-----------------------------------


  # reference classification for comparing predicted results
  #cellNames_subpop <- cbind(colnames(subpop1), colnames(subpop2)) #note: may edit
  # this to make 2 generic classes 1 and 2

  #keep all cells except for those used in the training set
  #recall that: subpop1_train_indx = subpop1cluster_indx[sample(1:ncol(mixedpop1),subsampling , replace = F)]

  cluster_select_indx_Round2 <- subpop1cluster_indx[- which(subpop1cluster_indx %in% subpop1_train_indx)]

  subremaining1_train_order <-  which(subremaining1_indx %in% subremaining1_train_indx)

  #check the SubSampling for the target population
  if (length(subremaining1_indx) - length(subremaining1_train_indx) > SubSampling) {
    cluster_compare_indx_Round2 <- sample(subremaining1_indx[-subremaining1_train_order], SubSampling, replace=F )
  } else {
    cluster_compare_indx_Round2 <- subremaining1_indx[-subremaining1_train_order] #keep everything in the predicted dat
  }

  temp_S2 <-SubPop1[genes_in_both, c(cluster_select_indx_Round2  , cluster_compare_indx_Round2 )]
  predictor_S2 <-t(temp_S2)

  #start prediction for estimating accuracy
  predict_clusters<-predict(cvfit, newx = predictor_S2,  type = "class", s = cvfit$lambda.min)
  #Done cross validation test to estimate accuracy------------------------------

  #to return the output as 5 lists
  return(list(predict_clusters, predictor_S2, sub_cvfit_out, dat_DE_fm_DE, cvfit,
              predictor_S1,fit.lda ))
}






