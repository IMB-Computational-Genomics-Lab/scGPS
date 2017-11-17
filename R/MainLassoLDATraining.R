#' Main training
#'
#' @description  Training LASSO and LDA
#' @param cluster_select_indx_S1 is the index from random sampling
#' genes are listed by the order of importance
#' @return a list containing: predict_clusters, predictor_S2, sub_cvfit_out, dat_DE_fm_DE, cvfit, predictor_S1,fit.lda
#' @author QN
#' @Example
#' prepare_2_subpop(mat1, mat2, geneNames)
#' \dontrun{
#'
#' prepare_2_subpop(subpop1 = NULL,
#'                  subpop2 = NULL,
#'                  genes = "")
#' }
#'

Lit_New_Lasso <-function(cluster_select_indx_S1=NULL) {
  #cluster_select_indx_S1 <-sample(cluster_select, subsampling, replace=F)
  #check if there is a very big cluster present in the dataset  C/2 = subsampling >length(total-C=cluster_compare)
  if (length(cluster_compare) > subsampling ) {
    cluster_compare_indx_S1 <- sample(cluster_compare, subsampling , replace=F)
  } else {
    cluster_compare_indx_S1 <- cluster_compare
  }

  c_compareID = 1:length(unique(my.clusters))
  c_compareID <- paste0(c_compareID[-which(c_compareID==c_selectID)], collapse="")
  M_or_DE_idx=DE_idx


  #prepare predictor matrix containing both clutering classes
  predictor_S1 <-ori_dat[M_or_DE_idx, c(cluster_select_indx_S1 , cluster_compare_indx_S1)]
  #generate categorical response
  #set all values to cluster select (character type)
  y_cat = rep(c_selectID,length(predictor_S1[1,]))
  #replace values for cluster compare
  #first get a new matrix
  ori_compare <-ori_dat[,cluster_compare]
  #get indexes for cells in predictor_S1 belong to cluster_compare class
  Sub_clustercompare_Indx_S1 <-which(colnames(predictor_S1) %in% colnames(ori_compare))
  #change value of the cluster id
  y_cat[Sub_clustercompare_Indx_S1] <-rep(c_compareID, length(Sub_clustercompare_Indx_S1))


  #fitting with cross validation to find the best LASSO model
  cvfit = cv.glmnet(t(predictor_S1), y_cat, family = "binomial", type.measure = "class")

  #fitting with lda, also with cross validation
  dataset <- t(predictor_S1) #(note predictor_S1 =t(gene_S1))
  dataset <- as.data.frame(dataset)

  Zero_col<-which(colSums(dataset)==0)
  duplicated_col <-which(duplicated(colnames(dataset))==TRUE)
  if(length(c(Zero_col, duplicated_col)) !=0){dataset <-dataset[,-c(Zero_col, duplicated_col)]}
  dataset$Cluster_class <- as.character(y_cat)

  trainControl <- trainControl(method="repeatedcv", number=10, repeats=3)
  metric <- "Accuracy"
  fit.lda <- train(Cluster_class ~., data=dataset, method="lda", metric=metric, trControl=trainControl, na.action=na.omit)

  #to extract coefficient Beta for a gene for an optimized lambda value
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

  #preparing for validation test to estimate accuracy
  #keep all cells except for those used in the training set
  cluster_select_indx_S2 <- cluster_select[-cluster_select_indx_S1]

  #check the subsampling for the target population
  if (length(cluster_compare[-cluster_compare_indx_S1]) > subsampling) {
    cluster_compare_indx_S2 <- sample(cluster_compare[-cluster_compare_indx_S1], subsampling, replace=F )
  } else {
    cluster_compare_indx_S2 <- cluster_compare #keep everything in the predicted dat
  }

  genes_S2 <-ori_dat[M_or_DE_idx, c(cluster_select_indx_S2 , cluster_compare_indx_S2)]
  predictor_S2 <-t(genes_S2)

  #start prediction for estimating accuracy
  predict_clusters<-predict(cvfit, newx = predictor_S2,  type = "class", s = cvfit$lambda.min)

  #to return the output as 5 lists
  return(list(predict_clusters, predictor_S2, sub_cvfit_out, dat_DE_fm_DE, cvfit, predictor_S1,fit.lda ))
}
