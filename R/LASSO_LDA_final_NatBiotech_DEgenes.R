#this script does generates all needed lists for 1 selected clusters!
#it compares within day and between day clusters
#it has both LASSO and LDA models
#it calculate accuracy and LASSO gene lists

library(glmnet)
library(caret)

argv <- commandArgs(TRUE)

c_selectID=argv[1]
dayID=argv[2]
dayID2=argv[3]

#c_selectID=1; dayID=2; dayID2=5

ori_dat <-readRDS(paste0("Exprs_DCVLnorm_unlog_minus1_pos_Day", dayID, ".RDS"))
my.clusters <-readRDS(paste0("my.clusters_0.45_day",dayID, ".RDS"))

#making sure only genes appearing in the two datasets are used to build cvfit
ori_dat_2 <-readRDS(paste0("../Exprs_DCVLnorm_unlog_minus1_pos_Day", dayID2, ".RDS"))
names <-rownames(ori_dat)
names <-gsub("_.*", '', names)
rownames(ori_dat) <-names

#get cluster IDs
n_clusters <-length(unique(my.clusters))
#extract clusters

cluster_select <- which(my.clusters==as.numeric(c_selectID))

cluster_compare <-  which(my.clusters !=as.numeric(c_selectID))

#this cellNames_cluster is to calculate the accuracy
cellNames_cluster <- cbind(colnames(ori_dat), my.clusters)

DE_result <- read.table(paste0('DEseq_Cluster', c_selectID, '_vs_OtherClusters_Day', dayID,'.txt_filtered_pAdjusted_sorted.txt'), header=T)

DEgenes <-DE_result$id

DEgenes <-gsub("_.*", '', DEgenes)

#select top 500, the DE_result is an already sorted list
if(length(DEgenes) >500){DEgenes <- DEgenes[1:500]}

#making sure only genes shared with oridat2 are kept
names <-gsub("_.*", '', rownames(ori_dat_2))
rownames(ori_dat_2) <- names
DEgenes_filtered <- DEgenes[which(DEgenes %in% names)]
DE_idx <-which(rownames(ori_dat) %in% DEgenes_filtered)

#taking a subsampling size of 50% of the cluster_select out for training
subsampling =round(length(cluster_select)/2)

#---------------------------------------------
#after this step we have: subsampling, cluster_select, cluster_compare, DE_idx, ori_dat, ori_dat2
#---------------------------------------------

##########################################################################################

#this function build cvfit model and evaluate prediction accuracy
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

#this function takes cvfit as the input and predict cells in another day or the same day
#this function predict the percent cells in the next cluster

predict_day <-function(dayID_same_or_nextday, predictor_S1, cvfit ){
	dayIDtarget <- dayID_same_or_nextday
	my.clusters <-readRDS(paste0("my.clusters_0.45_day",dayIDtarget, ".RDS"))
	ori_dat_2 <- readRDS(paste0("../Exprs_DCVLnorm_unlog_minus1_pos_Day", dayIDtarget, ".RDS"))
  names <-rownames(ori_dat_2)
  names <-gsub("_.*", '', names)
  rownames(ori_dat_2) <-names

  list_predict_clusters_LASSO <- list()
  list_predict_clusters_LDA <- list()

  for (clust in unique(my.clusters)){
		c_selectID_2 <-clust
		cluster_select <- which(my.clusters == as.numeric(c_selectID_2)) #select cells
    dataset <- predictor_S1 #predictor_S1 is the dataset used for the training phase (already transposed)
		gene_cvfit <- colnames(dataset)
		cvfitGenes_idx <-which(rownames(ori_dat_2) %in% gene_cvfit) #gene idx

		#predictor_S2_temp <-ori_dat_2[c(cvfitGenes_idx),cluster_select]

    #for getting genes
		cvfit_best <- cvfit
		gene_cvfit <- cvfit_best$glmnet.fit$beta@Dimnames[[1]]
		cvfitGenes_idx <-which(rownames(ori_dat_2) %in% gene_cvfit)

		to_add <-length(gene_cvfit) - length(cvfitGenes_idx)
		to_add_idx <- c()
		if(to_add >0) {
			to_add_idx <- sample(cvfitGenes_idx, to_add, replace=F)} else if (to_add <0) {
			cvfitGenes_idx <-sample(cvfitGenes_idx, length(gene_cvfit),replace=F)}

		predictor_S2_temp <-ori_dat_2[c(to_add_idx, cvfitGenes_idx),cluster_select]

		#predict LASSO:
		predict_clusters_LASSO <- predict(cvfit_best, newx = t(predictor_S2_temp),  type = "class", s = cvfit_best$lambda.min)
	  LASSO_result <- as.data.frame(table(predict_clusters_LASSO )) #convert table() to 2x2 dataframe, it will always have 2 variable names: $name, Freq
	  LASSO_cluster_idx <- which(LASSO_result[,1]== c_selectID)
	  predict_clusters_LASSO <-as.numeric(LASSO_result[LASSO_cluster_idx,2])/sum(as.numeric(LASSO_result[,2]))*100
	 #print out the prediction output
	  report_LASSO <- paste0('LASSO From cluster', c_selectID, ' Day ', dayID, ' to predict cluster ',c_selectID_2, 'Day ', dayIDtarget)
	  predict_clusters_LASSO <- list(predict_clusters_LASSO, report_LASSO)
    list_predict_clusters_LASSO <-c(list_predict_clusters_LASSO, predict_clusters_LASSO )
		#cluster_select_predict <- subset(compare, (as.numeric(compare$my.clusters) == c_selectID &  compare$`1`==c_selectID) | (compare$my.clusters!=c_selectID &  compare$`1`!=c_selectID))

		#predict LDA:
		newdataset <-t(predictor_S2_temp)
		newdataset <-as.data.frame(newdataset)

	  predict_clusters_LDA <- predict(fit.lda, newdataset)
	  LDA_result <- as.data.frame(table(predict_clusters_LDA )) #convert table() to 2x2 dataframe, it will always have 2 variable names: $name, Freq
	  LDA_cluster_idx <- which(LDA_result[,1] == c_selectID)
	  predict_clusters_LDA <-as.numeric(LDA_result[LDA_cluster_idx,2])/sum(as.numeric(LDA_result[,2]))*100

	  #write prediction output
	  report_LDA <- paste0('LDA From cluster', c_selectID, ' Day ', dayID, ' to predict cluster ',c_selectID_2, 'Day ', dayIDtarget)
	  predict_clusters_LDA <- list(predict_clusters_LDA, report_LDA)
    list_predict_clusters_LDA <-c(list_predict_clusters_LDA, predict_clusters_LDA)
		}
	return(list(predict_percent_LASSO =	list_predict_clusters_LASSO, predict_clusters_LDA = list_predict_clusters_LDA))
}

#Lists to be saved. Note that only list_predict_sameDay and list_predict_DifferentDay are day-specific; all others are shared for the same models LASSO or LDA
list_SigGenes = list() #for extracting genes post LASSO
list_Deviance = list() #for plotting deviance explained
list_cvFit = list() #for saving lasso model
list_lda = list() #for saving lda model

list_acc_inacc = list() #for plotting accuracy of a model
list_predict_sameDay = list() #for percent predicted by LASSO
list_predict_DifferentDay = list() #for percent predicted by LDA

#Run the loop below and save one big object containing all 7 lists above as the output
for (i in 1:100){

#c_compareID = 1:length(unique(my.clusters)) #not using this?
#c_compareID <- paste0(c_compareID[-which(c_compareID==c_selectID)], collapse="") #not using this?
#Call the big function here
#to randomise the input first

        cluster_select_indx_S1 <-sample(cluster_select, subsampling, replace=F)
        predict_marker<- Lit_New_Lasso(cluster_select_indx_S1 = cluster_select_indx_S1)

        predictor_S2_name <-predict_marker[[2]]
        predict_label <-predict_marker[[1]]

        #from here compare to original clusering classes to check for accuracy (just for the training set)
        #it uses the initial cellNames_cluster
        predict_index <-which(cellNames_cluster[,1] %in% row.names(predictor_S2_name))
        original_cluster <- cellNames_cluster[predict_index,]
        original_cluster <-original_cluster[order(original_cluster[,2], decreasing = T),]
        original_cluster <-as.data.frame(original_cluster)
        predict_clusters <-as.data.frame(predict_label)
        predict_clusters$cellnames <-row.names(predict_clusters)

        compare <-merge(original_cluster, predict_clusters, by.x='V1', by.y='cellnames')
        cluster_select_predict <- subset(compare, (as.numeric(compare$my.clusters) == c_selectID &  compare$`1`==c_selectID) | (compare$my.clusters!=c_selectID &  compare$`1`!=c_selectID)) #this is the accurate prediction
        accurate <-dim(cluster_select_predict)[1]
        inaccurate <-dim(compare)[1] - dim(cluster_select_predict)[1]

        list_acc_inacc[[i]] <- list(accurate, inaccurate)
        list_SigGenes[[i]] <-predict_marker[[3]]
        list_Deviance[[i]] <-predict_marker[[4]]
        list_cvFit[[i]]  <-predict_marker[[5]]
        list_lda[[i]] = predict_marker[[7]]

        #More new prediction
        cvfit=predict_marker[[5]]
        fit.lda =predict_marker[[7]]
        #predict lasso
        list_predict_sameDay[[i]] <- predict_day(dayID_same_or_nextday = dayID, predictor_S1 = predict_marker[[6]], cvfit= predict_marker[[5]])
        #predict LDA
        list_predict_DifferentDay[[i]] <- predict_day(dayID_same_or_nextday = dayID2, predictor_S1 = predict_marker[[6]], cvfit= predict_marker[[5]])

        }
###################################################################################################
list_all <-list(list_acc_inacc =list_acc_inacc,  list_SigGenes = list_SigGenes, list_Deviance = list_Deviance,  list_cvFit = list_cvFit,  list_predict_sameDay = list_predict_sameDay, list_predict_DifferentDay = list_predict_DifferentDay)
saveRDS(list_all, paste0('LASSO_LDA_DEgenes_cluster', c_selectID,'fromDay',dayID,'_toDay_', dayID2, '.RDS'))





