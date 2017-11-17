# ----------------------------------------------
#Summary LASSO QC results 
# ----------------------------------------------

##1: Summary Deviance Table 
##Note: no data for day0 C4 and day 2 C2  (For DE genes); no data for day0 C3 and day0 C4, 
##Note: for day 30 models, use day1=30, day2=30
#loading the list of 300 data-frames from 100 times bootstrap 
library(glmnet)
library(dplyr)
day1=0
day2=2
#dayID2=30
c_selectID_1=3

#LSOLDA_dat <-readRDS(paste0('LASSO_LDA_DEgenes_cluster',c_selectID_1,'fromDay',day1,'_toDay_', day2,'.RDS'))
#for Knownmarkers
LSOLDA_dat <-readRDS(paste0('LASSO_LDA_KnownMarker_Cluster',c_selectID_1,'fromDay',day1,'_toDay_', day2,'.RDS'))
#special case 
#LSOLDA_dat <-readRDS("LASSO_LDA_Day0_NatBiotech_Marker_SpecialCase_C3D0.PBS")

###STOP here to check the file is read before moving on############

names(LSOLDA_dat)

acc_inacc <- LSOLDA_dat$list_acc_inacc
deviDat <- LSOLDA_dat$list_Deviance
cvfit100 <-LSOLDA_dat$list_cvFit
genesSig <-LSOLDA_dat$list_SigGenes

#Save percent accuracy and deviance (deviance of the iteration with the highest deviance) in one object 
pcAcc <-as.vector(unlist(lapply(acc_inacc, function(x){x[[1]]/(x[[1]] +x[[2]])*100}))) 
maxAccurate_idx <-which(pcAcc == max(pcAcc,na.rm = T)) 

pcAcc[[maxAccurate_idx]] 

deviDat_max <-deviDat[[maxAccurate_idx]]
genesSig[[maxAccurate_idx]]
deviDat_max 


QC_all <-list("percentAcc" =pcAcc, "Deviance" = deviDat_max, "LASSO_genes" = genesSig[[maxAccurate_idx]])
str(QC_all)
#paste0("QCall_DEgenes_forCluster",c_selectID_1, "_fromDay_",day1,"_toDay_",day2,".RDS")

paste0("QCall_KnownMarkers_forCluster",c_selectID_1, "_fromDay_",day1,"_toDay_",day2,".RDS")

#change this to _KnownMarkers_
saveRDS(QC_all, paste0("QCall_KnownMarkers_forCluster",c_selectID_1, "_fromDay_",day1,"_toDay_",day2,".RDS"))

#saveRDS(QC_all, paste0("QCall_DEgenes_forCluster",c_selectID_1, "_fromDay_",day1,"_toDay_",day2,".RDS"))


# ----------------------------------------------
##Extract LASSO LDA prediction results (no need to rerun the above code)
##Note: no data for day0 C4 and day 2 C2  (For DE genes); no data for day0 C3 and day0 C4, day5 C4
##Note: for day 30 models, use day1=30, day2=30
# ----------------------------------------------

#A Function to find good prediction only during the boostrap
positive_pred <-c()
for (i in 1:100){
	for(j in c(1,3,5,7)){
		if(is.na(pred_sameDay[i][[1]][1][[1]][j] >0)){
			print(i)
			positive_pred <-c(positive_pred,i)

	}

} 
}

pred_sameDay <-pred_sameDay[-positive_pred] 
#Done function 
#-

#############################Extract LASSO LDA prediction results for the same day################################
library(glmnet)
library(dplyr)
day1=15
day2=30
#dayID2=30
c_selectID_1=2

LSOLDA_dat <-readRDS(paste0('LASSO_LDA_DEgenes_cluster',c_selectID_1,'fromDay',day1,'_toDay_', day2,'.RDS'))
#for Knownmarkers
#LSOLDA_dat <-readRDS(paste0('LASSO_LDA_KnownMarker_Cluster',c_selectID_1,'fromDay',day1,'_toDay_', day2,'.RDS'))
#special case 
#LSOLDA_dat <-readRDS("LASSO_LDA_Day0_NatBiotech_Marker_SpecialCase_C3D0.PBS")

###STOP here to check the file is read before moving on############

names(LSOLDA_dat)

pred_sameDay <- LSOLDA_dat$list_predict_sameDay 
##run the clean up function if needed##


pred_sameDay_tranformed <- as.vector(unlist(pred_sameDay))

toremove <- grep("cluster", pred_sameDay_tranformed)

percent_sameDay <- pred_sameDay_tranformed[-toremove]

length(percent_sameDay)
pred_sameDay_tranformed[1:20]
#the number here needed to be changed (to the number of predictions without NULL values)

pred_sameDay_mtrx <-matrix(percent_sameDay, nrow=4, ncol=100, byrow=F)

#the number here needed to be changed 
row_names <- pred_sameDay_tranformed[toremove[c(1:4)]]

row_names <-gsub(" ","_", row_names )

pred_sameDay_mtrx <-as.data.frame(pred_sameDay_mtrx)

pred_sameDay_mtrx$names <-row_names

pred_sameDay_mtrx

paste0('LASSO_LDA_DEgenes_PREDICTION_Cluster',c_selectID_1,'withinDay',day1, '.txt')

write.table(pred_sameDay_mtrx,paste0('LASSO_LDA_DEgenes_PREDICTION_Cluster',c_selectID_1,'withinDay',day1, '.txt'), sep="\t", col.names=F, row.names=F, quote=F)

dat <-read.table(paste0('LASSO_LDA_DEgenes_PREDICTION_Cluster',c_selectID_1,'withinDay',day1, '.txt'))

################################To predict different days##########################
#-------------------------------
#note: pred_DifferentDay_top <-pred_DifferentDay[1:30]

#A Function to find good prediction only during the boostrap
positive_pred <-c()
for (i in 1:100){
	for(j in c(1,3,5,7)){
		if(is.na(pred_DifferentDay[i][[1]][1][[1]][j] >0)){
			print(i)
			positive_pred <-c(positive_pred,i)

	}

} 
}

pred_DifferentDay <-pred_DifferentDay[-positive_pred] 
#Done function 
#-------------------------------


pred_DifferentDay <- LSOLDA_dat$list_predict_DifferentDay

#for cluster 3 day 2 use the first 30 elements of the list 


pred_DifferentDay_tranformed <- as.vector(unlist(pred_DifferentDay))

toremove <- grep("cluster", pred_DifferentDay_tranformed)

percent_DifferentDay <- pred_DifferentDay_tranformed[-toremove]

length(percent_DifferentDay)
pred_DifferentDay_tranformed[1:20] 

#the number here needed to be changed 

pred_DifferentDay_mtrx <-matrix(percent_DifferentDay, nrow=4, ncol=100, byrow=F)

#the number here needed to be changed 
row_names2 <- pred_DifferentDay_tranformed[toremove[c(1:4)]]
row_names2 <-gsub(" ","_", row_names2 )


pred_DifferentDay_mtrx  <-as.data.frame(pred_DifferentDay_mtrx )

pred_DifferentDay_mtrx$names <-row_names2

pred_DifferentDay_mtrx


paste0('LASSO_LDA_DEgenes_PREDICTION_Cluster',c_selectID_1,'fromDay',day1,'_toDay_', day2,'.txt')

write.table(pred_DifferentDay_mtrx,paste0('LASSO_LDA_DEgenes_PREDICTION_DifferentDay_Cluster',c_selectID_1,'fromDay',day1,'_toDay_', day2,'.txt'), sep="\t", col.names=F, row.names=F, quote=F)

dat2 <-read.table(paste0('LASSO_LDA_DEgenes_PREDICTION_DifferentDay_Cluster',c_selectID_1,'fromDay',day1,'_toDay_', day2,'.txt'))

##########################plotting see the rStudio version in hte labtop#####################

#-----To summary prediction accuracy------#

list_DE_files <- list.files(pattern="QCall_DEgenes_*")

read.files <-lapply(list_DE_files, function(x){readRDS(x)})

percentAcc_mtrx <-matrix(NA, nrow=100,ncol=length(read.files))

ShortNames <-gsub("QCall_DEgenes_forCluster", "S", list_DE_files)
ShortNames <-gsub("_fromDay_", ":D", ShortNames)
ShortNames <-gsub("_toDay_.*", "", ShortNames)
colnames(percentAcc_mtrx) <-ShortNames

for (i in 1:length(read.files)){
	percentAcc_mtrx[,i] <- read.files[[i]]$percentAcc
}



list_KnownMarkers_files <- list.files(pattern="QCall_KnownMarkers_*")

read.files_KnownMarkers <-lapply(list_KnownMarkers_files, function(x){readRDS(x)})

percentAcc_Kmarker_mtrx <-matrix(NA, nrow=100,ncol=length(list_KnownMarkers_files))

ShortNamesKmarkers <-gsub("QCall_KnownMarkers_forCluster", "S", list_KnownMarkers_files)
ShortNamesKmarkers <-gsub("_fromDay_", ":D", ShortNamesKmarkers)
ShortNamesKmarkers <-gsub("_toDay_.*", "", ShortNamesKmarkers)
colnames(percentAcc_Kmarker_mtrx) <-ShortNamesKmarkers

for (i in 1:length(list_KnownMarkers_files)){
	percentAcc_Kmarker_mtrx[,i] <- read.files_KnownMarkers[[i]]$percentAcc
}

NewNames <-c("D0:S1", "D15:S1", "D2:S1", "D30:S1", "D5:S1","D0:S2","D15:S2","D2:S2","D30:S2","D5:S2","D0:S3","D2:S3","D5:S3", "D5:S4")

colnames(percentAcc_Kmarker_mtrx) <-NewNames
colnames(percentAcc_mtrx) <-NewNames


percentAcc_mtrx <-as.data.frame(percentAcc_mtrx)
percentAcc_mtrx$Types <-"DEgenes"

percentAcc_Kmarker_mtrx <-as.data.frame(percentAcc_Kmarker_mtrx)
percentAcc_Kmarker_mtrx$Types <-"KnownMarkers"

merged_dat <-rbind(percentAcc_Kmarker_mtrx, percentAcc_mtrx)



write.table(merged_dat, "Summary_accuracy_all_KnownMarker_DEgenes.txt", col.names=T, row.names=F, quote=F, sep="\t")



#-----To summary prediction percent and stdev------#



list_withinDay_files <- list.files(pattern="LASSO_LDA_DEgenes_PREDICTION_Cluster*")

read.files <-lapply(list_withinDay_files, function(x){read.table(x)})

percentAcc_mtrx <-list()
for (i in 1:length(read.files)){
	temp <-as.matrix(read.files[[i]])
	t <- matrix(as.numeric(temp[, 1:ncol(temp) - 1]), ncol=ncol(temp) - 1)
	mean_temp <-rowMeans(t)
	stdev_temp <-apply(t, 1, sd)
	names_temp <-temp[, ncol(temp)]
	names_temp <-gsub("LASSO_From_", "", names_temp) 
	names_temp <-gsub("_to_predict_cluster_", "C", names_temp) 
	names_temp <-gsub("_From_cluster", "C", names_temp) 
	combined_temp <-cbind(mean_temp, stdev_temp,names_temp )
	percentAcc_mtrx[[i]] <- combined_temp
}


percentAcc_mtrx

#For prediction between different days

list_DiffDay_files <- list.files(pattern="LASSO_LDA_DEgenes_PREDICTION_DifferentDay_Cluster*")

read.files <-lapply(list_DiffDay_files, function(x){read.table(x)})

percentAcc_mtrx <-list()
for (i in 1:length(read.files)){
	temp <-as.matrix(read.files[[i]])
	t <- matrix(as.numeric(temp[, 1:ncol(temp) - 1]), ncol=ncol(temp) - 1)
	mean_temp <-rowMeans(t)
	stdev_temp <-apply(t, 1, sd)
	names_temp <-temp[, ncol(temp)]
	names_temp <-gsub("LASSO_From_", "", names_temp) 
	names_temp <-gsub("_to_predict_cluster_", "C", names_temp) 
	names_temp <-gsub("_From_cluster", "C", names_temp) 
	combined_temp <-cbind(mean_temp, stdev_temp,names_temp )
	percentAcc_mtrx[[i]] <- combined_temp
}








