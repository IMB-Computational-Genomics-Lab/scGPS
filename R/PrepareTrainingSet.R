#prepare subpopulations, the training subpop is subpop1
#inputs are two subpopulations (e.g. same day or different days) and a list of genes (the ordered list from most to least important)
#it generates one training dataset and a reference classification matrix
#the input matrix has gene names in rows and cell names in columns (make sure gene names in the same format)

#Prepare an object here, making sure it has 2 subpop slots

prepare_2_subpop <-function(subpop1, subpop2, genes){
  # making sure the gene list contains fewer than 500 genes
  # select top 500, the DE_result is an already sorted list
  if(length(genes) >500){genes <- genes[1:500]}

  names1 <-rownames(subpop1)
  subpop1_selected_genes <- names1[which(names1 %in% genes)]
  names2 <-rownames(subpop2)
  selected_genes_in_both <-names2[which(names2 %in% subpop1_selected_genes)]

  subpop1_train <- subpop1[which(names1 %in% subpop1_selected_genes)]
  subpop2_target <- subpop2[which(names2 %in% subpop1_selected_genes)]

  # reference classification for comparing predicted results
  cellNames_subpop <- cbind(colnames(subpop1), colnames(subpop2)) #note: may edit this to maje 2 genetic classes 1 and 2
}


#---------------------------------------------
#after this step we have: subsampling, cluster_select, cluster_compare, DE_idx, ori_dat, ori_dat2
#---------------------------------------------


