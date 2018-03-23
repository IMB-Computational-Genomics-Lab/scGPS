---
title: "scGPS introduction"
author: "Quan Nguyen"
date: "2018-03-23"

output: 
  html_document:
    keep_md: true

---



####**Installation instruction**


```r

#Prior to installing scGPS you need to install the SummarizedExperiment bioconductor package as the following
#source("https://bioconductor.org/biocLite.R")
#biocLite("SummarizedExperiment")

#R/3.4.1 or above is required 

#To install scGPS from github (Depending on the configuration of the local computer or HPC, possible custom C++ compilation may be required - see installation trouble-shootings below) 
devtools::install_github("IMB-Computational-Genomics-Lab/scGPS")

#for C++ compilation trouble-shooting, manual download and installation can be done from github
git clone https://github.com/IMB-Computational-Genomics-Lab/scGPS

#then check in scGPS/src if any of the precompiled (e.g. those with *.so and *.o) files exist and delete them before recompiling

#create a Makevars file in the  scGPS/src with one line: PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

#then with the scGPS as the R working directory, manually recompile scGPS in R using devtools
devools::document() 

#and then you can load the package: 
devtools::load_all()

```

####**A simple workflow of the scGPS: given a mixed population with known subpopulations, estimate transition scores between these subpopulation**


```r
devtools::load_all()

#load mixed population 1 (loaded from sample1 dataset, named it as day2)
day2 <- sample1
mixedpop1 <-NewscGPS_SME(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo,
                     CellMetadata = day2$dat2_clusters)
#load mixed population 2 (loaded from sample2 dataset, named it as day5)
day5 <- sample2
mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo,
                     CellMetadata = day5$dat5_clusters)
#load gene list (this can be any lists of user selected genes)
genes <-GeneList
genes <-genes$Merged_unique

#select a subpopulation
c_selectID <- 1

#run the test bootstrap
suppressWarnings(LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list()))
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>        Df       %Dev    Lambda
#>   [1,]  0 -2.563e-15 2.674e-01
#>   [2,]  1  3.515e-02 2.437e-01
#>   [3,]  2  7.249e-02 2.220e-01
#>   [4,]  2  1.076e-01 2.023e-01
#>   [5,]  2  1.380e-01 1.843e-01
#>   [6,]  2  1.647e-01 1.679e-01
#>   [7,]  2  1.884e-01 1.530e-01
#>   [8,]  3  2.134e-01 1.394e-01
#>   [9,]  3  2.358e-01 1.270e-01
#>  [10,]  4  2.570e-01 1.158e-01
#>  [11,]  5  2.795e-01 1.055e-01
#>  [12,]  7  3.006e-01 9.610e-02
#>  [13,]  7  3.202e-01 8.757e-02
#>  [14,]  7  3.378e-01 7.979e-02
#>  [15,]  7  3.536e-01 7.270e-02
#>  [16,] 10  3.695e-01 6.624e-02
#>  [17,] 10  3.865e-01 6.036e-02
#>  [18,] 11  4.020e-01 5.499e-02
#>  [19,] 14  4.168e-01 5.011e-02
#>  [20,] 16  4.328e-01 4.566e-02
#>  [21,] 21  4.530e-01 4.160e-02
#>  [22,] 23  4.734e-01 3.791e-02
#>  [23,] 23  4.919e-01 3.454e-02
#>  [24,] 24  5.093e-01 3.147e-02
#>  [25,] 24  5.249e-01 2.867e-02
#>  [26,] 29  5.414e-01 2.613e-02
#>  [27,] 31  5.596e-01 2.381e-02
#>  [28,] 37  5.791e-01 2.169e-02
#>  [29,] 38  5.990e-01 1.976e-02
#>  [30,] 40  6.181e-01 1.801e-02
#>  [31,] 41  6.358e-01 1.641e-02
#>  [32,] 41  6.520e-01 1.495e-02
#>  [33,] 43  6.669e-01 1.362e-02
#>  [34,] 44  6.817e-01 1.241e-02
#>  [35,] 48  6.959e-01 1.131e-02
#>  [36,] 48  7.097e-01 1.030e-02
#>  [37,] 50  7.227e-01 9.389e-03
#>  [38,] 50  7.354e-01 8.555e-03
#>  [39,] 50  7.473e-01 7.795e-03
#>  [40,] 50  7.586e-01 7.103e-03
#>  [41,] 53  7.694e-01 6.472e-03
#>  [42,] 56  7.804e-01 5.897e-03
#>  [43,] 57  7.912e-01 5.373e-03
#>  [44,] 56  8.013e-01 4.896e-03
#>  [45,] 56  8.107e-01 4.461e-03
#>  [46,] 57  8.195e-01 4.064e-03
#>  [47,] 58  8.284e-01 3.703e-03
#>  [48,] 61  8.369e-01 3.374e-03
#>  [49,] 62  8.452e-01 3.075e-03
#>  [50,] 67  8.538e-01 2.801e-03
#>  [51,] 68  8.623e-01 2.553e-03
#>  [52,] 68  8.705e-01 2.326e-03
#>  [53,] 71  8.785e-01 2.119e-03
#>  [54,] 72  8.865e-01 1.931e-03
#>  [55,] 72  8.944e-01 1.759e-03
#>  [56,] 72  9.020e-01 1.603e-03
#>  [57,] 72  9.093e-01 1.461e-03
#>  [58,] 73  9.162e-01 1.331e-03
#>  [59,] 73  9.228e-01 1.213e-03
#>  [60,] 73  9.290e-01 1.105e-03
#>  [61,] 73  9.349e-01 1.007e-03
#>  [62,] 73  9.403e-01 9.174e-04
#>  [63,] 74  9.454e-01 8.359e-04
#>  [64,] 74  9.502e-01 7.616e-04
#>  [65,] 74  9.546e-01 6.939e-04
#>  [66,] 74  9.588e-01 6.323e-04
#>  [67,] 74  9.625e-01 5.761e-04
#>  [68,] 74  9.658e-01 5.249e-04
#>  [69,] 75  9.689e-01 4.783e-04
#>  [70,] 75  9.717e-01 4.358e-04
#>  [71,] 75  9.743e-01 3.971e-04
#>  [72,] 75  9.766e-01 3.618e-04
#>  [73,] 75  9.788e-01 3.297e-04
#>  [74,] 75  9.807e-01 3.004e-04
#>  [75,] 75  9.824e-01 2.737e-04
#>  [76,] 75  9.840e-01 2.494e-04
#>  [77,] 75  9.855e-01 2.272e-04
#>  [78,] 75  9.868e-01 2.070e-04
#>  [79,] 75  9.880e-01 1.887e-04
#>  [80,] 75  9.890e-01 1.719e-04
#>  [81,] 75  9.900e-01 1.566e-04
#>  [82,] 75  9.909e-01 1.427e-04
#>  [83,] 75  9.917e-01 1.300e-04
#>  [84,] 75  9.925e-01 1.185e-04
#>  [85,] 75  9.931e-01 1.080e-04
#>  [86,] 75  9.937e-01 9.836e-05
#>  [87,] 75  9.943e-01 8.963e-05
#>  [88,] 75  9.948e-01 8.166e-05
#>  [89,] 75  9.953e-01 7.441e-05
#>  [90,] 75  9.957e-01 6.780e-05
#>  [91,] 75  9.961e-01 6.178e-05
#>  [92,] 75  9.964e-01 5.629e-05
#>  [93,] 75  9.967e-01 5.129e-05
#>  [94,] 75  9.970e-01 4.673e-05
#>  [95,] 75  9.973e-01 4.258e-05
#>  [96,] 75  9.975e-01 3.880e-05
#>  [97,] 75  9.977e-01 3.535e-05
#>  [98,] 76  9.979e-01 3.221e-05
#>  [99,] 76  9.981e-01 2.935e-05
#> [100,] 75  9.983e-01 2.674e-05
#> [1] "done 1"
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>        Df       %Dev    Lambda
#>   [1,]  0 -2.563e-15 2.505e-01
#>   [2,]  3  3.763e-02 2.282e-01
#>   [3,]  3  7.982e-02 2.079e-01
#>   [4,]  3  1.159e-01 1.895e-01
#>   [5,]  3  1.470e-01 1.726e-01
#>   [6,]  3  1.741e-01 1.573e-01
#>   [7,]  5  2.046e-01 1.433e-01
#>   [8,]  6  2.330e-01 1.306e-01
#>   [9,]  6  2.590e-01 1.190e-01
#>  [10,]  6  2.818e-01 1.084e-01
#>  [11,]  8  3.052e-01 9.879e-02
#>  [12,]  9  3.297e-01 9.001e-02
#>  [13,]  9  3.519e-01 8.201e-02
#>  [14,]  9  3.716e-01 7.473e-02
#>  [15,]  9  3.892e-01 6.809e-02
#>  [16,]  9  4.050e-01 6.204e-02
#>  [17,] 11  4.201e-01 5.653e-02
#>  [18,] 14  4.361e-01 5.151e-02
#>  [19,] 15  4.519e-01 4.693e-02
#>  [20,] 18  4.690e-01 4.276e-02
#>  [21,] 21  4.865e-01 3.896e-02
#>  [22,] 25  5.037e-01 3.550e-02
#>  [23,] 26  5.223e-01 3.235e-02
#>  [24,] 26  5.397e-01 2.947e-02
#>  [25,] 30  5.592e-01 2.686e-02
#>  [26,] 33  5.792e-01 2.447e-02
#>  [27,] 35  5.984e-01 2.230e-02
#>  [28,] 37  6.161e-01 2.032e-02
#>  [29,] 40  6.331e-01 1.851e-02
#>  [30,] 47  6.518e-01 1.687e-02
#>  [31,] 48  6.696e-01 1.537e-02
#>  [32,] 49  6.859e-01 1.400e-02
#>  [33,] 52  7.016e-01 1.276e-02
#>  [34,] 51  7.161e-01 1.163e-02
#>  [35,] 53  7.298e-01 1.059e-02
#>  [36,] 54  7.428e-01 9.651e-03
#>  [37,] 56  7.550e-01 8.794e-03
#>  [38,] 57  7.665e-01 8.013e-03
#>  [39,] 58  7.774e-01 7.301e-03
#>  [40,] 58  7.876e-01 6.652e-03
#>  [41,] 60  7.972e-01 6.061e-03
#>  [42,] 61  8.068e-01 5.523e-03
#>  [43,] 60  8.159e-01 5.032e-03
#>  [44,] 60  8.244e-01 4.585e-03
#>  [45,] 62  8.326e-01 4.178e-03
#>  [46,] 63  8.406e-01 3.807e-03
#>  [47,] 65  8.486e-01 3.469e-03
#>  [48,] 66  8.566e-01 3.160e-03
#>  [49,] 68  8.649e-01 2.880e-03
#>  [50,] 67  8.731e-01 2.624e-03
#>  [51,] 67  8.808e-01 2.391e-03
#>  [52,] 67  8.882e-01 2.178e-03
#>  [53,] 70  8.954e-01 1.985e-03
#>  [54,] 70  9.024e-01 1.808e-03
#>  [55,] 71  9.092e-01 1.648e-03
#>  [56,] 71  9.158e-01 1.501e-03
#>  [57,] 72  9.221e-01 1.368e-03
#>  [58,] 73  9.281e-01 1.247e-03
#>  [59,] 73  9.340e-01 1.136e-03
#>  [60,] 75  9.395e-01 1.035e-03
#>  [61,] 74  9.448e-01 9.430e-04
#>  [62,] 74  9.496e-01 8.592e-04
#>  [63,] 73  9.541e-01 7.829e-04
#>  [64,] 72  9.581e-01 7.133e-04
#>  [65,] 72  9.618e-01 6.499e-04
#>  [66,] 71  9.652e-01 5.922e-04
#>  [67,] 71  9.683e-01 5.396e-04
#>  [68,] 72  9.712e-01 4.917e-04
#>  [69,] 72  9.738e-01 4.480e-04
#>  [70,] 73  9.762e-01 4.082e-04
#>  [71,] 73  9.783e-01 3.719e-04
#>  [72,] 73  9.802e-01 3.389e-04
#>  [73,] 73  9.820e-01 3.088e-04
#>  [74,] 73  9.836e-01 2.813e-04
#>  [75,] 73  9.851e-01 2.563e-04
#>  [76,] 74  9.864e-01 2.336e-04
#>  [77,] 74  9.877e-01 2.128e-04
#>  [78,] 74  9.888e-01 1.939e-04
#>  [79,] 73  9.898e-01 1.767e-04
#>  [80,] 73  9.907e-01 1.610e-04
#>  [81,] 74  9.915e-01 1.467e-04
#>  [82,] 74  9.923e-01 1.337e-04
#>  [83,] 74  9.930e-01 1.218e-04
#>  [84,] 74  9.936e-01 1.110e-04
#>  [85,] 74  9.942e-01 1.011e-04
#>  [86,] 74  9.947e-01 9.213e-05
#>  [87,] 74  9.952e-01 8.394e-05
#>  [88,] 74  9.956e-01 7.649e-05
#>  [89,] 74  9.960e-01 6.969e-05
#>  [90,] 74  9.963e-01 6.350e-05
#>  [91,] 74  9.966e-01 5.786e-05
#>  [92,] 74  9.969e-01 5.272e-05
#>  [93,] 74  9.972e-01 4.804e-05
#>  [94,] 75  9.975e-01 4.377e-05
#>  [95,] 74  9.977e-01 3.988e-05
#>  [96,] 74  9.979e-01 3.634e-05
#>  [97,] 74  9.981e-01 3.311e-05
#>  [98,] 74  9.982e-01 3.017e-05
#>  [99,] 74  9.984e-01 2.749e-05
#> [100,] 74  9.985e-01 2.505e-05
#> [1] "done 2"

#display some results 
names(LSOLDA_dat)
#> [1] "Accuracy"     "LassoGenes"   "Deviance"     "LassoFit"    
#> [5] "LDAFit"       "predictor_S1" "LassoPredict" "LDAPredict"
LSOLDA_dat$LassoPredict
#> [[1]]
#> [[1]][[1]]
#> [1] "LASSO for subpop2 in target mixedpop2"
#> 
#> [[1]][[2]]
#> [1] 93.04029
#> 
#> [[1]][[3]]
#> [1] "LASSO for subpop3 in target mixedpop2"
#> 
#> [[1]][[4]]
#> [1] 10.97561
#> 
#> [[1]][[5]]
#> [1] "LASSO for subpop1 in target mixedpop2"
#> 
#> [[1]][[6]]
#> [1] 37.15847
#> 
#> [[1]][[7]]
#> [1] "LASSO for subpop4 in target mixedpop2"
#> 
#> [[1]][[8]]
#> [1] 56.12245
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] "LASSO for subpop2 in target mixedpop2"
#> 
#> [[2]][[2]]
#> [1] 58.97436
#> 
#> [[2]][[3]]
#> [1] "LASSO for subpop3 in target mixedpop2"
#> 
#> [[2]][[4]]
#> [1] 31.70732
#> 
#> [[2]][[5]]
#> [1] "LASSO for subpop1 in target mixedpop2"
#> 
#> [[2]][[6]]
#> [1] 50.27322
#> 
#> [[2]][[7]]
#> [1] "LASSO for subpop4 in target mixedpop2"
#> 
#> [[2]][[8]]
#> [1] 34.69388
LSOLDA_dat$LDAPredict
#> [[1]]
#> [[1]][[1]]
#> [1] "LDA for subpop 2 in target mixedpop2"
#> 
#> [[1]][[2]]
#> [1] 46.52015
#> 
#> [[1]][[3]]
#> [1] "LDA for subpop 3 in target mixedpop2"
#> 
#> [[1]][[4]]
#> [1] 7.317073
#> 
#> [[1]][[5]]
#> [1] "LDA for subpop 1 in target mixedpop2"
#> 
#> [[1]][[6]]
#> [1] 12.02186
#> 
#> [[1]][[7]]
#> [1] "LDA for subpop 4 in target mixedpop2"
#> 
#> [[1]][[8]]
#> [1] 21.42857
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] "LDA for subpop 2 in target mixedpop2"
#> 
#> [[2]][[2]]
#> [1] 40.29304
#> 
#> [[2]][[3]]
#> [1] "LDA for subpop 3 in target mixedpop2"
#> 
#> [[2]][[4]]
#> [1] 8.943089
#> 
#> [[2]][[5]]
#> [1] "LDA for subpop 1 in target mixedpop2"
#> 
#> [[2]][[6]]
#> [1] 4.918033
#> 
#> [[2]][[7]]
#> [1] "LDA for subpop 4 in target mixedpop2"
#> 
#> [[2]][[8]]
#> [1] 8.163265

#summary results LDA
summary_prediction_lda(LSOLDA_dat=LSOLDA_dat, nPredSubpop = 4)
#>                 V1               V2                                names
#> 1 46.5201465201465 40.2930402930403 LDA for subpop 2 in target mixedpop2
#> 2 7.31707317073171 8.94308943089431 LDA for subpop 3 in target mixedpop2
#> 3 12.0218579234973 4.91803278688525 LDA for subpop 1 in target mixedpop2
#> 4 21.4285714285714 8.16326530612245 LDA for subpop 4 in target mixedpop2

#summary results Lasso
summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop = 4)
#>                 V1               V2                                 names
#> 1  93.040293040293  58.974358974359 LASSO for subpop2 in target mixedpop2
#> 2 10.9756097560976 31.7073170731707 LASSO for subpop3 in target mixedpop2
#> 3 37.1584699453552 50.2732240437158 LASSO for subpop1 in target mixedpop2
#> 4 56.1224489795918 34.6938775510204 LASSO for subpop4 in target mixedpop2

#summary deviance 
summary_deviance(LSOLDA_dat)
#> $allDeviance
#> [1] "0.4919" "0.5397"
#> 
#> $DeviMax
#>          Dfd   Deviance        DEgenes
#> 1          0 -2.563e-15 genes_cluster1
#> 2          3     0.1741 genes_cluster1
#> 3          5     0.2046 genes_cluster1
#> 4          6     0.2818 genes_cluster1
#> 5          8     0.3052 genes_cluster1
#> 6          9      0.405 genes_cluster1
#> 7         11     0.4201 genes_cluster1
#> 8         14     0.4361 genes_cluster1
#> 9         15     0.4519 genes_cluster1
#> 10        18      0.469 genes_cluster1
#> 11        21     0.4865 genes_cluster1
#> 12        25     0.5037 genes_cluster1
#> 13        26     0.5397 genes_cluster1
#> 14 remaining          1        DEgenes
#> 
#> $LassoGenesMax
#>                                   1                   name
#> (Intercept)             0.703693814            (Intercept)
#> NPPB_ENSG00000120937    1.906902170   NPPB_ENSG00000120937
#> TPM3_ENSG00000143549    0.064487180   TPM3_ENSG00000143549
#> CXCR4_ENSG00000121966  -0.105050497  CXCR4_ENSG00000121966
#> TTN_ENSG00000155657     1.149149641    TTN_ENSG00000155657
#> FN1_ENSG00000115414    -0.044602262    FN1_ENSG00000115414
#> CLDN1_ENSG00000163347   0.344129424  CLDN1_ENSG00000163347
#> PDGFRA_ENSG00000134853  0.190914737 PDGFRA_ENSG00000134853
#> FOXC1_ENSG00000054598   0.080904513  FOXC1_ENSG00000054598
#> POU5F1_ENSG00000204531 -0.024385118 POU5F1_ENSG00000204531
#> GJA1_ENSG00000152661   -0.176438370   GJA1_ENSG00000152661
#> T_ENSG00000164458       0.241162180      T_ENSG00000164458
#> MYL7_ENSG00000106631    0.103592680   MYL7_ENSG00000106631
#> SNAI2_ENSG00000019549   0.255299986  SNAI2_ENSG00000019549
#> SOX17_ENSG00000164736  -0.007471097  SOX17_ENSG00000164736
#> HEY1_ENSG00000164683   -0.267241112   HEY1_ENSG00000164683
#> MYC_ENSG00000136997     0.047296408    MYC_ENSG00000136997
#> ZBTB16_ENSG00000109906  0.946598990 ZBTB16_ENSG00000109906
#> NODAL_ENSG00000156574   0.043598868  NODAL_ENSG00000156574
#> COL2A1_ENSG00000139219 -0.069607726 COL2A1_ENSG00000139219
#> FOXA1_ENSG00000129514  -0.395261449  FOXA1_ENSG00000129514
#> TPM1_ENSG00000140416   -0.059937858   TPM1_ENSG00000140416
#> MESP1_ENSG00000166823   0.057196801  MESP1_ENSG00000166823
#> FOXF1_ENSG00000103241   0.804699956  FOXF1_ENSG00000103241
#> FOXA2_ENSG00000125798  -0.441158964  FOXA2_ENSG00000125798
#> SNAI1_ENSG00000124216   0.274275947  SNAI1_ENSG00000124216
#> FOXA3_ENSG00000170608  -0.114814811  FOXA3_ENSG00000170608
```

####**A complete workflow of the scGPS: given an unknown mixed population, find clusters and estimate relationship between clusters** 



```r
#given a single cell expression matrix, without clustering information 
day5 <- sample2
cellnames <- colnames(day5$dat5_counts)
cluster <-day5$dat5_clusters
cellnames <-data.frame("Cluster"=cluster, "cellBarcodes" = cellnames)

mixedpop2 <-NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo, CellMetadata = cellnames ) 

#let's find the clusters
CORE_cluster <- CORE_scGPS(mixedpop2, remove_outlier = c(0))
#> [1] "Calculating distance matrix"
#> [1] "Performing hierarchical clustering"
#> [1] "Finding clustering information"
#> [1] "No more outliers detected after 1 filtering round"
#> [1] "writing clustering result for run 1"
#> [1] "writing clustering result for run 2"
#> [1] "writing clustering result for run 3"
#> [1] "writing clustering result for run 4"
#> [1] "writing clustering result for run 5"
#> [1] "writing clustering result for run 6"
#> [1] "writing clustering result for run 7"
#> [1] "writing clustering result for run 8"
#> [1] "writing clustering result for run 9"
#> [1] "writing clustering result for run 10"
#> [1] "writing clustering result for run 11"
#> [1] "writing clustering result for run 12"
#> [1] "writing clustering result for run 13"
#> [1] "writing clustering result for run 14"
#> [1] "writing clustering result for run 15"
#> [1] "writing clustering result for run 16"
#> [1] "writing clustering result for run 17"
#> [1] "writing clustering result for run 18"
#> [1] "writing clustering result for run 19"
#> [1] "writing clustering result for run 20"
#> [1] "writing clustering result for run 21"
#> [1] "writing clustering result for run 22"
#> [1] "writing clustering result for run 23"
#> [1] "writing clustering result for run 24"
#> [1] "writing clustering result for run 25"
#> [1] "writing clustering result for run 26"
#> [1] "writing clustering result for run 27"
#> [1] "writing clustering result for run 28"
#> [1] "writing clustering result for run 29"
#> [1] "writing clustering result for run 30"
#> [1] "writing clustering result for run 31"
#> [1] "writing clustering result for run 32"
#> [1] "writing clustering result for run 33"
#> [1] "writing clustering result for run 34"
#> [1] "writing clustering result for run 35"
#> [1] "writing clustering result for run 36"
#> [1] "writing clustering result for run 37"
#> [1] "writing clustering result for run 38"
#> [1] "writing clustering result for run 39"
#> [1] "writing clustering result for run 40"
#> [1] "Done clustering, moving to stability calculation..."
#> [1] "Done calculating stability..."
#> [1] "Start finding optimal clustering..."
#> [1] "Done finding optimal clustering..."
#let's plot all clusters
plot_CORE(CORE_cluster$tree, CORE_cluster$Cluster)
```

![](vignette_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

  
  #####Let's plot just the optimal clustering result (with colored dendrogram)

```r
optimal_index = which(CORE_cluster$optimalClust$KeyStats$Height == CORE_cluster$optimalClust$OptimalRes)

plot_optimal_CORE(original_tree= CORE_cluster$tree, optimal_cluster = unlist(CORE_cluster$Cluster[optimal_index]), shift = -200)
#> [1] "Ordering and assigning labels..."
#> [1] 2
#> [1] 204 424  NA  NA
#> [1] 3
#> [1] 204 424 536  NA
#> [1] 4
#> [1] 204 424 536 808
#> [1] "Plotting the colored dendrogram now...."
```

![](vignette_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
#> [1] "Plotting the bar underneath now...."
#> [1] "Users are required to check cluster labels...."
```
  
  #####Let's compare with other dimensional reduction methods 

```r
library(cidr)
t <- CIDR_scGPS(expression.matrix=assay(mixedpop2))
#> [1] "building cidr object..."
#> [1] "determine dropout candidates..."
#> [1] "determine the imputation weighting threshold..."
#> [1] "computes the _CIDR_ dissimilarity matrix..."
#> [1] "PCA plot with proportion of variance explained..."
#> [1] "find the number of PC..."
#> [1] "perform clustering..."
p2 <-plotReduced_scGPS(t, color_fac = factor(colData(mixedpop2)[,1]),palletes =1:length(unique(colData(mixedpop2)[,1])))
```

![](vignette_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
#may need to turn off the R graphic device dev.off() before plotting the following
p2
```

![](vignette_files/figure-html/unnamed-chunk-5-2.png)<!-- -->
  
  #####Find gene markers and annotate clusters 

```r

#load gene list (this can be any lists of user selected genes)
genes <-GeneList
genes <-genes$Merged_unique

#the gene list can also be generated objectively by differential expression analysis
#if the cluster information is in the mixedpop2 object, run this (if not run the CORE
#as described below) 
DEgenes <- findMarkers_scGPS(expression_matrix=assay(mixedpop2), cluster = mixedpop2@colData$Cluster)
#> [1] "Start estimate dispersions for cluster 2..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 2..."
#> [1] "Done nbinom test for cluster 2 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
#> [1] "Start estimate dispersions for cluster 3..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 3..."
#> [1] "Done nbinom test for cluster 3 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
#> [1] "Start estimate dispersions for cluster 1..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 1..."
#> [1] "Done nbinom test for cluster 1 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
#> [1] "Start estimate dispersions for cluster 4..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 4..."
#> [1] "Done nbinom test for cluster 4 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
names(DEgenes)
#> [1] "DE_Subpop2vsRemaining" "DE_Subpop3vsRemaining" "DE_Subpop1vsRemaining"
#> [4] "DE_Subpop4vsRemaining"

#you can annotate the identified clusters 
DEgeneList_3vsOthers <- DEgenes$DE_Subpop3vsRemaining$id
#format to ensembl genes 
DEgeneList_3vsOthers <-gsub("_.*", "", DEgeneList_3vsOthers )
head(DEgeneList_3vsOthers)
#> [1] "TTR"       "SERPINE2"  "APOA1"     "H2AFZ"     "APOA2"     "LINC01356"
#the following command saves the file "PathwayEnrichment.xlsx" to the working dir
enrichment_test <- annotate_scGPS(DEgeneList_3vsOthers, pvalueCutoff=0.05, gene_symbol=TRUE,output_filename = "PathwayEnrichment.xlsx", output_path = NULL )
#> 'select()' returned 1:many mapping between keys and columns
#> Warning in bitr(DEgeneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb
#> = "org.Hs.eg.db"): 13.12% of input gene IDs are fail to map...
#> [1] "Original gene number in geneList"
#> [1] 17402
#> [1] "Number of genes successfully converted"
#> [1] 15116
#note can conveniently plot the enrichment outputs by running the followings
dotplot(enrichment_test, showCategory=15)
```

![](vignette_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r

#users need to check the format of the gene input to make sure they are consistent to 
#the gene names in the expression matrix 

#add the CORE cluster information into the scGPS object 
Optimal_index <- which( CORE_cluster$optimalClust$KeyStats$Height == CORE_cluster$optimalClust$OptimalRes)
mixedpop2@colData$Cluster <- unlist(CORE_cluster$Cluster[[Optimal_index]])


#select a subpopulation
c_selectID <- 1

#run the test bootstrap with nboots = 2 runs
suppressWarnings(LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop2, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list()))
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>         Df       %Dev    Lambda
#>   [1,]   0 -2.883e-15 2.444e-01
#>   [2,]   2  3.887e-02 2.227e-01
#>   [3,]   3  7.612e-02 2.029e-01
#>   [4,]   4  1.092e-01 1.849e-01
#>   [5,]   5  1.421e-01 1.685e-01
#>   [6,]   6  1.715e-01 1.535e-01
#>   [7,]   6  1.979e-01 1.399e-01
#>   [8,]   7  2.214e-01 1.274e-01
#>   [9,]   7  2.423e-01 1.161e-01
#>  [10,]   7  2.608e-01 1.058e-01
#>  [11,]   9  2.782e-01 9.640e-02
#>  [12,]   9  2.958e-01 8.784e-02
#>  [13,]  10  3.125e-01 8.003e-02
#>  [14,]  10  3.278e-01 7.292e-02
#>  [15,]  12  3.420e-01 6.645e-02
#>  [16,]  13  3.549e-01 6.054e-02
#>  [17,]  15  3.696e-01 5.516e-02
#>  [18,]  18  3.856e-01 5.026e-02
#>  [19,]  21  4.009e-01 4.580e-02
#>  [20,]  22  4.155e-01 4.173e-02
#>  [21,]  22  4.296e-01 3.802e-02
#>  [22,]  23  4.429e-01 3.464e-02
#>  [23,]  23  4.550e-01 3.157e-02
#>  [24,]  25  4.666e-01 2.876e-02
#>  [25,]  29  4.780e-01 2.621e-02
#>  [26,]  32  4.904e-01 2.388e-02
#>  [27,]  37  5.033e-01 2.176e-02
#>  [28,]  41  5.159e-01 1.982e-02
#>  [29,]  43  5.286e-01 1.806e-02
#>  [30,]  48  5.408e-01 1.646e-02
#>  [31,]  54  5.530e-01 1.500e-02
#>  [32,]  57  5.650e-01 1.366e-02
#>  [33,]  60  5.763e-01 1.245e-02
#>  [34,]  63  5.868e-01 1.134e-02
#>  [35,]  67  5.968e-01 1.034e-02
#>  [36,]  71  6.062e-01 9.418e-03
#>  [37,]  73  6.152e-01 8.582e-03
#>  [38,]  76  6.237e-01 7.819e-03
#>  [39,]  78  6.322e-01 7.125e-03
#>  [40,]  81  6.405e-01 6.492e-03
#>  [41,]  82  6.482e-01 5.915e-03
#>  [42,]  84  6.553e-01 5.390e-03
#>  [43,]  84  6.625e-01 4.911e-03
#>  [44,]  87  6.693e-01 4.475e-03
#>  [45,]  88  6.764e-01 4.077e-03
#>  [46,]  92  6.839e-01 3.715e-03
#>  [47,]  95  6.909e-01 3.385e-03
#>  [48,]  95  6.975e-01 3.084e-03
#>  [49,]  96  7.039e-01 2.810e-03
#>  [50,]  98  7.107e-01 2.560e-03
#>  [51,]  99  7.181e-01 2.333e-03
#>  [52,] 100  7.256e-01 2.126e-03
#>  [53,] 101  7.330e-01 1.937e-03
#>  [54,] 103  7.402e-01 1.765e-03
#>  [55,] 104  7.469e-01 1.608e-03
#>  [56,] 105  7.534e-01 1.465e-03
#>  [57,] 105  7.596e-01 1.335e-03
#>  [58,] 105  7.654e-01 1.216e-03
#>  [59,] 105  7.708e-01 1.108e-03
#>  [60,] 105  7.759e-01 1.010e-03
#>  [61,] 105  7.807e-01 9.202e-04
#>  [62,] 106  7.852e-01 8.384e-04
#>  [63,] 107  7.895e-01 7.640e-04
#>  [64,] 105  7.936e-01 6.961e-04
#>  [65,] 106  7.975e-01 6.343e-04
#>  [66,] 107  8.012e-01 5.779e-04
#>  [67,] 107  8.047e-01 5.266e-04
#>  [68,] 109  8.079e-01 4.798e-04
#>  [69,] 109  8.112e-01 4.372e-04
#>  [70,] 109  8.144e-01 3.983e-04
#>  [71,] 109  8.172e-01 3.629e-04
#>  [72,] 111  8.200e-01 3.307e-04
#>  [73,] 112  8.226e-01 3.013e-04
#>  [74,] 112  8.250e-01 2.746e-04
#>  [75,] 113  8.271e-01 2.502e-04
#>  [76,] 113  8.291e-01 2.279e-04
#>  [77,] 113  8.311e-01 2.077e-04
#>  [78,] 113  8.329e-01 1.892e-04
#>  [79,] 112  8.346e-01 1.724e-04
#>  [80,] 112  8.362e-01 1.571e-04
#>  [81,] 112  8.379e-01 1.432e-04
#>  [82,] 111  8.396e-01 1.304e-04
#>  [83,] 110  8.413e-01 1.188e-04
#>  [84,] 110  8.428e-01 1.083e-04
#>  [85,] 107  8.443e-01 9.867e-05
#>  [86,] 107  8.454e-01 8.990e-05
#>  [87,] 106  8.466e-01 8.192e-05
#>  [88,] 107  8.479e-01 7.464e-05
#>  [89,] 107  8.492e-01 6.801e-05
#>  [90,] 109  8.503e-01 6.197e-05
#>  [91,] 109  8.515e-01 5.646e-05
#>  [92,] 109  8.526e-01 5.145e-05
#>  [93,] 109  8.533e-01 4.688e-05
#>  [94,] 111  8.544e-01 4.271e-05
#>  [95,] 111  8.549e-01 3.892e-05
#>  [96,] 111  8.559e-01 3.546e-05
#>  [97,] 111  8.562e-01 3.231e-05
#>  [98,] 111  8.571e-01 2.944e-05
#>  [99,] 110  8.573e-01 2.682e-05
#> [100,] 111  8.579e-01 2.444e-05
#> [1] "done 1"
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>         Df       %Dev    Lambda
#>   [1,]   0 -2.883e-15 2.506e-01
#>   [2,]   1  3.096e-02 2.284e-01
#>   [3,]   3  6.789e-02 2.081e-01
#>   [4,]   4  1.030e-01 1.896e-01
#>   [5,]   4  1.367e-01 1.728e-01
#>   [6,]   4  1.660e-01 1.574e-01
#>   [7,]   4  1.917e-01 1.434e-01
#>   [8,]   5  2.148e-01 1.307e-01
#>   [9,]   5  2.371e-01 1.191e-01
#>  [10,]   6  2.567e-01 1.085e-01
#>  [11,]   9  2.759e-01 9.886e-02
#>  [12,]  10  2.940e-01 9.008e-02
#>  [13,]  10  3.111e-01 8.207e-02
#>  [14,]  11  3.262e-01 7.478e-02
#>  [15,]  12  3.396e-01 6.814e-02
#>  [16,]  12  3.520e-01 6.209e-02
#>  [17,]  13  3.635e-01 5.657e-02
#>  [18,]  14  3.740e-01 5.155e-02
#>  [19,]  15  3.844e-01 4.697e-02
#>  [20,]  17  3.941e-01 4.279e-02
#>  [21,]  19  4.044e-01 3.899e-02
#>  [22,]  24  4.157e-01 3.553e-02
#>  [23,]  24  4.273e-01 3.237e-02
#>  [24,]  26  4.387e-01 2.950e-02
#>  [25,]  28  4.501e-01 2.688e-02
#>  [26,]  32  4.615e-01 2.449e-02
#>  [27,]  35  4.730e-01 2.231e-02
#>  [28,]  36  4.835e-01 2.033e-02
#>  [29,]  39  4.934e-01 1.852e-02
#>  [30,]  42  5.024e-01 1.688e-02
#>  [31,]  43  5.110e-01 1.538e-02
#>  [32,]  49  5.195e-01 1.401e-02
#>  [33,]  54  5.295e-01 1.277e-02
#>  [34,]  57  5.392e-01 1.163e-02
#>  [35,]  62  5.495e-01 1.060e-02
#>  [36,]  64  5.602e-01 9.659e-03
#>  [37,]  65  5.705e-01 8.801e-03
#>  [38,]  67  5.801e-01 8.019e-03
#>  [39,]  70  5.890e-01 7.306e-03
#>  [40,]  72  5.972e-01 6.657e-03
#>  [41,]  75  6.049e-01 6.066e-03
#>  [42,]  73  6.118e-01 5.527e-03
#>  [43,]  76  6.182e-01 5.036e-03
#>  [44,]  79  6.241e-01 4.589e-03
#>  [45,]  81  6.298e-01 4.181e-03
#>  [46,]  86  6.355e-01 3.810e-03
#>  [47,]  88  6.410e-01 3.471e-03
#>  [48,]  91  6.461e-01 3.163e-03
#>  [49,]  93  6.510e-01 2.882e-03
#>  [50,]  94  6.556e-01 2.626e-03
#>  [51,]  94  6.598e-01 2.393e-03
#>  [52,]  97  6.635e-01 2.180e-03
#>  [53,]  97  6.670e-01 1.986e-03
#>  [54,]  96  6.702e-01 1.810e-03
#>  [55,]  98  6.732e-01 1.649e-03
#>  [56,]  99  6.759e-01 1.503e-03
#>  [57,] 100  6.784e-01 1.369e-03
#>  [58,] 104  6.808e-01 1.247e-03
#>  [59,] 105  6.831e-01 1.137e-03
#>  [60,] 107  6.853e-01 1.036e-03
#>  [61,] 108  6.874e-01 9.437e-04
#>  [62,] 109  6.895e-01 8.598e-04
#>  [63,] 109  6.914e-01 7.834e-04
#>  [64,] 111  6.931e-01 7.138e-04
#>  [65,] 111  6.946e-01 6.504e-04
#>  [66,] 111  6.961e-01 5.926e-04
#>  [67,] 110  6.973e-01 5.400e-04
#>  [68,] 110  6.984e-01 4.920e-04
#>  [69,] 109  6.994e-01 4.483e-04
#>  [70,] 108  7.003e-01 4.085e-04
#>  [71,] 110  7.010e-01 3.722e-04
#>  [72,] 112  7.017e-01 3.391e-04
#>  [73,] 112  7.023e-01 3.090e-04
#>  [74,] 114  7.028e-01 2.816e-04
#>  [75,] 114  7.033e-01 2.565e-04
#>  [76,] 114  7.037e-01 2.338e-04
#>  [77,] 114  7.041e-01 2.130e-04
#>  [78,] 115  7.044e-01 1.941e-04
#>  [79,] 115  7.047e-01 1.768e-04
#>  [80,] 115  7.049e-01 1.611e-04
#>  [81,] 115  7.052e-01 1.468e-04
#>  [82,] 115  7.054e-01 1.338e-04
#>  [83,] 115  7.055e-01 1.219e-04
#>  [84,] 115  7.057e-01 1.111e-04
#>  [85,] 115  7.058e-01 1.012e-04
#>  [86,] 115  7.059e-01 9.220e-05
#>  [87,] 115  7.060e-01 8.401e-05
#>  [88,] 115  7.061e-01 7.654e-05
#>  [89,] 115  7.061e-01 6.974e-05
#>  [90,] 115  7.062e-01 6.355e-05
#>  [91,] 115  7.063e-01 5.790e-05
#>  [92,] 115  7.063e-01 5.276e-05
#>  [93,] 115  7.063e-01 4.807e-05
#>  [94,] 115  7.064e-01 4.380e-05
#>  [95,] 115  7.064e-01 3.991e-05
#>  [96,] 115  7.064e-01 3.636e-05
#>  [97,] 116  7.065e-01 3.313e-05
#>  [98,] 116  7.065e-01 3.019e-05
#>  [99,] 117  7.065e-01 2.751e-05
#> [100,] 117  7.065e-01 2.506e-05
#> [1] "done 2"

#display summary results 
#summary results LDA
row_cluster <-length(unique(mixedpop2@colData$Cluster))
summary_prediction_lda(LSOLDA_dat=LSOLDA_dat, nPredSubpop = row_cluster )
#>                 V1               V2                                names
#> 1 81.3725490196078 81.8627450980392 LDA for subpop 1 in target mixedpop2
#> 2 17.4285714285714 15.4285714285714 LDA for subpop 2 in target mixedpop2
#> 3 13.4715025906736 13.9896373056995 LDA for subpop 3 in target mixedpop2
#> 4               25           21.875 LDA for subpop 4 in target mixedpop2

#summary results Lasso
summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop = row_cluster)
#>                 V1               V2                                 names
#> 1 81.8627450980392 81.1274509803922 LASSO for subpop1 in target mixedpop2
#> 2 15.1428571428571 12.8571428571429 LASSO for subpop2 in target mixedpop2
#> 3  10.880829015544 11.9170984455959 LASSO for subpop3 in target mixedpop2
#> 4            18.75           21.875 LASSO for subpop4 in target mixedpop2

#summary deviance 
summary_deviance(LSOLDA_dat)
#> $allDeviance
#> [1] "0.6062" "0.374" 
#> 
#> $DeviMax
#>          Dfd   Deviance        DEgenes
#> 1          0 -2.883e-15 genes_cluster1
#> 2          2    0.03887 genes_cluster1
#> 3          3    0.07612 genes_cluster1
#> 4          4     0.1092 genes_cluster1
#> 5          5     0.1421 genes_cluster1
#> 6          6     0.1979 genes_cluster1
#> 7          7     0.2608 genes_cluster1
#> 8          9     0.2958 genes_cluster1
#> 9         10     0.3278 genes_cluster1
#> 10        12      0.342 genes_cluster1
#> 11        13     0.3549 genes_cluster1
#> 12        15     0.3696 genes_cluster1
#> 13        18     0.3856 genes_cluster1
#> 14        21     0.4009 genes_cluster1
#> 15        22     0.4296 genes_cluster1
#> 16        23      0.455 genes_cluster1
#> 17        25     0.4666 genes_cluster1
#> 18        29      0.478 genes_cluster1
#> 19        32     0.4904 genes_cluster1
#> 20        37     0.5033 genes_cluster1
#> 21        41     0.5159 genes_cluster1
#> 22        43     0.5286 genes_cluster1
#> 23        48     0.5408 genes_cluster1
#> 24        54      0.553 genes_cluster1
#> 25        57      0.565 genes_cluster1
#> 26        60     0.5763 genes_cluster1
#> 27        63     0.5868 genes_cluster1
#> 28        67     0.5968 genes_cluster1
#> 29        71     0.6062 genes_cluster1
#> 30 remaining          1        DEgenes
#> 
#> $LassoGenesMax
#>                                     1                    name
#> (Intercept)             -4.630438e-01             (Intercept)
#> NPPA_ENSG00000175206    -4.125500e-01    NPPA_ENSG00000175206
#> NPPB_ENSG00000120937     3.443263e-05    NPPB_ENSG00000120937
#> FCN3_ENSG00000142748     5.777290e-02    FCN3_ENSG00000142748
#> TAL1_ENSG00000162367     1.112969e-01    TAL1_ENSG00000162367
#> TPM3_ENSG00000143549    -2.111942e-01    TPM3_ENSG00000143549
#> RGS4_ENSG00000117152     3.626504e-02    RGS4_ENSG00000117152
#> NR5A2_ENSG00000116833    6.976268e-01   NR5A2_ENSG00000116833
#> TNNI1_ENSG00000159173    7.518817e-03   TNNI1_ENSG00000159173
#> CD34_ENSG00000174059    -1.295037e-02    CD34_ENSG00000174059
#> CXCR4_ENSG00000121966   -3.114906e-01   CXCR4_ENSG00000121966
#> MYO3B_ENSG00000071909    3.491400e-01   MYO3B_ENSG00000071909
#> FN1_ENSG00000115414      1.627898e-02     FN1_ENSG00000115414
#> EOMES_ENSG00000163508   -4.489676e-01   EOMES_ENSG00000163508
#> TNNC1_ENSG00000114854   -7.221320e-02   TNNC1_ENSG00000114854
#> HESX1_ENSG00000163666   -2.362801e-02   HESX1_ENSG00000163666
#> SST_ENSG00000157005     -2.418495e+00     SST_ENSG00000157005
#> MSX1_ENSG00000163132    -3.549699e-02    MSX1_ENSG00000163132
#> PDGFRA_ENSG00000134853   1.715565e-01  PDGFRA_ENSG00000134853
#> ODAM_ENSG00000109205     6.509572e-02    ODAM_ENSG00000109205
#> AFP_ENSG00000081051     -1.302778e+00     AFP_ENSG00000081051
#> CXCL5_ENSG00000163735   -7.453958e-01   CXCL5_ENSG00000163735
#> HAND2_ENSG00000164107    2.810149e-01   HAND2_ENSG00000164107
#> IRX2_ENSG00000170561     1.257410e+00    IRX2_ENSG00000170561
#> IRX1_ENSG00000170549    -3.219266e-01    IRX1_ENSG00000170549
#> IL6ST_ENSG00000134352   -1.725337e-01   IL6ST_ENSG00000134352
#> MEF2C_ENSG00000081189   -9.763946e-02   MEF2C_ENSG00000081189
#> HAND1_ENSG00000113196    4.749478e-02   HAND1_ENSG00000113196
#> POU5F1_ENSG00000204531  -1.543094e-01  POU5F1_ENSG00000204531
#> SRF_ENSG00000112658      4.702787e-01     SRF_ENSG00000112658
#> GJA1_ENSG00000152661    -2.305091e-01    GJA1_ENSG00000152661
#> T_ENSG00000164458       -1.153808e+00       T_ENSG00000164458
#> TBX20_ENSG00000164532    1.893405e-01   TBX20_ENSG00000164532
#> GATA4_ENSG00000136574    7.958437e-02   GATA4_ENSG00000136574
#> SNAI2_ENSG00000019549    6.275557e-04   SNAI2_ENSG00000019549
#> SOX17_ENSG00000164736   -6.226182e-01   SOX17_ENSG00000164736
#> HEY1_ENSG00000164683     2.727922e-01    HEY1_ENSG00000164683
#> SDC2_ENSG00000169439    -1.870796e-01    SDC2_ENSG00000169439
#> MYC_ENSG00000136997      8.570504e-01     MYC_ENSG00000136997
#> THY1_ENSG00000154096     9.837578e-02    THY1_ENSG00000154096
#> VIM_ENSG00000026025      9.650017e-03     VIM_ENSG00000026025
#> NODAL_ENSG00000156574    1.568184e-01   NODAL_ENSG00000156574
#> HHEX_ENSG00000152804    -6.586621e-02    HHEX_ENSG00000152804
#> CACNA1C_ENSG00000151067  2.328150e-01 CACNA1C_ENSG00000151067
#> CD9_ENSG00000010278     -2.449195e-01     CD9_ENSG00000010278
#> NANOG_ENSG00000111704    1.687831e-01   NANOG_ENSG00000111704
#> ATP2A2_ENSG00000174437   1.758317e-01  ATP2A2_ENSG00000174437
#> TBX3_ENSG00000135111    -1.746431e-01    TBX3_ENSG00000135111
#> CDX2_ENSG00000165556    -6.426846e-03    CDX2_ENSG00000165556
#> KLF5_ENSG00000102554    -2.093684e+00    KLF5_ENSG00000102554
#> PAPLN_ENSG00000100767   -1.627414e-01   PAPLN_ENSG00000100767
#> GSC_ENSG00000133937     -1.102509e-01     GSC_ENSG00000133937
#> ACTC1_ENSG00000159251    1.717451e-02   ACTC1_ENSG00000159251
#> MEIS2_ENSG00000134138    6.218361e-02   MEIS2_ENSG00000134138
#> MESP1_ENSG00000166823    7.269706e-01   MESP1_ENSG00000166823
#> MESP2_ENSG00000188095    9.588729e-02   MESP2_ENSG00000188095
#> NR2F2_ENSG00000185551    4.783700e-01   NR2F2_ENSG00000185551
#> CDH5_ENSG00000179776    -2.673434e-01    CDH5_ENSG00000179776
#> FOXF1_ENSG00000103241   -2.367287e+00   FOXF1_ENSG00000103241
#> FOXC2_ENSG00000176692   -1.564603e-01   FOXC2_ENSG00000176692
#> NOS2_ENSG00000007171     2.392834e-02    NOS2_ENSG00000007171
#> VTN_ENSG00000109072     -2.571607e-02     VTN_ENSG00000109072
#> HNF1B_ENSG00000275410   -3.970306e-01   HNF1B_ENSG00000275410
#> MYL4_ENSG00000198336     3.121863e-02    MYL4_ENSG00000198336
#> GATA6_ENSG00000141448   -2.946932e-02   GATA6_ENSG00000141448
#> FOXA2_ENSG00000125798   -1.187547e+00   FOXA2_ENSG00000125798
#> HNF4A_ENSG00000101076   -1.200698e-01   HNF4A_ENSG00000101076
#> SNAI1_ENSG00000124216   -5.598797e-01   SNAI1_ENSG00000124216
#> LYL1_ENSG00000104903    -2.004452e-02    LYL1_ENSG00000104903
#> PLVAP_ENSG00000130300   -1.769311e-01   PLVAP_ENSG00000130300
#> FOXA3_ENSG00000170608   -2.928374e-01   FOXA3_ENSG00000170608
#> TNNT1_ENSG00000105048   -8.181756e-02   TNNT1_ENSG00000105048
```

