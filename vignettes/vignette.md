---
title: "scGPS introduction"
author: "Quan Nguyen"
date: "2018-04-08"

output: 
  html_document:
    keep_md: TRUE
    theme: united
    highlight: tango

---



#Installation instruction


```r
# Prior to installing scGPS you need to install the SummarizedExperiment
# bioconductor package as the following
# source('https://bioconductor.org/biocLite.R') biocLite('SummarizedExperiment')

# R/3.4.1 or above is required

# To install scGPS from github (Depending on the configuration of the local
# computer or HPC, possible custom C++ compilation may be required - see
# installation trouble-shootings below)
devtools::install_github("IMB-Computational-Genomics-Lab/scGPS")

# for C++ compilation trouble-shooting, manual download and installation can be
# done from github

# git clone https://github.com/IMB-Computational-Genomics-Lab/scGPS

# then check in scGPS/src if any of the precompiled (e.g.  those with *.so and
# *.o) files exist and delete them before recompiling

# create a Makevars file in the scGPS/src with one line: PKG_LIBS =
# $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

# then with the scGPS as the R working directory, manually recompile scGPS in R
# using devtools
devools::document()

```

#A simple workflow of the scGPS: given a mixed population with known subpopulations, estimate transition scores between these subpopulation


```r
devtools::load_all()

# load mixed population 1 (loaded from sample1 dataset, named it as day2)
day2 <- sample1
mixedpop1 <- NewscGPS_SME(ExpressionMatrix = day2$dat2_counts, GeneMetadata = day2$dat2geneInfo, 
    CellMetadata = day2$dat2_clusters)
# load mixed population 2 (loaded from sample2 dataset, named it as day5)
day5 <- sample2
mixedpop2 <- NewscGPS_SME(ExpressionMatrix = day5$dat5_counts, GeneMetadata = day5$dat5geneInfo, 
    CellMetadata = day5$dat5_clusters)
# load gene list (this can be any lists of user selected genes)
genes <- GeneList
genes <- genes$Merged_unique

# select a subpopulation
c_selectID <- 1

# run the test bootstrap
suppressWarnings(LSOLDA_dat <- bootstrap_scGPS(nboots = 2, mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, genes = genes, c_selectID, listData = list()))
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>        Df       %Dev    Lambda
#>   [1,]  0 -2.563e-15 2.802e-01
#>   [2,]  1  3.866e-02 2.553e-01
#>   [3,]  2  7.328e-02 2.326e-01
#>   [4,]  3  1.166e-01 2.119e-01
#>   [5,]  4  1.587e-01 1.931e-01
#>   [6,]  4  2.009e-01 1.760e-01
#>   [7,]  4  2.378e-01 1.603e-01
#>   [8,]  5  2.718e-01 1.461e-01
#>   [9,]  6  3.046e-01 1.331e-01
#>  [10,]  7  3.347e-01 1.213e-01
#>  [11,]  7  3.617e-01 1.105e-01
#>  [12,]  8  3.865e-01 1.007e-01
#>  [13,]  9  4.097e-01 9.174e-02
#>  [14,]  9  4.313e-01 8.359e-02
#>  [15,]  9  4.506e-01 7.617e-02
#>  [16,]  9  4.678e-01 6.940e-02
#>  [17,]  9  4.832e-01 6.323e-02
#>  [18,]  9  4.969e-01 5.762e-02
#>  [19,]  9  5.092e-01 5.250e-02
#>  [20,] 12  5.212e-01 4.783e-02
#>  [21,] 14  5.355e-01 4.358e-02
#>  [22,] 16  5.502e-01 3.971e-02
#>  [23,] 19  5.643e-01 3.618e-02
#>  [24,] 19  5.778e-01 3.297e-02
#>  [25,] 23  5.919e-01 3.004e-02
#>  [26,] 27  6.059e-01 2.737e-02
#>  [27,] 29  6.197e-01 2.494e-02
#>  [28,] 31  6.338e-01 2.273e-02
#>  [29,] 33  6.483e-01 2.071e-02
#>  [30,] 36  6.628e-01 1.887e-02
#>  [31,] 40  6.781e-01 1.719e-02
#>  [32,] 41  6.937e-01 1.566e-02
#>  [33,] 41  7.084e-01 1.427e-02
#>  [34,] 42  7.220e-01 1.300e-02
#>  [35,] 44  7.349e-01 1.185e-02
#>  [36,] 45  7.476e-01 1.080e-02
#>  [37,] 47  7.594e-01 9.837e-03
#>  [38,] 49  7.706e-01 8.963e-03
#>  [39,] 49  7.813e-01 8.167e-03
#>  [40,] 52  7.919e-01 7.441e-03
#>  [41,] 51  8.018e-01 6.780e-03
#>  [42,] 52  8.111e-01 6.178e-03
#>  [43,] 53  8.198e-01 5.629e-03
#>  [44,] 53  8.281e-01 5.129e-03
#>  [45,] 54  8.359e-01 4.673e-03
#>  [46,] 55  8.435e-01 4.258e-03
#>  [47,] 58  8.521e-01 3.880e-03
#>  [48,] 58  8.608e-01 3.535e-03
#>  [49,] 59  8.695e-01 3.221e-03
#>  [50,] 62  8.781e-01 2.935e-03
#>  [51,] 62  8.866e-01 2.674e-03
#>  [52,] 63  8.949e-01 2.437e-03
#>  [53,] 62  9.030e-01 2.220e-03
#>  [54,] 62  9.106e-01 2.023e-03
#>  [55,] 62  9.179e-01 1.843e-03
#>  [56,] 63  9.247e-01 1.680e-03
#>  [57,] 63  9.311e-01 1.530e-03
#>  [58,] 64  9.372e-01 1.394e-03
#>  [59,] 63  9.426e-01 1.271e-03
#>  [60,] 64  9.478e-01 1.158e-03
#>  [61,] 64  9.525e-01 1.055e-03
#>  [62,] 63  9.568e-01 9.611e-04
#>  [63,] 63  9.607e-01 8.757e-04
#>  [64,] 63  9.642e-01 7.979e-04
#>  [65,] 62  9.675e-01 7.270e-04
#>  [66,] 62  9.704e-01 6.624e-04
#>  [67,] 62  9.730e-01 6.036e-04
#>  [68,] 61  9.755e-01 5.500e-04
#>  [69,] 62  9.777e-01 5.011e-04
#>  [70,] 62  9.797e-01 4.566e-04
#>  [71,] 62  9.815e-01 4.160e-04
#>  [72,] 62  9.832e-01 3.791e-04
#>  [73,] 63  9.847e-01 3.454e-04
#>  [74,] 63  9.860e-01 3.147e-04
#>  [75,] 63  9.873e-01 2.868e-04
#>  [76,] 64  9.884e-01 2.613e-04
#>  [77,] 64  9.895e-01 2.381e-04
#>  [78,] 64  9.904e-01 2.169e-04
#>  [79,] 64  9.913e-01 1.977e-04
#>  [80,] 64  9.920e-01 1.801e-04
#>  [81,] 65  9.927e-01 1.641e-04
#>  [82,] 65  9.934e-01 1.495e-04
#>  [83,] 65  9.940e-01 1.362e-04
#>  [84,] 65  9.945e-01 1.241e-04
#>  [85,] 65  9.950e-01 1.131e-04
#>  [86,] 65  9.954e-01 1.031e-04
#>  [87,] 65  9.958e-01 9.390e-05
#>  [88,] 65  9.962e-01 8.556e-05
#>  [89,] 65  9.965e-01 7.796e-05
#>  [90,] 65  9.969e-01 7.103e-05
#>  [91,] 65  9.971e-01 6.472e-05
#>  [92,] 66  9.974e-01 5.897e-05
#>  [93,] 66  9.976e-01 5.373e-05
#>  [94,] 66  9.978e-01 4.896e-05
#>  [95,] 66  9.980e-01 4.461e-05
#>  [96,] 66  9.982e-01 4.065e-05
#>  [97,] 66  9.983e-01 3.704e-05
#>  [98,] 66  9.985e-01 3.375e-05
#>  [99,] 66  9.986e-01 3.075e-05
#> [100,] 66  9.987e-01 2.802e-05
#> [1] "done 1"
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>        Df       %Dev    Lambda
#>   [1,]  0 -2.563e-15 2.711e-01
#>   [2,]  2  4.118e-02 2.470e-01
#>   [3,]  2  8.227e-02 2.251e-01
#>   [4,]  2  1.175e-01 2.051e-01
#>   [5,]  2  1.479e-01 1.869e-01
#>   [6,]  3  1.753e-01 1.703e-01
#>   [7,]  4  2.025e-01 1.551e-01
#>   [8,]  5  2.317e-01 1.413e-01
#>   [9,]  5  2.578e-01 1.288e-01
#>  [10,]  6  2.810e-01 1.173e-01
#>  [11,]  6  3.029e-01 1.069e-01
#>  [12,]  8  3.230e-01 9.743e-02
#>  [13,]  8  3.448e-01 8.877e-02
#>  [14,]  9  3.644e-01 8.088e-02
#>  [15,]  9  3.819e-01 7.370e-02
#>  [16,] 10  3.986e-01 6.715e-02
#>  [17,] 11  4.141e-01 6.119e-02
#>  [18,] 11  4.284e-01 5.575e-02
#>  [19,] 11  4.412e-01 5.080e-02
#>  [20,] 14  4.543e-01 4.629e-02
#>  [21,] 15  4.681e-01 4.217e-02
#>  [22,] 18  4.838e-01 3.843e-02
#>  [23,] 20  5.016e-01 3.501e-02
#>  [24,] 20  5.201e-01 3.190e-02
#>  [25,] 21  5.370e-01 2.907e-02
#>  [26,] 23  5.539e-01 2.649e-02
#>  [27,] 30  5.713e-01 2.413e-02
#>  [28,] 34  5.889e-01 2.199e-02
#>  [29,] 34  6.057e-01 2.004e-02
#>  [30,] 38  6.214e-01 1.826e-02
#>  [31,] 40  6.379e-01 1.663e-02
#>  [32,] 42  6.533e-01 1.516e-02
#>  [33,] 43  6.676e-01 1.381e-02
#>  [34,] 46  6.813e-01 1.258e-02
#>  [35,] 50  6.951e-01 1.147e-02
#>  [36,] 50  7.087e-01 1.045e-02
#>  [37,] 51  7.216e-01 9.519e-03
#>  [38,] 51  7.337e-01 8.673e-03
#>  [39,] 53  7.452e-01 7.902e-03
#>  [40,] 53  7.560e-01 7.200e-03
#>  [41,] 56  7.667e-01 6.561e-03
#>  [42,] 57  7.780e-01 5.978e-03
#>  [43,] 58  7.887e-01 5.447e-03
#>  [44,] 60  7.995e-01 4.963e-03
#>  [45,] 61  8.104e-01 4.522e-03
#>  [46,] 61  8.209e-01 4.120e-03
#>  [47,] 62  8.309e-01 3.754e-03
#>  [48,] 61  8.407e-01 3.421e-03
#>  [49,] 61  8.502e-01 3.117e-03
#>  [50,] 63  8.596e-01 2.840e-03
#>  [51,] 64  8.687e-01 2.588e-03
#>  [52,] 64  8.776e-01 2.358e-03
#>  [53,] 65  8.861e-01 2.148e-03
#>  [54,] 65  8.947e-01 1.958e-03
#>  [55,] 67  9.029e-01 1.784e-03
#>  [56,] 67  9.111e-01 1.625e-03
#>  [57,] 67  9.186e-01 1.481e-03
#>  [58,] 69  9.257e-01 1.349e-03
#>  [59,] 68  9.323e-01 1.229e-03
#>  [60,] 68  9.383e-01 1.120e-03
#>  [61,] 68  9.439e-01 1.021e-03
#>  [62,] 68  9.489e-01 9.300e-04
#>  [63,] 66  9.534e-01 8.474e-04
#>  [64,] 66  9.576e-01 7.721e-04
#>  [65,] 66  9.615e-01 7.035e-04
#>  [66,] 66  9.650e-01 6.410e-04
#>  [67,] 66  9.681e-01 5.841e-04
#>  [68,] 66  9.710e-01 5.322e-04
#>  [69,] 66  9.736e-01 4.849e-04
#>  [70,] 66  9.760e-01 4.418e-04
#>  [71,] 66  9.782e-01 4.026e-04
#>  [72,] 65  9.801e-01 3.668e-04
#>  [73,] 65  9.819e-01 3.342e-04
#>  [74,] 65  9.835e-01 3.045e-04
#>  [75,] 66  9.850e-01 2.775e-04
#>  [76,] 67  9.863e-01 2.528e-04
#>  [77,] 66  9.876e-01 2.304e-04
#>  [78,] 66  9.887e-01 2.099e-04
#>  [79,] 66  9.897e-01 1.913e-04
#>  [80,] 67  9.906e-01 1.743e-04
#>  [81,] 67  9.914e-01 1.588e-04
#>  [82,] 67  9.922e-01 1.447e-04
#>  [83,] 67  9.929e-01 1.318e-04
#>  [84,] 67  9.935e-01 1.201e-04
#>  [85,] 67  9.941e-01 1.094e-04
#>  [86,] 67  9.946e-01 9.972e-05
#>  [87,] 67  9.951e-01 9.086e-05
#>  [88,] 67  9.955e-01 8.279e-05
#>  [89,] 67  9.959e-01 7.543e-05
#>  [90,] 67  9.963e-01 6.873e-05
#>  [91,] 67  9.966e-01 6.263e-05
#>  [92,] 67  9.969e-01 5.706e-05
#>  [93,] 67  9.972e-01 5.199e-05
#>  [94,] 67  9.974e-01 4.737e-05
#>  [95,] 67  9.977e-01 4.317e-05
#>  [96,] 67  9.979e-01 3.933e-05
#>  [97,] 67  9.980e-01 3.584e-05
#>  [98,] 67  9.982e-01 3.265e-05
#>  [99,] 67  9.984e-01 2.975e-05
#> [100,] 69  9.985e-01 2.711e-05
#> [1] "done 2"

# display some results
names(LSOLDA_dat)
#> [1] "Accuracy"     "LassoGenes"   "Deviance"     "LassoFit"    
#> [5] "LDAFit"       "predictor_S1" "LassoPredict" "LDAPredict"
LSOLDA_dat$LassoPredict
#> [[1]]
#> [[1]][[1]]
#> [1] "LASSO for subpop1 in target mixedpop2"
#> 
#> [[1]][[2]]
#> [1] 70.05348
#> 
#> [[1]][[3]]
#> [1] "LASSO for subpop2 in target mixedpop2"
#> 
#> [[1]][[4]]
#> [1] 96.42857
#> 
#> [[1]][[5]]
#> [1] "LASSO for subpop3 in target mixedpop2"
#> 
#> [[1]][[6]]
#> [1] 46.61654
#> 
#> [[1]][[7]]
#> [1] "LASSO for subpop4 in target mixedpop2"
#> 
#> [[1]][[8]]
#> [1] 67.5
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] "LASSO for subpop1 in target mixedpop2"
#> 
#> [[2]][[2]]
#> [1] 67.91444
#> 
#> [[2]][[3]]
#> [1] "LASSO for subpop2 in target mixedpop2"
#> 
#> [[2]][[4]]
#> [1] 93.57143
#> 
#> [[2]][[5]]
#> [1] "LASSO for subpop3 in target mixedpop2"
#> 
#> [[2]][[6]]
#> [1] 45.86466
#> 
#> [[2]][[7]]
#> [1] "LASSO for subpop4 in target mixedpop2"
#> 
#> [[2]][[8]]
#> [1] 60
LSOLDA_dat$LDAPredict
#> [[1]]
#> [[1]][[1]]
#> [1] "LDA for subpop 1 in target mixedpop2"
#> 
#> [[1]][[2]]
#> [1] 16.04278
#> 
#> [[1]][[3]]
#> [1] "LDA for subpop 2 in target mixedpop2"
#> 
#> [[1]][[4]]
#> [1] 43.57143
#> 
#> [[1]][[5]]
#> [1] "LDA for subpop 3 in target mixedpop2"
#> 
#> [[1]][[6]]
#> [1] 9.774436
#> 
#> [[1]][[7]]
#> [1] "LDA for subpop 4 in target mixedpop2"
#> 
#> [[1]][[8]]
#> [1] 27.5
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] "LDA for subpop 1 in target mixedpop2"
#> 
#> [[2]][[2]]
#> [1] 98.39572
#> 
#> [[2]][[3]]
#> [1] "LDA for subpop 2 in target mixedpop2"
#> 
#> [[2]][[4]]
#> [1] 77.85714
#> 
#> [[2]][[5]]
#> [1] "LDA for subpop 3 in target mixedpop2"
#> 
#> [[2]][[6]]
#> [1] 97.74436
#> 
#> [[2]][[7]]
#> [1] "LDA for subpop 4 in target mixedpop2"
#> 
#> [[2]][[8]]
#> [1] 82.5

# summary results LDA
summary_prediction_lda(LSOLDA_dat = LSOLDA_dat, nPredSubpop = 4)
#>                 V1               V2                                names
#> 1 16.0427807486631 98.3957219251337 LDA for subpop 1 in target mixedpop2
#> 2 43.5714285714286 77.8571428571429 LDA for subpop 2 in target mixedpop2
#> 3 9.77443609022556 97.7443609022556 LDA for subpop 3 in target mixedpop2
#> 4             27.5             82.5 LDA for subpop 4 in target mixedpop2

# summary results Lasso
summary_prediction_lasso(LSOLDA_dat = LSOLDA_dat, nPredSubpop = 4)
#>                 V1               V2                                 names
#> 1 70.0534759358289 67.9144385026738 LASSO for subpop1 in target mixedpop2
#> 2 96.4285714285714 93.5714285714286 LASSO for subpop2 in target mixedpop2
#> 3 46.6165413533835 45.8646616541353 LASSO for subpop3 in target mixedpop2
#> 4             67.5               60 LASSO for subpop4 in target mixedpop2

# summary deviance
```

#A complete workflow of the scGPS: given an unknown mixed population, find clusters and estimate relationship between clusters


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

##Let's plot just the optimal clustering result (with colored dendrogram)
  

```r
optimal_index = which(CORE_cluster$optimalClust$KeyStats$Height == CORE_cluster$optimalClust$OptimalRes)

plot_optimal_CORE(original_tree= CORE_cluster$tree, optimal_cluster = unlist(CORE_cluster$Cluster[optimal_index]), shift = -1500)
#> [1] "Ordering and assigning labels..."
#> [1] 2
#> [1] 128 270  NA
#> [1] 3
#> [1] 128 270 393
#> [1] "Plotting the colored dendrogram now...."
```

![](vignette_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
#> [1] "Plotting the bar underneath now...."
#> [1] "Users are required to check cluster labels...."
```
  
##Let's compare with other dimensional reduction methods 


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
  
##Find gene markers and annotate clusters


```r

#load gene list (this can be any lists of user selected genes)
genes <-GeneList
genes <-genes$Merged_unique

#the gene list can also be generated objectively by differential expression analysis
#if the cluster information is in the mixedpop2 object, run this (if not run the CORE
#as described below) 
DEgenes <- findMarkers_scGPS(expression_matrix=assay(mixedpop2), cluster = colData(mixedpop2)[,1])
#> [1] "Start estimate dispersions for cluster 1..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 1..."
#> [1] "Done nbinom test for cluster 1 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
#> [1] "Start estimate dispersions for cluster 2..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 2..."
#> [1] "Done nbinom test for cluster 2 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
#> [1] "Start estimate dispersions for cluster 3..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 3..."
#> [1] "Done nbinom test for cluster 3 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
#> [1] "Start estimate dispersions for cluster 4..."
#> [1] "Done estimate dispersions. Start nbinom test for cluster 4..."
#> [1] "Done nbinom test for cluster 4 ..."
#> [1] "Adjust foldchange by subtracting basemean to 1..."
names(DEgenes)
#> [1] "DE_Subpop1vsRemaining" "DE_Subpop2vsRemaining" "DE_Subpop3vsRemaining"
#> [4] "DE_Subpop4vsRemaining"

#you can annotate the identified clusters 
DEgeneList_3vsOthers <- DEgenes$DE_Subpop3vsRemaining$id
#format to ensembl genes 
DEgeneList_3vsOthers <-gsub("_.*", "", DEgeneList_3vsOthers )
head(DEgeneList_3vsOthers)
#> [1] "TTR"       "SERPINE2"  "APOA1"     "APOA2"     "H2AFZ"     "LINC01356"
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
colData(mixedpop2)[,1] <- unlist(CORE_cluster$Cluster[[Optimal_index]])
```

##Start the prediction


```r

#select a subpopulation
c_selectID <- 1

#run the test bootstrap with nboots = 2 runs
LSOLDA_dat <- bootstrap_scGPS(nboots = 2,mixedpop1 = mixedpop2, mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list())
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>        Df       %Dev    Lambda
#>   [1,]  0 -1.922e-15 2.653e-01
#>   [2,]  1  3.460e-02 2.418e-01
#>   [3,]  1  6.389e-02 2.203e-01
#>   [4,]  4  9.254e-02 2.007e-01
#>   [5,]  4  1.248e-01 1.829e-01
#>   [6,]  4  1.528e-01 1.666e-01
#>   [7,]  5  1.781e-01 1.518e-01
#>   [8,]  5  2.024e-01 1.383e-01
#>   [9,]  5  2.238e-01 1.261e-01
#>  [10,]  5  2.427e-01 1.149e-01
#>  [11,]  5  2.595e-01 1.046e-01
#>  [12,]  5  2.744e-01 9.535e-02
#>  [13,]  5  2.876e-01 8.688e-02
#>  [14,]  5  2.993e-01 7.916e-02
#>  [15,]  7  3.117e-01 7.213e-02
#>  [16,]  8  3.237e-01 6.572e-02
#>  [17,] 12  3.368e-01 5.988e-02
#>  [18,] 12  3.500e-01 5.456e-02
#>  [19,] 13  3.640e-01 4.972e-02
#>  [20,] 17  3.781e-01 4.530e-02
#>  [21,] 23  3.925e-01 4.128e-02
#>  [22,] 26  4.128e-01 3.761e-02
#>  [23,] 28  4.311e-01 3.427e-02
#>  [24,] 29  4.479e-01 3.122e-02
#>  [25,] 33  4.656e-01 2.845e-02
#>  [26,] 34  4.847e-01 2.592e-02
#>  [27,] 41  5.039e-01 2.362e-02
#>  [28,] 44  5.235e-01 2.152e-02
#>  [29,] 48  5.423e-01 1.961e-02
#>  [30,] 50  5.614e-01 1.787e-02
#>  [31,] 51  5.792e-01 1.628e-02
#>  [32,] 52  5.957e-01 1.483e-02
#>  [33,] 54  6.117e-01 1.352e-02
#>  [34,] 56  6.280e-01 1.232e-02
#>  [35,] 58  6.438e-01 1.122e-02
#>  [36,] 62  6.593e-01 1.022e-02
#>  [37,] 63  6.744e-01 9.316e-03
#>  [38,] 65  6.889e-01 8.488e-03
#>  [39,] 66  7.027e-01 7.734e-03
#>  [40,] 67  7.158e-01 7.047e-03
#>  [41,] 69  7.288e-01 6.421e-03
#>  [42,] 72  7.415e-01 5.851e-03
#>  [43,] 73  7.543e-01 5.331e-03
#>  [44,] 76  7.669e-01 4.857e-03
#>  [45,] 78  7.791e-01 4.426e-03
#>  [46,] 80  7.911e-01 4.033e-03
#>  [47,] 80  8.024e-01 3.674e-03
#>  [48,] 80  8.132e-01 3.348e-03
#>  [49,] 80  8.233e-01 3.051e-03
#>  [50,] 81  8.330e-01 2.780e-03
#>  [51,] 81  8.424e-01 2.533e-03
#>  [52,] 82  8.512e-01 2.308e-03
#>  [53,] 84  8.597e-01 2.103e-03
#>  [54,] 84  8.679e-01 1.916e-03
#>  [55,] 85  8.757e-01 1.746e-03
#>  [56,] 86  8.833e-01 1.591e-03
#>  [57,] 88  8.909e-01 1.449e-03
#>  [58,] 88  8.980e-01 1.321e-03
#>  [59,] 88  9.048e-01 1.203e-03
#>  [60,] 89  9.114e-01 1.096e-03
#>  [61,] 90  9.177e-01 9.989e-04
#>  [62,] 89  9.238e-01 9.102e-04
#>  [63,] 88  9.296e-01 8.293e-04
#>  [64,] 87  9.351e-01 7.557e-04
#>  [65,] 87  9.403e-01 6.885e-04
#>  [66,] 85  9.451e-01 6.274e-04
#>  [67,] 84  9.496e-01 5.716e-04
#>  [68,] 85  9.538e-01 5.208e-04
#>  [69,] 85  9.578e-01 4.746e-04
#>  [70,] 85  9.615e-01 4.324e-04
#>  [71,] 84  9.649e-01 3.940e-04
#>  [72,] 84  9.680e-01 3.590e-04
#>  [73,] 83  9.708e-01 3.271e-04
#>  [74,] 83  9.735e-01 2.980e-04
#>  [75,] 83  9.759e-01 2.716e-04
#>  [76,] 83  9.780e-01 2.474e-04
#>  [77,] 84  9.800e-01 2.255e-04
#>  [78,] 84  9.818e-01 2.054e-04
#>  [79,] 83  9.834e-01 1.872e-04
#>  [80,] 84  9.849e-01 1.706e-04
#>  [81,] 85  9.863e-01 1.554e-04
#>  [82,] 85  9.875e-01 1.416e-04
#>  [83,] 85  9.886e-01 1.290e-04
#>  [84,] 85  9.896e-01 1.176e-04
#>  [85,] 85  9.905e-01 1.071e-04
#>  [86,] 85  9.914e-01 9.760e-05
#>  [87,] 84  9.921e-01 8.893e-05
#>  [88,] 84  9.928e-01 8.103e-05
#>  [89,] 84  9.935e-01 7.383e-05
#>  [90,] 84  9.940e-01 6.727e-05
#>  [91,] 85  9.946e-01 6.129e-05
#>  [92,] 84  9.950e-01 5.585e-05
#>  [93,] 83  9.955e-01 5.089e-05
#>  [94,] 83  9.959e-01 4.637e-05
#>  [95,] 83  9.962e-01 4.225e-05
#>  [96,] 84  9.965e-01 3.849e-05
#>  [97,] 84  9.968e-01 3.507e-05
#>  [98,] 83  9.971e-01 3.196e-05
#>  [99,] 84  9.974e-01 2.912e-05
#> [100,] 84  9.976e-01 2.653e-05
#> [1] "done 1"
#> 
#> Call:  glmnet(x = t(predictor_S1), y = y_cat, family = "binomial") 
#> 
#>         Df       %Dev    Lambda
#>   [1,]   0 -1.922e-15 2.687e-01
#>   [2,]   2  3.739e-02 2.448e-01
#>   [3,]   2  7.869e-02 2.230e-01
#>   [4,]   2  1.139e-01 2.032e-01
#>   [5,]   4  1.464e-01 1.852e-01
#>   [6,]   4  1.760e-01 1.687e-01
#>   [7,]   4  2.017e-01 1.537e-01
#>   [8,]   5  2.249e-01 1.401e-01
#>   [9,]   5  2.458e-01 1.276e-01
#>  [10,]   5  2.641e-01 1.163e-01
#>  [11,]   5  2.801e-01 1.060e-01
#>  [12,]   6  2.946e-01 9.655e-02
#>  [13,]   6  3.078e-01 8.797e-02
#>  [14,]   6  3.193e-01 8.016e-02
#>  [15,]   7  3.309e-01 7.304e-02
#>  [16,]   9  3.444e-01 6.655e-02
#>  [17,]  10  3.582e-01 6.064e-02
#>  [18,]  10  3.707e-01 5.525e-02
#>  [19,]  11  3.819e-01 5.034e-02
#>  [20,]  14  3.957e-01 4.587e-02
#>  [21,]  17  4.100e-01 4.179e-02
#>  [22,]  21  4.237e-01 3.808e-02
#>  [23,]  23  4.375e-01 3.470e-02
#>  [24,]  26  4.525e-01 3.162e-02
#>  [25,]  27  4.661e-01 2.881e-02
#>  [26,]  29  4.787e-01 2.625e-02
#>  [27,]  31  4.915e-01 2.392e-02
#>  [28,]  34  5.042e-01 2.179e-02
#>  [29,]  36  5.167e-01 1.986e-02
#>  [30,]  38  5.289e-01 1.809e-02
#>  [31,]  45  5.408e-01 1.648e-02
#>  [32,]  44  5.520e-01 1.502e-02
#>  [33,]  50  5.632e-01 1.369e-02
#>  [34,]  51  5.748e-01 1.247e-02
#>  [35,]  55  5.867e-01 1.136e-02
#>  [36,]  58  5.979e-01 1.035e-02
#>  [37,]  62  6.093e-01 9.433e-03
#>  [38,]  64  6.202e-01 8.595e-03
#>  [39,]  64  6.301e-01 7.831e-03
#>  [40,]  64  6.390e-01 7.136e-03
#>  [41,]  67  6.471e-01 6.502e-03
#>  [42,]  71  6.551e-01 5.924e-03
#>  [43,]  73  6.631e-01 5.398e-03
#>  [44,]  73  6.704e-01 4.918e-03
#>  [45,]  77  6.775e-01 4.481e-03
#>  [46,]  79  6.847e-01 4.083e-03
#>  [47,]  83  6.921e-01 3.721e-03
#>  [48,]  83  7.000e-01 3.390e-03
#>  [49,]  84  7.074e-01 3.089e-03
#>  [50,]  84  7.143e-01 2.814e-03
#>  [51,]  87  7.208e-01 2.564e-03
#>  [52,]  91  7.271e-01 2.337e-03
#>  [53,]  94  7.332e-01 2.129e-03
#>  [54,]  95  7.390e-01 1.940e-03
#>  [55,]  94  7.443e-01 1.768e-03
#>  [56,]  94  7.491e-01 1.611e-03
#>  [57,]  94  7.534e-01 1.467e-03
#>  [58,]  94  7.573e-01 1.337e-03
#>  [59,]  94  7.609e-01 1.218e-03
#>  [60,]  97  7.646e-01 1.110e-03
#>  [61,]  97  7.687e-01 1.011e-03
#>  [62,]  97  7.728e-01 9.216e-04
#>  [63,]  96  7.768e-01 8.397e-04
#>  [64,]  97  7.805e-01 7.651e-04
#>  [65,]  99  7.844e-01 6.972e-04
#>  [66,] 100  7.884e-01 6.352e-04
#>  [67,] 101  7.927e-01 5.788e-04
#>  [68,] 101  7.972e-01 5.274e-04
#>  [69,] 100  8.020e-01 4.805e-04
#>  [70,]  98  8.070e-01 4.378e-04
#>  [71,]  98  8.123e-01 3.989e-04
#>  [72,]  97  8.178e-01 3.635e-04
#>  [73,] 100  8.243e-01 3.312e-04
#>  [74,] 101  8.316e-01 3.018e-04
#>  [75,] 100  8.392e-01 2.750e-04
#>  [76,]  99  8.464e-01 2.505e-04
#>  [77,]  99  8.549e-01 2.283e-04
#>  [78,]  99  8.644e-01 2.080e-04
#>  [79,]  99  8.775e-01 1.895e-04
#>  [80,]  99  8.888e-01 1.727e-04
#>  [81,]  98  8.993e-01 1.574e-04
#>  [82,]  98  9.089e-01 1.434e-04
#>  [83,]  97  9.175e-01 1.306e-04
#>  [84,]  98  9.253e-01 1.190e-04
#>  [85,]  99  9.329e-01 1.085e-04
#>  [86,]  98  9.396e-01 9.882e-05
#>  [87,]  98  9.453e-01 9.004e-05
#>  [88,]  98  9.505e-01 8.204e-05
#>  [89,]  98  9.551e-01 7.475e-05
#>  [90,]  98  9.593e-01 6.811e-05
#>  [91,]  98  9.630e-01 6.206e-05
#>  [92,]  98  9.663e-01 5.655e-05
#>  [93,]  97  9.695e-01 5.153e-05
#>  [94,]  97  9.721e-01 4.695e-05
#>  [95,]  97  9.744e-01 4.278e-05
#>  [96,]  98  9.764e-01 3.898e-05
#>  [97,]  99  9.782e-01 3.551e-05
#>  [98,]  99  9.800e-01 3.236e-05
#>  [99,]  98  9.817e-01 2.948e-05
#> [100,]  98  9.832e-01 2.687e-05
#> [1] "done 2"
```

##Display summary results for the prediction


```r
#summary results LDA
row_cluster <-length(unique(colData(mixedpop2)[,1]))
summary_prediction_lda(LSOLDA_dat=LSOLDA_dat, nPredSubpop = row_cluster )
#>                 V1               V2                                names
#> 1        78.515625          76.5625 LDA for subpop 1 in target mixedpop2
#> 2 17.2093023255814 15.8139534883721 LDA for subpop 2 in target mixedpop2
#> 3 17.2413793103448 13.7931034482759 LDA for subpop 3 in target mixedpop2

#summary results Lasso
summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop = row_cluster)
#>                 V1               V2                                 names
#> 1        78.515625           78.125 LASSO for subpop1 in target mixedpop2
#> 2 18.1395348837209 19.0697674418605 LASSO for subpop2 in target mixedpop2
#> 3 27.5862068965517 24.1379310344828 LASSO for subpop3 in target mixedpop2

#summary deviance 
summary_deviance(LSOLDA_dat)
#> $allDeviance
#> [1] "0.2876" "0.3078"
#> 
#> $DeviMax
#>         Dfd   Deviance        DEgenes
#> 1         0 -1.922e-15 genes_cluster1
#> 2         2     0.1139 genes_cluster1
#> 3         4     0.2017 genes_cluster1
#> 4         5     0.2801 genes_cluster1
#> 5         6     0.3078 genes_cluster1
#> 6 remaining          1        DEgenes
#> 
#> $LassoGenesMax
#>                                   1                   name
#> (Intercept)            -1.090622240            (Intercept)
#> TNNI1_ENSG00000159173   0.008930123  TNNI1_ENSG00000159173
#> VIM_ENSG00000026025     0.015868582    VIM_ENSG00000026025
#> HHEX_ENSG00000152804   -0.065894344   HHEX_ENSG00000152804
#> TPM1_ENSG00000140416    0.010002491   TPM1_ENSG00000140416
#> TMEM88_ENSG00000167874  0.033797085 TMEM88_ENSG00000167874
#> MYL4_ENSG00000198336    0.003407876   MYL4_ENSG00000198336
```
##R Setting Information

```r
sessionInfo()
#> R version 3.4.2 (2017-09-28)
#> Platform: x86_64-apple-darwin15.6.0 (64-bit)
#> Running under: OS X El Capitan 10.11.6
#> 
#> Matrix products: default
#> BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] scGPS_0.1.0                xlsx_0.5.7                
#>  [3] xlsxjars_0.6.1             rJava_0.9-9               
#>  [5] org.Hs.eg.db_3.5.0         AnnotationDbi_1.40.0      
#>  [7] clusterProfiler_3.6.0      ReactomePA_1.22.0         
#>  [9] DOSE_3.4.0                 bindrcpp_0.2              
#> [11] DESeq_1.30.0               locfit_1.5-9.1            
#> [13] cowplot_0.9.1              cidr_0.1.5                
#> [15] RColorBrewer_1.1-2         reshape2_1.4.3            
#> [17] rmarkdown_1.8              dynamicTreeCut_1.63-1     
#> [19] dplyr_0.7.4                caret_6.0-79              
#> [21] ggplot2_2.2.1              lattice_0.20-35           
#> [23] glmnet_2.0-16              foreach_1.4.4             
#> [25] Matrix_1.2-12              SingleCellExperiment_1.0.0
#> [27] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
#> [29] matrixStats_0.52.2         GenomicRanges_1.30.1      
#> [31] GenomeInfoDb_1.14.0        IRanges_2.12.0            
#> [33] S4Vectors_0.16.0           Biobase_2.38.0            
#> [35] BiocGenerics_0.24.0       
#> 
#> loaded via a namespace (and not attached):
#>   [1] backports_1.1.2           fastmatch_1.1-0          
#>   [3] RcppEigen_0.3.3.3.1       igraph_1.1.2             
#>   [5] plyr_1.8.4                lazyeval_0.2.1           
#>   [7] splines_3.4.2             BiocParallel_1.12.0      
#>   [9] digest_0.6.12             GOSemSim_2.4.0           
#>  [11] htmltools_0.3.6           GO.db_3.5.0              
#>  [13] viridis_0.4.0             checkmate_1.8.5          
#>  [15] magrittr_1.5              memoise_1.1.0            
#>  [17] cluster_2.0.6             sfsmisc_1.1-2            
#>  [19] recipes_0.1.2             fastcluster_1.1.24       
#>  [21] annotate_1.56.1           gower_0.1.2              
#>  [23] RcppParallel_4.4.0        dimRed_0.1.0             
#>  [25] colorspace_1.3-2          rappdirs_0.3.1           
#>  [27] blob_1.1.0                crayon_1.3.4             
#>  [29] RCurl_1.95-4.8            RcppArmadillo_0.8.400.0.0
#>  [31] graph_1.56.0              roxygen2_6.0.1           
#>  [33] genefilter_1.60.0         bindr_0.1                
#>  [35] survival_2.41-3           iterators_1.0.9          
#>  [37] glue_1.2.0                DRR_0.0.2                
#>  [39] gtable_0.2.0              ipred_0.9-6              
#>  [41] zlibbioc_1.24.0           XVector_0.18.0           
#>  [43] graphite_1.24.0           kernlab_0.9-25           
#>  [45] ddalpha_1.3.1.1           prabclus_2.2-6           
#>  [47] DEoptimR_1.0-8            scales_0.5.0             
#>  [49] mvtnorm_1.0-6             DBI_0.7                  
#>  [51] Rcpp_0.12.16              viridisLite_0.2.0        
#>  [53] xtable_1.8-2              reactome.db_1.62.0       
#>  [55] bit_1.1-12                foreign_0.8-69           
#>  [57] mclust_5.4                lava_1.5.1               
#>  [59] prodlim_1.6.1             httr_1.3.1               
#>  [61] fgsea_1.4.0               fpc_2.1-10               
#>  [63] clusterCrit_1.2.7         modeltools_0.2-21        
#>  [65] pkgconfig_2.0.1           XML_3.98-1.9             
#>  [67] flexmix_2.3-14            nnet_7.3-12              
#>  [69] tidyselect_0.2.4          labeling_0.3             
#>  [71] rlang_0.2.0               munsell_0.4.3            
#>  [73] tools_3.4.2               RSQLite_2.0              
#>  [75] ade4_1.7-8                devtools_1.13.4          
#>  [77] broom_0.4.4               evaluate_0.10.1          
#>  [79] stringr_1.3.0             yaml_2.1.18              
#>  [81] bit64_0.9-7               ModelMetrics_1.1.0       
#>  [83] knitr_1.20                robustbase_0.92-8        
#>  [85] purrr_0.2.4               dendextend_1.6.0         
#>  [87] nlme_3.1-131              whisker_0.3-2            
#>  [89] RcppRoll_0.2.2            DO.db_2.9                
#>  [91] xml2_1.1.1                compiler_3.4.2           
#>  [93] rstudioapi_0.7            e1071_1.6-8              
#>  [95] testthat_1.0.2            geneplotter_1.56.0       
#>  [97] tibble_1.4.2              stringi_1.1.7            
#>  [99] desc_1.1.1                trimcluster_0.1-2        
#> [101] commonmark_1.4            psych_1.8.3.3            
#> [103] pillar_1.2.1              data.table_1.10.4-3      
#> [105] bitops_1.0-6              qvalue_2.10.0            
#> [107] R6_2.2.2                  gridExtra_2.3            
#> [109] codetools_0.2-15          MASS_7.3-47              
#> [111] assertthat_0.2.0          CVST_0.2-1               
#> [113] rprojroot_1.3-2           minpack.lm_1.2-1         
#> [115] withr_2.1.2               mnormt_1.5-5             
#> [117] GenomeInfoDbData_0.99.1   diptest_0.75-7           
#> [119] grid_3.4.2                rpart_4.1-11             
#> [121] timeDate_3043.102         tidyr_0.8.0              
#> [123] NbClust_3.0               class_7.3-14             
#> [125] rvcheck_0.0.9             lubridate_1.7.3
#render("/Users/quan.nguyen/Documents/Powell_group_MacQuan/AllCodes/scGPS/vignettes/vignette.Rmd","all")
```

