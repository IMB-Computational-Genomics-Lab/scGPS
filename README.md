# _scGPS_ - Single Cell Global fate Potential of Subpopulations 
<p align="center">
	<img src="man/figures/scGPSlogo.png" width="200px">
</p>

The _scGPS_ package website is available at: https://imb-computational-genomics-lab.github.io/scGPS/index.html 

## _scGPS_ general description
A complete  single cell RNA analysis framework to answer the common question: how to decompose a mixed population into clusters (SCORE) and to analyse the relationship between clusters (scGPS). 

## _scGPS_ workflow

_scGPS_ takes scRNA expression dataset(s) for one or more unknown sample(s) to find subpopulations and relationship between these subpopulations. The input dataset(s) contains mixed, heterogeous cells. _scGPS_ first uses _SCORE_ (or _CORE_ V2.0) to identify homogenous subpopulations. scGPS contains a number of functions to verify the subpopulations identified by _SCORE_ (e.g. functions to compare with results from PCA, tSNE and the imputation method CIDR). scGPS also has options to find gene markers that distinguish a subpopulation from the remaining cells. In the second stage, _scGPS_ applies a machine learning procedure to select optimal gene predictors and to build prediction models that can estimate between-subpopulation transition scores, which are the probability of cells from one subpopulation that can likely transition to the other subpopulation.

 
<p align="center">
	<img src="man/figures/packagePlan.png" width="400px"> <br>
Figure 1. scGPS workflow. Yellow boxes show inputs, and green boxes show main scGPS analysis.  
</p>





