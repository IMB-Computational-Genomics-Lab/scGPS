# scGPS
A complete  single cell RNA analysis frame work to answer the question: how to decompose a mixed population into clusters (CORE) and to analyse the relationship between clusters (scGPS). 

## Workflow

scGPS takes scRNA dataset(s) for one more more unknown sample(s) to find subpopulations and relationship between these subpopulations. The input dataset(s) contains mixed, heterogeous cells. scGPS first uses CORE (V2.0) to identify homogenous subpopulations. It contains a number of functions to verify the identified subpopulations (e.g. by comparing with results from PCA, tSNE and the imputation method CIDR). If has options to find gene markers that distinguish a subpopulation from the remaining cells. In the second stage, scGPS applies a machine learning procedure to select optimal gene predictors and to build a prediction model that can estimate transition score, which is the probability of cells from one subpopulation that can likely transition to the other subpopulation.
 
![Alt text](./packagePlan.png?raw=true "scGPS")

