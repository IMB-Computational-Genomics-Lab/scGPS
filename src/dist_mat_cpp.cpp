//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double distvec(NumericVector x, NumericVector y){
  double d = sqrt( sum( pow(x - y, 2) ) );
  return d;
}

// [[Rcpp::export]]
NumericMatrix calcDist(NumericMatrix x){
  int outrows = x.nrow();
  int outcols = x.nrow();
  NumericMatrix out(outrows,outcols);

  for (int i = 0 ; i < outrows - 1; i++){
    for (int j = i + 1  ; j < outcols ; j ++) {
      NumericVector v1 = x.row(i);
      NumericVector v2 = x.row(j-1);
      double d = distvec(v1, v2);
      out(j-1,i) = d;
      out(i,j-1)= d;
    }
  }
  return (out) ;
}


//' Compute Euclidean distance matrix by rows
//'
//' Used in consmx function
//'
//' @param x A numeric matrix.
// [[Rcpp::export]]
arma::mat calcDistArma(const arma::mat & x) {
  unsigned int outrows = x.n_rows, i = 0, j = 0;
  double d;
  mat out = zeros<mat>(outrows, outrows);

  for (i = 0; i < outrows - 1; i++) {
    arma::rowvec v1 = x.row(i);
    for (j = i + 1; j < outrows; j++) {
      d = sqrt(sum(pow(v1 - x.row(j), 2.0)));
      out(j, i) = d;
      out(i, j) = d;
    }
  }
  return out;
}

//Test a large matrix-----------------------------------------------------------
//path_exprs <- '/Users/quan.nguyen/Documents/Powell_group_MacQuan/CardioDiff/DataShared_JP_NP/MoreInfo/'
//Exprs<-readRDS(paste0(path_exprs,'Exprs_DCVLnorm_unlog_minus1_pos.RDS'))
//Exprs_test <-as.matrix(Exprs[1:17000,1:1000])
//system.time(calcDistArma(Exprs_test))
//user;  elapsed; 332.263; 348.759
//system.time(dist(Exprs_test))
//user;  elapsed; 889.812; 934.161
//Done test a large matrix------------------------------------------------------
