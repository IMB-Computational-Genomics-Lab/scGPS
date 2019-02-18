#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Calculate mean
//'
//' @param x integer.
//' @return a scalar value
//' @examples
//' mean_cpp(c(1:10^6))
//'
// [[Rcpp::export]]
double mean_cpp(NumericVector x){
  int n = x.size(); // Size of vector
  double sum = 0; // Sum value
  // For loop, note cpp index shift to 0
  for(int i = 0; i < n; i++){
    // Shorthand for sum = sum + x[i]
    sum += x[i];
  }
  return sum/n; // Obtain and return the Mean
}


//' Calculate variance
//'
//' @param x a vector of gene expression.
//' @param bias degree of freedom
//' @return a variance value
//' @examples
//'var_cpp(c(1:10^6))
//'
// [[Rcpp::export]]

double var_cpp(NumericVector x, bool bias = true){
  // Calculate the mean using C++ function
  double mean = mean_cpp(x);
  double sum = 0;
  int n = x.size();
  for(int i = 0; i < n; i++){
    sum += pow(x[i] - mean, 2.0); // Square
  }
  return sum/(n-bias); // Return variance
}


//' Transpose a matrix
//'
//' @param X  an R matrix (expression matrix)
//' @return a transposed matrix
//' @examples
//' mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=100)
//' tp_mat <- tp_cpp(mat_test)
// [[Rcpp::export]]
arma::mat tp_cpp(const arma::mat X) {
  return arma::trans(X);
}


//' Subset a matrix
//'
//' @param m1in an R matrix (expression matrix)
//' @param rowidx_in a numeric vector of rows to keep
//' @param colidx_in a numeric vector of columns to keep
//' @return a subsetted matrix
//' @examples
//' mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=100)
//' subset_mat <- subset_cpp(mat_test, rowidx_in=c(1:10), colidx_in=c(100:500))
//' dim(subset_mat)
// [[Rcpp::export]]

arma::mat subset_cpp(NumericMatrix m1in, NumericVector rowidx_in, NumericVector colidx_in){
  mat m1 = Rcpp::as<mat>(m1in);
  uvec rowidx = as<uvec>(rowidx_in) - 1;
  uvec colidx = as<uvec>(colidx_in) - 1;
  mat s1 = m1.submat(rowidx, colidx);
  return(s1);
}


//' Principal component analysis
//'
//' @description This function provides significant speed gain if the input matrix
//' is big
//' @param X  an R matrix (expression matrix), rows are genes, columns are cells
//' @return a list with three list pca lists
//' @examples
//' mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=1000)
//' #library(microbenchmark)
//' #microbenchmark(PrinComp_cpp(mat_test), prcomp(mat_test), times=3)
//'
// [[Rcpp::export]]
List PrinComp_cpp(const arma::mat X) {
  arma::mat coeff;
  arma::mat score;
  arma::vec latent;
  arma::princomp(coeff, score, latent, X);
  return List::create(Named("coefficients") = coeff,
                      Named("scores")       =  score,
                      Named("eigenValues")  = latent);

}

/*** R

#mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=1000)
#microbenchmark(PrinComp_cpp(mat_test), prcomp(mat_test), times=3)

*/
