#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Calculate mean
//'
//' @param N  integer.
//' @param thin integer
//' @export
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
//' @export
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
//' @export
//'
// [[Rcpp::export]]
arma::mat tp_cpp(const arma::mat X) {
  return arma::trans(X);
}


//' Subset a matrix
//'
//' @param X an R matrix (expression matrix)
//' @export
//'
// [[Rcpp::export]]

arma::mat subset_cpp(NumericMatrix m1in, NumericVector rowidx_in,NumericVector colidx_in){
  mat m1 = Rcpp::as<mat>(m1in);
  uvec rowidx = as<uvec>(rowidx_in) - 1;
  uvec colidx = as<uvec>(colidx_in) - 1;
  mat s1 = m1.submat(rowidx, colidx);
  return(s1);
}



