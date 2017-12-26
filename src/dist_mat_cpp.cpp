#include <Rcpp.h>
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




