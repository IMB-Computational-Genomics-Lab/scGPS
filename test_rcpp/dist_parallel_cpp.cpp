#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]


double dist1 (NumericVector x, NumericVector y){
  int n = y.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += pow(x(i)-y(i),2.0);
  }
  total = sqrt(total);
  return total;
}

struct EclDistance : public Worker {
  // input matrix to read from
  const RMatrix<double> mat;
  // output matrix to write to
  RMatrix<double> rmat;
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  EclDistance(const NumericMatrix mat, NumericMatrix rmat)
    : mat(mat), rmat(rmat) {}
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end -1; i++) {
      for (std::size_t j = i+1; j < end ; j++) {
        // rows we will operate on
        RMatrix<double>::Row row1 = mat.row(i);
        RMatrix<double>::Row row2= mat.row(j);


        double total = dist1(row1.begin(), row1.end(), avg.begin());
        rmat(j-1,i) = total;
        rmat(i,j-1)= total;
      }
    }
  }
};

  // Now that we have the JsDistance function object we can pass it to parallelFor,
  // specifying an iteration range based on the number of rows in the input matrix:

  // [[Rcpp::export]]
  NumericMatrix rcpp_parallel_distance(NumericMatrix mat) {

    // allocate the matrix we will return
    NumericMatrix rmat(mat.nrow(), mat.nrow());

    // create the worker
    EclDistance EclDistanceFunction(mat, rmat);

    // call it with parallelFor
    parallelFor(0, mat.nrow(), EclDistanceFunction); //Quan: parallel for each row (not between rows)

    return rmat;
  }



#include <RcppArmadillo.h>

  using namespace arma;

  //' Compute Euclidean distance matrix by rows
  //'
  //' Used in consmx function
  //'
  //' @param x A numeric matrix.
  // [[Rcpp::export]]
  arma::mat ED1(const arma::mat & x) {
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
  *//


