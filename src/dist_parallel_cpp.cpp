#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace RcppParallel;
// [[Rcpp::depends(RcppParallel)]]

// generic function for accumulating sum square distance
template <typename InputIterator1, typename InputIterator2>
inline double EuclDistWhole(InputIterator1 begin1, InputIterator1 end1,
                       InputIterator2 begin2) {
  // distance value for two cells (two rows)
  double rval=0;

  // set iterators to beginning of ranges (for the entire row)
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;

  // accumulate, take the value and increment the iterator
  while (it1 != end1) {
  double d1 = *it1++;
  double d2 = *it2++;
  rval += std::pow((d1-d2),2);
  }
  return std::sqrt(rval);
}

struct EclDistance : public Worker {

  const RMatrix<double> mat;
  RMatrix<double> rmat;

  EclDistance(const NumericMatrix mat, NumericMatrix rmat)
    : mat(mat), rmat(rmat) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) { // i from begin to end not 0 to end
      for (std::size_t j = 0; j < i ; j++) { // j from 0 to i not 0 to end

        RMatrix<double>::Row row1 = mat.row(i);
        RMatrix<double>::Row row2= mat.row(j);

        double total = EuclDistWhole(row1.begin(),row1.end(), row2.begin());
        rmat(i,j) = total;
        rmat(j,i)= total;
      }
    }
  }
};

//' Function to calculate Eucledean distance matrix with paralleled C++
//'
//' @param X an R matrix (expression matrix), with cells in rows and genes in columns
//' @export
//'
// [[Rcpp::export]]
NumericMatrix rcpp_parallel_distance(NumericMatrix mat) {

  NumericMatrix rmat(mat.nrow(), mat.nrow());

  // create the worker
  EclDistance EclDistanceFunc(mat, rmat);

  // call it with parallelFor
  parallelFor(0, mat.nrow(), EclDistanceFunc);

  return rmat;
}


// Non parallel version-------------------------------------------------------

//' Function to calculate Eucledean distance matrix without parallelisation
//'
//' @param X an R matrix (expression matrix), with cells in rows and genes in columns
//' @export
//'
// [[Rcpp::export]]

NumericMatrix rcpp_Eucl_distance_NotPar(NumericMatrix mat) {

  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow(), mat.nrow());

  for (int i = 0; i < rmat.nrow(); i++) {
    for (int j = 0; j < i; j++) {

      // rows we will operate on
      NumericMatrix::Row row1 = mat.row(i);
      NumericMatrix::Row row2 = mat.row(j);

      // calculate divergences
      double total = EuclDistWhole(row1.begin(), row1.end(), row2.begin());
      rmat(i,j) = total;
      rmat(j,i) = total;
      // write to output matrix
    }
  }
  return rmat;
}

//----------------------------------------------------------------------------

// sourceCpp("/Users/quan.nguyen/Documents/Powell_group_MacQuan/AllCodes/scGPS/src/dist_parallel_cpp.cpp")

/*** R
  #mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=10000)
  #test_cpp_par <- rcpp_parallel_distance(mat_test)
  #test_cpp_notPar <-rcpp_Eucl_distance_NotPar(mat_test)
  #test_Rdist <- dist((mat_test))
  #test_Rdist <-as.matrix(test_Rdist)
  #all(test_Rdist == test_cpp_par)
  #microbenchmark(rcpp_parallel_distance(mat_test), dist(mat_test), times=3)

  */
