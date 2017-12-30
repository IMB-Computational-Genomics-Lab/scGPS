#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;
// [[Rcpp::depends(RcppParallel)]]

// function for accumulating sum square distance
template <typename InputIterator1>
inline double EuclDist(InputIterator1 begin1, InputIterator1 end1) {

  double rval=0;

  InputIterator1 it1 = begin1;

  while (it1 != end1) {

    double d1 = *it1++;

    rval += d1;
  }
  return std::sqrt(rval);
}


// helper function for squaring the substraction of two numbers
inline double SumSquare(double val1, double val2) {
  return std::pow((val1 - val2),2);
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

        std::vector<double> SSQ(row1.length());

        std::transform(row1.begin(), row1.end(), // input range 1
                       row2.begin(),             // input range 2
                       SSQ.begin(),              // output range
                       SumSquare);                 // function to apply

        double total = EuclDist(SSQ.begin(), SSQ.end());
        rmat(i,j) = total;
        rmat(j,i)= total;
      }
    }
  }
};


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

// generic function for distance
template <typename InputIterator1, typename InputIterator2>
inline double EuclDist(InputIterator1 begin1, InputIterator1 end1,
                       InputIterator2 begin2) {
  // value to return
  double rval=0;

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;

  // for each input item
  while (it1 != end1) {

    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    rval += std::pow((d1-d2),2);
  }
  double test = std::sqrt(rval);
  return test;
}



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
      double total = EuclDist(row1.begin(), row1.end(), row2.begin());
      rmat(i,j) = total;
      rmat(j,i) = total;
      // write to output matrix
    }
  }
  return rmat;
}

//----------------------------------------------------------------------------
// A shorter function for distance--------------------------------------------


/*** R
#sourceCpp("/Users/quan.nguyen/Documents/Powell_group_MacQuan/AllCodes/scGPS/test_rcpp/dist_parallel_cpp_v2.cpp")
mat_test <-matrix(rnorm(100000,10,1),nrow=10000)
  test_cpp_par <- rcpp_parallel_distance(mat_test)
  test_Rdist <- dist((mat_test))
  test_Rdist <-as.matrix(test_Rdist)
  all(test_Rdist == test_cpp_par)
#rcpp_Eucl_distance_NotPar(mat_test)
  */
