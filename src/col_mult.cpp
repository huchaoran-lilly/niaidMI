#include <Rcpp.h>
using namespace Rcpp;

// hi Multiplies each col of a matrix my the mating element in the v
// [[Rcpp::export]]
NumericMatrix col_mult( NumericMatrix M , NumericVector v ){
  NumericMatrix out(M) ;
  for (int j = 0; j < M.ncol(); j++) {
    for (int i = 0; i < M.nrow(); i++) {
      out(i,j) = M(i,j) * v[j];
    }
  }
  return out ;
}
