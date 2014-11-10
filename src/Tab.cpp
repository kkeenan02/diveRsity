#include <Rcpp.h>
using namespace Rcpp;

// cpp table function

// [[Rcpp::export]]
IntegerVector Tab(CharacterVector x) {
  return table(x);
}
