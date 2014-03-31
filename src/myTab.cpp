#include <Rcpp.h>
using namespace Rcpp;

// Rcpp table function is faster than base R
// Working! (3x speedup)

// [[Rcpp::export]]
NumericVector myTab(CharacterVector x){
  IntegerVector tab = table(x);
  double n = sum(tab);
  NumericVector out(tab.size());
  for(int i = 0; i < tab.size(); ++i){
    out[i] = tab[i] / n;
  }
  return out;
}