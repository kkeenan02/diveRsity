#include <Rcpp.h>
using namespace Rcpp;

// Calculate af from genos



// [[Rcpp::export]]
double bsHetCalc(CharacterMatrix af) {
  
  CharacterVector all = wrap(na_omit(af)) ;
  
  double n = all.size() + 0.0 ;
  
  IntegerVector al = table(all) ;
  
  NumericVector al2 = pow(pow(al, 1.0) / n, 2.0) ;
  
  return 1.0 - sum(al2) ;
}