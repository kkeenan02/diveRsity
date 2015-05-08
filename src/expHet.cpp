#include <Rcpp.h>
using namespace Rcpp;

// Calculate expected heterozygosity from allele frequency data
// Kevin Keenan, 2015

// [[Rcpp::export]]
NumericVector expHet(NumericMatrix af) {
  
  int n = af.ncol() ;
  NumericVector out(n) ;
  
  for(int i = 0; i < n; i++){
    out[i] = 1 - sum(pow(af(_,i), 2.0)) ;
  }
  
  return out ;
}