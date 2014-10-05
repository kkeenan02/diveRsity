#include <Rcpp.h>
using namespace Rcpp;

// Function for calculating harmonic sample sizes
// R version
//   diffCalcHarm <- function(idt, pw){
//     ps <- apply(pw, 2, function(x){
//      return((0.5*((idt[x[1]]^-1) + idt[x[2]]^-1))^-1)
//    })
//    return(ps)
//  }
// [[Rcpp::export]]
NumericVector diffCalcHarm(NumericVector idt, NumericMatrix pw){
  int n = pw.ncol();
  NumericVector out(n);
  for(int i = 0; i < n; ++i){
    out[i] = pow(0.5*(pow(idt[pw(0,i)], -1.0)+pow(idt[pw(1,i)], -1.0)), -1.0);
  }
  return out;
}