#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericMatrix genos2mat(NumericMatrix mat, IntegerVector ip, NumericVector na){
  int n = ip.size();
  NumericMatrix out(mat.nrow(), mat.ncol());
  LogicalVector tst = is_na(ip);
  for(int i = 0; i < n; ++i){
    if(tst[i] == false){
      out(i,ip[i]) = 0.5;
    } else {
      out(i,_) = na;
    }
  }
  return out;
}
