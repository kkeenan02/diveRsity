#include <Rcpp.h>
using namespace Rcpp;

// pairwise heterozygosity calculation

// R-code
//pwDivCalc <- function(x, pw, npops){
//  ht <- matrix(ncol = npops, nrow = npops)
//  hs <- matrix(ncol = npops, nrow = npops)
//  for(i in 1:ncol(pw)){
//    gamma <- sum(sqrt(abs(x[,pw[1,i]] * x[,pw[2,i]])))^-1 
//    f <- gamma * sqrt(x[,pw[1,i]] * x[,pw[2,i]])
//    ht[pw[1,i],pw[2,i]] <- 1 - sum(((f + x[,pw[1,i]])/2)^2)
//    ht[pw[2,i],pw[1,i]] <- 1 - sum(((f + x[,pw[2,i]])/2)^2)
//    hs[pw[1,i],pw[2,i]] <- 1 - sum((f^2 + x[,pw[1,i]]^2)/2)
//    hs[pw[2,i],pw[1,i]] <- 1 - sum((f^2 + x[,pw[2,i]]^2)/2)
//  }
//  ht[is.nan(ht)] <- 0
//  hs[is.nan(hs)] <- 0
//  list(ht = ht, 
//       hs = hs)
//}


// [[Rcpp::export]]
List pwHt(NumericMatrix af, IntegerMatrix pw) {
  int np = af.ncol();
  int n = pw.ncol();
  //int na = af.nrow();
  NumericMatrix ht(np, np);
  NumericMatrix hs(np, np);
  for(int i = 0; i < n; ++i){
    int p1 = pw(0,i);
    int p2 = pw(1,i);
    NumericVector ppreg = sqrt((af(_,p1) * af(_,p2)));
    double preg = sum(ppreg);
    double g = pow(preg, -1.0);
    NumericVector f = g * ppreg;
    // calculate ht's
    NumericVector af1 = (f + af(_,p1))/2.0;
    NumericVector af2 = (f + af(_,p2))/2.0;
    NumericVector aff1 = pow(af1, 2.0);
    NumericVector aff2 = pow(af2, 2.0);
    ht(p2, p1) = 1.0 - sum(aff1);
    ht(p1, p2) = 1.0 - sum(aff2);
    // calulate hs's
    NumericVector af1sq = (pow(f, 2.0) + pow(af(_,p1), 2.0))/2.0;
    NumericVector af2sq = (pow(f, 2.0) + pow(af(_,p2), 2.0))/2.0;
    hs(p2,p1) = 1.0 - sum(af1sq);
    hs(p1,p2) = 1.0 - sum(af2sq);
  }
  return List::create(
    _["ht"] = ht,
    _["hs"] = hs
  );
}
