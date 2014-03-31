#include <Rcpp.h>
using namespace Rcpp;

// global het function.

// R version
//  # locus global stats (replace with Rcpp function)
//  varFunc <- function(af, nHarm = NULL){
//    ht <- 1 - sum(rowMeans(af, na.rm = TRUE)^2)
//    hs <- 1 - sum(rowMeans(af^2, na.rm = TRUE))
//    if(!is.null(nHarm)){
//      hsEst <- hs * ((2*nHarm)/((2*nHarm)-1))
//      htEst <- ht + (hsEst/(2*nHarm*ncol(af)))
//      list(ht = ht, hs = hs, htEst = htEst, hsEst = hsEst)
//    } else {
//      list(ht = ht, hs = hs)
//    }
//  }

// C++ version

// [[Rcpp::export]]
List varFunc(NumericMatrix af, double sHarm){
  int np = af.ncol();
  int na = af.nrow();
  NumericVector afsum(na);
  NumericVector afsumSq(na);
  for(int i = 0; i < na; ++i){
    afsum[i] = sum(af(i,_))/np;
    afsumSq[i] = sum(pow(af(i,_), 2.0))/np;
  }
  double ht1 = 1 - sum(pow(afsum, 2.0));
  double hs1 = 1 - sum(afsumSq);
  double hs = hs1 * ((2*sHarm)/((2*sHarm)-1));
  double ht = ht1 + (hs/(2*sHarm*np));
  return List::create(
    _["htEst"] = ht,
    _["hsEst"] = hs
  );
}