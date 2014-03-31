#include <Rcpp.h>
using namespace Rcpp;

// pwHetCalc function.

// [[Rcpp::export]]

List pwHCalc(NumericMatrix af, NumericVector sHarm, IntegerMatrix pw){
  int n = pw.ncol();
  NumericVector ht(n);
  NumericVector hs(n);
  
  for(int i = 0; i < n; ++i){
    ht[i] = 1.0 - sum(pow(((af(_, pw(0,i)) + af(_,pw(1,i)))/2.0), 2.0));
    hs[i] = 1.0 - sum((pow(af(_, pw(0,i)), 2.0) + pow(af(_, pw(1,i)), 2.0))/2.0);
    hs[i] = hs[i] * ((2.0 * sHarm[i]) / (2.0 * sHarm[i] - 1));
    ht[i] = ht[i] + (hs[i]/(4 * sHarm[i]));
  }
  return List::create(
    _["hsEst"] = hs,
    _["htEst"] = ht
  );
}
