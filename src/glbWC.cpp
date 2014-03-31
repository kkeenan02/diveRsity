#include <Rcpp.h>
using namespace Rcpp;


// C++ version of glbWC function

// [[Rcpp::export]]
List glbWCcpp(IntegerVector hsum, NumericMatrix af, NumericVector indtyp){
  // fix hsum
  Function rownames("rownames");
  Function names("names");
  CharacterVector hsnms = names(hsum);
  CharacterVector nms = rownames(af);
  NumericVector hs(nms.size());
  IntegerVector idx = match(hsnms, nms);
  int n(idx.size());
  for(int i = 0; i < n; ++i){
    int j = idx[i]-1;
    hs[j] = hsum[i];
  }
  // include fixes for unscored loci
  LogicalVector tf = !is_na(indtyp);
  NumericVector indtypfix = indtyp[tf];
  double np = indtyp.size();
  double r = indtypfix.size();
  double na = af.nrow();
  double indT = sum(indtypfix);
  // double sumSq = sum(pow(indtypfix, 2.0));
  double nbar = indT/r;
  NumericVector hbar(na);
  for(int i = 0; i < na; ++i){
    hbar[i] = hs[i]/indT;
  }
  // sd is a C++ function
  double C2 = pow(sd(indtypfix) / (sum(indtypfix)/r), 2.0);
  double nc = nbar * (1.0 - (C2/r));
  NumericMatrix p(na, np);
  NumericVector pbar(na);
  // calculate p
  for(int i = 0; i < np; ++i){
    p(_,i) = af(_,i) * (indtyp[i]*2.0);
  }
  // calculate pbar
  for(int i = 0; i < na; ++i){
    pbar[i] = sum(p(i,_) / (2.0 * indT));
  }
  // calculate pp
  NumericMatrix pp(na, np);
  for(int i = 0; i < np; ++i){
    pp(_,i) = pow(af(_,i) - pbar, 2.0) * indtyp[i];
  }
  NumericVector s2(na);
  for(int i = 0; i < na; ++i){
    s2[i] = sum(pp(i,_))/((np - 1.0) * nbar);
  }
  NumericVector A = pbar * (1.0 - pbar) - (((np - 1.0)/np) * s2);
  double a = nbar*(sum(s2) - (sum(A) - (sum(hbar)/4.0))/(nbar - 1.0)) / nc;
  double b = sum((nbar/(nbar-1.0)) * (A - ((2.0 * nbar - 1.0)/(4.0*nbar))*hbar));
  double c = sum(0.5 * hbar);
  return List::create(
    _["a"] = a,
    _["b"] = b,
    _["c"] = c
  );
}