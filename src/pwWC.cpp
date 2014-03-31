#include <Rcpp.h>
using namespace Rcpp;

// C++ version of glbWC function

// [[Rcpp::export]]
List pwWCcpp(List hsum1, NumericMatrix af1, NumericVector indtyp1,
            IntegerMatrix pw){
  int n = pw.ncol();
  NumericVector a(n);
  NumericVector b(n);
  NumericVector c(n);
  Function rownames("rownames");
  Function names("names");
  CharacterVector nms = rownames(af1);
  for(int i = 0; i < n; ++i){
    NumericVector indtyp(2);
    indtyp[0] = indtyp1[pw(0,i)];
    indtyp[1] = indtyp1[pw(1,i)];
    LogicalVector tf = is_na(indtyp);
    if(tf[0] == true || tf[1] == true){
      a[i] = NA_REAL;
      b[i] = NA_REAL;
      c[i] = NA_REAL;
    } else {
      NumericVector hsum11 = hsum1[i];
      CharacterVector hsnms = names(hsum1[i]);
      NumericMatrix af(af1.nrow(), 2);
      af(_,0) = af1(_,pw(0,i));
      af(_,1) = af1(_,pw(1,i));
      NumericVector hsum(nms.size());
      IntegerVector idx = match(hsnms, nms);
      int n1(idx.size());
      for(int j = 0; j < n1; ++j){
        int k = idx[j]-1;
        hsum[k] = hsum11[j];
      }
      double r = af.ncol();
      double na = af.nrow();
      double indT = sum(indtyp);
      double nbar = indT/r;
      NumericVector hbar(na);
      for(int j = 0; j < na; ++j){
        hbar[j] = hsum[j]/(indT + 0.0);
      }
      // sd is a C++ function
      double C2 = pow(sd(indtyp) / (sum(indtyp)/r), 2.0);
      double nc = nbar * (1.0 - (C2/r));
      NumericMatrix p(na, 2);
      NumericVector pbar(na);
      // calculate p
      for(int j = 0; j < r; ++j){
        p(_,j) = af(_,j) * (indtyp[j]*2.0);
      }
      // calculate pbar
      for(int j = 0; j < na; ++j){
        pbar[j] = sum(p(j,_) / (2.0 * indT));
      }
      // calculate pp
      NumericMatrix pp(na, 2);
      for(int j = 0; j < r; ++j){
        pp(_,j) = pow(af(_,j) - pbar, 2.0) * indtyp[j];
      }
      NumericVector s2(na);
      for(int j = 0; j < na; ++j){
        s2[j] = sum(pp(j,_))/((r - 1.0) * nbar);
      }
      NumericVector A = pbar * (1.0 - pbar) - (((r - 1.0)/r) * s2);
      a[i] = nbar*(sum(s2) - (sum(A) - (sum(hbar)/4.0))/(nbar - 1.0)) / nc;
      b[i] = sum((nbar/(nbar-1.0)) * (A - ((2.0 * nbar - 1.0)/(4.0*nbar))*hbar));
      c[i] = sum(0.5 * hbar); 
    }
  }
  return List::create(
    _["a"] = a,
    _["b"] = b,
    _["c"] = c
  );
}