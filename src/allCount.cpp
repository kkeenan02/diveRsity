#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
int allCount(CharacterMatrix x) {
  CharacterVector al1 = unique(x(_,0)) ;
  CharacterVector al2 = unique(x(_,1)) ;
  CharacterVector out( al1.size() + al2.size() ) ;
  for(int i = 0; i < al1.size(); i++){
    out[i] = al1[i] ;
  }
  int n = al1.size() ;
  int imax = out.size() ;
  for(int i = n; i < imax; i++){
    int idx = i - al1.size() ;
    out[i] = al2[idx] ;
  }
  CharacterVector unq = unique(out) ;
  LogicalVector idx = is_na(unq) ;
  CharacterVector newUnq = unq[!idx] ;
  return newUnq.size() ;
}
