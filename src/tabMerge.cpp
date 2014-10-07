#include <Rcpp.h>
using namespace Rcpp;

// A c++ implementation of the tabMerge function

// [[Rcpp::export]]
NumericVector tabMerge(List ip) {
    
    RCPP_UNORDERED_MAP<std::string,double> out ;
    
    int n = ip.size() ;
    
    for(int i = 0; i < n; i++){
      NumericVector x = ip[i] ;
      CharacterVector names = x.attr("names") ;
      int m = x.size() ;
      
      for(int j = 0; j < m; j++){
        String name = names[j] ;
        out[ name ] += x[j] ;
      }
    }
  return wrap(out) ;
}
