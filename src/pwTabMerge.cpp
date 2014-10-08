#include <Rcpp.h>
using namespace Rcpp;

// A c++ implementation of the tabMerge function

// [[Rcpp::export]]
List pwTabMerge(List hsum, NumericMatrix pw) {
    
    int n = pw.ncol() ;
    
    List out(n) ;
    
    for(int k = 0; k < n; k++){
      
      RCPP_UNORDERED_MAP<std::string,double> preout ;
      
      List ip(2) ;
      ip[0] = hsum[pw(0,k)] ;
      ip[1] = hsum[pw(1,k)] ;
      
      for(int i = 0; i < 2; i++){
        NumericVector x = ip[i] ;
        CharacterVector names = x.attr("names") ;
        
        int m = x.size() ;
        
        for(int j = 0; j < m; j++){
          String name = names[j] ;
          preout[ name ] += x[j] ;
        }
      }
      out[k] = preout ;
    }
  return wrap(out) ;
}