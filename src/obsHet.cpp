#include <Rcpp.h>
using namespace Rcpp;

// Observed heterozygosity calculation

// Kevin Keenan, 2015


// [[Rcpp::export]]
double obsHet(CharacterMatrix in_mat) {
  
  //CharacterVector new_mat = in_mat(_,1) ;
  
  int n = in_mat.nrow() ;
  double total = 0 ;
  double ns = 0 ;
  
  for(int i = 0; i < n; i++){
    
    CharacterVector tst = in_mat(i,_) ;
    LogicalVector id = is_na(tst) ;
    
    if(tst[0] != tst[1]){
      total += 1 ;
    }
    
    if(!id[0]) {
      ns += 1 ;
    }
    
  }
  
  return total/ns ;
}