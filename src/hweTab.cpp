#include <Rcpp.h>
using namespace Rcpp;

// try a table option

// [[Rcpp::export]]
NumericMatrix hweTab(CharacterMatrix af) {
  
  // Subset the matrix in to allele components
  CharacterVector al1(af.nrow()) ;
  CharacterVector al2(af.nrow()) ;
  
  int n0 = af.nrow() ;
  
  for(int i = 0; i < n0; i++){
    al1[i] = af(i,0) ;
    al2[i] = af(i,1) ;
  }
  
  // Remove any NAs
  CharacterVector al1_na = na_omit(al1) ;
  CharacterVector al2_na = na_omit(al2) ;
  
  // Combine the vectors to fine unique alleles
  CharacterVector all_al( al1.size() + al2.size() ) ;
  
  for(int i = 0; i < al1.size(); i++){
    all_al[i] = al1[i] ;
  }
  
  int n1 = al1.size() ;
  int imax = all_al.size() ;
  
  for(int i = n1; i < imax; i++){
    int idx = i - al1.size() ;
    all_al[i] = al2[idx] ;
  }
  
  // get unique alleles
  CharacterVector prelev = na_omit(all_al) ;
  CharacterVector lev = unique(prelev) ;
  
  NumericMatrix out( lev.size(), lev.size() ) ;
  //IntegerMatrix tout( lev.size(), lev.size() ) ;
  rownames(out) = lev ;
  colnames(out) = lev ;
  //rownames(tout) = lev ;
  //colnames(tout) = lev ;
  
  // create a named index vector
  IntegerVector suber = seq_len(lev.size() ) - 1;
  suber.attr("names") = lev ;
  
  int n2 = al1_na.size() ;
  
  for(int i = 0; i < n2; i++){
    String id1 = al1_na[i] ;
    String id2 = al2_na[i] ;
    int idx1 = suber[id1] ;
    int idx2 = suber[id2] ;
    if(id1 == id2){
      out(idx1, idx2) += 1.0 ;
    } else{
      out(idx1, idx2) += 1.0 ;
      out(idx2, idx1) += 1.0 ; 
    }
  }
  
  //NumericVector op = pow((out + tout), 0.5);
  
  
  return out ;
  
}
