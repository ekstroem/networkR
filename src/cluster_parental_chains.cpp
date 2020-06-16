// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_set>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

//
// This function replaces 
//
//
IntegerVector replace_nazero2(IntegerVector x, int replacement) {
  int N = x.size() ;

  IntegerVector out(N+1);
  
  for(int i=0; i<N; i++){
    if( IntegerVector::is_na( x[i] ) )
      out[i] = replacement-1;
    else if (x[i] == 0) {
      out[i] = replacement-1;
    } else {
      out[i] = x[i]-1;
    }
  }
  out[N]=N;
  return(out);
}


//' Construct parental chain id vector from pedigree trio information
//'
//' @description Create a vector of length n, giving the id of
//' parental chains.  If the pedigree is totally connected, then everyone
//' will end up in tree 1, otherwise the tree numbers represent the
//' disconnected subfamilies.  Singleton subjects each have unique
//' family numbers.
//'
//' No check is done to ensure that the id, fid, and mid actually refer to proper family structure.
//' References to ids in the fid and mid arguments that are not part of the id vector are considered founders.
//' @param id Numeric vector of ids
//' @param fid Numeric vector of ids of the father. This should be NA or 0 for a founder.
//' @param mid Numeric vector of ids of the mother. This should be NA or 0 for a founder.
//' @return Returns an integer vector giving the family index
//' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
//' @keywords manip
//' @examples
//'
//' id <- 1:11
//' fid <- c(0,0,1,0,0,4,0,0,3,7,7)
//' mid <- c(0,0,2,0,0,5,0,0,6,6,8)
//' make_parental_chain(id, fid, mid)
//'
//' @export 
// [[Rcpp::export]]
IntegerVector make_parental_chain(const NumericVector& id, const NumericVector& fid, const NumericVector& mid) {

  unsigned int N = id.size();

  // Sanity checks
  if (fid.size() != N)
    stop("id, fid, and mid must have the same length");
  if (mid.size() != N)
    stop("id, fid, and mid must have the same length");

      
  IntegerVector fatherid = replace_nazero2(match(fid, id), N+1);
  IntegerVector motherid = replace_nazero2(match(mid, id), N+1);
  IntegerVector family = seq_len(N);
  //  IntegerVector newid = seq_len(N+1);
 
  for (unsigned int i=0; i<N; i++) {
    for (unsigned int j=i; j<N; j++) {

          if (fatherid[j]<N+1 & motherid[j]<N+1) {

	    family[fatherid[j]] = std::min(family[motherid[j]], family[fatherid[j]]);
	    family[motherid[j]] = family[fatherid[j]];
	    
	  }	  
    }
  }


  // Insert check here?


  IntegerVector famnames = sort_unique(family);
  IntegerVector idx = seq_len(famnames.size())-1;
  family = family[family<N];
  
  return(match(family, famnames));



  //  return(family);

}

/*

int unique_set_native(Rcpp::CharacterVector x) {
  std::unordered_set<SEXP> tab(x.begin(), x.end());
  return tab.size();
}

*/
