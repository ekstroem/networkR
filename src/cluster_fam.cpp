// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_set>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

//
// This function replaces subtracts 1 from all indices after it has replaced NAs and 0s with replacement.
//
//
IntegerVector replace_nazero(IntegerVector x, int replacement) {
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


//' Construct family id vector from pedigree trio information
//'
//' @description Create a vector of length n, giving the family id of
//' each subject.  If the pedigree is totally connected, then everyone
//' will end up in tree 1, otherwise the tree numbers represent the
//' disconnected subfamilies.  Singleton subjects each have unique
//' family numbers.
//'
//' No check is done to ensure that the id, fid, and mid actually refer to proper family structure.
//' References to ids in the fid and mid arguments that are not part of the id vector are considered founders and are thus replaced by NA or 0s after being used to group full and half-sibs.
//' @param id Numeric vector of ids
//' @param fid Numeric vector of ids of the father. This should be NA or 0 for a founder.
//' @param mid Numeric vector of ids of the mother. This should be NA or 0 for a founder.
//' @return Returns an integer vector giving the family index of each individual
//' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
//' @keywords manip
//' @examples
//'
//' id <- 1:11
//' fid <- c(NA, NA, 1, 1, NA, 23, 45, 5, 5, 7, NA)
//' mid <- c(NA, NA, 2, 2, 65, NA, 46, 6, 6, 6, 0)
//' make_family_id(id, fid, mid)
//'
//' @export 
// [[Rcpp::export]]
IntegerVector make_family_id(const NumericVector& id, const NumericVector& fid, const NumericVector& mid) {

  unsigned int N = id.size();

  // Sanity checks
  if (fid.size() != N)
    stop("id, fid, and mid must have the same length");
  if (mid.size() != N)
    stop("id, fid, and mid must have the same length");

      
  IntegerVector family = seq_len(N+1)-1;
  IntegerVector newid = seq_len(N+1);

  // Combine full and half siblings into same families before starting
  // This will ensure that some trios are grouped together
  // Using the fact tha family ids are sorted at this point
  for (unsigned int i=0; i<N-1; i++) {
    
    for (unsigned int j=i+1; j<N; j++) {
      if (fid[i] == fid[j] & fid[i]>0 & !NumericVector::is_na( fid[i] )) { // If same father with an actual id then same family
	family[j] = family[i];
      }
      if (mid[i] == mid[j] & mid[i]>0 & !NumericVector::is_na( mid[i] )) { // If same mother with an actual id then same family      
	family[j] = family[i];
      }
    }
  }

  IntegerVector fatherid = replace_nazero(match(fid, id), N+1);
  IntegerVector motherid = replace_nazero(match(mid, id), N+1);

  // We need to run through all individuals over and onver since the ordering might be messed up in the input vectors
  for (unsigned int i=0; i<N; i++) {

    // First we set the familyid of individual j to be the smallest of the parents and the id itself
    for (unsigned int j=0; j<N; j++) {
      newid[j] = std::min(std::min(family[j], family[motherid[j]]), family[fatherid[j]]);
      //      Rcout << "  Changed newid for position " << j << " to " << newid[j] << std::endl;
    }


    // Could be sped up
    for (unsigned int j=0; j<N; j++) {
      // Fix the mothers and fathers
      //      Rcout << "    Change newid for position " << motherid[j] << " from " << newid[motherid[j]] ;      
      newid[motherid[j]] = std::min(newid[j], newid[motherid[j]]);
      //      Rcout << " to " << newid[motherid[j]] << std::endl;
      newid[N] = N;
      newid[fatherid[j]] = std::min(newid[j], newid[fatherid[j]]);
      //      Rcout << "    Changed newid for position " << fatherid[j] << " to " << newid[fatherid[j]] << std::endl;      
      newid[N] = N;      
    }
    //    Rf_PrintValue(newid);
    // Rf_PrintValue(family);
    if (is_true(all((newid==family)==TRUE)))
      break;
    else if (i < N)
      std::copy( newid.begin(), newid.end(), family.begin() ) ;
    //      family = newid;
  }

  // Insert check here?
  IntegerVector famnames = sort_unique(family);
  IntegerVector idx = seq_len(famnames.size())-1;
  family = family[family<N];
  
  return(match(family, famnames));
}

int unique_set_native(Rcpp::CharacterVector x) {
  std::unordered_set<SEXP> tab(x.begin(), x.end());
  return tab.size();
}
