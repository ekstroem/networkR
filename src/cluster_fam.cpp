#include <Rcpp.h>
using namespace Rcpp;

IntegerVector replace_na(IntegerVector x, int replacement) {
  long N = x.size() ;

  IntegerVector out(N+1);
  
  for( long i=0; i<N; i++){
    if( IntegerVector::is_na( x[i] ) )
      out[i] = replacement-1;
    else {
      out[i] = x[i]-1;
    }
  }
  out[N]=N;
  return(out);
}


//' Cluster families
//'
//' @description Computes a vector of groupings in families based on id, father id, and mother id. No check is done to ensure that the id, fid, and mid actually refere to a proper family structure. References to ids in the fid and mid arguments that are not part of the id vector are considered founders.
//' @param id Numeric vector of ids
//' @param fid Numeric vector of ids of the father
//' @param mid Numeric vector of ids of the mother
//' @return Returns an integer vector giving the family index
//' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
//' @keywords manip
//' @examples
//'
//' id <- 1:11
//' fid <- c(NA, NA, 1, 1, NA, 23, 45, 5, 5, 7, NA)
//' mid <- c(NA, NA, 2, 2, 65, NA, 46, 6, 6, 6, NA)
//' cluster_families(id, fid, mid)
//'
//' @export 
// [[Rcpp::export]]
IntegerVector cluster_families(NumericVector id, NumericVector fid, NumericVector mid) {

  long N = id.size();

  IntegerVector fatherid = replace_na(match(fid, id), N+1);
  IntegerVector motherid = replace_na(match(mid, id), N+1);
  IntegerVector family = seq_len(N+1)-1;
  IntegerVector newid = seq_len(N);
  
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      newid[j] = std::min(std::min(family[j], family[motherid[j]]), family[fatherid[j]]);
      // Fix the mother and fathes
      newid[motherid[j]] = newid[j];
      newid[N] = N;
      newid[fatherid[j]] = newid[j];
      newid[N] = N;      
    }
    // Rf_PrintValue(newid);
    if (is_true(all((newid==family)==TRUE)))
      break;
    else if (i < N) 
      family = newid;
  }

  // Insert check here?

  IntegerVector famnames = unique(family);
  std::sort(famnames.begin(), famnames.end());
 
  return(match(family, famnames));
}
