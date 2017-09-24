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
//' @description Computes the correlation matrix distance between two correlation matrices
//' @param id Numeric vector of ids
//' @param fid Numeric vector of ids of the father
//' @param mid Numeric vector of ids of the father
//' @return Returns a vector giving the xxxxx
//' zero for equal correlation matrices and unity if they differ to a maximum extent.
//' @author Claus Ekstrom \email{claus@@rprimer.dk}
//' @keywords manip
//' @examples
//'
//' x <- 1:5
//'
// [[Rcpp::export]]
IntegerVector cluster_families(NumericVector id, NumericVector fid, NumericVector mid) {

  long N = id.size();

  IntegerVector fatherid = replace_na(match(fid, id), N+1);
  IntegerVector motherid = replace_na(match(mid, id), N+1);
  IntegerVector family = seq_len(N+1)-1;

  for (int i=0; i<N+1; i++) {
    for (int j=0; j<N; j++) {
      family[j] = std::min(std::min(family[j], family[motherid[j]]), family[fatherid[j]]);
      // Fix the mother and fathes
      family[motherid[j]] = family[j];
      family[N] = N;
      family[fatherid[j]] = family[j];
      family[N] = N;      
    }
  }

  return(family);
}
