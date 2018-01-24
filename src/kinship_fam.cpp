// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_set>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;




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
//' ## cluster_families(id, fid, mid)
//'
//' @export 
// [[Rcpp::export]]
arma::sp_mat kinship(IntegerVector id, IntegerVector fid, IntegerVector mid) {
  long N = id.size();

  arma::sp_mat res(N,N);

  // Get the clusters
  // 1) For each cluster

  
  
  //IntegerVector clusters = make_family_id(id, fid, mid);
  //  int nclusters = unique_set_native(as<CharacterVector>(clusters));

  //  for (int cluster=1; cluster<=nclusters; cluster++) {
    // 2 Extract the submatrix
    // 3 Computer the kinship coefficient
    //   arma::sp_mat oyoy  = compute_kinship(id[nclusters==cluster], fid[nclusters==cluster], mid[nclusters==cluster]);
    // 4 Insert the kinship coefficient
  //  }

  return(res);
}


