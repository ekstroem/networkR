// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <unordered_set>
// #include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;


template <int RTYPE>
IntegerVector order_impl(const Vector<RTYPE>& x, bool desc) {
    auto n = x.size();
    IntegerVector idx = no_init(n);
    std::iota(idx.begin(), idx.end(), static_cast<size_t>(1));
    if (desc) {
        auto comparator = [&x](size_t a, size_t b){ return x[a - 1] > x[b - 1]; };
        std::stable_sort(idx.begin(), idx.end(), comparator);
    } else {
        auto comparator = [&x](size_t a, size_t b){ return x[a - 1] < x[b - 1]; };
        std::stable_sort(idx.begin(), idx.end(), comparator);
        // simulate na.last
        size_t nas = 0;
        for (size_t i = 0; i < n; ++i, ++nas)
            if (!Vector<RTYPE>::is_na(x[idx[i] - 1])) break;
        std::rotate(idx.begin(), idx.begin() + nas, idx.end());
    }
    return idx;
}





//' Cluster network
//'
//' @description Identifies cluster in a symmetric, weighted network based on an edgelist. The edgelist should contain positive integers for the node numbers. Sequential node numbers will speed up computations.
//' @param from Integer vector of nodes where the link is from
//' @param to Integer vector of nodes where the link is to
//' @param weight Numeric vector of link weights
//' @param threshold A double indicating the minimum link weight that should be considered. Only weights equal to or greater than the threshold will be considered in the clustering
//' @return Returns an integer vector giving the cluster index for each node
//' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
//' @keywords manip
//' @examples
//'
//' aaa <- structure(c(1, 6, 5, 4, 5, 11, 1, 5, 3, 13, 16, 15, 18, 6, 7, 8, 10, 12, 14, 15, 15, 16, 17, 17, 18, 20), .Dim = c(13L, 2L))
//' cluster_network(aaa[,1], aaa[,2], weight=rep(1, nrow(aaa))) 
//'
//' @export 
// [[Rcpp::export]]
IntegerVector cluster_network(const IntegerVector& from, const IntegerVector& to, const NumericVector& weight, double threshold=0) {


  // Start by checking lengths - they should be identical
  if (from.size() != to.size())
    stop("from, to, and weight must all have the same length");
  if (from.size() != weight.size())
    stop("from, to, and weight must all have the same length");

  
  // Remove weights under threshold ?!?
  
  long Nedges = from.size();

  // Swap directions
  /*  
  long keep;
  for (long i=0; i<Nedges; i++) {
    // Skip if weight is under threshold
    if (from[i] > to[i]) {
      keep = from[i];
      from[i] = to[i];
      to[i] = keep;
    }
  }
  */
  
  //  long N = max(200); 

  long N = std::max(max(to), max(from));
  IntegerVector order = order_impl(from, false)-1;

  IntegerVector clusterid = seq_len(N);
  bool change;
  for (long j=0; j<N; j++) {
    change = false;
    
    for (long i=0; i<Nedges; i++) {
      // Skip if weight is under threshold
      if (weight[order[i]]<threshold)
	continue;
      
      //      Rcout << i << " for " << from[order[i]]-1 << "  " << to[order[i]]-1 << arma::endl;
      
      if (clusterid[from[order[i]]-1] < clusterid[to[order[i]]-1]) {
	//	Rcout << "   Comp : " << clusterid[to[order[i]]-1]  << " is set to " <<  clusterid[from[order[i]]-1] << arma::endl;
	clusterid[to[order[i]]-1] = clusterid[from[order[i]]-1];
	change=true;
      } else if (clusterid[from[order[i]]-1] > clusterid[to[order[i]]-1]) {
	clusterid[from[order[i]]-1] = clusterid[to[order[i]]-1];
	change=true;

      }
      
    }
    
    if (!change) {
      break;
    }
	
  }

  // Now we should just clean up the cluster names
  IntegerVector clusternames = sort_unique(clusterid);
  return(match(clusterid, clusternames));
}

