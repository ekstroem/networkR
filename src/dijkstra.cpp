# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;

using namespace Rcpp;

#include <queue>


#define pii pair<int,int>
#define F first
#define S second
#define mp make_pair
#define pb emplace_back

bool vis[100001];
int dis[100001];
vector<pii> a[100001];



class prioritize {
public: bool operator ()(pii &p1 , pii &p2) {
  return p1.S > p2.S;
}
};

void Dijkstra_cpp(int s, int n) {
  
  for (int i = 0; i <= n; i++) {
    vis[i] = false;
    dis[i] = INT_MAX;
  }
  priority_queue<pii, vector<pii>, prioritize> pq;
  pq.push(mp(s, dis[s] = 0));
  while (!pq.empty()) {
    pii cur = pq.top(); pq.pop();
    int cv = cur.F, cw = cur.S;
    if (vis[cv]) continue;
    vis[cv] = true;
    for (pii x : a[cv]) {
      if (!vis[x.F] && (cw + x.S) < dis[x.F]) {
	pq.push(mp(x.F, dis[x.F] = cw + x.S));
      }
    }
  }
}



typedef pair<double, double> ii;   //The Only Macros I use...usually
typedef vector<double> vi;
vector<vector< ii > > graph;
vi dist;
 
void dijkstra2(int s)
{
  //  int n;
    dist.assign( graph.size() , 100000 ); dist[s] = 0;
    priority_queue< ii, vector<ii>, greater<ii> > pq; pq.push({0, s});
 
    while(!pq.empty())
    {
        ii v = pq.top(); pq.pop();
 
        for(int i = 0 ; v.first < dist[v.second] && i < graph[v.second].size(); i++)
        {
            ii u = graph[v.second][i];
 
            if(dist[v.second] + u.second < dist[u.first])
                pq.push({dist[u.first] = dist[v.second] + u.second, u.first});
        }
 
    }
}




struct edge { int to, length; };
    
int DIJkstra(const vector< vector<edge> > &graph, int source, int target) {
    vector<int> min_distance( graph.size(), INT_MAX );
    min_distance[ source ] = 0;
    set< pair<int,int> > active_vertices;
    active_vertices.insert( {0,source} );
        
    while (!active_vertices.empty()) {
        int where = active_vertices.begin()->second;
        if (where == target) return min_distance[where];
        active_vertices.erase( active_vertices.begin() );
        for (auto ed : graph[where]) 
            if (min_distance[ed.to] > min_distance[where] + ed.length) {
                active_vertices.erase( { min_distance[ed.to], ed.to } );
                min_distance[ed.to] = min_distance[where] + ed.length;
                active_vertices.insert( { min_distance[ed.to], ed.to } );
            }
    }
    return INT_MAX;
}




IntegerVector distances44(const IntegerVector& from, const IntegerVector& to, const IntegerVector& weight) {

  long nVertices = 1000;
  long nEdges = from.size();

  for (int i=0; i<nEdges; i++) {
    a[from[i]].push_back(make_pair(to[i],weight[i]));
    //    a[to[i]].push_back(make_pair(from[i],weight[i]));  // If symmetric
  }
  int s=1;

  for (int j=1; j<=nVertices; j++) {
  Dijkstra_cpp(j, nVertices);
  /* for (int i = 1; i <= nVertices; i++) {
    /*    if (dis[i] != INT_MAX) {
      Rcout << dis[i] << " ";
    } else {
      Rcout << "-1 ";
    }
  }
  */
   
  }
  return wrap(0);
}


//' ' ' ' ' ' '' '  ' @export
//' ' ' ' ' ' ' ' '' ' ' [[Rcpp::export]]
IntegerVector Distances(const IntegerVector& from, const IntegerVector& to, const NumericVector& weight) {

  long nVertices = 25000;
  long nEdges = from.size();

  for (int i=0; i<nEdges; i++) {
    a[from[i]].push_back(make_pair(to[i],weight[i]));
    //    a[to[i]].push_back(make_pair(from[i],weight[i]));  // If symmetric
  }
  int s=1;

  for (int j=1; j<=nVertices; j++) {
  Dijkstra_cpp(j, nVertices);
  for (int i = 1; i <= nVertices; i++) {
    /*    if (dis[i] != INT_MAX) {
      Rcout << dis[i] << " ";
    } else {
      Rcout << "-1 ";
    }
    */
  }
  }
  return wrap(0);
}
