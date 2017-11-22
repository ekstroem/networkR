#' Hyperlink-induced topic search
#'
#' Hyperlink-induced topic search (HITS) is a link analysis algorithm
#' that is also known as hubs and authorities. It rates nodes by
#' comparing arrows pointing in and out of nodes in an assymetrical graph.
#' 
#' @param adjmatrix an adjacency matrix
#' @param maxiter non-negative integer
#' @return Returns a list of authorities and hubs
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @keywords manip
#' @examples
#'
#' from <- c("A", "A", "A", "B", "C")
#' to <- c("B", "C", "D", "D", "E")
#' hits(adjacency(from, to))
#'
#' @importFrom Matrix sparseMatrix
#' @export hits
hits <- function(adjmatrix, maxiter=25L) {

    ## 
    
    N <- nrow(adjmatrix)

    h <- rep(1, N)
    for (i in 1:maxiter) {
        a <- as.vector(h %*% adjmatrix)
        a <- a / sqrt(sum(a^2))
        h <- as.vector(adjmatrix %*% a)
        h <- h / sqrt(sum(h^2))        
    }

    list(authorities=a, hubs=h)
}

