#' Hyperlink-induced topic search
#'
#' Hyperlink-induced topic search (HITS) is a link analysis algorithm
#' that is also known as hubs and authorities. It rates nodes by
#' comparing arrows pointing in and out of nodes in an asymmetrical graph.
#'
#' Hubs are nodes with a lot of arrows pointing out while authorities are node with a lot of arrows pointing in.
#' 
#' @param adjmatrix an adjacency matrix
#' @param maxiter non-negative integer
#' @param tol positive numeric value to be used as tolerance threshold for convergence
#' @return Returns a list with three elements: authorities (a vector) of  and hubs (a vector), and number of iterations used.
#' @references Kleinberg, Jon (1999). "Authoritative sources in a hyperlinked environment" (PDF). Journal of the ACM. 46 (5): 604–632. doi:10.1145/324133.324140 
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
#' @keywords manip
#' @examples
#'
#' from <- c("A", "A", "A", "B", "C")
#' to <- c("B", "C", "D", "D", "E")
#' hits(adjacency(from, to))
#'
#' @importFrom Matrix sparseMatrix
#' @export hits
hits <- function(adjmatrix, maxiter=100L, tol=1e-05) {

    ## Sanity checks
    maxiter <- as.integer(maxiter)
    if (maxiter<=0)
        stop("maxiter should be a positive integer value")

    if (tol<=0)
        stop("tolerance (tol) should be a positive numeric value")

    if (NROW(adjmatrix) != NCOL(adjmatrix))
        stop("adjmatrix must be a square matrix")

    if (any(adjmatrix<0))
        stop("adjacency matrix cannot have negative weights")


    N <- NROW(adjmatrix)

    h <- rep(1, N)

    oldh <- rep(2, N)
    olda <- rep(2, N)
    
    for (i in 1:maxiter) {
        a <- as.vector(h %*% adjmatrix)
        a <- a / sqrt(sum(a^2))
        h <- as.vector(adjmatrix %*% a)
        h <- h / sqrt(sum(h^2))
        
        ## Check that 
        if ( (sqrt(sum(a-olda)^2) < sqrt(N)*tol) && (sqrt(sum(h-oldh)^2) < sqrt(N)*tol) )
            break
        olda <- a
        oldh <- h
    }

    list(authorities=a,
         hubs=h,
         iterations=i)
}

