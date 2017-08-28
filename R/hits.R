#' Create adjacency matrix
#' 
#' Split a matrix into block diagonal sub matrices according to clusters and
#' combine the lower triangular parts into a vector
#' 
#' @param adjmatrix an adjacency matrix
#' @param maxiter non-negative integer
#' @return Returns a sparse adjacency matrix
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @keywords manip
#' @examples
#'
#' from <- c("A", "A", "A", "B", "C")
#' to <- c("B", "C", "D", "D", "E")
#' adjacency(from, to) %>% hits()
#'
#' @importFrom Matrix sparseMatrix
#' @export hits
hits <- function(adjmatrix, maxiter=25L) {

    ## 
    
    N <- nrow(adjmatrix)

    h <- rep(1, N)
    for (i in 1:maxiter) {
        a <- crossprod(adjmatrix, h)
        a <- a / sqrt(sum(a^2))
        h <- adjmatrix %*% a
        h <- h / sqrt(sum(h^2))        
    }

    list(authorities=a, hubs=h)
}

