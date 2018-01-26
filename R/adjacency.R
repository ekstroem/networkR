#' Create adjacency matrix
#' 
#' Create an adjacency matrix from a set of nodes and edges.
#' 
#' @param from a vector of nodes where the edges originate
#' @param to a vector of nodes where the edges point to
#' @param weight a numeric vector of weights 
#' @param directed logical. Are the edges directed (TRUE, the default) or bidirected(FALSE).
#' @return Returns a sparse adjacency matrix
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
#' @keywords manip
#' @examples
#'
#' from <- c("A", "A", "A", "B", "C")
#' to <- c("B", "C", "D", "D", "E")
#' adjacency(from, to)
#'
#' from <- c("A", "A", "A", "B", "C")
#' to <- c("B", "C", "D", "D", "E")
#' weights <- c(1, .5, 1, .7, 1)
#' adjacency(from, to, weights)
#'
#' @importFrom fastmatch fmatch
#' @importFrom Matrix sparseMatrix
#' @export adjacency
adjacency <- function(from, to, weight=1, directed=TRUE) {

    if (length(from) != length(to))
        stop("The from and to vectors must have the same length")
    
    if (length(weight)>1 & (length(weight) != length(from)))
        stop("The weight vectors must have the same length as from or be a scalar")
    
    entries <- unique(c(from, to))
    N <- length(entries)
    
    rows <- fastmatch::fmatch(from, entries)
    cols <- fastmatch::fmatch(to, entries)
    
    Matrix::sparseMatrix(i = rows, j = cols, x=weight, dims=c(N,N),
                         symmetric = !directed)
}
