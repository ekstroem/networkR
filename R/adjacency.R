#' Create adjacency matrix
#' 
#' Split a matrix into block diagonal sub matrices according to clusters and
#' combine the lower triangular parts into a vector
#' 
#' @param from a vector
#' @param to a vector
#' @param type numeric 
#' @param directed logical. Are the edges directed (TRUE, the default) or bidirected(FALSE).
#' @return Returns a sparse adjacency matrix
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @keywords manip
#' @examples
#'
#' from <- 1:5
#'
#' @importFrom fastmatch fmatch
#' @importFrom Matrix sparseMatrix
#' @export adjacency
adjacency <- function(from, to, type=1, directed=TRUE) {

    if (length(from) != length(to))
        stop("XXXXInput must be a square matrix")
    
    entries <- unique(c(from, to))
    N <- length(entries)

    rows <- fastmatch::fmatch(from, entries)
    cols <- fastmatch::fmatch(to, entries)

    Matrix::sparseMatrix(i = rows, j = cols, x=type, 
             symmetric = !directed)
}
