#' Validate pedigree trio information consistency
#' 
#' Simple tests to chech the consistency of the pedigree trio family data. 
#' Currently the following checks are undertaken:
#' 1) that no duplicates ids are found;
#' 2) that the primary id is not missing for anyone;
#' 3) that founders have both the father and mother id missing;
#' 4) that individuals are not both classified as male (fathers and mothers);
#' 
#' There are no checks of persons being both mother and father, nor being its
#' own parent and incest checks are not performed. In other words, the
#' \code{obj} is assumed to be sane, but possibly immoral.
#' 
#' @param id Numeric. The id of the individual. These values should be unique
#' @param fid Numeric. The father id. NA or 0 are used for missing.
#' @param mid Numeric. The mother id. NA or 0 are used for missing.
#' @param sex An optional numeric vector with the sex of the individual. Only four values should be present 1 (male), 2 (female), 0 or NA (missing)
#' @return Throws an error if an inconsistency is found. Otherwise returns TRUE.
#' @author Claus Thorn Ekstrøm, \email{ekstrom@@sund.ku.dk}
#' @keywords manip
#' @examples
#' 
#' library("data.table")
#' id <- 1:12
#' fid <- c(NA,  0, 1, 1, NA, 23, 45, 5, 5, 7, 10, 10)
#' mid <- c(NA, NA, 2, 2,  0, 56, 46, 6, 6, 6, 9, 11)
#'
#' validate_trio_consistency(id, fid, mid)
#' 
#' @import data.table
#' @export
validate_trio_consistency <- function(id, fid, mid, sex=NULL) {

    if (! class(id)[1] %in% c("numeric", "integer"))
        stop("id should be a numeric vector")

    if (! class(fid)[1] %in% c("numeric", "integer"))
        stop("fid should be a numeric vector")
    
    if (! class(mid)[1] %in% c("numeric", "integer"))
        stop("mid should be a numeric vector")

    if (is.null(sex)) {
        DT <- data.table(id, fid, mid)
    } else {
        if (! class(sex)[1] %in% c("numeric", "integer"))
            stop("sex should be a numeric vector")

        ## Should check that NA, 0, 1, 2 are the only possibilities
        if (any(!(sex %in% c(NA, 0, 1, 2))))
            stop("sex should be one of the following values NA, 0, 1, 2")
            
        DT <- data.table(id, fid, mid, sex)
    }
    
    ## Replace all NA with 0's
    for (j in seq_len(ncol(DT)))
        set(DT,which(is.na(DT[[j]])),j,0)

    ## Do the following checks
    ## 1) No duplicate ids
    if (any(duplicated(DT[id!=0,list(id)])))
        stop("duplicated ids found")

    ## 1) No zeroes for id
    if (DT[, list(zeroids=sum(id == 0))] > 0) {
        stop("the id variable should not contain 0 or NAs")
    }

    ## 2) Check that founders are proper founders (ie have both variables missing)
    if (NROW(DT[,list(number=.N), by=list(fid>0, mid>0)][fid!=mid,])>0) {
        stop("some individuals are only part founders (have exactly one of fid or mid missing)")
    }

    ## 3) Check that an id is not present both as father and mother
    if (NROW(DT[,list(number=.N), by=list(fid, mid)][mid==fid & fid>0,])>0) {
        
        stop("some individuals are classified both as fathers and mothers")
    }

    if (!is.null(sex)) {
        ## Males
        if (NROW(DT[sex==1, list(id)][DT[mid>0,list(id=mid)], on="id", nomatch=0])>0)
            stop("some individuals with sex=1 are found as mothers of others")
        if (NROW(DT[sex==2, list(id)][DT[fid>0,list(id=fid)], on="id", nomatch=0])>0)
            stop("some individuals with sex=2 are found as fathers of others")         
    }

    return(TRUE)
}
