
#' Generate variables (or lists) of siblings from a file of ids of persons and
#' their father and mother.
#' 
#' The function generates for each person lists of maternal half-sibs, paternal
#' half-sibs and full sibs. Optionally these are expanded to separate columns
#' in a data.table.
#' 
#' There are no checks of persons being both mother and father, nor being its
#' own parent and incest checks are not performed. In other words, the
#' \code{obj} is assumed to be sane, but possibly immoral.
#' 
#' @param obj A 3-column structure with column names \code{id}, \code{pid}
#' (paternal id) and \code{mid} (maternal id).
#' @param ns Integer. The maximal no of sibs of each type to include in the
#' result if sibling ids are required in separate columns.
#' @param expand.vars Logical. Should the sibling ids be returned in separate
#' columns. If \code{FALSE} they will be returned i three columns of lists.
#' @return A \code{data.table} with the columns of the \code{obj} and columns
#' for \code{ns} maternal, paternal and full sibs, named \code{ms1},
#' \code{ms2}, \ldots{} \code{ps1}, \code{ps2}, \ldots{} \code{fs1},
#' \code{fs2}.
#' 
#' If \code{expand.vars=FALSE} there will instead be three columns of lists
#' named \code{msibs}, \code{psibs} and \code{fsibs}.
#' @author Claus Thorn Ekstrøm, \email{ekstrom@@sund.ku.dk}, Bendix Carstensen,
#' \email{b@@bxc.dk}
#' @keywords manip
#' @examples
#' 
#' library( data.table )
#' id <- 1:12
#' pid <- c(NA,  1, 1, 1, NA, 23, 45, 5, 5, 7, 12, NA)
#' mid <- c(NA, NA, 2, 2, 12, NA, 46, 6, 6, 6, NA, 12)
#' indd <- data.table( id, mid, pid )
#' indata <- copy( indd )
#' indata
#' 
#' str( xx <- mksib( indata ) )
#' xx
#' 
#' zz <- mksib( indata, 2, e=FALSE )
#' zz
#'
#' @import data.table
#' @export
mksib <- function( obj, ns=3, expand.vars=TRUE ) {
    indata <- copy( obj )

    ## Check that the ...

    ## Setting these variables to prevent "no visible binding for global variable"

    id <- NULL
    pid <- NULL
    mid <- NULL
    fsibs <- NULL
    msibs <- NULL
    psibs <- NULL
    
    ## missing parent id causes very large fsibships
    ## missing replaced by bogus ids beyond maximal id
    mxid <- max( c(indata$pid,
                   indata$mid), na.rm=TRUE )
    pna <- sum( is.na(indata[,"pid"]) )
    mna <- sum( is.na(indata[,"mid"]) )
    indata[is.na(pid),pid := mxid + 1:pna]
    indata[is.na(mid),mid := mxid + pna + 1:mna]
                                        # this is the core
    setDT(indata)[,msibs := list(list(id)), by = "mid"][         # lists of maternal children
       ,msibs := mapply(setdiff, msibs, id)][      # but not your own sib
       ,psibs := list(list(id)), by = "pid"][         # lists of paternal children
       ,psibs := mapply(setdiff, psibs, id)][      # but not your own sib
       ,fsibs := mapply(intersect, msibs, psibs)][ # full sibs
       ,msibs := mapply(setdiff, msibs, fsibs)][   # remove full from maternal
       ,psibs := mapply(setdiff, psibs, fsibs)]    # remove full from paternal
    if( expand.vars ) {  
                                        # expand the lists to variables
        for( i in 1:ns )
            indata[,paste("ms",i,sep="")] <- sapply(indata$msibs, function(x) x[i])
        for( i in 1:ns )
            indata[,paste("ps",i,sep="")] <- sapply(indata$psibs, function(x) x[i])
        for( i in 1:ns )
            indata[,paste("fs",i,sep="")] <- sapply(indata$fsibs , function(x) x[i])
    }
    ## bogus ids reset to NA
    indata[pid>mxid,pid:=NA]
    indata[mid>mxid,mid:=NA]
    ## remove unwanted columns in the table
    indata[,c("msib","psib") := NULL]
    if( expand.vars ) indata[,c("msibs","psibs","fsibs") := NULL]
    ## return object
    indata
}
