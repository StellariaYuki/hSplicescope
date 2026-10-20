#' Replace NA values with mean of reference sample 
#'
#' @param module exon splicing matrix
#'
#' @return module exon splicing matrix with NA replaced
#' @export
#'
#' @examples
#' na.replace(data)
na.replace <- function(data) {
    
    rep <- function(x) {
        idx <- which(is.na(x))
        x[idx] <- brainspan.Smodulemean[idx]
        x
    }
    data.narm <- t(apply(as.matrix(data), 1, rep))
    
    data.narm
}
