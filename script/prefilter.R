#' Filter module exons by ID
#'
#' @param data Data splicing matrix
#'
#' @return Module exon splicing matrix
#' @export
#'
#' @examples
#' prefilter(data)
prefilter <- function(data) {
    
    sidx <- match(rownames(brainspan.Smodule), data$event_id, nomatch = 0)
    sidx <- sidx[sidx > 0]
    rownames(data) <- data[, 1]
    data.new <- data[sidx, -c(1:2)]
    # rownames(data.new) <- rownames(module)
    
    
    data.new
}
