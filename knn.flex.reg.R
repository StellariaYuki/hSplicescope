#' KNN predictor
#'
#' @param dis distance vector output by function get.dis
#' @param k  parameter for neighbor number
#'
#' @return Maturation stage
#' @export
#'
#' @examples
#' knn.flex.reg(dis, k)
knn.flex.reg <- function(dis, k) {
    
    label <- Stage
    
    class.beta <- order(dis)[1:k]
    class.res <- Stage[class.beta]
    # names(class.res) <- colnames(distance.matrix)
    class.res
}
