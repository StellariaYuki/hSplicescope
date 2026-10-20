#' Calculate distance between reference sample and predicated samples with exons in each module
#'
#' @param beta regression object
#'
#' @return vector of distance to reference samples
#' @export
#'
#' @examples
#' get.dis.module(lm.obj, module)
get.dis.module <- function(lm.obj, module) {
  # y_hat <- lm.obj$fitted.values
  dis <- as.matrix(dist(t(cbind(lm.obj$fitted.values, module))))[1, -1]
  dis <- dis/sum(dis)
  dis
}