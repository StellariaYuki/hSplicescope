#' Calculate distance between reference sample and predicated samples with specifc RBP binding exons
#'
#' @param beta regression object
#'
#' @return vector of distance to reference samples
#' @export
#'
#' @examples
#' get.dis.module(lm.obj, brainspan.RBP)
get.dis.rbp <- function(lm.obj, brainspan.RBP) {
  dis <- as.matrix(dist(t(cbind(lm.obj$fitted.values, brainspan.RBP))))[1, -1]
  dis <- dis/sum(dis)
  dis
}