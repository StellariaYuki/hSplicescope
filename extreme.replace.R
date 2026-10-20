#' Replace 0 and 1 in splicing matrix with approximate value
#'
#' @param splicing matrix
#'
#' @return splicing matrix with 0 and 1 replaced
#' @export
#'
#' @examples
#' extreme.replace(data)
extreme.replace <- function(data) {
    rep.extreme <- function(x) {
        n <- length(x)
        x[x == 0 | x == 1] <- (x[x == 0 | x == 1] * (n - 1) + 0.5)/n
        x
    }
    data.nonextreme <- apply(as.matrix(data), 2, rep.extreme)
    data.nonextreme
}
