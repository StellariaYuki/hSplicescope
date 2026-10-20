#' predict neuronal maturation with exons in each module using beta regression methods
#'
#' @param data processed splicing matrix with 'NA' and extreme value removed
#'
#' @return data.frame with Predicated maturation stage and Predication confidence score in each module
#' @export
#'
#' @examples
#' module.beta.predict(data, module)
module.beta.predict <- function(data, module) {
  library(betareg)
  res <- data.frame(matrix("", ncol = 2, nrow = ncol(data)), stringsAsFactors = F)
  rownames(res) <- colnames(data)
  colnames(res) <- c("Maturation", "ConfidenceScore")
  
  for (i in 1:ncol(data)) {
    if(max(data[,i])-min(data[,i])>0){
      beta.res <- betareg(data[, i] ~ ., as.data.frame(module))
      beta.Sdis <- get.dis.module(beta.res, module)
      knn.res <- knn.flex.reg(beta.Sdis, 1)
      confi.score <- 1 - min(beta.Sdis)/max(beta.Sdis)
      beta.class <- unlist(sapply(knn.res, function(x) return(Score[[x]])))
      res[i, 1] <- beta.class
      res[i, 2] <- confi.score
    }else{
      res[i, 1] <- NA
      res[i, 2] <- NA
    }
  }
  
  res
  
}