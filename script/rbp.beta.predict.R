#' predict neuronal maturation with all module exons using beta regression methods with RBP specif exons
#'
#' @param data processed splicing matrix with 'NA' and extreme value removed
#'
#' @return data.frame with Predicated maturation stage and Predication confidence score
#' @export
#'
#' @examples
#' rbp.beta.predict(data, brainspan.RBP)
rbp.beta.predict <- function(data, brainspan.RBP) {
  library(betareg)
  res <- data.frame(matrix("", ncol = 2, nrow = ncol(data)), stringsAsFactors = F)
  rownames(res) <- colnames(data)
  colnames(res) <- c("Maturation", "ConfidenceScore")
  
  for (i in 1:ncol(data)) {
    if(max(data[,i])-min(data[,i])>0){
      beta.res <- betareg(data[, i] ~ ., as.data.frame(brainspan.RBP))
      beta.Sdis <- get.dis.rbp(beta.res, brainspan.RBP)
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
