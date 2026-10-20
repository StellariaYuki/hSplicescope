#' predict neuronal maturation with all module exons using beta regression methods
#'
#' @param data processed splicing matrix with 'NA' and extreme value removed
#'
#' @return data.frame with Predicated maturation stage and Predication confidence score
#' @export
#'
#' @examples
#' beta.predict(processed_splicing_matrix)
beta.predict <- function(data) {
    library(betareg)
    res <- data.frame(matrix("", ncol = 2, nrow = ncol(data)), stringsAsFactors = F)
    rownames(res) <- colnames(data)
    colnames(res) <- c("Maturation", "ConfidenceScore")
    
    ##if the row number of input is less than brainspan.Smodule2, filt it based on input rownames
    sidx <- match(rownames(data), rownames(brainspan.Smodule2), nomatch = 0)
    sidx <- sidx[sidx > 0]
    brainspan.Smodule2.filt <<- brainspan.Smodule2[sidx, c(1:9)]
    
    
    for (i in 1:ncol(data)) {
        
        beta.res <- betareg(data[, i] ~ ., as.data.frame(brainspan.Smodule2.filt)) ##devcortex.Smodule2->brainspan.Smodule2
        
        ##summary(beta.res)
        beta.Sdis <- get.dis(beta.res)
        knn.res <- knn.flex.reg(beta.Sdis, 1)
        confi.score <- 1 - min(beta.Sdis)/max(beta.Sdis)
        beta.class <- unlist(sapply(knn.res, function(x) return(Score[[x]]))) 
        res[i, 1] <- beta.class
        res[i, 2] <- confi.score
    }
    
    res <- as.data.frame(lapply(res, as.numeric))
    rownames(res) <- colnames(data)
    res
    
}
