.hsplicescope.supported.rbps <- function() {
  c("QKI", "U2AF2", "U2AF1", "PTBP1", "RBFOX2", "SRSF1", "KHSRP", "TIA1", "IGF2BP1")
}

.hsplicescope.normalize.rbps <- function(rbps, available.rbps) {
  if (is.null(rbps)) {
    return(available.rbps)
  }
  rbps <- trimws(as.character(rbps))
  rbps <- trimws(unlist(strsplit(rbps, ",", fixed = TRUE)))
  rbps <- rbps[rbps != ""]
  if (length(rbps) == 0 || any(tolower(rbps) %in% c("a", "all"))) {
    return(available.rbps)
  }
  rbp.map <- setNames(available.rbps, toupper(available.rbps))
  normalized <- toupper(rbps)
  invalid <- rbps[!normalized %in% names(rbp.map)]
  if (length(invalid) > 0) {
    stop(
      "Invalid RBP value(s): ", paste(invalid, collapse = ", "),
      "\nValid values: ", paste(available.rbps, collapse = ", "),
      ", or a/all for all supported RBPs.\n"
    )
  }
  unname(rbp.map[normalized])
}

#' Pipeline function for neuronal prediction with RBP sepcific exons
#'
#' @param data input splicing matrix with NA and extreme value replaced
#' @param rbps RBP names to run. Use NULL, "a", or "all" for all supported RBPs.
#' @return data.frame with Predicated maturation stage and Predication confidence score
#' @export
#'
#' @examples
#' rbp.predict(data)
rbp.predict <- function(data, rbps = NULL) {
  RBP <- 12:20
  rbp.names <- colnames(RBPall.bn)[RBP]
  canonical.rbps <- .hsplicescope.supported.rbps()
  if (all(canonical.rbps %in% rbp.names)) {
    rbp.names <- canonical.rbps
  }
  selected.rbps <- .hsplicescope.normalize.rbps(rbps, rbp.names)

  res <- data.frame(matrix("", ncol = 0, nrow = ncol(data)), stringsAsFactors = F)
  rownames(res) <- colnames(data)
  for (rbp in selected.rbps) {
    idx <- RBP[match(rbp, colnames(RBPall.bn)[RBP])]
    RBP.bn <- subset(RBPall.bn, RBPall.bn[, idx] != 0)
    RBP.bn <- RBP.bn[order(RBP.bn[, idx]), ]
    RBP.tissue.Smodule <- data[match(RBP.bn$event_id, rownames(data), nomatch = 0), , drop = FALSE]
    brainspan.RBP.Smodule2 <- brainspan.Smodule2[match(RBP.bn$event_id, rownames(brainspan.Smodule2), 
                                                       nomatch = 0), , drop = FALSE]
    # rbp.tissue.Smodule2 <- na.replace(rbp.tissue.Smodule)
    tmp.res <- rbp.beta.predict(RBP.tissue.Smodule, brainspan.RBP.Smodule2)
    colnames(tmp.res) <- c(paste(colnames(RBPall.bn)[idx], "_Maturation", sep = ""), paste(colnames(RBPall.bn)[idx], "_ConfidenceScore", sep = ""))
    res <- cbind(res, tmp.res)
    
  }
  if (ncol(res) > 0) {
    res <- as.data.frame(lapply(res, as.numeric))
  }
  rownames(res) <- colnames(data) 
  newcolnames <- unlist(lapply(selected.rbps, function(rbp) {
    c(paste(rbp, "_Maturation", sep = ""), paste(rbp, "_ConfidenceScore", sep = ""))
  }))
  res <- res[, match(newcolnames, colnames(res)), drop = FALSE]
  res
}
