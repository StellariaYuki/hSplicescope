.hsplicescope.module.ranges <- function() {
  list("1" = 1:1183, "2" = 1184:2095, "3" = 2096:2991)
}

.hsplicescope.normalize.modules <- function(modules, available.modules) {
  if (is.null(modules)) {
    return(available.modules)
  }
  modules <- trimws(as.character(modules))
  modules <- trimws(unlist(strsplit(modules, ",", fixed = TRUE)))
  modules <- modules[modules != ""]
  if (length(modules) == 0 || any(tolower(modules) %in% c("a", "all"))) {
    return(available.modules)
  }
  if (any(tolower(modules) %in% c("none", "off", "false", "0"))) {
    return(character(0))
  }
  normalized <- sub("^M", "", toupper(modules))
  invalid <- modules[!normalized %in% available.modules]
  if (length(invalid) > 0) {
    stop(
      "Invalid module value(s): ", paste(invalid, collapse = ", "),
      "\nValid values: 1, 2, 3, M1, M2, M3, a/all, or none.\n"
    )
  }
  unique(normalized)
}

#' Pipeline function for neuronal prediction with exons in each module
#'
#' @param data input splicing matrix with NA and extreme value replaced
#' @param modules Modules to run. Use NULL, "a", or "all" for all modules; use
#'   "none" to skip module predictions.
#' @return data.frame with Predicated maturation stage and Predication confidence score
#' @export
#'
#' @examples
#' module.predict(data)
module.predict <- function(data, modules = NULL) {
  
  res <- data.frame(matrix("", ncol = 0, nrow = ncol(data)), stringsAsFactors = F)
  rownames(res) <- colnames(data)
  
  module.ranges <- .hsplicescope.module.ranges()
  selected.modules <- .hsplicescope.normalize.modules(modules, names(module.ranges))
  for (module.id in selected.modules){
    j <- module.ranges[[module.id]]
    
    brainspan.Smodule2.module <- brainspan.Smodule2[j, , drop = FALSE]
    data.module <- data[match(rownames(brainspan.Smodule2)[j], rownames(data), nomatch = 0), , drop = FALSE]

    tmp.res <- module.beta.predict(data.module, brainspan.Smodule2.module)
    colnames(tmp.res) <- c(paste("M", module.id, "_Maturation", sep = ""), paste("M", module.id, "_ConfidenceScore", sep = ""))
    res <- cbind(res, tmp.res)
    
  }
  if (ncol(res) > 0) {
    res <- as.data.frame(lapply(res, as.numeric))
  }
  rownames(res) <- colnames(data)
  #res <- res[, match(newcolnames, colnames(res))]
  res
}
