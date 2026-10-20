#' Pipeline function for neuronal prediction
#'
#' @param test.data input splicing matrix
#' @param outfile  output file for prediction result
#' @param plot default is TRUE to plot PCA for all the data
#' @param rbps optional RBP names to run. Use "a" or "all" for all supported RBPs.
#' @param modules optional module IDs to run. Use "a" or "all" for all modules.
#'
#' @return maturation prediction and confidence with PCA result
#' @export
#'
#' @examples
#' maturation.predict(data, res)
maturation.predict <- function(test.data, outfile, plot = T, verbose = verbose, label=1, rbps = NULL, modules = NULL) {
    if (verbose) {
        cat("Processing data......\n")
    }
    test.Smodule <- prefilter(test.data)
    na_count <- colSums(is.na(test.Smodule))
    write.table(na_count, file = "data_NA_count.txt", quote = F, 
                sep = "\t", col.names = T, row.names = T)
    if (verbose) {
        cat("if the number of NA count is more than 1500, be care of the outcome\n")
    }
    test.Smodule2 <- na.replace(test.Smodule)
    test.Smodule3 <- extreme.replace(test.Smodule2)
    if (verbose) {
        cat("Predicting maturation......\n")
    }
    
    test.betareg <- beta.predict(test.Smodule3)
    test.module.betareg <- module.predict(test.Smodule3, modules = modules)
    test.rbp.betareg <- NULL
    if (!is.null(rbps)) {
        test.rbp.betareg <- rbp.predict(test.Smodule3, rbps = rbps)
    }
    if (verbose) {
        cat("Output results......\n")
    }
    pca.plot.files <- character(0)
    extra.pca.tables <- list()
    if (!is.null(outfile)) {
        write.table(cbind(test.betareg, test.module.betareg), file = "result_prediction.txt", quote = F, 
            sep = "\t", col.names = T, row.names = T) ##rm test.rbp.betareg
        if (!is.null(test.rbp.betareg)) {
            write.table(cbind(test.betareg, test.rbp.betareg), file = "result_RBP_prediction.txt", quote = F, 
            sep = "\t", col.names = T, row.names = T)
        }
    }
    test.pca <- as.data.frame(predict(brainspan.pca, t(test.Smodule2))[, 1:2])
    if (plot) {
        
        test.plot <- data.frame(matrix("", ncol = ncol(brainspan.plot), nrow = ncol(test.Smodule3)), stringsAsFactors = F)
        colnames(test.plot) <- colnames(brainspan.plot)
        test.plot$Group <- "UserDefined"
        test.plot$Name <- colnames(test.Smodule3)
        test.plot$TrueAge <- NA
    #     test.plot$True.Label <- NA
        test.plot$PC1 <- test.pca$PC1
        test.plot$PC2 <- test.pca$PC2
        test.plot$Maturation <- test.betareg$Maturation
        #user.plot <- rbind(data.plot, test.plot)
        pcaplot(brainspan.plot, test.plot, plot, label)

        module.ranges <- .hsplicescope.module.ranges()
        selected.modules <- .hsplicescope.normalize.modules(modules, names(module.ranges))
        for (module.id in selected.modules) {
            plot.data <- .hsplicescope.get.plot.data(paste("M", module.id, "_plot", sep = ""))
            if (!is.null(plot.data)) {
                maturation.col <- paste("M", module.id, "_Maturation", sep = "")
                module.maturation <- NULL
                if (maturation.col %in% colnames(test.module.betareg)) {
                    module.maturation <- test.module.betareg[, maturation.col]
                    names(module.maturation) <- rownames(test.module.betareg)
                }
                pca.result <- .hsplicescope.make.subset.pca.plot(
                    test.Smodule2,
                    brainspan.Smodule2[module.ranges[[module.id]], , drop = FALSE],
                    plot.data,
                    paste("M", module.id, "_PCA", sep = ""),
                    label,
                    module.maturation
                )
                pca.plot.files <- c(pca.plot.files, pca.result$files)
                if (!is.null(pca.result$pca)) {
                    colnames(pca.result$pca$test) <- c(paste("M", module.id, "_PC1", sep = ""), paste("M", module.id, "_PC2", sep = ""))
                    colnames(pca.result$pca$reference) <- c(paste("M", module.id, "_PC1", sep = ""), paste("M", module.id, "_PC2", sep = ""))
                    extra.pca.tables[[paste("M", module.id, sep = "")]] <- pca.result$pca
                }
            }
        }

        if (!is.null(rbps)) {
            RBP <- 12:20
            selected.rbps <- .hsplicescope.normalize.rbps(rbps, colnames(RBPall.bn)[RBP])
            for (rbp in selected.rbps) {
                plot.data <- .hsplicescope.get.plot.data(paste(rbp, "_plot", sep = ""))
                if (!is.null(plot.data)) {
                    idx <- RBP[match(rbp, colnames(RBPall.bn)[RBP])]
                    RBP.bn <- subset(RBPall.bn, RBPall.bn[, idx] != 0)
                    RBP.bn <- RBP.bn[order(RBP.bn[, idx]), ]
                    reference.idx <- match(RBP.bn$event_id, rownames(brainspan.Smodule2), nomatch = 0)
                    reference.idx <- reference.idx[reference.idx > 0]
                    maturation.col <- paste(rbp, "_Maturation", sep = "")
                    rbp.maturation <- NULL
                    if (!is.null(test.rbp.betareg) && maturation.col %in% colnames(test.rbp.betareg)) {
                        rbp.maturation <- test.rbp.betareg[, maturation.col]
                        names(rbp.maturation) <- rownames(test.rbp.betareg)
                    }
                    pca.result <- .hsplicescope.make.subset.pca.plot(
                        test.Smodule2,
                        brainspan.Smodule2[reference.idx, , drop = FALSE],
                        plot.data,
                        paste(rbp, "_PCA", sep = ""),
                        label,
                        rbp.maturation
                    )
                    pca.plot.files <- c(pca.plot.files, pca.result$files)
                    if (!is.null(pca.result$pca)) {
                        colnames(pca.result$pca$test) <- c(paste(rbp, "_PC1", sep = ""), paste(rbp, "_PC2", sep = ""))
                        colnames(pca.result$pca$reference) <- c(paste(rbp, "_PC1", sep = ""), paste(rbp, "_PC2", sep = ""))
                        extra.pca.tables[[rbp]] <- pca.result$pca
                    }
                }
            }
        }
    }
    if (!is.null(outfile)) {
        result.pca <- rbind(test.pca, brainspan.pca$x[, 1:2])
        colnames(result.pca) <- c("PC1", "PC2")
        for (pca.table in extra.pca.tables) {
            pca.add <- data.frame(matrix(NA, nrow = nrow(result.pca), ncol = ncol(pca.table$test)))
            rownames(pca.add) <- rownames(result.pca)
            colnames(pca.add) <- colnames(pca.table$test)

            matched.rows <- match(rownames(result.pca), rownames(pca.table$test))
            matched.idx <- which(!is.na(matched.rows))
            if (length(matched.idx) > 0) {
                pca.add[matched.idx, ] <- pca.table$test[matched.rows[matched.idx], , drop = FALSE]
            }

            user.n <- nrow(test.pca)
            result.reference.rows <- seq(user.n + 1, nrow(result.pca))
            if (length(result.reference.rows) == nrow(pca.table$reference)) {
                pca.add[result.reference.rows, ] <- pca.table$reference
            }

            result.pca <- cbind(result.pca, pca.add)
        }
        write.table(result.pca, file = "result_pca.txt", quote = F, sep = "\t", 
            col.names = T, row.names = T)
    }
    # if (verbose) {
    #     cat("Generating HTML file......\n")
    # }
    # generate.html(test.betareg) ##remove test.rbp.betareg
    # if (verbose) {
    #     cat("Generating zip file......\n")
    # }
    if (!is.null(outfile)) {
        zip.files <- c("dataPCA.png", "dataPCA.pdf", "result_prediction.txt", "result_pca.txt")
        if (label == 1) {
            zip.files <- c("dataPCAwithlabel.pdf", zip.files)
        }
        if (!is.null(test.rbp.betareg)) {
            zip.files <- c(zip.files, "result_RBP_prediction.txt")
        }
        zip.files <- c(zip.files, pca.plot.files)
        system(paste("zip -qm", outfile, paste(zip.files, collapse = " ")))  ## remove index.html --yuki
    }
    if (verbose) {
        cat("Done......\n")
    }
}
