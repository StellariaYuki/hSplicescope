#' Plot PCA scatter plot using ggplot2
#'
#' @param data.plot reference PCA plot data
#' @param test.plot user PCA plot data
#' @param plot default is TRUE to output the figure
#'
#' @param label output labelled plot when label is 1
#' @param prefix output file prefix
#' @return PCA pdf & png files
#' @export
#'
#' @examples
#' pcaplot(data.plot, test.plot)
#' 
pcaplot <- function(data.plot, test.plot, plot = T, label = 1, prefix = "dataPCA") {
    library(ggplot2)
    #data.plot$Name <- NA
    user.plot <- rbind(data.plot, test.plot)
    color.code <- c("#045a8d", "#43a2ca", "#48C5CC", "#1b9e77", "#a6d854", "#aab71d", "#ffe25e", "#fc8d59", "#d53e4f")
    shape.code <- c(16, 2, 0)
    size.code <- c(7, 4, 4)
    user.plot$Group <- as.factor(user.plot$Group)
    user.plot$Maturation <- as.factor(user.plot$Maturation)
    
    plot.res <- ggplot(data = user.plot, aes(x = PC1, y = PC2)) + theme(panel.grid.major = element_line(), 
        panel.grid.minor = element_line(), panel.background = element_blank(), axis.line.y = element_line(colour = "black"), 
        axis.line.x = element_line(colour = "black")) + 
    geom_point(aes(shape = Group, colour = Maturation, size = Group)) + scale_shape_manual(values = shape.code) + 
        scale_size_manual(values = size.code) + scale_colour_manual(values = color.code) + xlab("PC1") + 
        ylab("PC2") + theme(legend.key = element_blank())
   if (label==1) {
        reswithlabel <- plot.res + geom_text(aes(label = Name), size = 3, show.legend = F,  vjust = 2)
        ggsave(paste(prefix, "withlabel.pdf", sep = ""), plot = reswithlabel, width = 10, height = 6)
	ggsave(paste(prefix, ".pdf", sep = ""), plot = plot.res, width = 10, height = 6)
	ggsave(paste(prefix, ".png", sep = ""), plot=reswithlabel, width=10, height=6)
    }else{
 	ggsave(paste(prefix, ".png", sep = ""), plot = plot.res, width = 10, height = 6)
   	 ggsave(paste(prefix, ".pdf", sep = ""), plot = plot.res, width = 10, height = 6)
	}


}

.hsplicescope.get.plot.data <- function(plot.name) {
    plot.data <- get0(plot.name, mode = "list", inherits = TRUE)
    if (is.null(plot.data)) {
        warning("Skip PCA plot: object ", plot.name, " was not found in package data.")
        return(NULL)
    }
    plot.data
}

.hsplicescope.calculate.subset.pca <- function(test.data, reference.data) {
    shared.events <- intersect(rownames(reference.data), rownames(test.data))
    if (length(shared.events) < 2) {
        return(NULL)
    }

    reference.subset <- reference.data[shared.events, , drop = FALSE]
    test.subset <- test.data[shared.events, , drop = FALSE]
    pca.model <- prcomp(t(reference.subset), center = TRUE, scale. = FALSE)
    test.pca <- as.data.frame(predict(pca.model, t(test.subset))[, 1:2, drop = FALSE])

    list(
        test.pca = test.pca,
        test.sample.names = colnames(test.subset)
    )
}

.hsplicescope.make.subset.pca.plot <- function(test.data, reference.data, plot.data, prefix, label = 1, maturation = NULL) {
    pca.result <- .hsplicescope.calculate.subset.pca(test.data, reference.data)
    if (is.null(pca.result)) {
        warning("Skip PCA plot for ", prefix, ": fewer than two shared events.")
        return(list(files = character(0), pca = NULL))
    }
    test.pca <- pca.result$test.pca

    test.plot <- data.frame(matrix("", ncol = ncol(plot.data), nrow = length(pca.result$test.sample.names)), stringsAsFactors = F)
    colnames(test.plot) <- colnames(plot.data)
    test.plot$Group <- "UserDefined"
    test.plot$Name <- pca.result$test.sample.names
    test.plot$TrueAge <- NA
    test.plot$PC1 <- test.pca$PC1
    test.plot$PC2 <- test.pca$PC2
    if (is.null(maturation)) {
        test.plot$Maturation <- NA
    } else {
        test.plot$Maturation <- maturation[pca.result$test.sample.names]
    }

    pcaplot(plot.data, test.plot, plot = T, label = label, prefix = prefix)
    reference.pca <- as.data.frame(plot.data[, c("PC1", "PC2"), drop = FALSE])
    if (nrow(reference.pca) == ncol(reference.data)) {
        rownames(reference.pca) <- colnames(reference.data)
    }
    pca.table <- list(test = pca.result$test.pca, reference = reference.pca)
    if (label == 1) {
        return(list(
            files = c(paste(prefix, "withlabel.pdf", sep = ""), paste(prefix, ".pdf", sep = ""), paste(prefix, ".png", sep = "")),
            pca = pca.table
        ))
    }
    list(
        files = c(paste(prefix, ".pdf", sep = ""), paste(prefix, ".png", sep = "")),
        pca = pca.table
    )
}
