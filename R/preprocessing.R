#' pre-filtering low quality cells
#'
#' Pre-filtering low quality cells based on library complexity (i.e., the number of genes expressed in a library)
#'
#' @param object A sincera object
#' @param min.expression The threshold value of expressed genes
#' @param min.genes Cells with less than this number of expressed genes will be excluded from downstream analysis
#' @param do.plot If TRUE, plot the distribution of library complexity
#' @return An update sincera object with low quality cells pre-filtered
#'
setGeneric("filterLowQualityCells", function(object, min.expression=1, min.genes=500, do.plot=FALSE, ...) standardGeneric("filterLowQualityCells"))
#' @export
setMethod("filterLowQualityCells","sincera",
          function(object, min.expression=1, min.genes=500, do.plot=FALSE, ...) {

            # select cells with at least min.genes number of genes with expression at least min.expression
            ngenes <- colSums(getExpression(object)>=min.expression)

            cat("\nNumber of expressed genes per cell\n\n")

            print(tapply(ngenes, getCellMeta(object, name="GROUP"), summary))

            if (do.plot==TRUE) {


              if (FALSE) {
                g <- ggplot(data.frame(Group=factor(getCellMeta(object,name="GROUP")), Num_genes=ngenes), aes(x=Group, y=Num_genes))
                g <- g + geom_jitter(height=0) + geom_violin(aes(fill=Group),  alpha=0.3, scale="width",adjust=0.75)
                g <- g + ggtitle(label=paste("Number of expressed genes per cell ( min.exp=", min.expression," )", sep=""))
                g <- g + ylab("# expressed genes")
                g <- g + sincera_theme()
                print(g)
              }

              g <- ggplot(data.frame(Group=factor(getCellMeta(object,name="GROUP")), Num_genes=ngenes), aes(x=Num_genes))
              g <- g + geom_histogram(binwidth = 500, fill="grey",  col="black") + xlim(c(0,20000))
              g <- g + ggtitle("Histogram of library complexity of all cells")
              g <- g + sincera_theme()
              g1 <- g

              g <- ggplot(data.frame(Group=factor(getCellMeta(object,name="GROUP")), Num_genes=ngenes), aes(x=Num_genes))
              g <- g + facet_wrap(~Group, scales="free") + geom_histogram(binwidth = 500, fill="grey",  col="black") + xlim(c(0,20000))
              g <- g + ggtitle("Histogram of library complexity per sample")
              g <- g + sincera_theme()
              g2 <- g

              library(gridExtra)
              library(grid)
              grid.arrange(g1, g2, ncol=1)

            }

            cells.selected <- which(ngenes>=min.genes)
            cat("\n", length(cells.selected), " (", length(cells.selected)*100/length(ngenes), "%) cells passed the prefiltering.\n", sep="")
            object@data <- object@data[, cells.selected]

            return(object)
          }
)

#' Detecting and removing potentially contaminated cells
#'
#' Detecting contaminated cells based on the co-expression of highly selectively markers of two different cell types, such as epithelial and endothelial
#'
#' @param object A sincera object
#' @param do.rm If TRUE, remove detected contaminants
#' @param markers.1 The vector encoding the markers for one cell type
#' @param markers.2 The vector encoding the markers for the other cell type
#' @param min.expression.1 The threshold value of the expression of markers in markers.1
#' @param min.quantile.1 The threshold value of the expression quantile of markers in markers.1
#' @param min.expression.2 The threshold value of the expression of markers in markers.2
#' @param min.quantile.2 The threshold value of the expression quantile of markers in markers.2
#' @return An updated sincera object with detected contaminated cells removed
#'
setGeneric("filterContaminatedCells", function(object, do.rm=T, markers.1, markers.2, min.expression.1=1, min.quantile.1=0.99, min.expression.2=1, min.quantile.2=0.99, ...) standardGeneric("filterContaminatedCells"))
#' @export
setMethod("filterContaminatedCells","sincera",
          function(object, do.rm=T, markers.1, markers.2, min.expression.1=1, min.quantile.1=0.99, min.expression.2=1, min.quantile.2=0.99, ...) {

            cat("\nDetecting contaminated cells:\n")

            m1 <- which(markers.1 %in% getGenes(object))
            m2 <- which(markers.2 %in% getGenes(object))

            if (length(m1)==0 | length(m2)==0) {
              cat("No contaminated cells were detected.\n")
              return(object)
            }

            cellnames <- getCells(object)

            c1 <- cellnames
            for (m in markers.1[m1]) {
              m.values <- as.numeric(getExpression(object)[m, ])
              m.topquantile <- which(m.values >= quantile(m.values, probs=c(min.quantile.1)))
              m.expressed <- which(m.values >= min.expression.1)
              m.cells <- intersect(m.topquantile, m.expressed)
              c1 <- intersect(c1, cellnames[m.cells])
            }
            cat(length(c1), "cells co-highly-expressed markers of first type:", paste(c1, collapse=", ", sep=""),"\n")

            c2 <- cellnames
            for (m in markers.2[m2]) {
              m.values <- as.numeric(getExpression(object)[m, ])
              m.topquantile <- which(m.values >= quantile(m.values, probs=c(min.quantile.2)))
              m.expressed <- which(m.values >= min.expression.2)
              m.cells <- intersect(m.topquantile, m.expressed)
              c2 <- intersect(c2, cellnames[m.cells])
            }
            cat(length(c2), "cells co-highly-expressed markers of second type:", paste(c2, collapse=", ", sep=""),"\n")

            contaminated <- intersect(c1, c2)
            if (length(contaminated) > 0) {
              cat(length(contaminated),"cells satisfied the contamination criteria and were removed: ", paste(contaminated, collapse=", ", sep=""), "\n")
              cells.selected <- setdiff(cellnames, contaminated)
              if (length(cells.selected) < 2) {
                  stop("Only 1 cells passed the contamination detection.")
              }
              if (do.rm==TRUE) {
                object@data <- object@data[, cells.selected]
                cat(length(contaminated), "contaminants removed\n")
              } else {
                cat("NOTE: because do.rm==FALSE, the detected contaminants were not removed\n")
              }
            } else {
              cat("No contaminated cells were detected.\n")
            }

            return(object)
          }
)

#' Set the lower bound of gene expression
#'
#' Set the lower bound of gene expression
#'
#' @param object A sincera object
#' @param value The lower bound value
#' @return An updated sincera object with the manipulated expression values, use getExpression function to access the expression values
#'
setGeneric("expr.minimum", function(object, value=0.01, ...) standardGeneric("expr.minimum"))
#' @export
setMethod("expr.minimum","sincera",
          function(object, value=0.01, ...) {

              expr <- getExpression(object)

              n <- length(which(expr<value))
              n.all <- getCellNum(object)*getGeneNum(object)
              expr[expr<value] <- value
              object <- setExpression(object, value = expr)
              cat("The minimum expression value was set to ", value, "\n", sep="")
              cat(n, " (", 100*n/n.all, "%) expression values are affected.\n\n", sep="")

              return(object)
          }
)


#' Pre-filtering genes with low abundancy
#'
#' @param object A sincera object
#' @param pergroup If TRUE, the prefiltering will be based on genes' per sample abundancy; otherwise, it will be based on genes' abundancy in all cells
#' @param min.expression The threshold value of gene expression
#' @param min.cells The minimum number of cells
#' @param min.samples The selected genes must pass the abundancy criteria in at least \"min.samples\" samples
#' @return An update sincera object with prefiltered genes removed
#'
setGeneric("prefilterGenes", function(object, pergroup=TRUE, min.expression=5, min.cells=2, min.samples=1, ...) standardGeneric("prefilterGenes"))
#' @export
setMethod("prefilterGenes","sincera",
          function(object, pergroup=TRUE, min.expression=5, min.cells=2, min.samples=1, ...) {

              cat("\nPrefitering genes with low expression abundancy\n")

              cellsample <- getCellMeta(object, name="GROUP")

              if (pergroup==FALSE) {
                  cellsample <- rep("sample1", length(cellsample))
                  min.samples=1
              }

              samples <- sort(unique(as.character(cellsample)))

              stats <- data.frame(SYMBOL=getGenes(object))
              rownames(stats) <- stats$SYMBOL
              for (s in samples) {
                s.idx <- which(cellsample == s)
                stats[, s] <- rowSums(getExpression(object)[, s.idx]>=min.expression)
              }

              #stats <- stats[, -1]

              if (length(samples)>1) {
                nsamples <- rowSums(stats[, -1]>=min.cells)
              } else {
                nsamples <- rep(0, dim(stats)[1])
                nsamples[which(as.numeric(stats[, 2])>=min.cells)] <- 1
              }

              genes.selected <- which(nsamples >= min.samples)

              if (length(genes.selected)<2) {
                stop("Less than two genes passed the prefiltering. Please select appropriate prefiltering criteria.")
              }

              object@data <- object@data[genes.selected, ]

              cat(length(genes.selected)," (", length(genes.selected)*100/length(nsamples), "%) genes passed the prefiltering criteria.\n", sep="")

              return(object)
          }
)

#' Analyzing batch difference
#'
#' @param object A sincera object
#' @param analysis A vector containing the analyses to be performed, including
#'        "q" - plotting the quantiles of gene expression in individual cells
#'        "qq" - Q-Q plot per sample pair
#'        "MA" - MA plot per sample pair
#'        "ngene" - Number of expressed genes in cells of different batches
#'        "distribution" - The distribution of expression values in cells of different batches
#' @param min.expression The threshold of expressed genes
#' @return The updated sincera object with batch analyses results
#'
setGeneric("batch.analysis", function(object, analysis=c("q", "qq", "ma","ngenes", "distribution"), min.expression=1, ...) standardGeneric("batch.analysis"))
#' @export
setMethod("batch.analysis","sincera",
          function(object, analysis=c("q", "qq", "ma","ngenes", "distribution"), min.expression=1, ...) {

            cellsample <- getCellMeta(object, name="GROUP")
            #expr <- getExpression(object)

            samples <- sort(unique(as.character(cellsample)))
            sample.pairs <- combn(samples, 2)

            if ("q" %in% analysis) { # plot quantiles of cells
                getExprQuantiles(object)
                pause()
            }

            if ("qq" %in% analysis) { # plot per sample pair Q-Q plot

                ES <- getES(object)
                for (i in 1:dim(sample.pairs)[2]) {

                    s1.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[1,i])
                    s2.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[2,i])

                    s1.label <- as.character(sample.pairs[1,i])
                    s2.label <- as.character(sample.pairs[2,i])

                    batch.analysis.QQ(ES, s1.cols, s2.cols, do.log2=TRUE, fig.filename=NULL, s1.label, s2.label, fig.main="Q-Q Plot")

                }
                rm(ES)
                pause()
            }

            if ("ma" %in% analysis) { # plot per sample pair MA plot
                ES <- getES(object)
                for (i in 1:dim(sample.pairs)[2]) {

                    s1.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[1,i])
                    s2.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[2,i])

                    s1.label <- as.character(sample.pairs[1,i])
                    s2.label <- as.character(sample.pairs[2,i])

                    batch.analysis.MA(ES, s1.cols, s2.cols, fig.filename=NULL, fig.main=paste("MA of ", s1.label, " vs. ", s2.label, sep=""))

                }
                rm(ES)
                pause()
            }

            if ("ngenes" %in% analysis) { # plot number of expressed genes
                ngenes <- colSums(getExpression(object)>min.expression)
                viz <- data.frame(nGenes=ngenes, Group=factor(getCellMeta(object,name="GROUP")))
                g <- ggplot(viz, aes(x=Group, y=nGenes))
                g <- g + geom_jitter(height=0, alpha=0.3) + geom_violin(aes(fill=Group), alpha=0.6, scale="width",adjust=0.75)
                g <- g + ggtitle(paste("Number of expressed genes per cell (min.exp=", min.expression,")", sep=""))
                g <- g + sincera_theme()
                g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
                g <- g + theme(panel.border = element_blank())
                g <- g + theme(axis.line.y=element_line(), axis.line.x=element_line())
                g <- g + theme(text=element_text(size=12))
                g <- g + theme(axis.title.x=element_text(size=rel(1.2)), axis.title.y=element_text(size=rel(1.2)), strip.text.x=element_text(size=rel(2)))
                print(g)
                pause()
            }

            if ("distribution" %in% analysis) {
                viz <- getExpression(object)
                viz <- as.data.frame(t(viz))
                viz$Group <- getCellMeta(object, name="GROUP")
                viz <- melt(viz, id.vars="Group")
                viz <- viz[which(viz$value>min.expression), ]
                viz$Expression <- log(viz$value, 2)

                if (FALSE) {
                  g <- ggplot(viz, aes(x=Group, y=Expression, fill=Group))
                  g <- g + geom_violin()
                  g <- g + ggtitle(paste("Distribution of expression values ( >", min.expression, ")", sep=""))
                  g <- g + ylab("Expression (log2)")
                  g <- g + sincera_theme()
                  g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
                  g <- g + theme(panel.border = element_blank())
                  g <- g + theme(axis.line.y=element_line(), axis.line.x=element_line())
                  g <- g + theme(text=element_text(size=12))
                  g <- g + theme(axis.title.x=element_text(size=rel(1.2)), axis.title.y=element_text(size=rel(1.2)), strip.text.x=element_text(size=rel(2)))

                  print(g)
                  pause()
                }

                g <- ggplot(viz, aes(x=Group, y=Expression, fill=Group))
                g <- g + geom_boxplot()
                g <- g + ggtitle(paste("Distribution of expression values", sep=""))
                g <- g + ylab("Expression (log2)")
                g <- g + sincera_theme()
                g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
                g <- g + theme(panel.border = element_blank())
                g <- g + theme(axis.line.y=element_line(), axis.line.x=element_line())
                g <- g + theme(text=element_text(size=12))
                g <- g + theme(axis.title.x=element_text(size=rel(1.2)), axis.title.y=element_text(size=rel(1.2)), strip.text.x=element_text(size=rel(2)))

                print(g)
                pause()
            }

            return(object)
          }
)

#' Get and plot the quantiles of gene expression per cell
#'
#' Enable analysis of batch difference by calculation and plotting the quantiles of cellular expression profiles
#' Colors indicate cells from different batches
#'
#' @param object A sincera object
#' @param probs a vector of quantiles to be computed
#' @param lower.bound The minimum value for quantile calculation
#' @param upper.bound The maximum value for quantile calculation
#' @param do.plot If true, visualize the quantiles in individual cells with colors indicating different batches
#' @return A data frame containing the quantile values of individual cells
#'
getExprQuantiles <- function(object, probs=c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 1),
                             lower.bound=-1, upper.bound=Inf, do.plot=T) {
    expr <- getExpression(object)
    cellname <- getCells(object)
    cellsample <- getCellMeta(object, name="GROUP")

    if (lower.bound >= upper.bound) {
        stop("invalid lower.bound or upper.bound")
    }
    if (any(probs>1) | any(probs<0)) {
        stop("Elements in probs must be between 0 and 1.")
    }

    quants <- t(apply(expr, 2, function(y) quantile.helper(as.vector(y), lower.bound=lower.bound, upper.bound=upper.bound, probs=probs)))
    colnames(quants) <- paste("Quantile ", probs, split="", sep="")
    quants <- data.frame(Cell=cellname, Group=cellsample, quants)

    if (do.plot==T) {
        qnm <- melt(quants, id.vars=c("Cell","Group"))
        colnames(qnm) <- c("Cell","Group","Quantile","Expression")
        g <- ggplot(qnm, aes(x=Cell, y=Expression)) + facet_wrap(~Quantile, scale="free_y")
        g <- g + geom_point(aes(col=Group))
        g <- g + ggtitle("Expression Quantiles")
        g <- g + sincera_theme()
        g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
        #g <- g + theme(text=element_text(size=12))
        #g <- g + theme(axis.title.x=element_text(size=rel(1.2)), axis.title.y=element_text(size=rel(1.2)), strip.text.x=element_text(size=rel(2)))

        print(g)
    }

    return(quants)
}




#' Performing zscore scaling of gene expression
#'
#' Perform global or per-sample zscore scaling of gene expression for cell cluster identification and expression pattern visualization
#'
#' @param object A sincera object
#' @param pergroup If true, performing per cell group scaling; otherwise, take all cells as from the same sample and perform scaling.
#' @return An updated sincera object with scaled expression values, use getExpresssion function with scaled=TRUE to access the scaled expression
#'
setGeneric("normalization.zscore", function(object, pergroup=TRUE, ...) standardGeneric("normalization.zscore"))
#' @export
setMethod("normalization.zscore","sincera",
          function(object, pergroup=TRUE, ...) {
            es <- NULL
            if (class(object)=="sincera") {
              es <- getES(object)
            } else if (class(object)=="ExpressionSet") {
              es <- object
            }

            if (TRUE==pergroup) {
              es <- normalization.zscore.old(es, group.by="GROUP", groups=NULL, verbose=T)
            } else {
              zexprs <- scale(t(exprs(es)), center=T, scale=T)
              exprs(es) <- t(zexprs)
            }

            if (class(object)=="sincera") {
              object <- setExpression(object, value=exprs(es), scaled=TRUE)
            } else if (class(object)=="ExpressionSet") {
              object <- es
            }

            return(object)
          }
)



