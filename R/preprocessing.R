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

                
                for (i in 1:dim(sample.pairs)[2]) {

                    s1.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[1,i])
                    s2.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[2,i])

                    s1.label <- as.character(sample.pairs[1,i])
                    s2.label <- as.character(sample.pairs[2,i])

                    batch.analysis.QQ(getExpression(object), s1.cols, s2.cols, do.log2=TRUE, fig.filename=NULL, s1.label, s2.label, fig.main="Q-Q Plot")

                }
                
                pause()
            }

            if ("ma" %in% analysis) { # plot per sample pair MA plot
                #ES <- getES(object)
                for (i in 1:dim(sample.pairs)[2]) {

                    s1.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[1,i])
                    s2.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[2,i])

                    s1.label <- as.character(sample.pairs[1,i])
                    s2.label <- as.character(sample.pairs[2,i])

                    batch.analysis.MA(getExpression(object), s1.cols, s2.cols, fig.filename=NULL, fig.main=paste("MA of ", s1.label, " vs. ", s2.label, sep=""))

                }
                #rm(ES)
                pause()
            }
            
            if ("isccd" %in% analysis) {

              for (i in 1:dim(sample.pairs)[2]) {
                
                s1.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[1,i])
                s2.cols <- which(getCellMeta(object, name="GROUP") == sample.pairs[2,i])
                
                s1.label <- as.character(sample.pairs[1,i])
                s2.label <- as.character(sample.pairs[2,i])
                
                batch.analysis.ISCCD(getExpression(object), s1.cols, s2.cols, s1.label="s1", s2.label="s2", file.prefix=NULL, fig.res=150)
              }
              
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
  
    quantile.helper <- function(y, lower.bound=-1, upper.bound=Inf, probs) {
      y <- y[which(y>lower.bound)]
      y <- y[which(y<upper.bound)]
      if (length(y)>0) {
        return(quantile(y, probs))
      }
      return(NA)
    }
  
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
              es <- normalization.zscore.1(es, group.by="GROUP", groups=NULL, verbose=T)
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

#' Per group zscore normalization of gene expression profiles
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the sample information
#' @param groups (character) samples for measuring abundancy
#' @param verbose (logical) verbose the function output
#' @return an ExpressionSet object with per-sample zscore normalized expression values
normalization.zscore.1 <- function(ES, group.by="GROUP", groups=NULL, verbose=TRUE) {
  if (verbose) {
    cat("Sincera: row-by-row per-group zscore normalization ...")
  }
  if (is.null(groups)) {
    groups <- sort(unique(pData(ES)[,group.by]))
  }
  for (i in groups) {
    i.cells <- rownames(subset(pData(ES), pData(ES)[,group.by] %in% i))
    i.fpkm <- exprs(ES)[,i.cells]
    i.fpkm <- apply(i.fpkm, 1, function(y) z.transform.helper(y))
    i.fpkm <- t(i.fpkm)
    exprs(ES)[,i.cells] <- i.fpkm
  }
  if (verbose) {
    cat("done\n\n")
  }
  return(ES)
}

z.transform.helper <- function(x) {
  x <- as.numeric(x)
  x.mu <- mean(x)
  x.sd <- sd(x)
  if (x.sd == 0) {
    x <- rep(0, length(x))
  } else {
    x <- (x-x.mu)/x.sd
  }
  return(x)
}



#' Analyze the inter-sample correlation and distance
#'
#' @param x The expression matrix; rows are genes, cells are columns.
#' @param s1.cols (numeric) the column indices for the first set of cells
#' @param s2.cols (numeric) the column indices for the second set of cells
#' @param s1.label (character) the name for the first set of cells
#' @param s2.label (character) the name for the second set of cells
#' @param fig.prefix (character) the file prefix of the output figure; if file.prefix is null, the figure will be plotted to the screen
#' @return NULL
batch.analysis.ISCCD <- function(x, s1.cols, s2.cols, s1.label="s1", s2.label="s2", file.prefix=NULL, fig.res=150) {
  
  common <- intersect(s1.cols, s2.cols)
  if (length(common)>0) {
    stop("invalide s1.cols or s2.cols")
  }
  
  cor.mat <- cor(x)
  
  s1.df <- data.frame(intra.max.cor=rep(-1, length(s1.cols)), inter.max.cor=rep(-1, length(s1.cols)), inter.cor=NA, inter.distance=NA)
  
  for (i in 1:length(s1.cols)) {
    s1.df$intra.max.cor[i] <- max(cor.mat[s1.cols[i], s1.cols[-i]])
    s1.df$inter.max.cor[i] <- max(cor.mat[s1.cols[i], s2.cols])
    s1.df$inter.cor[i] <- s1.df$inter.max.cor[i]
    s1.df$inter.distance[i] <- s1.df$intra.max.cor[i] - s1.df$inter.max.cor[i]
  }
  
  s1.df <- cbind(CELL=colnames(x)[s1.cols], s1.df)
  
  s2.df <- data.frame(intra.max.cor=rep(-1, length(s2.cols)), inter.max.cor=rep(-1, length(s2.cols)), inter.cor=NA, inter.distance=NA)
  
  for (i in 1:length(s2.cols)) {
    s2.df$intra.max.cor[i] <- max(cor.mat[s2.cols[i], s2.cols[-i]])
    s2.df$inter.max.cor[i] <- max(cor.mat[s2.cols[i], s1.cols])
    s2.df$inter.cor[i] <- s2.df$inter.max.cor[i]
    s2.df$inter.distance[i] <- s2.df$intra.max.cor[i] - s2.df$inter.max.cor[i]
  }
  
  s2.df <- cbind(CELL=colnames(x)[s2.cols], s2.df)
  
  if (!is.null(file.prefix)) {
    write.table(s1.df, file=paste(file.prefix, ".", s1.label, ".txt", sep=""), sep="\t", col.names=T, row.names=F)
    write.table(s2.df, file=paste(file.prefix, ".", s2.label, ".txt", sep=""), sep="\t", col.names=T, row.names=F)
  } else {
    write.table(s1.df, file=paste("ISSCD.", s1.label, ".txt", sep=""), sep="\t", col.names=T, row.names=F)
    write.table(s2.df, file=paste("ISSCD.", s2.label, ".txt", sep=""), sep="\t", col.names=T, row.names=F)
  }
  
  
  
  if (!is.null(file.prefix)) {
    tiff(file=paste(file.prefix, ".tif", sep=""), width = 4, height =2, units = "in", res = fig.res, compression="lzw", bg = "white")
  }
  
  gs <- list()
  
  corr.data <- data.frame(group=c(s1.label, s2.label),
                          MEAN=c(mean(s1.df$inter.cor), mean(s2.df$inter.cor)),
                          SD=c(sd(s1.df$inter.cor), sd(s2.df$inter.cor)),
                          SE=c(sd(s1.df$inter.cor)/sqrt(length(s1.df$inter.cor)), sd(s2.df$inter.cor)/sqrt(length(s2.df$inter.cor))))
  
  ymax <- max(1, max(corr.data$MEAN[1]+corr.data$SE[1], corr.data$MEAN[2]+corr.data$SE[2]))
  ymin <- min(0, min(corr.data$MEAN[1]-corr.data$SE[1], corr.data$MEAN[2]-corr.data$SE[2]))
  
  
  corr.data$ORDER = factor(corr.data$group, as.character(corr.data$group))
  
  i.g <- ggplot(corr.data, aes(x=ORDER, y=MEAN), fill=group)
  i.g <- i.g + geom_bar(position=position_dodge(), stat="identity")
  i.g <- i.g + geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.2)#, position=position_dodge(.9))
  i.g <- i.g + labs(x="SAMPLE", y="CORRELATION", title="Inter-Sample \nCell Correlation")
  i.g <- i.g + scale_y_continuous(limits=c(ymin, ymax))
  i.g <- i.g + sincera_theme()
  
  
  gs[[1]] <- i.g
  
  
  dist.data <- data.frame(group=c(s1.label, s2.label),
                          MEAN=c(mean(s1.df$inter.distance), mean(s2.df$inter.distance)),
                          SD=c(sd(s1.df$inter.distance), sd(s2.df$inter.distance)),
                          SE=c(sd(s1.df$inter.distance)/sqrt(length(s1.df$inter.distance)), sd(s2.df$inter.distance)/sqrt(length(s2.df$inter.distance))))
  
  ymax <- max(1, max(dist.data$MEAN[1]+dist.data$SE[1], dist.data$MEAN[2]+dist.data$SE[2]))
  ymin <- min(0, min(dist.data$MEAN[1]-dist.data$SE[1], dist.data$MEAN[2]-dist.data$SE[2]))
  
  dist.data$ORDER = factor(dist.data$group, as.character(dist.data$group))
  
  i.g <- ggplot(dist.data, aes(x=ORDER, y=MEAN), fill=group)
  i.g <- i.g + geom_bar(position=position_dodge(), stat="identity")
  i.g <- i.g + geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.2)#, position=position_dodge(.9))
  i.g <- i.g + labs(x="SAMPLE", y="CORRELATION", title="Inter-Sample \nCell Distance")
  i.g <- i.g + scale_y_continuous(limits=c(ymin, ymax))
  i.g <- i.g + sincera_theme()
  
  
  gs[[2]] <- i.g
  
  
  multiplot(plotlist=gs, cols=2)
  
  if (!is.null(file.prefix)) {
    dev.off()
  }
}


#' plot the QQ plot for two sets of cells
#'
#' @param x The expression matrix; rows are genes, cells are columns
#' @param s1.cols (numeric) the column indices for the first set of cells
#' @param s2.cols (numeric) the column indices for the second set of cells
#' @param s1.label (character) the name for the first set of cells
#' @param s2.label (character) the name for the second set of cells
#' @param do.log2 (logical) log2 transformation of mean expression values of each gene
#' @param fig.filename (character) the file name of the output figure; if fig.filename is null, the figure will be plotted to the screen
#' @param fig.main (character) the title of the figure
#' @return NULL
batch.analysis.QQ <- function(x, s1.cols, s2.cols, do.log2=TRUE, fig.filename=NULL, s1.label="s1", s2.label="S2", fig.main="QQ plot", fig.res=150) {
  
  common <- intersect(s1.cols, s2.cols)
  if (length(common)>0) {
    stop("invalide s1.cols or s2.cols")
  }
  
  q1 <- apply(x[,s1.cols],1,mean)
  q2 <- apply(x[,s2.cols],1,mean)
  
  if (do.log2==TRUE) {
    if (any(q1==0) || any(q2==0)) {
      q1 <- q1 + 1
      q2 <- q2 + 1
    }
    q1 <- log(q1,2)
    q2 <- log(q2,2)
  }
  
  if (!is.null(fig.filename)) {
    tiff(file=fig.filename, width = 4, height = 4, units = "in", res = fig.res, compression="lzw", bg = "white")
  }
  
  mq=qqplot(q1, q2, main=fig.main, xlab=paste( s1.label, " Quantiles", sep=""), ylab=paste( s2.label, " Quantiles", sep=""), col="red")
  segments(min(c(q1,q2)),min(c(q1,q2)),max(c(q1,q2)), max(c(q1,q2)), col="black")
  box(lwd=5)
  
  if (!is.null(fig.filename)) {
    dev.off()
  }
}


#' plot the MA plot for two sets of cells
#'
#' @param x The expression matrix; rows are genes, cells are columns
#' @param s1.cols (numeric) the column indices for the first set of cells
#' @param s2.cols (numeric) the column indices for the second set of cells
#' @param s1.label (character) the name for the first set of cells
#' @param s2.label (character) the name for the second set of cells
#' @param do.log (logical) log2 transformation of mean expression values of each gene
#' @param fig.filename (character) the file name of the output figure; if fig.filename is null, the figure will be plotted to the screen
#' @param fig.main (character) the title of the figure
#' @return NULL
batch.analysis.MA <- function(x, s1.cols, s2.cols, fig.filename=NULL, fig.main="MA plot", fig.res=150) {
  
  common <- intersect(s1.cols, s2.cols)
  if (length(common)>0) {
    stop("invalide s1.cols or s2.cols")
  }
  
  s1.means <- apply(x[,s1.cols],1,mean)
  s2.means <- apply(x[,s2.cols],1,mean)
  
  if (any(s1.means==0) || any(s2.means==0)) {
    s1.means <- s1.means + 1
    s2.means <- s2.means + 1
  }
  
  A <- 0.5*(log(s1.means,2)+log(s2.means,2))
  M <- log(s1.means, 2) - log(s2.means,2)
  
  #A <- 0.5*(log(apply(x[,s1.cols],1,mean),2)+log(apply(x[,s2.cols],1,mean),2))
  #M <- log(apply(x[,s1.cols],1,mean),2)-log(apply(x[,s2.cols],1,mean),2)
  
  x.min <- floor(min(A))
  x.max <- ceiling(max(A))
  y.min <- floor(min(M))
  y.max <- ceiling(max(M))
  
  if (!is.null(fig.filename)) {
    tiff(file=fig.filename, width = 4, height = 4, units = "in", res = fig.res, compression="lzw", bg = "white")
  }
  
  ma <- plot(A, M, col="green", main=fig.main, xlim=c(x.min,x.max), ylim=c(y.min, y.max), xlab="A", ylab="M", pch=19, cex.lab=1, cex=0.5)
  segments(x.min,0,x.max,0, col = "black", lwd=1)
  box(lwd=2)
  
  if (!is.null(fig.filename)) {
    dev.off()
  }
}




#' log transformation
#'
#' @param x The expression matrix; rows are genes, cells are columns
#' @param log.base (numeric) the base of the logarithm
#' @param min.impute (numeric) a numeric value to impute the min expression values
#' @return The nomalized expression matrix
normalization.log <- function(x, log.base=2, min.impute = 0.01, increase=1) {
  if (log.base<=1 | !(log.base%%1==0)) {
    stop("log.base must be an integer greater than 1")
  }
  if (min.impute<=0) {
    stop("min.impute must be greater than 0")
  }
  if (increase<0) {
    stop("increase must be greater than or equal to 0")
  }
  x[x < min.impute] <- min.impute
  x <- log(x+increase, log.base)
  return(x)
}


#' square root transformation
#'
#' @param x The expression matrix; rows are genes, cells are columns
#' @return The normalized expression matrix
normalization.sqrt <- function(x) {
  return(sqrt(x))
}


#' Trimmed mean normalization
#'
#' @param x The expression matrix; rows are genes, cells are columns
#' @param up (numeric) upper quantile to be removed for calculation of the trimmed mean
#' @param lo (numeric) lower quantile to be removed for calculation of the trimmed mean
#' @param do.scaling (logical) scale the trimmed mean normalized expression values
#' @return The normalized expression matrix
normalization.trimmed.mean <- function(x, up=0.9, lo=0.1, do.scaling=FALSE) {
  #x <- exprs(ES)
  trimmed.means <- apply(x, 2, function(y) trimmed.mean.helper(y, up=up, lo=lo) )
  x <- x/trimmed.means
  if (do.scaling==TRUE) {
    mean.trimmed.means <- mean(trimmed.means)
   x <- x*mean.trimmed.means
  }
  return(x)
}
trimmed.mean.helper <- function(x, up=0.9, lo=0.1) {
  x <- as.numeric(x)
  bounds <- as.numeric(quantile(x, probs=c(lo, up)))
  if (bounds[1] < bounds[2]) {
    x <- x[x>bounds[1]]
    if (length(x) > 0) {
      x <- x[x<bounds[2]]
      if (length(x) > 0) {
        x.mu <- mean(x)
        if (x.mu==0) {
          return(NA)
        }
        return(x.mu)
      }
    }
  }
  return(NA)
}
