##########################################################################################
# NOTE: This file contains SINCERA functions in version a10242015.
#       Some of the functions are utilized in the current version.
#       We are working on upgrading those functions.
#       Once it is done, this file will be removed.
##########################################################################################

dir.delim <- "/"
TIMESTAMP.FORMAT <- '%m%d%Y_%Hh%Mm%Ss'

GENE.SYMBOL.LABEL = "SYMBOL"					    # for labeling the column encoding the official gene symbol
SAMPLE.LABEL = "SAMPLE"                        		# for labeling the column encoding the sample information
CLUSTER.LABEL = "CLUSTER"  					   		# for labeling the column encoding the results of cluster assignment
TF.LABEL = "TF"										# for labeling the column encoding the transcription factors/cofactors


COMMON.TRESHOLD=5                                   # in common gene metric, expression > COMMON.TRESHOLD will considered as expressed
COMMON.PERCENTAGE=0.8                               # in common gene metric, genes that express in >COMMON.PRECENTAGE cluster cells will be considered as a common gene shared by cluster cells.
UNIQUE.RATIO=2                                      # parameters for unique gene metric
UNIQUE.QUANTILE=0.85								# parameters for unique gene metric


EXPR.SPECIFICITY.PREFIX <- "specificity_"      		# for labeling the columns encoding the results of expression specificity calculation for each group
EXPR.ABUNDANCY.PREFIX <- "abundancy_"          		# for labeling the columns encoding the results of expression abundance calculation
PREFILTERING.PREFIX <- "use_for_analysis"	   		# for labeling the column encoding the result of prefiltering


DIFF.EXPR.PREFIX = "test_"                     		# for labeling the columns encoding the results of differential expression test for each group
CELLTYPE.ENRICHMENT.PREFIX = "use_for_celltype_"	# for labeling the columns encoding the gene list used for cell type enrichment analysis for each group
MARKER.PREFIX = "use_as_marker_"					# for labeling the columns encoding the biomarkers for each group


COMMON.PREFIX="common_"								# for labeling the columns encoding the results of common metric for each group
UNIQUE.PREFIX="unique_"								# for labeling the columns encoding the results of unique metric for each group
TEST.STATS.METRIC.PREFIX = "test_statistic_"		# for labeling the columns encoding the results of test statistic metric for each group
SYN.SIM.PREFIX="synthetic_similiary_"				# for labeling the columns encoding the results of synthetic profile similarity metric for each group
TRAINSET.PREFIX="use_for_train_"					# for labeling the columns encoding the training instances for each group
TRAIN.CLASS.PREFIX="train_class_"					# for labeling the columns encoding the class of group specific training instances
TESTSET.PREFIX="use_for_test_"						# for labeling the columns encoding the testing instances for each group
SIGNATURE.PREFIX="sig_pred_"						# for labeling the columns encoding the results of signature prediction


TG.PREFIX="use_for_tg_" 							# for labeling the columns encoding the group specific candidate regulatory targets
TF.PREFIX="use_for_tf_"								# for labeling the columns encoding the group specific candidate transcription factors/cofactors


#######################################################
#                Preprocessing                        #
#######################################################

#' Truncate expression values that are out of bounds
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param upper.bound (numeric) the upper bound of expression
#' @param lower.bound (numeric) the lower bound of expression
#' @param verbose (logical) verbose the function output
#' @return an ExpressionSet object containing the truncated expression values

exprs.truncate <- function(ES, upper.bound=NULL, lower.bound=NULL, verbose=TRUE) {

    # truncate all values above a user-defined upper.bound
    if (!is.null(upper.bound)) {
        if (verbose) {
            cat("Sincera: converting expression >", upper.bound, "to", upper.bound, "... ")
        }
        expr <- exprs(ES)
        n <- length(which(expr>upper.bound))
        n.all <- dim(exprs(ES))[1]*dim(exprs(ES))[2]
        expr[expr > upper.bound] <- upper.bound
        exprs(ES) <- expr
        if (verbose) {
            cat("done\n", n, " (", 100*n/n.all, "%) expression values are affected.\n\n", sep="")
        }
    }

    # truncate all values below a user-defined lower.bound
    if (!is.null(lower.bound)) {
        if (verbose) {
            cat("Sincera: converting expression <", lower.bound, "to", lower.bound, "... ")
        }
        expr <- exprs(ES)
        n <- length(which(expr<lower.bound))
        n.all <- dim(exprs(ES))[1]*dim(exprs(ES))[2]
        expr[expr < lower.bound] <- lower.bound
        exprs(ES) <- expr
        if (verbose) {
            cat("done\n", n, " (", 100*n/n.all, "%) expression values are affected.\n\n", sep="")
        }
    }
    return(ES)
}

#' Per group expression specificity calculation
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the sample information
#' @param groups (character) samples for measuring specificity
#' @param specificity.prefix (character) the prefix for labeling the columns encoding the results of expression specificity calculation for each group
#' @return an ExpressionSet object with calculated specificities encoded in fData attributes
#' @details if a gene has zero expression across all cells, its specificity is zero

exprs.specificity <- function(ES, group.by=SAMPLE.LABEL, groups=NULL, specificity.prefix=EXPR.SPECIFICITY.PREFIX) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    for (i in groups) { # per group
        i.cells <- rownames(subset(pData(ES), pData(ES)[,group.by] %in% i))
        i.name <- paste(specificity.prefix,i,sep="")
        if (length(i.cells) <=1) {
            i.specificity <- rep(1, dim(exprs(ES))[1])
            fData(ES)[,i.name] <- i.specificity
        } else {
            i.specificity <- apply(exprs(ES[,i.cells]), 1, function(y) specificity.helper(y))
            fData(ES)[,i.name] <- i.specificity
        }
    }

    return(ES)
}
specificity.helper <- function(x) {
    x <- as.numeric(x)
    if (!any(x!=0)) { # if all zeros, return zero
        return(0)
    }
    specificity <- 0
    ma <- max(x)
    x <- x/ma
    x <- 1-x
    specificity <- sum(x)/(length(x)-1)
    return(specificity)
}


#' Plotting the distribution of specificity for each group
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the sample information
#' @param groups (character) samples for plotting the distrbution of specificity
#' @param ref.idx (numeric) a vectorcontaining the row indices of the set of reference genes, whose specificities will be plotting using green vertical lines
#' @param criterion (numeric) a specificity criterion that will be plotted using a red vertical line
#' @param specificity.prefix (character) the prefix labeling the columns encoding the results of expression specificity calculation for each group
#' @param fig.filename (character) the filename of the output figure; if the fig.filename is null, the figure will be plotted to the screen
#' @param fig.res (numeric) resolution of the figure
#' @return NULL

plot.exprs.specificity <- function(ES, group.by=SAMPLE.LABEL, groups=NULL, ref.idx=NULL, criterion=0.7, specificity.prefix=EXPR.SPECIFICITY.PREFIX, fig.filename=NULL, fig.res=300) {

    if (criterion > 1 | criterion <0) {
        stop("Criterion should be within 0 and 1.")
    }

    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[, group.by]))
    }
    n <- length(groups)


    if (!is.null(fig.filename)) {
        #fig.filename <- "Specificity_Distribution.tif"
        tiff(file=fig.filename, width=3, height=2*n, unit="in", res=fig.res, pointsize=2, compression="lzw")
    }

    gs <- list()

    for (i in 1:n) {
        i.name <- paste(specificity.prefix, groups[i], sep="")
        i.specificity <- data.frame(V1=as.numeric(fData(ES)[,i.name]))
        colnames(i.specificity) <- i.name
        i.g <- ggplot(i.specificity, aes_string(x=i.name))
        i.g <- i.g + geom_density();
        i.g <- i.g + scale_x_continuous(breaks = seq(0,1, 0.1))

        if (length(ref.idx) > 0) {
            #ref.cols <- rep("green", length(ref.idx))
            i.g <- i.g + geom_vline(xintercept=i.specificity[ref.idx, i.name], colour="green", show_guide=FALSE)
        }

        i.g <- i.g + geom_vline(xintercept=criterion, colour="red", show_guide=FALSE)

        i.g <- i.g + labs(title=paste("Distribution of Specificity\n(Sample", groups[i], ")"))
        i.g <- i.g + sincera_theme()
        gs[[i]] <- i.g
    }

    multiplot(plotlist=gs, cols=1)

    if (!is.null(fig.filename)) {
        dev.off()
    }
}


#' Determine the threshold for specificity filter based on a set of reference genes
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the sample information
#' @param groups (character) samples involves in determining the threshold
#' @param ref.idx (numeric) a vectorcontaining the row indices of the set of reference genes
#' @param exp.t (numeric) min expression value for the reference genes
#' @param exp.p (numeric) min expression percentage per sample for the reference genes
#' @param tolerance (numeric) percentage of reference genes that will be allowed to pass the specificity filters
#' @param step (numeric) a numeric value between 0 and 1, specifying the granularity of
#' @param do.plot (logical) plotting the distribution of specificity, the specificities of the reference genes (green lines), and the selected threshold (red line, if exists)
#' @param verbose (logical) verbose the function output
#' @return a numeric value between 0 and 1 if a threshold exists, otherwise, -1

specificity.criterion.selection <- function(ES, group.by=SAMPLE.LABEL, groups=NULL, ref.idx=NULL, exp.t=5, exp.p = 0.95, tolerance=0.05, step=0.1, do.plot=FALSE, verbose=T) {

    if (!is.null(ref.idx)) {
        if (verbose) {
            cat("Sincera: determining a threshold for specificity filter..\n")
        }
        n.ref <- length(ref.idx)

        if (is.null(groups)) {
            groups <- sort(unique(pData(ES)[, group.by]))
        }
        n <- length(groups)

        # keep reference genes with expression >= exp.t in at least exp.p percentage of cells per sample
        temp <- matrix(0, nrow=length(ref.idx), ncol=n)
        cnames <- paste("exp.", groups, sep="")
        colnames(temp) <- cnames
        ref.exp <- data.frame(IDX=ref.idx, temp)

        if (n==1) {
            col.idx <- 1:dim(pData(ES))[1]
            ref.exp[which(rowSums(exprs(ES)[ref.idx, col.idx] > exp.t) > floor(length(col.idx)*exp.p)), paste("exp.", groups, sep="")] <- 1
            ref.idx <- ref.idx[which(ref.exp[, cnames]>=length(groups))]
        } else if (n>1) { # more than one sample
            for (g in 1:n) {
                col.idx <- which(pData(ES)[, group.by] %in% groups[g])
                ref.exp[which(rowSums(exprs(ES)[ref.idx, col.idx] > exp.t) > floor(length(col.idx)*exp.p)), paste("exp.", groups[g], sep="")] <- 1
            }
            ref.idx <- ref.idx[which(rowSums(ref.exp[, cnames])>=length(groups))]
        }

        if (length(ref.idx) > 0) {

            if (verbose) {
                cat("\t", length(ref.idx), "/", n.ref, " reference genes passed the abundancy criteria and are used for determining the specificity threshold\n", sep="")
            }

            # measure the specificity of ref genes
            ES <- exprs.specificity(ES, group.by=group.by, groups=groups, specificity.prefix = "specificity_")

            specificity.cols <- paste("specificity_", groups, sep="")

            criterion.candidates <- seq(0, 1, by=step)

            specificity <- fData(ES)[ref.idx, specificity.cols]

            for (i in 1:length(criterion.candidates)) {
                if (n==1) {
                    if (length(which(as.numeric(specificity) < criterion.candidates[i])) > floor(length(ref.idx)*(1-tolerance))) {
                        if (do.plot) {
                            plot.exprs.specificity(ES, group.by=group.by, groups=groups, specificity.prefix="specificity_", fig.filename="Determine Specificity Criterion.tiff", ref.idx=ref.idx, criterion=criterion.candidates[i])
                        }
                        if (verbose) {
                            cat("Threshold determined:", criterion.candidates[i],"\n")
                        }
                        return(criterion.candidates[i])
                    }
                } else if (n>1) {
                    if (length((which(rowSums(specificity > criterion.candidates[i]) < length(groups)))) > floor(length(ref.idx)*(1-tolerance))) {
                        if (do.plot) {
                            plot.exprs.specificity(ES, group.by=group.by, groups=groups, specificity.prefix="specificity_", fig.filename="Determine Specificity Criterion.tiff", ref.idx=ref.idx, criterion=criterion.candidates[i])
                        }
                        if (verbose) {
                            cat("Threshold determined:", criterion.candidates[i],"\n")
                        }
                        return(criterion.candidates[i])
                    }
                }
            }
        } else {
            stop("No reference genes are ubiquitously. Please change the settings for exp.t or exp.p.")
        }
    } else {
        stop("Please set the indices of reference genes through ref.idx.\n")
    }

    if (verbose) {
        cat("No threshold found\n")
    }
    return(-1)
}




#' Per group expression abundancy calculation
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the sample information
#' @param groups (vector) samples for measuring abundancy
#' @param abundancy.threshold (numeric) the min value of expression
#' @param abundancy.prefix (character) the prefix for labeling the columns encoding the results of expression abundancy calculation for each group
#' @return an ExpressionSet object with calculated abundancies encoded in fData attributes

exprs.abundancy <- function(ES, group.by=SAMPLE.LABEL, groups=NULL, abundancy.threshold=5, abundancy.prefix=EXPR.ABUNDANCY.PREFIX) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    for (i in groups) {
        i.cells <- rownames(subset(pData(ES), pData(ES)[,group.by] %in% i))
        i.abundancy <- apply(exprs(ES[,i.cells]), 1, function(y) abundancy.helper(y, abundancy.threshold))
        i.name <- paste(abundancy.prefix,i,sep="")
        fData(ES)[,i.name] <- i.abundancy
    }
    return(ES)
}
abundancy.helper <- function(x, threshold=5) {
    x <- as.numeric(x)
    return(length(which(x>=threshold)))
}

#' Calculating the expression quantiles per cell
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param probs (numeric) a set of quantiles to be calculated
#' @param lower.bound (numeric) expression values greater than the lower.bound in each cell will be used for the per cell quantile calculation
#' @param lower.bound (numeric) expression values less than the upper.bound in each cell will be used for the per cell quantile calculation
#' @param do.plot (logical) plot the per cell quantiles
#' @param fig.filename (character) file name of the output figure; if fig.filename is null, the figure will be plotted to the screen
#' @return an ExpressionSet object with calculated quantiles encoded in fData attributes

expr.quantiles <- function(ES, probs=c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 1), lower.bound=-1, upper.bound=Inf, do.plot=FALSE, color.by="SAMPLE", fig.filename=NULL, fig.res=300) {

    if (lower.bound >= upper.bound) {
        stop("invalid lower.bound or upper.bound")
    }
    if (any(probs>1) | any(probs<0)) {
        stop("Elements in probs must be between 0 and 1.")
    }

    quants <- t(apply(exprs(ES), 2, function(y) quantile.helper(as.vector(y), lower.bound=lower.bound, upper.bound=upper.bound, probs=probs)))
    colnames(quants) <- paste("quantile_", probs*100, split="", sep="")
    pData(ES) <- cbind(pData(ES), quants)

    if (do.plot == TRUE) {
        np <- dim(quants)[2]
        nx <- dim(quants)[1]

        if (!is.null(fig.filename)) {
            # bug: unable to start tiff() device (probably due to the image size)
            tiff(file=fig.filename, width=4, height=1.2*np, unit="in", res=fig.res, pointsize=2, compression="lzw")
        }

        gs <- list()

        for (i in 1:np) {
            i.quants <- data.frame(CELL=rownames(quants), Expression=as.numeric(quants[,i]))
            i.quants[,color.by] <- as.factor(pData(ES)[,color.by])
            i.g <- ggplot(i.quants, aes(x=CELL, y=Expression))
            if (is.null(color.by)) {
                i.g <- i.g + geom_point()
            } else {
                i.g <- i.g + geom_point(aes_string(color=color.by))
            }
            i.g <- i.g + labs(xlab="CELL", title=paste("Quantile ", probs[i]*100, "%", sep=""))
            i.g <- i.g + sincera_theme()
            i.g <- i.g + theme(axis.text.x=element_blank())
            gs[[i]] <- i.g
        }

        multiplot(plotlist=gs, cols=1)

        if (!is.null(fig.filename)) {
            dev.off()
        }

    }
    return(ES)
}
quantile.helper <- function(y, lower.bound=-1, upper.bound=Inf, probs) {
    y <- y[which(y>lower.bound)]
    y <- y[which(y<upper.bound)]
    if (length(y)>0) {
        return(quantile(y, probs))
    }
    return(NA)
}


#' Gene prefiltering based on expression specificity and abundancy
#'
prefiltering <- function(ES,
                         group.by=SAMPLE.LABEL,
                         groups = NULL,
                         specificity.threshold=0.7,
                         abundancy.threshold=5,
                         abundancy.count=2,
                         mode = "union",
                         do.batch.analysis = FALSE,
                         specificity.prefix = EXPR.SPECIFICITY.PREFIX,
                         abundancy.prefix = EXPR.ABUNDANCY.PREFIX,
                         prefiltering.prefix=PREFILTERING.PREFIX,
                         fig.res=150,
                         verbose=T,
                         export=T, export.components="fd")
{
    modes <- c("union", "intersect")
    m.id <- pmatch(mode, modes)
    if (is.na(m.id)) {
        stop("invalid mode. should be union or intersect.")
    }
    mode <- modes[m.id]

    if (verbose) {
        cat("Sincera: prefiltering start ... \n")
    }
    wd <- paste(getwd(), dir.delim, "sincera.prefiltering.", getTimestamp(), dir.delim, sep="")
    dir.create(wd)

    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }

    if (verbose) {
        cat("\tMeasuring per-group specificity\n")
    }
    ES <- exprs.specificity(ES, group.by=group.by, groups=groups, specificity.prefix = specificity.prefix)

    filename=paste(wd,"distribution_of_specificity.tif", sep="")
    plot.exprs.specificity(ES, group.by=group.by, groups=groups, ref.idx=NULL, criterion=specificity.threshold, specificity.prefix=specificity.prefix, fig.filename=filename, fig.res=fig.res)


    if (verbose) {
        cat("\tAnalyzing cellular expression quantiles\n")
    }
    filename=paste(wd,"cellular_epxression_quantiles.tif", sep="")
    expr.quantiles(ES, probs=c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 1), lower.bound=-1, upper.bound=Inf, do.plot=TRUE, color.by=group.by, fig.filename=filename, fig.res=fig.res)

    if (verbose) {
        cat("\tMeasuring per-group abundancy\n")
    }
    ES <- exprs.abundancy(ES, group.by=group.by, groups=groups, abundancy.threshold=abundancy.threshold, abundancy.prefix=abundancy.prefix)

    groups <- sort(unique(pData(ES)[,group.by]))
    s_colnames <- paste(specificity.prefix, groups, sep="")
    a_colnames <- paste(abundancy.prefix, groups, sep="")
    s_cols <- which(colnames(fData(ES)) %in% s_colnames)
    a_cols <- which(colnames(fData(ES)) %in% a_colnames)

    fData(ES)[, prefiltering.prefix] <- apply(fData(ES), 1, function(y) prefiltering.helper(y, s_cols, a_cols, s_t=specificity.threshold, a_c=abundancy.count, mode=mode))

    if (length(groups)==1) {
        cat("\tOnly one sample is detected. Ignore batch analysis configurations\n")
    } else {
        if (do.batch.analysis == TRUE) {

            cat("\tGenerating Q-Q, MA, and ISCCD figures before and after pre-filtering...")

            sample.pairs <- combn(groups, 2)

            for (i in 1:dim(sample.pairs)[2]) {

                s1.cols <- which(pData(ES)[,group.by] == sample.pairs[1,i])
                s2.cols <- which(pData(ES)[,group.by] == sample.pairs[2,i])

                s1.label <- as.character(sample.pairs[1,i])
                s2.label <- as.character(sample.pairs[2,i])

                filename=paste(wd,"QQ.", s1.label, "vs", s2.label,".before.tif", sep="")
                batch.analysis.QQ(ES, s1.cols, s2.cols, do.log2=TRUE, fig.filename=filename, s1.label, s2.label, fig.main=paste("Q-Q of ", s1.label, " vs. ", s2.label, " (before)", sep=""), fig.res=fig.res)
                filename=paste(wd,"QQ.", s1.label, "vs", s2.label,".after.tif", sep="")
                batch.analysis.QQ(ES[which(fData(ES)[,prefiltering.prefix]==1),], s1.cols, s2.cols, do.log2=TRUE, fig.filename=filename, s1.label, s2.label, fig.main=paste("Q-Q of ", s1.label, " vs. ", s2.label, " (after)", sep=""), fig.res=fig.res)

                filename=paste(wd,"MA.", s1.label, "vs", s2.label,".before.tif", sep="")
                batch.analysis.MA(ES, s1.cols, s2.cols, fig.main=paste("MA of ", s1.label, " vs. ", s2.label, " (before)", sep=""), fig.filename=filename, fig.res=fig.res)
                filename=paste(wd,"MA.", s1.label, "vs", s2.label,".after.tif", sep="")
                batch.analysis.MA(ES[which(fData(ES)[,prefiltering.prefix]==1),], s1.cols, s2.cols, fig.main=paste("MA of ", s1.label, " vs. ", s2.label, " (after)", sep=""), fig.filename=filename, fig.res=fig.res)

                file.prefix=paste(wd,"ISCCD.before", sep="")
                batch.analysis.ISCCD(ES, s1.cols, s2.cols, s1.label, s2.label, file.prefix=file.prefix, fig.res=fig.res)
                file.prefix=paste(wd,"ISCCD.after", sep="")
                batch.analysis.ISCCD(ES[which(fData(ES)[,prefiltering.prefix]==1),], s1.cols, s2.cols, s1.label, s2.label, file.prefix=file.prefix, fig.res=fig.res)
            }
        }

        cat(" done\n")
    }

    if (verbose) {
        cat("\t", length(which(fData(ES)[,prefiltering.prefix]==1)), " profiles passed both filters\n", sep="")
    }

    if (export) {
        if (verbose) {
            cat("\tExporting results to", wd, "\n")
        }
        exportES(ES, prefix=wd, suffix=".prefiltering.txt", components=export.components)
    }

    if (verbose) {
        cat("Sincera: prefiltering completed\n\n")
    }
    return(ES)
}
prefiltering.helper <- function(x, s_cols, a_cols, s_t=0.7, a_c=2, mode="intersect") {
    if (mode=="intersect") {
        section <- 1
        if (any(as.numeric(x[s_cols])<s_t) | any(as.numeric(x[a_cols])<a_c)) {
            section <- 0
        }
    } else if (mode=="union") {
        section <- 0
        if (any(as.numeric(x[s_cols])>=s_t) & any(as.numeric(x[a_cols])>=a_c)) {
            section <- 1
        }
    }
    return(section)
}


#' Per group zscore normalization of gene expression profiles
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the sample information
#' @param groups (character) samples for measuring abundancy
#' @param verbose (logical) verbose the function output
#' @return an ExpressionSet object with per-sample zscore normalized expression values
normalization.zscore.old <- function(ES, group.by=SAMPLE.LABEL, groups=NULL, verbose=TRUE) {
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


#' log transformation
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param log.base (numeric) the base of the logarithm
#' @param min.impute (numeric) a numeric value to impute the min expression values
#' @return an ExpressionSet object with log transformed expression values
normalization.log.old <- function(ES, log.base=2, min.impute = 0.01, increase=1) {
    if (log.base<=1 | !(log.base%%1==0)) {
        stop("log.base must be an integer greater than 1")
    }
    if (min.impute<=0) {
        stop("min.impute must be greater than 0")
    }
    if (increase<0) {
        stop("increase must be greater than or equal to 0")
    }
    exprs(ES)[exprs(ES) < min.impute] <- min.impute
    exprs(ES) <- log(exprs(ES)+increase, log.base)
    return(ES)
}

#' square root transformation
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @return an ExpressionSet object with square root transformed expression values
normalization.sqrt.old <- function(ES) {
    exprs(ES) <- sqrt(exprs(ES))
    return(ES)
}

#' log transformation
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param up (numeric) upper quantile to be removed for calculation of the trimmed mean
#' @param lo (numeric) lower quantile to be removed for calculation of the trimmed mean
#' @param do.scaling (logical) scale the trimmed mean normalized expression values
#' @return an ExpressionSet object with trimmed mean normalized expression values
normalization.trimmed.mean.old <- function(ES, up=0.9, lo=0.1, do.scaling=FALSE) {
    #x <- exprs(ES)
    trimmed.means <- apply(exprs(ES), 2, function(y) trimmed.mean.helper(y, up=up, lo=lo) )
    exprs(ES) <- exprs(ES)/trimmed.means
    if (do.scaling==TRUE) {
        mean.trimmed.means <- mean(trimmed.means)
        exprs(ES) <- exprs(ES)*mean.trimmed.means
    }
    return(ES)
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

#' plot the QQ plot for two sets of cells
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param s1.cols (numeric) the column indices for the first set of cells
#' @param s2.cols (numeric) the column indices for the second set of cells
#' @param s1.label (character) the name for the first set of cells
#' @param s2.label (character) the name for the second set of cells
#' @param do.log2 (logical) log2 transformation of mean expression values of each gene
#' @param fig.filename (character) the file name of the output figure; if fig.filename is null, the figure will be plotted to the screen
#' @param fig.main (character) the title of the figure
#' @return NULL
batch.analysis.QQ <- function(ES, s1.cols, s2.cols, do.log2=TRUE, fig.filename=NULL, s1.label="s1", s2.label="S2", fig.main="QQ plot", fig.res=150) {

    common <- intersect(s1.cols, s2.cols)
    if (length(common)>0) {
        stop("invalide s1.cols or s2.cols")
    }

    x <- exprs(ES)

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
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param s1.cols (numeric) the column indices for the first set of cells
#' @param s2.cols (numeric) the column indices for the second set of cells
#' @param s1.label (character) the name for the first set of cells
#' @param s2.label (character) the name for the second set of cells
#' @param do.log (logical) log2 transformation of mean expression values of each gene
#' @param fig.filename (character) the file name of the output figure; if fig.filename is null, the figure will be plotted to the screen
#' @param fig.main (character) the title of the figure
#' @return NULL
batch.analysis.MA <- function(ES, s1.cols, s2.cols, fig.filename=NULL, fig.main="MA plot", fig.res=150) {

    common <- intersect(s1.cols, s2.cols)
    if (length(common)>0) {
        stop("invalide s1.cols or s2.cols")
    }

    x <- exprs(ES)

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

#' Analyze the inter-sample correlation and distance
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param s1.cols (numeric) the column indices for the first set of cells
#' @param s2.cols (numeric) the column indices for the second set of cells
#' @param s1.label (character) the name for the first set of cells
#' @param s2.label (character) the name for the second set of cells
#' @param fig.prefix (character) the file prefix of the output figure; if file.prefix is null, the figure will be plotted to the screen
#' @return NULL
batch.analysis.ISCCD <- function(ES, s1.cols, s2.cols, s1.label="s1", s2.label="s2", file.prefix=NULL, fig.res=150) {

    common <- intersect(s1.cols, s2.cols)
    if (length(common)>0) {
        stop("invalide s1.cols or s2.cols")
    }

    x <- exprs(ES)

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


#######################################################
#                  Clustering                         #
#######################################################


#' Assign cells to clusters
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param cluster.lable (character) the name of the column that will be used to encode the result of cluster assignment
#' @param clustering.method (character) the clustering method: hc - hierarchical clustering, tight - tight clustering, consensus - consensus clustering
#' @param target (numeric) the targeted number of clusters for tight clustering; required for tight clustering
#' @param k.min (numeric) required for tight clustering
#' @param min.are.increase (numeric) required for consensus clustering
#' @param maxK (numeric) parameter for consensus clustering
#' @param reps (numeric) parameter for consensus clustering
#' @param pItem (numeric) parameter for consensus clustering
#' @param pFeature (numeric) parameter for consensus clustering
#' @param clusterAlg (character) the clustering method used in the consensus clustering
#' @param distance.method (character) the distance method used in the hierarchical and consensus clustering: pearson - (1-pearons'correlation)/2, spearman - (1-spearman's correlation)/2, euclidean - euclidean distance
#' @param num.singleton (numeric) when hierarchical clustering is used, the number of singleton clusters that will be allowed in the generated clustering scheme
#' @param h (numeric) parameter for hierarchical clustering; the distance threshold that will be used to generate cell clusters
#' @param k (numeric) parameter for hierarchical clustering; the number of clusters that will be generated
#' @param do.shift (logical) parameter for hierarchical clustering
#' @param do.plot (logical) parameter for hierarchical clustering; plot the hierarchical tree
#' @param verbose (logical)
#' @param export (logical) export the clustering
#' @param export.components (character) the components of an ExpressionSet object that will be exported
#' @return a list of two elements: cell.cluster - a vector encoding the result of cluster assignment; cell.order - when hc is used, the cells in the order as they appeared in the hierarchical tree

cluster.assignment.old <- function(ES, cluster.label=CLUSTER.LABEL,
                                   clustering.method="hc", #"tight", "consensus"
                                   distance.method="pearson", #"spearman","euclidean"
                                   linkage.method="average", # parameters for hierarchical clustering
                                   num.singleton=0, h=NULL, k=NULL, do.shift=TRUE, do.plot=TRUE, # parameters for hierarchical clustering
                                   target=1, k.min=25, # required parameters for tight clustering
                                   min.area.increase=0.2, maxK=10, reps=10, pItem=0.8, pFeature=1, clusterAlg="hc", # parameters for consensus clustering
                                   verbose=TRUE, export=TRUE, export.components="pd") {

    cmethods <- c("tight", "hc","consensus")
    dmethods <- c("pearson","spearman","euclidean")
    cmethod.idx <- pmatch(clustering.method, cmethods)
    if (is.na(cmethod.idx)) {
        stop("invalid clustering method")
    }
    dmethod.idx = pmatch(distance.method, dmethods)
    if (is.na(dmethod.idx)) {
        stop("invalid distance method")
    }
    clustering.method <- cmethods[cmethod.idx]
    distance.method <- dmethods[dmethod.idx]

    if (verbose) {
        cat("Sincera: cluster assignment start ... \n")
    }
    if (clustering.method=="tight") {
        if (!require(tightClust)) {
            stop("The package 'tightClust' is required for tight clustering.")
        }
    } else if (clustering.method=="consensus") {
        if (!require(ConsensusClusterPlus)) {
            stop("The package 'ConsensusClusterPlus' is required for consensus clustering.")
        }
    }

    wd <- paste(getwd(), dir.delim, "sincera.cluster.assignment.", getTimestamp(), dir.delim, sep="")
    dir.create(wd)

    fpkm <- exprs(ES)
    clusters <- NULL
    hc.obj <- NULL
    hc.k <- 1

    if (clustering.method=="tight") {

        if (verbose) {
            cat("\tUsing tight clustering to find ", as.numeric(target), " cell clusters \n", sep="")
        }

        ret <- NULL
        ret <- tight.clust(t(fpkm), target=target, k.min=k.min)
        #if (do.plot==TRUE) {
        #    plot(ret)
        #}
        clusters <- ret$cluster
        names(clusters) <- rownames(pData(ES))

        pData(ES)[,cluster.label] <- clusters
        cell_order <- rownames(pData(ES))[order(pData(ES)[,cluster.label])]
        ES <- cluster.ordering(ES, col.order=cell_order, row.order=NULL, verbose=FALSE)
        cell_order <- rownames(pData(ES))

    } else if (clustering.method=="consensus") {

        if (verbose) {
            cat("Using consensus clustering", sep="")
        }

        ret <- NULL
        ret <- ConsensusClusterPlus(fpkm, maxK=maxK, reps=reps, pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg, innerLinkage=linkage.method,finalLinkage=linkage.method, distance=distance.method)

        kk <- length(ret)
        areaK <- c()
        breaks <- 100


        for (i in 2:kk) {
            v <- triangle(ret[[i]][["ml"]], mode=1)
            h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
            h$counts = cumsum(h$counts)/sum(h$counts)

            thisArea <- 0
            for (bi in 1:(length(h$breaks)-1)){
                thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
                bi = bi + 1
            }
            areaK = c(areaK,thisArea)
        }
        deltaK=areaK[1] #initial auc at k=2
        for(i in 2:(length(areaK))){
            #proportional increase relative to prior K.
            deltaK = c(deltaK,(areaK[i] - areaK[i-1])/areaK[i-1])
        }

        i <- 3
        while(deltaK[i] > min.area.increase && i<kk) {
            i <- i+1
        }
        #plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF curve",main="Delta area",type="b")

        clusters <- as.character(ret[[i]][["consensusClass"]])
        names(clusters) <- rownames(pData(ES))


        pData(ES)[,cluster.label] <- clusters
        cell_order <- rownames(pData(ES))[order(pData(ES)[,cluster.label])]
        ES <- cluster.ordering(ES, col.order=cell_order, row.order=NULL, verbose=FALSE)
        cell_order <- rownames(pData(ES))

    } else if (clustering.method=="hc") {

        if (do.shift) {
            fpkm <- apply(fpkm, 2, function(y) y-mean(y))
        }

        if (distance.method=="euclidean") {
            dd <- dist(t(fpkm), method=distance.method)
        } else  {
            dd <- as.dist((1-cor(fpkm, method=distance.method))/2)
        }

        hc <- hclust(dd, method=linkage.method)

        cell_order <- hc$labels[hc$order]

        if (is.null(k)) {
            if (is.null(h) || h<=0) {
                kk <- 1
                for (i in 2:dim(pData(ES))[1]) {
                    clusters <- cutree(hc, k=i)
                    clustersizes <- as.data.frame(table(clusters))
                    singleton.clusters <- which(clustersizes$Freq < 2)
                    if (length(singleton.clusters) <= num.singleton) {
                        kk <- i
                    } else {
                        break;
                    }
                }
                clusters <- cutree(hc, k=kk)
            } else {
                clusters <- cutree(hc, h=h)
                kk <- length(unique(clusters))
            }
        } else {
            if (k>0) {
                kk <- k
                clusters <- cutree(hc, k=kk)
            }
        }

        names(clusters) <- rownames(pData(ES))

        hc.obj <- hc
        hc.k <- kk

        if (do.plot == TRUE) {

            cut <- kk
            dendr <- dendro_data(hc, type = "rectangle")
            clust <- cutree(hc, k = cut)               # find 'cut' clusters
            clust.df <- data.frame(label = names(clust), cluster = clust)

            # Split dendrogram into upper grey section and lower coloured section
            height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
            cut.height <- mean(c(height[cut], height[cut-1]))
            dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend &
                                              dendr$segments$y > cut.height, 1, 2)
            dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)

            # Number the clusters
            dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
            change <- which(dendr$segments$cluster == 1)
            for (i in 1:cut) dendr$segments$cluster[change[i]] = i + 1
            dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1,
                                              ifelse(dendr$segments$cluster == 0, NA, dendr$segments$cluster))
            dendr$segments$cluster <- na.locf(dendr$segments$cluster)

            # Consistent numbering between segment$cluster and label$cluster
            clust.df$label <- factor(clust.df$label, levels = levels(dendr$labels$label))
            clust.df <- arrange(clust.df, label)
            clust.df$cluster <- factor((clust.df$cluster), levels = unique(clust.df$cluster), labels = (1:cut) + 1)
            dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")

            # Positions for cluster labels
            n.rle <- rle(dendr$segments$cluster)
            N <- cumsum(n.rle$lengths)
            N <- N[seq(1, length(N), 2)] + 1
            N.df <- dendr$segments[N, ]
            N.df$cluster <- N.df$cluster - 1

            # Plot the dendrogram
            g <- ggplot() +
                geom_segment(data = segment(dendr),
                             aes(x=x, y=y, xend=xend, yend=yend, size=factor(line), colour=factor(cluster)),
                             lineend = "square", show.legend = TRUE) +
                scale_colour_manual(values = c("grey60", rainbow(cut)), guide=guide_legend(title="Cluster")) +
                scale_size_manual(values = c(.1, 1), guide=FALSE) +
                #geom_text(data = N.df, aes(x = x, y = y, label = cluster,  colour = factor(cluster + 1)), hjust = 1.5, show_guide = FALSE) +
                geom_text(data = label(dendr), aes(x, y, label = label, colour = factor(cluster)), hjust = -0.2, size = 2, show.legend = FALSE) +
                scale_y_reverse(expand = c(0.2, 0)) +
                labs(x = NULL, y = NULL) +
                coord_flip() +
                theme(axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      panel.background = element_rect(fill = "white"),
                      panel.grid = element_blank())

            ggsave(g, file=paste(wd, "clusters.pdf", sep=""), scale=2)

        }

        if (verbose) {
            cat("\tCells are divided into", kk, "clusters\n")
            cat("\tCluster membership is encoded in pData(ES)$", cluster.label, "\n",sep="")
        }
    }

    if (export) {
        if (verbose) {
            cat("\tExporting results to", wd, "\n")
        }
        exportES(ES, prefix=wd, suffix=".clusterassignment.txt", components=export.components)
    }

    if (verbose) {
        cat("Sincera: cluster assignment completed\n\n")
    }
    #return(list(ES=ES, cell.cluster=y$CLUSTER, cell.order=cell_order))
    return(list(cell.cluster=clusters, cell.order=cell_order, hc.obj=hc.obj, hc.k=hc.k))
}

# this function is borrowed from the package: ConsensusClusterPlus
triangle = function(m,mode=1){
    #mode=1 for CDF, vector of lower triangle.
    #mode==3 for full matrix.
    #mode==2 for calcICL; nonredundant half matrix coun
    #mode!=1 for summary
    n=dim(m)[1]
    nm = matrix(0,ncol=n,nrow=n)
    fm = m

    nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half

    fm = t(nm)+nm
    diag(fm) = diag(m)

    nm=fm
    nm[upper.tri(nm)] = NA
    diag(nm) = NA
    vm = m[lower.tri(nm)]

    if(mode==1){
        return(vm) #vector
    }else if(mode==3){
        return(fm) #return full matrix
    }else if(mode == 2){
        return(nm) #returns lower triangle and no diagonal. no double counts.
    }
}



#' Ordering rows or columns according to a specified order
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param col.order (numeric) the names of columns in the new order
#' @param row.order (numeric) the names of rows in the new order
#' @param verbose (logical)
#' @return an ExpressionSet object with re-ordered data

cluster.ordering <- function(ES, col.order=NULL, row.order=NULL, verbose=TRUE) {

    if (!is.null(col.order)) {
        if (verbose) {
            cat("Sincera: ordering columns... ")
        }
        exprs.m <- exprs(ES)[,col.order]
        cells <- pData(ES)[col.order,]
        exprs(ES) <- exprs.m
        pData(ES) <- cells
        if (verbose) {
            cat("done\n")
        }
    }
    if (!is.null(row.order)) {
        if (verbose) {
            cat("Sincera: ordering rows... ")
        }
        exprs.m <- exprs(ES)[row.order,]
        genes <- fData(ES)[row.order,]
        exprs(ES) <- exprs.m
        fData(ES) <- genes
        if (verbose) {
            cat("done\n")
        }
    }
    return(ES)
}


#' Export the clustering resuls for visualization using external tools, e.g., GENE-E
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param signature (data.frame) a data.frame object containing cell type signature, shall contain at least two columns SYMBOL and TYPE
#' @param cell.order (numeric) the names of columns in the new order
#' @param verbose (logical)
#' @return an ExpressionSet object with re-ordered data

expression.export4viz <- function(ES, signature=NULL, cell.attrs=NULL, file.prefix="sincera.cluster.export", verbose=T) {
    if (verbose) {
        cat("Sincera: exporting expression for visualization ... \n")
    }
    wd <- paste(getwd(), dir.delim, "sincera.expression.export.", getTimestamp(), dir.delim, sep="")
    dir.create(wd)

    filename <- paste(wd, file.prefix, ".txt", sep="")

    if (is.null(signature)) {
        if (is.null(cell.attrs)) {
            header=t(data.frame(CID=rownames(pData(ES))))
        } else {
            header=t(data.frame(CID=rownames(pData(ES)), pData(ES)[,cell.attrs]))
        }
        expressions <- data.frame(SYMBOL=rownames(ES), exprs(ES))
    } else {
        signature <- signature[which(signature$SYMBOL %in% fData(ES)[,GENE.SYMBOL.LABEL]), ]
        signature <- signature[order(signature$TYPE, signature$SYMBOL), ]
        expressions <- data.frame(SYMBOL=fData(ES)[, GENE.SYMBOL.LABEL], exprs(ES))
        expressions <- merge(signature, expressions, by.x="SYMBOL", by.y="SYMBOL", sort=FALSE)
        expressions <- expressions[,which(!(colnames(expressions) %in% "TYPE"))]
        if (is.null(cell.attrs)) {
            header=t(data.frame(CID=rownames(pData(ES))))
        } else {
            header=t(data.frame(CID=rownames(pData(ES)), pData(ES)[,cell.attrs]))
        }
    }

    write.table(header, file=filename, col.names=F, row.names=T, sep="\t")
    write.table(expressions, file=filename, append=T, col.names=F, row.names=F, sep="\t")

    if (verbose) {
        cat("done\n\n")
    }
}



#' Permutation Analysis for determining significance of cluster assignments
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param n (numeric) the number of permutations
#' @param distance.method (character) the distance method: pearson or spearman - (1-correlation)/2; euclidean - euclidean distance
#' @param log.base (numeric) the base of logorithm; if log.base <=1 or log.base is NULL, no log transformations will be applied
#' @param verbose (logical)
#' @return NULL

cluster.permutation.analysis.old <- function(ES, group.by=CLUSTER.LABEL, n=20, distance.method="euclidean", log.base=2, verbose=TRUE) {

    if (!(n > 1)) {
        stop("invalid number of permutations. n should be an integer and greater than 1")
    }
    dmethods <- c("pearson", "euclidean")
    dmethod.id <- pmatch(distance.method, dmethods)
    if (is.na(dmethod.id)) {
        stop("invalid distance method")
    }
    distance.method <- dmethods[dmethod.id]
    if (is.null(log.base)) {
        log.base=1
    }
    if (log.base>1 & any(exprs(ES)<=0)) {
        log.base=1
        stop("There are zero or negative expression values. Please set log.base to NULL\n")
    }

    cat("Sincera: permutation analysis of cluster assignment\n")
    cat("\tgenerating", n, "random assignments for obtaining a background distribution\n")

    # order expression profile
    cell_order <- rownames(pData(ES))[order(pData(ES)[,group.by])]
    ES <- cluster.ordering(ES, col.order=cell_order, row.order=NULL, verbose=FALSE)

    # cluster size info
    cs <- as.data.frame(table(pData(ES)[, group.by]))
    colnames(cs) <- c(group.by, "SIZE")

    # vector to store quality scores
    qs <- rep(NA, n+1)

    # observed assignment
    oa <- 1:dim(pData(ES))[1]

    # quality of observed assignment
    qs[1] <- pa.helper(ES, group.by=group.by, log.base=log.base, cs, oa,  distance.method=distance.method)

    # generate random assignments and evaluate their quality score
    for (i in 2:(n+1)) {
        # obtains a random permutation
        ra <- sample(oa, replace=FALSE)
        qs[i] <- pa.helper(ES, group.by=group.by, log.base=log.base, cs, ra, distance.method=distance.method)
    }

    qs.sd <- sd(qs[-1])
    qs.mu <- mean(qs[-1])

    cat("\tThe quality scores of random assignments has a mean=", qs.mu, " and standard deviation=", qs.sd, "\n", sep="")
    cat("\tThe quality score of the input cluster assignment is", qs[1],"\n")

    # using approximated normal distribution to compute p-value
    p.value <- NA
    if (qs[1]>qs.mu) {
        p.value <- 1-pnorm(qs[1], mean=qs.mu, sd=qs.sd)
    } else {
        p.value <- pnorm(qs[1], mean=qs.mu, sd=qs.sd)
    }
    cat("\tThe p-value of the quality of the input cluster assignment is ", p.value, "\n")
    cat("sincera: permutation analysis completed\n\n")
}
pa.helper <- function(ES, group.by=CLUSTER.LABEL, log.base=2, cs, a, distance.method="euclidean") {
    s <- 0
    idx <- 1
    x.mu <- matrix(0, nrow=dim(exprs(ES))[1], ncol=length(cs[, group.by]))
    # obtain cluster centroids
    for (i in 1:length(cs[, group.by])) {
        i.x.mu <- apply(exprs(ES)[, a[idx:(idx+cs$SIZE[i]-1)]], 1, mean)
        i.s <- sum(apply(exprs(ES)[, a[idx:(idx+cs$SIZE[i]-1)]], 2, function(z) pa.distance.helper(mu=as.numeric(i.x.mu), y=as.numeric(z), log.base=log.base, distance.method=distance.method)))
        s <- s + i.s
        idx <- idx+cs$SIZE[i]
    }

    return(s)
}
pa.distance.helper <- function(mu, y, log.base=2, distance.method="euclidean") {
    d <- NULL
    if (log.base > 1) {
        mu <- log(mu, log.base)
        y <- log(y, log.base)
    }
    if (distance.method == "pearson" | distance.method=="spearman") {
        d <- (1-cor(mu, y, method="pearson"))/2
    } else if (distance.method == "euclidean") {
        d <- sqrt(sum((y-mu)^2))
    }
    return(d)
}


#######################################################
#              Differential Expression                #
#######################################################

#' Detecting Cluster Specific Differentially Expressed Genes
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for detecting cluster specific differentially expressed genes
#' @param method (character) differential test method: welch - Welch's t-test, wilcoxon - Wilcoxon rank sum test, samseq - SAMseq
#' @param do.fdr (logical) performing B&H fdr when method is welch or wilcoxon
#' @param samseq.fdr (numeric) parameter for SAMseq
#' @param samseq.nperms (numeric) parameter for SAMseq
#' @param samseq.nresamp (numeric) parameter for SAMseq
#' @param diffexpr.prefix (character) the prefix to label the columns encoding the results of differential test
#' @param verbose (logical)
#' @param export (logical) wheterh to export the ExpressionSet
#' @param export.components (character) the components of an ExpressionSet object that will be exported
#' @return an ExpressionSet object containing the results of differential tests in the attributes of fData

diff.test <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, method="welch", do.fdr=FALSE, samseq.fdr=0.2, samseq.nperms=10, samseq.nresamp=20, diffexpr.prefix=DIFF.EXPR.PREFIX, verbose=T, export=T, export.components="fd") {

    tmethods <- c("welch","wilcoxon", "samseq")
    tm.id <- pmatch(method, tmethods)
    if (is.na(tm.id)) {
        stop("invalid test method")
    }
    method <- tmethods[tm.id]

    if (method=="samseq") {
        if (!require(samr)) {
            stop("The R package 'samr' is required for SAMseq based differential expression test")
        }
    }
    if (verbose) {
        cat("Sincera: differential expression test ... \n")
    }
    wd <- paste(getwd(), dir.delim, "sincera.diffexpr.", getTimestamp(), dir.delim, sep="")
    dir.create(wd)

    if (method=="welch") {
        ES <- diff.test.welch(ES, group.by=group.by, groups=groups, do.fdr=do.fdr, diffexpr.prefix=diffexpr.prefix, verbose=verbose)
    } else if (method=="wilcoxon") {
        ES <- diff.test.wilcoxon(ES, group.by=group.by, groups=groups, do.fdr=do.fdr, diffexpr.prefix=diffexpr.prefix, verbose=verbose)
    } else if (method=="samseq") {
        ES <- diff.test.samseq(ES, group.by=group.by, groups=groups, samseq.fdr=samseq.fdr, samseq.nperms=sameseq.nperms, samseq.nresamp=samseq.nresamp, diffexpr.prefix=diffexpr.prefix, verbose=verbose)
    }

    if (export) {
        if (verbose) {
            cat("\tExporting results to", wd, "\n")
        }
        exportES(ES, prefix=wd, suffix=".difftest.txt", components=export.components)
    }

    if (verbose) {
        cat("Sincera: differential test completed\n\n")
    }
    return(ES)
}

diff.test.welch <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, do.fdr=FALSE, diffexpr.prefix=DIFF.EXPR.PREFIX, verbose=T) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    cells <- rownames(pData(ES))
    for (i in groups) {
        if (verbose) {
            cat("\tDifferential expression testing for group", i, "using Welch's t-test ...")
        }
        i.cells <- rownames(subset(pData(ES), pData(ES)[, group.by] %in% i))
        i.cells.o <- setdiff(cells, i.cells)
        if (length(i.cells) > 1 & length(i.cells.o)>1) {
            i.cells.idx <- which(colnames(exprs(ES)) %in% i.cells)
            i.cells.o.idx <- which(colnames(exprs(ES)) %in% i.cells.o)


            i.test <- apply(exprs(ES), 1, function(y) welch.test.helper(a=y, idx.i=i.cells.idx, idx.o=i.cells.o.idx))
            i.name <- paste(diffexpr.prefix,i,sep="")
            fData(ES)[,i.name] <- i.test

            if (do.fdr == TRUE) {
                i.fdr.name <- paste(diffexpr.prefix, "fdr_", i, sep="")
                fData(ES)[,i.fdr.name] <- p.adjust(i.test, method="fdr")
            }
        }
        if (verbose) {
            cat(" done\n")
        }
    }
    return(ES)
}
welch.test.helper <- function(a, idx.i, idx.o) {
    a <- as.numeric(a)
    a_p <- 1
    if (sd(a[idx.i]) == 0 && sd(a[idx.o]) == 0 ) {
        a_p <- 1
    } else {
        a_p <- t.test(a[idx.i], a[idx.o], alternative="greater", paired = FALSE, var.equal = FALSE)$p.value
    }
    return(a_p)
}

diff.test.wilcoxon <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, do.fdr=FALSE, diffexpr.prefix=DIFF.EXPR.PREFIX, verbose=T) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    cells <- rownames(pData(ES))
    for (i in groups) {
        if (verbose) {
            cat("\tDifferential expression testing for group", i, "using Wilcoxon rank sum test ...")
        }
        i.cells <- rownames(subset(pData(ES), pData(ES)[, group.by] %in% i))
        i.cells.o <- setdiff(cells, i.cells)
        if (length(i.cells)>1 & length(i.cells.o)>1) {
            i.cells.idx <- which(colnames(exprs(ES)) %in% i.cells)
            i.cells.o.idx <- which(colnames(exprs(ES)) %in% i.cells.o)

            i.test <- apply(exprs(ES), 1, function(y) wilcoxon.test.helper(a=y, idx.i=i.cells.idx, idx.o=i.cells.o.idx))
            i.name <- paste(diffexpr.prefix,i,sep="")
            fData(ES)[,i.name] <- i.test

            if (do.fdr == TRUE) {
                i.fdr.name <- paste(diffexpr.prefix, "fdr_", i, sep="")
                fData(ES)[,i.fdr.name] <- p.adjust(i.test, method="fdr")
            }
        }
        if (verbose) {
            cat(" done\n")
        }
    }
    return(ES)
}
wilcoxon.test.helper <- function(a, idx.i, idx.o) {
    a <- as.numeric(a)
    a_p <- 1
    if (sd(a[idx.i]) == 0 && sd(a[idx.o]) == 0 ) {
        a_p <- 1
    } else {
        a_p <- wilcox.test(a[idx.i], a[idx.o], alternative="greater")$p.value
    }
    return(a_p)
}

diff.test.samseq <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, samseq.fdr=0.2, samseq.nperms=10, samseq.nresamp=20, diffexpr.prefix=DIFF.EXPR.PREFIX, verbose=T) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    cells <- rownames(pData(ES))
    for (i in groups) {
        if (verbose) {
            cat("\nDifferential expression testing for group", i, "using SAMseq ...\n")
        }
        i.cells <- rownames(subset(pData(ES), pData(ES)[, group.by] %in% i))
        i.cells.o <- setdiff(cells, i.cells)
        if (length(i.cells)>1 & length(i.cells.o)>1) {
            i.cells.idx <- which(colnames(exprs(ES)) %in% i.cells)
            i.cells.o.idx <- which(colnames(exprs(ES)) %in% i.cells.o)

            x <- exprs(ES)
            x <- ceiling(x)
            y <- rep(-1, length(cells))
            y[i.cells.idx] <- 2
            y[i.cells.o.idx] <- 1

            samfit <- SAMseq(x, y, resp.type = "Two class unpaired", nperms = 10, random.seed = NULL, nresamp = 20, fdr.output = samseq.fdr)

            i.tt <- samfit$samr.obj$tt
            i.tt.name <- paste(diffexpr.prefix,i,sep="")
            fData(ES)[,i.tt.name] <- i.tt

            i.fc <- samfit$samr.obj$foldchange
            i.fc.name <- paste(diffexpr.prefix,"fc_", i,sep="")
            fData(ES)[,i.fc.name] <- i.fc

            if (samfit$siggenes.table$ngenes.up>0) {
                i.up.name <- paste(diffexpr.prefix,"genesup_", i,sep="")
                fData(ES)[,i.up.name] <- 0
                fData(ES)[as.numeric(samfit$siggenes.table$genes.up[,2]), i.up.name] <- 1
            }
            if (samfit$siggenes.table$ngenes.lo>0) {
                i.lo.name <- paste(diffexpr.prefix,"geneslo_", i,sep="")
                fData(ES)[,i.lo.name] <- 0
                fData(ES)[as.numeric(samfit$siggenes.table$genes.lo[,2]), i.lo.name] <- 1
            }
        }
        if (verbose) {
            cat(" done\n")
        }
    }
    return(ES)
}

#' Retrieve Differentially Expressed Genes
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for selecting cluster specific differentially expressed genes
#' @param diffexpr.prefix (character) the prefix to label the columns encoding the results of differential test
#' @param threshold (numeric) selecting diff genes <threshold
#' @return a data frame containing the retreived differentially expressed genes

get.diff.genes <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, diffexpr.prefix=DIFF.EXPR.PREFIX, threshold=0.05) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }

    groups <- sort(unique(pData(ES)[,group.by]))
    diff.genes <- data.frame(CLUSTER=NULL, GID=NULL, SYMBOL=NULL, PVALUE=NULL)
    for (i in groups) {
        i.diff.col <- paste(diffexpr.prefix, i, sep="")
        i.diff.rows <- which(fData(ES)[,i.diff.col] < threshold)
        i.diff.genes <- data.frame(CLUSTER=as.character(i), GID=rownames(fData(ES))[i.diff.rows], SYMBOL=fData(ES)[i.diff.rows,GENE.SYMBOL.LABEL], PVALUE=fData(ES)[i.diff.rows, i.diff.col])
        diff.genes <- rbind(diff.genes, i.diff.genes)
    }

    diff.genes <- diff.genes[order(diff.genes$CLUSTER, diff.genes$PVALUE, decreasing=FALSE),]

    cat("\nThe Number of Diff Genes in Each Cluster\n", sep="")
    tt <- table(diff.genes$CLUSTER)
    print(tt)
    cat("\n")

    return(diff.genes)
}

#######################################################
#          Cell Type Enrichment Analysis              #
#######################################################



#' Contruct the knowledge base for cell type enrichment analysis
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param KB (list) the knowledge base prepared by celltype.enrichment.initKB()
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for cell type enrichment
#' @param species (character) the species of the genome: MUSMU - mouse, HOMSA - human
#' @param id.type (character) the type of ids of the genes in the genome: ENSEMBL - ENSEMBL gene id, SYMBOL -  Entrez Gene Symbol, EG - Entrez Gene Id
#' @param celltype.enrichment.prefix (character) the prefix of columns encoding the cluster-specific gene list for cell type enrichment analysis
#' @param verbose (logical)
#' @return NULL
#'
celltype.enrichment.old <- function(ES, KB=NULL, group.by="CLUSTER", groups=NULL, species="MUSMU", id.type="SYMBOL", celltype.enrichment.prefix=CELLTYPE.ENRICHMENT.PREFIX, verbose=TRUE) {

    if (verbose) {
        cat("Sincera: cell type enrichment analysis ... \n")
    }
    wd <- paste(getwd(), dir.delim, "sincera.celltype.enrichment.", getTimestamp(), dir.delim, sep="")
    dir.create(wd)

    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }

    x <- KB$celltype.gene.association
    types.count <- KB$celltype.genome.count
    genome.type <- paste(species,".", id.type, sep="")

    ret <- list()

    for (i in groups) {
        if (verbose) {
            cat("\tEnriching cell types for group", i, "...")
        }
        i.de.name <- paste(celltype.enrichment.prefix, i, sep="")
        i.de <- rownames(fData(ES))[which(fData(ES)[, i.de.name]==1)]

        if (!(genome.type == "MUSMU.ENSEMBL")) { # if not MUSMU.ENSEMBL, map to MUSMU.ENSEMBL
            i.de <- unique(KB$genome[which(KB$genome[, genome.type] %in% i.de),"MUSMU.ENSEMBL"])
        }

        i.x.idx <- NULL
        if (length(i.de) > 0) {
            i.x.idx <- which(x$GENE.ID %in% i.de)
        }

        if (length(i.x.idx) > 0 ) {
            i.x <- x[i.x.idx,]
            i.types <- unique(i.x$ANNOTATION.ID)
            i.types.count <- data.frame(TYPE=i.types, CL.1=NA, CL.0=NA, ALL.1=NA, ALL.0=NA, Fisher.PV=NA, Hits_SYMBOL=NA, Hits_MUSMU_ENSEMBL=NA)
            h <- 1

            for (j in 1:dim(i.types.count)[1]) {
                j.type <- i.types.count$TYPE[j]
                i.j.idx <- which(i.x$ANNOTATION.ID %in% j.type)
                i.j.de <- unique(i.x$GENE.ID[i.j.idx])
                j.idx <- which(types.count$TYPE %in% j.type)
                i.types.count$CL.1[j] <- length(i.j.de)
                i.types.count$ALL.1[j] <- types.count$COUNT[j.idx]
                i.types.count$CL.0[j] <- length(i.de)-length(i.j.idx)
                i.types.count$ALL.0[j] <- length(unique(x$GENE.ID)) - types.count$COUNT[j.idx]
                f <- fisher.test(matrix(c(i.types.count$CL.1[j],i.types.count$CL.0[j],i.types.count$ALL.1[j],i.types.count$ALL.0[j]), nrow=2), alternative="greater")
                i.types.count$Fisher.PV[j] <- f$p.value
                if (genome.type =="MUSMU.ENSEMBL") {
                    i.j.de.originalid <- i.j.de
                } else {
                    i.j.de.originalid <- unique(KB$genome[which(KB$genome[,"MUSMU.ENSEMBL"] %in% i.j.de),genome.type])
                }
                if (!is.null(fData(ES)[, GENE.SYMBOL.LABEL])) { #%%
                    i.j.de.symbols <- unique(fData(ES)[which(rownames(fData(ES)) %in% i.j.de.originalid), GENE.SYMBOL.LABEL]) #%%
                    i.types.count$Hits_MUSMU_ENSEMBL[j] <- paste(i.j.de, collapse=",")
                }
                i.types.count$Hits_SYMBOL[j] <- paste(i.j.de.symbols, collapse=",")
                h <- h+1
            }
        }
        if (verbose) {
            cat("done\n")
        }
        if (verbose) {
            cat("\tExporting enrichment results of group", i, "to", wd,"... ")
        }

        i.types.count <- i.types.count[order(i.types.count$Fisher.PV), ]
        ret[[as.character(i)]] <- i.types.count

        write.table(i.types.count, file=paste(wd, i, "-celltype-enrichment.txt", sep=""), sep="\t", col.name=T, row.name=F)
        if (verbose) {
            cat("done\n")
        }
    }

    if (verbose) {
        cat("Sincera: cell type enrichment completed\n\n")
    }

    return(ret)
}





#######################################################
#                 Signature Prediction                #
#######################################################


#' Cell type specific signature prediction
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param trainset.prefix (character) prefix for labeling the columns encoding the training instances for each group
#' @param train.class.prefix (character) prefix for labeling the columns encoding the class of group specific training instances
#' @param marker.prefix (character) prefix for labeling the columns encoding the biomarkers for each group
#' @param testset.prefix (character) prefix for labeling the columns encoding the testing instances for each group
#' @param common.prefix (character) prefix for labeling the columns encoding the results of common metric for each group
#' @param common.threshold (numeric) in common gene metric, expression > COMMON.TRESHOLD will considered as expressed
#' @param common.percentage (numeric) in common gene metric, genes that express in >COMMON.PRECENTAGE cluster cells will be considered as a common
#' @param unique.prefix (character) prefix for labeling the columns encoding the results of unique metric for each group
#' @param unique.ratio (numeric) parameter for unique gene metric
#' @param unique.quantile (numeric) parameters for unique gene metric
#' @param test.statistic.metric.prefix (character) prefix for labeling the columns encoding the results of test statistic metric for each group
#' @param diff.expr.prefix (character) prefix for labeling the columns encoding the results of differential test
#' @param diff.expr.threshold (numeric) genes with differential expression p-value<0.05 will be selected as candidates for signature
#' @param syn.sim.prefix (character) prefix for labeling the columns encoding the results of synthetic profile similarity metric for each group
#' @param signature.prefix (character) prefix for labeling the columns encoding the results of signature prediction
#' @param gene.symbol.label (character) the name of the column encoding gene symbols
#' @param verbose (logical)
#' @param export (logical) wheterh to export the ExpressionSet
#' @param export.components (character) the components of an ExpressionSet object that will be exported
#' @return an ExpressionSet object containing the results of signature prediction in the attributes of fData
#'
signature.analysis <- function(ES, group.by = CLUSTER.LABEL, groups=c(1,2,3,5,7,8,9),
                               # training set
                               trainset.prefix = TRAINSET.PREFIX,
                               train.class.prefix = TRAIN.CLASS.PREFIX,
                               marker.prefix=MARKER.PREFIX,
                               # testing set
                               testset.prefix=TESTSET.PREFIX,
                               # common gene metric
                               common.prefix = COMMON.PREFIX,
                               common.threshold=COMMON.TRESHOLD, common.percentage=COMMON.PERCENTAGE,
                               # unique gene metric
                               unique.prefix = UNIQUE.PREFIX,
                               unique.ratio=UNIQUE.RATIO, unique.quantile=UNIQUE.QUANTILE,
                               # test statistic metric
                               test.statistic.metric.prefix = TEST.STATS.METRIC.PREFIX, log.base=2,
                               diff.expr.prefix=DIFF.EXPR.PREFIX,
                               diff.expr.threshold = 0.05,
                               # synthetic profile similarity metric
                               syn.sim.prefix = SYN.SIM.PREFIX,
                               signature.prefix = SIGNATURE.PREFIX,
                               gene.symbol.label = GENE.SYMBOL.LABEL,
                               verbose=TRUE,
                               export=TRUE, export.components="fd", dir.prefix=NULL ) {

    if (verbose) {
        cat("Sincera: gene signature prediction ... \n")
    }

    wd <- NULL
    if (export == TRUE) {
        if (is.null(dir.prefix)) {
            wd <- paste(getwd(), dir.delim, "sincera.signature.", getTimestamp(), dir.delim, sep="")
        } else {
            wd <- paste(dir.prefix, dir.delim, sep="")
        }

        dir.create(wd)
    }

    # calculate metrics
    ES <- exprs.test.statistic.metric(ES, group.by=group.by, groups=groups, log.base=log.base, test.statistic.metric.prefix=test.statistic.metric.prefix, diff.expr.prefix=diff.expr.prefix)

    ES <- exprs.common(ES, group.by=group.by, groups=groups, common.threshold=common.threshold, common.percentage=common.percentage, common.prefix=common.prefix)

    ES <- exprs.unique(ES, group.by=group.by, groups=groups, unique.ratio=unique.ratio, unique.quantile=unique.quantile, unique.prefix=unique.prefix)

    ES <- exprs.synthetic.similarity(ES, group.by=group.by, groups=groups, syn.sim.prefix=syn.sim.prefix)

    # construct training sets
    ES <- sigpred.training.sets(ES, group.by=group.by, groups=groups,
                                trainset.prefix=trainset.prefix,
                                train.class.prefix=train.class.prefix,
                                test.statistic.metric.prefix = test.statistic.metric.prefix,
                                common.prefix = common.prefix,
                                unique.prefix = unique.prefix,
                                syn.sim.prefix = syn.sim.prefix,
                                marker.prefix = marker.prefix
    )
    # construct testing sets
    ES <- sigpred.testing.sets(ES, group.by=group.by, groups=groups, testset.prefix=testset.prefix, diff.expr.prefix=diff.expr.prefix, diff.expr.threshold=diff.expr.threshold)

    # performing signature prediction
    ES <- signature.prediction(ES, group.by=group.by, groups=groups,
                               trainset.prefix=trainset.prefix,
                               train.class.prefix=train.class.prefix,
                               testset.prefix=testset.prefix,
                               test.statistic.metric.prefix = test.statistic.metric.prefix,
                               common.prefix = common.prefix,
                               unique.prefix = unique.prefix,
                               syn.sim.prefix = syn.sim.prefix,
                               signature.prefix = signature.prefix,
                               gene.symbol.label = gene.symbol.label,
                               verbose=verbose,
                               export=export,
                               dir.prefix = wd
    )
    if (export) {
        if (verbose) {
            cat("\tExporting results to", wd, "\n")
        }
        exportES(ES, prefix=wd, suffix=".sigpred.txt", components=export.components)
    }

    if (verbose) {
        cat("Sincera: gene signature prediction completed\n\n")
    }
    return(ES)

}



#' Identifying group specific common genes
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param common.prefix (character) prefix for labeling the columns encoding the results of common metric for each group
#' @param common.threshold (numeric) in common gene metric, expression > COMMON.TRESHOLD will considered as expressed
#' @param common.percentage (numeric) in common gene metric, genes that express in >COMMON.PRECENTAGE cluster cells will be considered as a common
#' @return an ExpressionSet object containing the results of common gene metrics in the attributes of fData
#'
exprs.common <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, common.threshold=5, common.percentage=0.8, common.prefix=COMMON.PREFIX) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    for (i in groups) {
        i.cells <- rownames(subset(pData(ES), pData(ES)[, group.by] %in% i))
        i.cells.idx <- which(colnames(exprs(ES)) %in% i.cells)
        i.common <- apply(exprs(ES), 1, function(y) common.helper(y, i.cells.idx, common.threshold, common.percentage))
        i.name <- paste(common.prefix, i, sep="")
        fData(ES)[,i.name]<-i.common
    }
    return(ES)
}
common.helper <- function(a, idx.i, common.threshold=5, common.percentage=0.8) {
    a <- as.numeric(a)
    a_common_threshold <- ceiling(length(idx.i)*common.percentage)
    a_common <- 0
    if (length(which(a[idx.i] >= common.threshold)) >= a_common_threshold) {
        a_common <- 1
    }
    return(a_common)
}


#' Identifying group specific unique genes
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param unique.prefix (character) prefix for labeling the columns encoding the results of unique metric for each group
#' @param unique.ratio (numeric) parameter for unique gene metric
#' @param unique.quantile (numeric) parameters for unique gene metric
#' @return an ExpressionSet object containing the results of unique gene metric in the attributes of fData
#'
exprs.unique <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, unique.ratio=2, unique.quantile=0.85, unique.prefix=UNIQUE.PREFIX) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    cells <- rownames(pData(ES))
    for (i in groups) {
        i.cells <- rownames(subset(pData(ES), pData(ES)[, group.by] %in% i))
        i.cells.o <- setdiff(cells, i.cells)
        i.cells.idx <- which(colnames(exprs(ES)) %in% i.cells)
        i.cells.o.idx <- which(colnames(exprs(ES)) %in% i.cells.o)
        i.unique <- apply(exprs(ES), 1, function(y) unique.helper(y, i.cells.idx, i.cells.o.idx, unique.ratio, unique.quantile))
        i.name <- paste(unique.prefix, i, sep="")
        fData(ES)[,i.name]<-i.unique
    }
    return(ES)
}
unique.helper <- function(a, idx.i, idx.o, unique.ratio=2, unique.quantile=0.85) {
    a <- as.numeric(a)
    a_mu <- mean(a[idx.i])
    a_qn <- quantile(a[idx.o], unique.quantile)[[1]]
    a_unique <- 0

    if (a_mu/a_qn >= unique.ratio) {
        a_unique <- 1
    }

    return(a_unique)
}


#' Normalized and smoothened results of group specific differential expression test
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param log.base (numeric) the base of logarithm
#' @param test.statistic.metric.prefix (character) prefix for labeling the columns encoding the results of test statistic metric for each group
#' @param diff.expr.prefix (character) prefix for labeling the columns encoding the results of differential test
#' @return an ExpressionSet object containing the results of test statistic metrics in the attributes of fData
#'
exprs.test.statistic.metric <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, log.base=2, test.statistic.metric.prefix=TEST.STATS.METRIC.PREFIX, diff.expr.prefix=DIFF.EXPR.PREFIX) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    #cells <- rownames(pData(ES))
    for (i in groups) {
        i.name <- paste(test.statistic.metric.prefix, i, sep="")
        i.diffexpr.test.name <- paste(diff.expr.prefix, i, sep="")

        fData(ES)[,i.name] <- fData(ES)[, i.diffexpr.test.name]
        fData(ES)[,i.name] <- -log(fData(ES)[,i.name], log.base)
        i.min <- min(fData(ES)[,i.name])
        i.max <- max(fData(ES)[,i.name])
        fData(ES)[,i.name] <- fData(ES)[,i.name]/i.max
    }
    return(ES)
}


#' Calculating the similarity between each gene profile and group specific synthetic reference expression profile
#' synthetic profile is created by setting fpkm(cluster cells)=1, fpkm(non cluster cells)=0
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param syn.sim.prefix (character) prefix for labeling the columns encoding the results of synthetic profile similarity metric for each group
#' @return an ExpressionSet object containing the results of synthetic profile similarity metrics in the attributes of fData
#'
exprs.synthetic.similarity <- function(ES, group.by=CLUSTER.LABEL, groups=NULL, syn.sim.prefix=SYN.SIM.PREFIX) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    cells <- rownames(pData(ES))
    for (i in groups) {
        i.cells <- rownames(subset(pData(ES), pData(ES)[, group.by] %in% i))
        i.cells.o <- setdiff(cells, i.cells)
        i.cells.idx <- which(colnames(exprs(ES)) %in% i.cells)
        i.cells.o.idx <- which(colnames(exprs(ES)) %in% i.cells.o)
        i.synthetic.profile <-  rep(0, dim(exprs(ES))[2])
        i.synthetic.profile[i.cells.idx] <- 1
        i.ss <- apply(exprs(ES), 1, function(y) (1+cor(i.synthetic.profile, as.numeric(y)))/2)
        i.name <- paste(syn.sim.prefix, i, sep="")
        fData(ES)[,i.name]<-i.ss
    }
    return(ES)
}

#' Constructing group specific training sets for signature prediction
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param trainset.prefix (character) prefix for labeling the columns encoding the training instances for each group
#' @param train.class.prefix (character) prefix for labeling the columns encoding the class of group specific training instances
#' @param marker.prefix (character) prefix for labeling the columns encoding the biomarkers for each group
#' @param common.prefix (character) prefix for labeling the columns encoding the results of common metric for each group
#' @param unique.prefix (character) prefix for labeling the columns encoding the results of unique metric for each group
#' @param test.statistic.metric.prefix (character) prefix for labeling the columns encoding the results of test statistic metric for each group
#' @param syn.sim.prefix (character) prefix for labeling the columns encoding the results of synthetic profile similarity metric for each group
#' @return an ExpressionSet object containing the constructed training sets in the attributes of fData
#'
sigpred.training.sets <- function(ES, group.by=CLUSTER.LABEL, groups=NULL,
                                  trainset.prefix=TRAINSET.PREFIX,
                                  train.class.prefix=TRAIN.CLASS.PREFIX,
                                  test.statistic.metric.prefix = TEST.STATS.METRIC.PREFIX,
                                  common.prefix = COMMON.PREFIX,
                                  unique.prefix = UNIQUE.PREFIX,
                                  syn.sim.prefix = SYN.SIM.PREFIX,
                                  marker.prefix = MARKER.PREFIX
) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }

    for (i in groups) {
        i.train.name <- paste(trainset.prefix, i, sep="")
        fData(ES)[,i.train.name] <- 0
        i.train.class <- paste(train.class.prefix, i, sep="")
        fData(ES)[,i.train.class] <- NA

        i.p.candidates <- c()
        i.n.candidates <- c()

        i.test.name <- paste(test.statistic.metric.prefix, i, sep="")
        i.common.name <- paste(common.prefix, i, sep="")
        i.unique.name <- paste(unique.prefix, i, sep="")
        i.ss.name <- paste(syn.sim.prefix, i, sep="")

        fd <- fData(ES)
        fd <- fd[order(fd[,i.test.name]),]

        i.marker.name <- paste(marker.prefix, i, sep="")
        #i.p.symbols <- markers$SYMBOL[which(markers$CLUSTER %in% i)]
        i.p.candidates <- rownames(fd)[which(fd[, i.marker.name] == 1)]

        lambda_n =length(i.p.candidates)

        fd <- fData(ES)
        #fd <- fd[order(fd[,i.test.name], decreasing=T),]
        fd <- fd[order(fd[,i.test.name], decreasing=F),]
        fd <- subset(fd, fd[,i.common.name]==0 & fd[,i.unique.name]==0)
        i.n.candidates <- rownames(fd)[1:lambda_n]

        if (lambda_n > 0) {

            i.p.idx <- which(rownames(fData(ES)) %in% i.p.candidates)
            i.n.idx <- which(rownames(fData(ES)) %in% i.n.candidates)

            fData(ES)[c(i.p.idx, i.n.idx), i.train.name] <- 1
            fData(ES)[i.p.idx, i.train.class] <- 1
            fData(ES)[i.n.idx, i.train.class] <- 0
        }
    }
    return(ES)
}


#' Constructing group specific testing sets for signature prediction
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param testset.prefix (character) prefix for labeling the columns encoding the testing instances for each group
#' @param diff.expr.prefix (character) prefix for labeling the columns encoding the results of differential test
#' @param diff.expr.threshold (numeric) genes with differential expression p-value<0.05 will be selected as candidates for signature
#' @return an ExpressionSet object containing the constructed testing sets in the attributes of fData
#'
sigpred.testing.sets <- function(ES, group.by=CLUSTER.LABEL,  groups=NULL, testset.prefix=TESTSET.PREFIX, diff.expr.prefix=DIFF.EXPR.PREFIX, diff.expr.threshold=0.05) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    for (i in groups) {
        i.name <- paste(testset.prefix, i, sep="")
        fData(ES)[,i.name] <- 0
        i.test.name <- paste(diff.expr.prefix, i, sep="")
        i.test.idx <- which(fData(ES)[,i.test.name] < diff.expr.threshold)
        if (length(i.test.idx) > 0) {
            fData(ES)[i.test.idx, i.name] <- 1
        }
    }
    return(ES)
}



#' logistic regression based signature prediction and cross cluster validation
#' Cell type specific signature prediction
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param groups (character) the clusters for signature prediction
#' @param trainset.prefix (character) prefix for labeling the columns encoding the training instances for each group
#' @param train.class.prefix (character) prefix for labeling the columns encoding the class of group specific training instances
#' @param testset.prefix (character) prefix for labeling the columns encoding the testing instances for each group
#' @param common.prefix (character) prefix for labeling the columns encoding the results of common metric for each group
#' @param unique.prefix (character) prefix for labeling the columns encoding the results of unique metric for each group
#' @param test.statistic.metric.prefix (character) prefix for labeling the columns encoding the results of test statistic metric for each group
#' @param syn.sim.prefix (character) prefix for labeling the columns encoding the results of synthetic profile similarity metric for each group
#' @param signature.prefix (character) prefix for labeling the columns encoding the results of signature prediction
#' @param gene.symbol.label (character) the name of the column encoding gene symbols
#' @param verbose (logical)
#' @param export (logical) wheterh to export the ExpressionSet
#' @param export.components (character) the components of an ExpressionSet object that will be exported
#' @return an ExpressionSet object containing the results of signature prediction in the attributes of fData
#'
signature.prediction.old <- function(ES, group.by="CLUSTER", groups=NULL,
                                     trainset.prefix=TRAINSET.PREFIX,
                                     train.class.prefix=TRAIN.CLASS.PREFIX,
                                     testset.prefix=TESTSET.PREFIX,
                                     test.statistic.metric.prefix = TEST.STATS.METRIC.PREFIX,
                                     common.prefix = COMMON.PREFIX,
                                     unique.prefix = UNIQUE.PREFIX,
                                     syn.sim.prefix = SYN.SIM.PREFIX,
                                     signature.prefix = SIGNATURE.PREFIX,
                                     gene.symbol.label = GENE.SYMBOL.LABEL,
                                     verbose=TRUE,
                                     export=TRUE,
                                     dir.prefix = ""
) {
    if (is.null(groups)) {
        groups <- sort(unique(pData(ES)[,group.by]))
    }
    predictions <- list()
    for (i in groups) {

        i.train.name <- paste(trainset.prefix, i, sep="")

        i.test.name <- paste(test.statistic.metric.prefix, i, sep="")
        i.common.name <- paste(common.prefix, i, sep="")
        i.unique.name <- paste(unique.prefix, i, sep="")
        i.ss.name <- paste(syn.sim.prefix, i, sep="")
        i.train.class <- paste(train.class.prefix, i, sep="")

        # train
        i.train.data <- subset(fData(ES)[, c(i.common.name, i.unique.name, i.test.name, i.ss.name, i.train.class)], fData(ES)[, i.train.name] == 1)

        if (verbose == T) {
            cat("\tConstructing model for cluster", i,":")
        }

        formula_for_train <- paste(i.train.class, "~", i.test.name, "+", i.ss.name, sep="")

        i.train.p <- i.train.data[which(i.train.data[,i.train.class]==1),]
        i.train.p.sd <- apply(as.matrix(i.train.p), 2, sd)

        # removing dominant factors

        factors <- c(i.common.name, i.unique.name, i.test.name, i.ss.name)

        if (i.train.p.sd[i.common.name] == 0) {
            if (verbose) {
                cat("ignore common gene metric ")
            }
        } else {
            formula_for_train <- paste(formula_for_train, "+", i.common.name, sep="")
        }

        if (i.train.p.sd[i.unique.name] == 0) {
            if(verbose) {
                cat(" ignore unique gene metric ")
            }
        } else {
            formula_for_train <- paste(formula_for_train, "+", i.unique.name, sep="")
        }

        i.glm <- glm(as.formula(formula_for_train),
                     data=i.train.data,
                     family = "binomial")
        if (verbose) {
            cat(" done\n")
            cat("\tpredicting signature genes for group", i, " ")
        }


        # prediction
        i.testing.name <- paste(testset.prefix, i, sep="")
        i.test.data <- subset(fData(ES)[, c(gene.symbol.label, i.common.name, i.unique.name, i.test.name, i.ss.name)], fData(ES)[,i.testing.name]==1)
        i.prediction <- predict(i.glm, i.test.data)
        i.pred.name <- paste(signature.prefix, i, sep="")
        #predictions[i.pred.name] <- data.frame(GENE=i.test.data$SYMBOL, P=as.numeric(i.prediction))

        i.prediction <- as.numeric(i.prediction)
        i.prediction <- (i.prediction-min(i.prediction))/(max(i.prediction)-min(i.prediction))

        fData(ES)[,i.pred.name] <- NA
        fData(ES)[rownames(i.test.data), i.pred.name] <- i.prediction
        if (export) {
            write.table(data.frame(GENE=rownames(i.test.data), SYMBOL=i.test.data[,gene.symbol.label], PREDICTION=as.numeric(i.prediction)), file=paste(dir.prefix , "signature_prediction_", i, ".txt", sep=""), sep="\t", col.names=T, row.names=F) #%%
        }
        if (verbose) {
            cat("done\n")
        }
    }

    #print(glms)
    return(ES)
}



#######################################################
#             Driving Force Analysis                  #
#######################################################











#######################################################
#                AUXILIARY FUNCTIONS                  #
#######################################################



# Exporting specified Expression Set components
exportES <- function(ES, prefix="", suffix="", components=c("fd", "pd", "exprs")) {
    if ("fd" %in% components) {
        fd.df <- fData(ES)
        fd.df <- cbind(RID=rownames(fd.df), fd.df)
        write.table(fd.df, file=paste(prefix, "fd", suffix, sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
    }
    if("pd" %in% components) {
        pd.df <- pData(ES)
        pd.df <- cbind(RID=rownames(pd.df), pd.df)
        write.table(pd.df, file=paste(prefix, "pd", suffix, sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
    }
    if("exprs" %in% components) {
        exprs.m <- as.data.frame(exprs(ES))
        exprs.m <- cbind(RID=rownames(exprs.m), exprs.m)
        write.table(exprs.m, file=paste(prefix, "exprs", suffix, sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
    }
}




multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # Make the panel
    plotCols = cols                          # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }

}


#######################################################
#                Plotting FUNCTIONS                  #
#######################################################

plotProfiles <- function(ES, genes, color.by=CLUSTER.LABEL,
                         fig.filename=NULL, fig.res=150, fig.w=5, fig.h=2,
                         xlabel="E16.5 Lung Single Cells", ylabel="Expression", point.size=3, font.size=12) {


    genes <- genes[which(genes %in% fData(ES)[,GENE.SYMBOL.LABEL])]

    if (length(genes) <=0) {
        stop("no matched genes found.")
    }

    n.genes <- length(genes)

    n.cells <- dim(pData(ES))[1]

    if (!is.null(fig.filename)) {
        # bug: unable to start tiff() device (probably due to the image size)
        tiff(file=fig.filename, width=fig.w, height=fig.h*n.genes, unit="in", res=fig.res, compression="lzw")
    }

    gs <- list()

    for (i in 1:n.genes) {
        gene.idx <- which(fData(ES)[,GENE.SYMBOL.LABEL] %in% genes[i])
        i.expr <- NULL
        i.expr <- data.frame(CID=1:dim(pData(ES))[1], CELL=rownames(pData(ES)), EXPRESSION=as.numeric(exprs(ES)[gene.idx,]))
        if (!is.null(color.by)) {
            i.expr[, color.by] <- as.factor(pData(ES)[,color.by])
        }
        i.g <- ggplot(i.expr, aes(x=CID, y=EXPRESSION))
        if (is.null(color.by)) {
            i.g <- i.g + geom_point(size=point.size)
        } else {
            i.g <- i.g + geom_point(aes_string(colour=color.by), size=point.size)
            i.g <- i.g + guides(colour=guide_legend(ncol=ceiling(length(unique(i.expr[,color.by]))/5)))
        }

        i.g <- i.g + labs(xlab=xlabel, title=as.character(fData(ES)[gene.idx, GENE.SYMBOL.LABEL]))
        i.g <- i.g + ylab(ylabel)
        i.g <- i.g + sincera_theme()
        i.g <- i.g + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
        i.g <- i.g + theme(text= element_text (size=font.size ))
        i.g <- i.g + theme(axis.title.y=element_text(size=font.size), axis.text.y=element_text(size=font.size))
        i.g <- i.g + theme(plot.title = element_text(size = rel(1.5), face="bold.italic"))

        gs[[i]] <- i.g
    }

    multiplot(plotlist=gs, cols=1)

    if (!is.null(fig.filename)) {
        dev.off()
    }

}
