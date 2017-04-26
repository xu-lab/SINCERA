#' Predicting cell type signature
#'
#' By default, the function predicts cell type signature by using percentage filter and fold-induction filter to refine cluster specific differentially expressed genes
#' When the marker information for some cell types is available, one can choose to predict signature using a logistic regression model integrating four metrics
#'
#' @param object A sincera object
#' @param groups The cell groups for signature prediction
#' @param use.logireg If TRUE, use logistic regression model based method to predict signature for cell groups (cell types) with marker information available;
#'                    use setCellTypeMarkers to add marker information to sincera, and use getCellTypeMarkers to check added cell type information
#' @param sigs.n A vector containing the number of signatures to be predicted for each cell group; if NULL, use genes with top 20 prediction score as the signature
#' @param diff.method Take effect when use.logireg==FALSE. The method of differentiatial expression
#' @param use.fdr Take effect when use.logireg==FALSE. If TRUE, use BH to adjust pvalues of differential expression
#' @param diff.thresh Take effect when use.logireg==FALSE. The threshold of pvalue or fdr.
#' @param avg.thresh A cell type signature gene must have at least avg.thresh expression in the defined cell type
#' @param min.expression The threshold for determining expressed genes
#' @param pct.thresh A cell type signature gene must express in at least pct.thresh percentage of the cells of the defined cell type
#' @param fc.thresh A cell type signature gene must have at least fc.thresh fold of average expression induction in its defined cell type cells when compared to its average expression in all the other cells
#' @return An updated sincera object, use getSigGenes to access the predicted signature genes'
#'
setGeneric("signature.prediction", function(object, groups=NULL, use.logireg=FALSE, sigs.n=NULL, diff.method="welch", use.fdr=FALSE, diff.thresh=0.05,
                                            avg.thresh=1, min.expression=1, pct.thresh=0.8, fc.thresh=1.5, ...) standardGeneric("signature.prediction"))
#' @export
setMethod("signature.prediction","sincera",
          function(object, groups=NULL, use.logireg=FALSE, sigs.n=NULL, diff.method="welch", use.fdr=FALSE, diff.thresh=0.05,
                   avg.thresh=1, min.expression=1, pct.thresh=0.8, fc.thresh=1.5, ...) {

            if (is.null(groups)) {
              groups <- sort(unique(getCellMeta(object, name="GROUP")))
            }

            clusters <- groups

            if (!use.logireg) { # use per cluster statistics to refine diff genes for signature prediction

                cat("Use percentage and fold_induction filters to refine differentially expressed genes for signature prediction")

                gs <- GeneStats(dp=getExpression(object), ident=getCellMeta(object, name="GROUP"),
                                groups=clusters, min.expression=min.expression, min.gsize=2, min.avg=1,
                                stats=c(diff.method, "ident.avg", "ident.pct","ident.fc"))

                if (use.fdr) {
                    sigs <- data.frame(GROUP=NULL, SYMBOL=NULL, P=NULL, FDR=NULL, AVG=NULL, PCT=NULL, FC=NULL)
                } else {
                    sigs <- data.frame(GROUP=NULL, SYMBOL=NULL, P=NULL, AVG=NULL, PCT=NULL, FC=NULL)
                }


                for (g in clusters) {
                    diff.col <- paste(diff.method, ".", g, sep="")
                    diff.fdr.col <- paste(diff.method, ".", g, ".fdr", sep="")
                    fc.col <- paste("ident.fc.", g, sep="")
                    avg.col <- paste("ident.avg.", g, sep="")
                    pct.col <- paste("ident.pct.", g, sep="")
                    sig.col <- paste("sig.", g, sep="")

                    g.fc.idx <- which(gs[, fc.col] >= fc.thresh)
                    g.avg.idx <- which(gs[, avg.col] >= avg.thresh)
                    g.pct.idx <- which(gs[, pct.col] >= pct.thresh)

                    g.sig.idx <- c()
                    g.sig.idx <- intersect(g.fc.idx, g.avg.idx)
                    g.sig.idx <- intersect(g.sig.idx, g.pct.idx)

                    if (use.fdr) {
                        gs[, diff.fdr.col] <- NA
                        gs[g.sig.idx, diff.fdr.col] <- p.adjust(as.numeric(gs[g.sig.idx, diff.col]), method="BH")
                        g.diff.idx <- which(gs[, diff.fdr.col]<diff.thresh)
                    } else {
                        g.diff.idx <- which(gs[, diff.col] < diff.thresh)
                    }

                    g.sig.idx <- intersect(g.sig.idx, g.diff.idx)

                    gs[, sig.col] <- 0
                    if (length(g.sig.idx) > 0) {
                        gs[g.sig.idx, sig.col] <- 1
                        if (use.fdr) {
                            g.sig <- data.frame(GROUP=g, SYMBOL=rownames(gs)[g.sig.idx],
                                                P=gs[g.sig.idx, diff.col],
                                                FDR=gs[g.sig.idx, diff.fdr.col],
                                                AVG=gs[g.sig.idx, avg.col],
                                                PCT=gs[g.sig.idx, pct.col],
                                                FC=gs[g.sig.idx, fc.col])
                        } else {
                            g.sig <- data.frame(GROUP=g, SYMBOL=rownames(gs)[g.sig.idx],
                                                P=gs[g.sig.idx, diff.col],
                                                AVG=gs[g.sig.idx, avg.col],
                                                PCT=gs[g.sig.idx, pct.col],
                                                FC=gs[g.sig.idx, fc.col])
                        }
                        sigs <- rbind(sigs, g.sig)
                    }

                }

                sigs <- sigs[order(sigs$GROUP, sigs$P, -sigs$FC), ]

                object <- setSigGenes(object, value=sigs)


            } else { # use logistic regression model to predict signature for clusters with markers

                es <- getES(object)
                es <- es[getGenesForClustering(object), ]

                markers <- getCellTypeMarkers(object)

                types <- c()
                if (dim(markers)[1]>0) {
                  types <- sort(unique(markers$TYPE))
                }

                celltypes <- getCellType(object)

                clusters.withMarkers <- c()

                for (i in 1:length(types)) {
                  i.type <- types[i]
                  i.cluster <- names(celltypes)[which(celltypes == i.type)]
                  i.name <- paste("use_as_marker_", i.cluster, sep="")
                  i.markers <- markers$SYMBOL[which(markers$TYPE==i.type)]
                  fData(es)[, i.name] <- 0
                  fData(es)[which(rownames(fData(es)) %in% i.markers), i.name] <- 1
                  clusters.withMarkers <- c(clusters.withMarkers, i.cluster)
                }
                clusters.withMarkers <- sort(clusters.withMarkers)

                if (length(clusters.withMarkers)>0) {

                    cat("Use logistic regression and markers to predict signature for cell group(s):", paste(clusters.withMarkers, split=" "))

                    difftest <- getDiffTest(object)
                    es <- es[rownames(difftest), ]
                    fData(es)[, colnames(difftest)] <- difftest

                    wd <- NULL

                    group.by = "GROUP"
                    groups=clusters.withMarkers
                    # training set
                    trainset.prefix = TRAINSET.PREFIX
                    train.class.prefix = TRAIN.CLASS.PREFIX
                    marker.prefix=MARKER.PREFIX
                    # testing set
                    testset.prefix=TESTSET.PREFIX
                    # common gene metric
                    common.prefix = COMMON.PREFIX
                    common.threshold=COMMON.TRESHOLD
                    common.percentage=COMMON.PERCENTAGE
                    # unique gene metric
                    unique.prefix = UNIQUE.PREFIX
                    unique.ratio=UNIQUE.RATIO
                    unique.quantile=UNIQUE.QUANTILE
                    # test statistic metric
                    test.statistic.metric.prefix = TEST.STATS.METRIC.PREFIX
                    log.base=2
                    diff.expr.prefix="welch."
                    diff.expr.threshold = 1
                    # synthetic profile similarity metric
                    syn.sim.prefix = SYN.SIM.PREFIX
                    signature.prefix = SIGNATURE.PREFIX
                    gene.symbol.label = GENE.SYMBOL.LABEL
                    verbose=TRUE
                    export=FALSE
                    export.components="fd"
                    dir.prefix=NULL

                    # calculate metrics
                    es <- exprs.test.statistic.metric(es, group.by=group.by, groups=groups, log.base=log.base, test.statistic.metric.prefix=test.statistic.metric.prefix, diff.expr.prefix=diff.expr.prefix)

                    es <- exprs.common(es, group.by=group.by, groups=groups, common.threshold=common.threshold, common.percentage=common.percentage, common.prefix=common.prefix)

                    es <- exprs.unique(es, group.by=group.by, groups=groups, unique.ratio=unique.ratio, unique.quantile=unique.quantile, unique.prefix=unique.prefix)

                    es <- exprs.synthetic.similarity(es, group.by=group.by, groups=groups, syn.sim.prefix=syn.sim.prefix)

                    # construct training sets
                    es <- sigpred.training.sets(es, group.by=group.by, groups=groups,
                                                trainset.prefix=trainset.prefix,
                                                train.class.prefix=train.class.prefix,
                                                test.statistic.metric.prefix = test.statistic.metric.prefix,
                                                common.prefix = common.prefix,
                                                unique.prefix = unique.prefix,
                                                syn.sim.prefix = syn.sim.prefix,
                                                marker.prefix = marker.prefix)

                    # construct testing sets
                    es <- sigpred.testing.sets(es, group.by=group.by, groups=groups, testset.prefix=testset.prefix, diff.expr.prefix=diff.expr.prefix, diff.expr.threshold=diff.expr.threshold)

                    # performing signature prediction
                    es <- signature.prediction.old(es, group.by=group.by, groups=groups,
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
                                               dir.prefix = wd)

                    sig.cols <- paste("sig_pred_", clusters.withMarkers, sep="")
                    if (is.null(sigs.n)) sigs.n <- rep(20, length(clusters.withMarkers))
                    sigs <- data.frame(GROUP=NULL, SYMBOL=NULL, SCORE=NULL, stringsAsFactors = FALSE)
                    for (i in 1:length(sig.cols)) {
                        tmp <- as.numeric(fData(es)[, sig.cols[i]])
                        if (sigs.n[i] > length(tmp)) sigs.n[i] <- length(tmp)
                        i.idx <- order(tmp, decreasing=T)[1:sigs.n[i]]
                        i.idx <- i.idx[which(!is.na(tmp[i.idx]))]
                        if (length(i.idx)>0) {
                            i.sigs <- data.frame(GROUP=clusters.withMarkers[i], SYMBOL=rownames(fData(es))[i.idx], SCORE=as.numeric(fData(es)[i.idx, sig.cols[i]]), stringsAsFactors = FALSE)
                            sigs <- rbind(sigs, i.sigs)
                        }
                    }

                    if (dim(sigs)[1]>0) sigs <- sigs[order(sigs$GROUP, -sigs$SCORE), ]

                    object <- setSigGenes(object, value=sigs)

                }
            }

            return(object)
          }
)



#' Validating cell type signature using random subsampling mehod
#'
#' @param object A sincera object
#' @param siggenes A data frame containing the . The first column of the data frame contains cell cluster information, while the second column contains the genes
#' @param pct The percentage of cells to be used to construct testing set
#' @param repeats The number of subsamplings
#' @return An updated sincera object, use getSigValidation to access the validation results
#'
setGeneric("signature.validation", function(object, siggenes=NULL, pct=0.2, repeats=3, ...) standardGeneric("signature.validation"))
#' @export
setMethod("signature.validation","sincera",
          function(object, siggenes=NULL, pct=0.2, repeats=3, ...) {

            if (is.null(siggenes)) {

            } else {

              cellcluster <- getCellMeta(object, name="GROUP")
              all.groups <- names(table(cellcluster))
              cluster.cells.n <- as.numeric(table(cellcluster))
              names(cluster.cells.n) <- all.groups

              sig.groups <- sort(as.character(unique(siggenes$GROUP)))
              sig.groups <- intersect(sig.groups, all.groups)

              if (length(sig.groups)==0) {
                stop("Please check cluster definition in the siggenes.")
              }

              cat("Sincera: performing repeated random subsampling validation of signature prediction for groups:", paste(sig.groups, split=" "), "\n")

              nsiggrp <- length(sig.groups)
              nallgrp <- length(all.groups)

              cluster.cells.idx <- list()
              j <- 1
              for (i in all.groups) {
                  idx.i <- which(cellcluster==i)
                  cluster.cells.idx[[j]] <- idx.i
                  j <- j+1
              }
              names(cluster.cells.idx) <- all.groups



              cells.idx <- 1:getCellNum(object)

              tr = NULL
              tr <- data.frame(TID=1:(nsiggrp*nallgrp*repeats), I=NA, J=NA, P=NA, IN=NA, JN=NA, R=NA, IC=NA, JC=NA, TA=NA, TBA=NA, TF1=NA, PA=NA, PBA=NA, PF1=NA)
              z <- 1

              ES <- getES(object)

              for (i in 1:nsiggrp) { # evaluate the signature of cluster i


                  i.sigs <- as.character(siggenes[which(siggenes$GROUP == sig.groups[i]), "SYMBOL"])

                  if (length(i.sigs)>0) {

                    i.sigs.idx <- which(rownames(fData(ES)) %in% i.sigs)

                    for (j in 1:nallgrp) { # using cluster j

                        p = pct
                        #for (p in pct) { # percentage of cells from both clusters for testing
                            for (r in 1:repeats) {

                                if (sig.groups[i] != all.groups[j]) {

                                    i.n <- max(floor(cluster.cells.n[sig.groups[i]] * p),1)
                                    j.n <- max(floor(cluster.cells.n[all.groups[j]] * p),1)
                                    i.cells.test.idx <- sample(cluster.cells.idx[[sig.groups[i]]], i.n, replace=FALSE)
                                    j.cells.test.idx <- sample(cluster.cells.idx[[all.groups[j]]], j.n, replace=FALSE)


                                    i.trainset <- data.frame(t(exprs(ES)[i.sigs.idx, setdiff(cluster.cells.idx[[sig.groups[i]]], i.cells.test.idx)]), CLUSTER=paste("C", sig.groups[i], sep=""))
                                    j.trainset <- data.frame(t(exprs(ES)[i.sigs.idx, setdiff(cluster.cells.idx[[all.groups[j]]], j.cells.test.idx)]), CLUSTER=paste("C", all.groups[j], sep=""))

                                    i.testset <- data.frame(t(exprs(ES)[i.sigs.idx, i.cells.test.idx]), CLUSTER=paste("C", sig.groups[i], sep=""))
                                    j.testset <- data.frame(t(exprs(ES)[i.sigs.idx, j.cells.test.idx]), CLUSTER=paste("C", all.groups[j], sep=""))

                                    trainset = NULL
                                    testset = NULL
                                    trainset <- rbind(i.trainset, j.trainset)
                                    testset <- rbind(i.testset, j.testset)

                                    category.idx <- dim(trainset)[2]

                                    library(e1071)

                                    mm = NULL
                                    mm <- svm(CLUSTER ~ ., data=trainset)

                                    prediction = NULL
                                    prediction <- predict(mm, trainset[,-category.idx])

                                    tp = NA
                                    fn = NA
                                    tn = NA
                                    fp = NA
                                    ta = NA
                                    tba = NA
                                    tf1 = NA

                                    # training accuracy
                                    tab = NULL
                                    tab <- table(pred = prediction, true = trainset[,category.idx])
                                    tp <- tab[1,1]
                                    fn <- tab[2,1]
                                    fp <- tab[1,2]
                                    tn <- tab[2,2]
                                    ta <- (tp + tn)/(tp+fp+fn+tn)
                                    tba <- (tp/(tp+fn) + tn/(tn+fp))/2
                                    tf1 <- 2*tp / (2*tp + fp + fn)

                                    tp = NA
                                    fn = NA
                                    tn = NA
                                    fp = NA
                                    pa = NA
                                    pba = NA
                                    pf1 = NA

                                    # prediction accuracy
                                    prediction = NULL
                                    prediction <- predict(mm, testset[, -category.idx])

                                    tab = NULL
                                    tab <- table(pred = prediction, true = testset[,category.idx])
                                    tp <- tab[1,1]
                                    fn <- tab[2,1]
                                    fp <- tab[1,2]
                                    tn <- tab[2,2]
                                    pa <- (tp + tn)/(tp+fp+fn+tn)
                                    pba <- (tp/(tp+fn) + tn/(tn+fp))/2
                                    pf1 <- 2*tp / (2*tp + fp + fn)

                                    tr$I[z] <- sig.groups[i]
                                    tr$J[z] <- all.groups[j]
                                    tr$P[z] <- p
                                    tr$IN[z] <- i.n
                                    tr$JN[z] <- j.n
                                    tr$R[z] <- r
                                    tr$IC[z] <- paste(as.character(rownames(pData(ES))[i.cells.test.idx]), collapse=", ")
                                    tr$JC[z] <- paste(as.character(rownames(pData(ES))[j.cells.test.idx]), collapse=", ")
                                    tr$TA[z] <- ta
                                    tr$TBA[z] <- tba
                                    tr$TF1[z] <- tf1
                                    tr$PA[z] <- pa
                                    tr$PBA[z] <- pba
                                    tr$PF1[z] <- pf1


                                    z <- z+1

                                }
                           # }
                        }
                    }
                  }

                  # print(i)
              }

              tr <- tr[which(!is.na(tr$I)),]


              #write.table(tr, file=paste(wd, "cross-validation-log.txt", sep=""), sep="\t", col.names=T, row.names=F)

              tr.mean <- aggregate(tr[,"PA"], list(tr$I, tr$J, tr$P), mean)
              tr.sd <- aggregate(tr[,"PA"], list(tr$I, tr$J, tr$P), sd)

              colnames(tr.mean) <- c("I","J","P", "PA.mean")
              colnames(tr.sd) <- c("I","J","P", "PA.sd")

              tr.mean <- tr.mean[order(tr.mean$I),]
              tr.sd <- tr.sd[order(tr.sd$I),]

              tr.summary <- cbind(tr.mean, PA.sd=tr.sd$PA.sd, PA.se=tr.sd$PA.sd/sqrt(repeats))

              #
              tr.summary.2 <-tr.summary
              colnames(tr.summary.2) <- c("cluster.i","cluster.j","sampling.p","accuracy.mean","accuracy.sd", "accuracy.se")
              #write.table(tr.summary.2, file=paste(wd, "cross-validation-summary.txt", sep=""), sep="\t", col.names=T, row.names=F)

              np <- nsiggrp
              gs <- list()

              for (i in 1:np) {

                  tr.i <- tr.summary[which(tr.summary$I == sig.groups[i]),]
                  tr.i <- tr.i[which(!(tr.i$I==tr.i$J)),]
                  tr.i$Group <- paste(tr.i$J, sep="")

                  g <- ggplot(tr.i, aes(x=Group, y=PA.mean)) +
                      geom_bar(position=position_dodge(), stat="identity", width=.6) +
                      geom_errorbar(aes(ymin=PA.mean-PA.se, ymax=PA.mean+PA.se), width=.2) +
                      ylab("Accuracy") +
                      ggtitle(paste("Classificiation accurary of predicted signature of cell group ", sig.groups[i], " (", (1-pct)*100, "% training; ", repeats, " repeats)", sep="")) +
                      sincera_theme()
                  gs[[i]] <- g

                  print(g)
                  pause()

              }

              #tiff(file=paste(wd, "accuracy.tiff",sep=""), width=4, height=2*np, unit="in", res=150, pointsize=2, compression="lzw")
              #multiplot(plotlist=gs, cols=1)
              #dev.off()



              cat("Sincera: done\n")



              return(object)
            }
        }

)
