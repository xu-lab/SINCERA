

#' Selecting a threshold value for specificity
#'
#' @param object A sincera object
#' @param genes The set of housekeeping genes used to determine the specificity threshold
#' @export
specificity.thresholdSelection <- function(object, genes) {

    ES <- getES(object)
    ref.idx <- which(rownames(fData(ES)) %in% genes)
    specificity.threshold=specificity.criterion.selection(ES, group.by="SAMPLE", groups=NULL, ref.idx=ref.idx, exp.t=5, exp.p = 0.95, tolerance=0.05, step=0.01)

    return(specificity.threshold)
}

#' Selecting genes for cell cluster identification
#'
#' Currently support the specificity metric (Guo et al., PLoS Comp Bio 2015); other metrics will be supported soon
#'
#' @param object A sincera object
#' @param method The selection method, possible values include specificity
#' @param pergroup If TRUE, calculate metrics for each cell groups (samples)
#' @param min.sample The selected genes must pass the selection criteria in at least min.sample samples
#' @param specifity.thresh The specificity threshold
#' @param do.plot If TRUE, plot figures of selection
#' @return The updated sincera object with the selected genes in the genes.forclustering slot
#'
setGeneric("cluster.geneSelection", function(object, method="specificity",
                                             pergroup=TRUE, min.samples=2,
                                             specifity.thresh=0.7,
                                             do.plot=T,
                                             ...) standardGeneric("cluster.geneSelection"))
#' @export
setMethod("cluster.geneSelection","sincera",
          function(object, method="specificity", 
                   pergroup=TRUE, min.samples=2,
                   specifity.thresh=0.7,
                   do.plot=T,
                   ...) {

            genes.selected <- as.character(getGenes(object))
            n <- length(genes.selected)
            cellsample <- getCellMeta(object, name="GROUP")

            if (TRUE==pergroup) {
              samples <- sort(unique(as.character(cellsample)))
            } else {
              samples <- c("sample1")
              cellsample <- rep("sample1", getCellNum(object))
              min.samples <- 1
            }

            if (method=="specificity") {

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


              cat("\nUse expression specificity to select genes for cluster identification\n")

              stats <- data.frame(SYMBOL=genes.selected)
              rownames(stats) <- stats$SYMBOL
              for (s in samples) {
                s.idx <- which(cellsample == s)
                if (length(s.idx)==1) {
                  stats[, s] <- rep(1, dim(stats)[1]) # singleton cluster has no effect on selection
                } else {
                  stats[, s] <- apply(exprs(object@data)[, s.idx], 1, function(x) specificity.helper(x))
                }
              }


              if (length(samples)>1) {
                nsamples <- rowSums(stats[, -1] >= specifity.thresh)
              } else {
                nsamples <- rep(0, length(genes.selected))
                nsamples[which(as.numeric(stats[, 2]) >= specifity.thresh)] <- 1
              }

              genes.selected <- which(nsamples >= min.samples)
              n <- length(genes.selected)
              if (n<2) {
                stop("Too few number of genes passed the selection criteria.")
              }

              if (TRUE==do.plot) {
                cat("\nPlot distributions of expression specificity\n")
                stats.melt <- melt(stats, i.vars="SYMBOL")
                colnames(stats.melt) <- c("SYMBOL", "Group", "Specificity")
                g <- ggplot(stats.melt, aes(x=Specificity))
                g <- g + geom_histogram(col="grey")
                g <- g + facet_wrap(~Group)
                g <- g + geom_vline(xintercept=specifity.thresh, col="red")
                g <- g + ggtitle(label="Distribution of expression specificity")
                g <- g + sincera_theme()
                print(g)
              }

            } else if (method=="cv-mean") {
              #TODO: check CV-avg calculation
            }

            object <- setGenesForClustering(object, value=getGenes(object)[genes.selected])

            cat(n, " genes selected for identifying cell clusters\n", sep="")
            cat("Use getGenesForClustering() to view the selected genes\n\n")


            return(object)
          }
)




#' Partition cells into disjoint clusters
#'
#' @param object A sincera object
#' @param feature.type The feature space used for clustering: "gene" - means gene expression space, "pca" - means reduced dimension space using PCA
#' @param rds.dims A numberic vector specified the reduced dimensions used for clustering. Only take effect when the feature.type is "pca" or "tsne"
#' @param update.cellgroup The clustering result will be saved to the "CLUSTER" meta data. If TRUE, the GROUP meta data will be updated as well.
#' @param clustering.method The cluster identification algorithm, possible values include pam, kmeans, hc, graph, tight, consensus
#' @param verbose If TRUE, print the verbose messages
#' @return The update sincera object with clustering results in the "CLUSTER" meta data, use getCellMeata with name="CLUSTER" to assess the result
#' @details
#' The default clustering method is hc - hierarchical clustering with Pearson's correlation based distance and average linkage. 
#' Possible parameters include:
#'   h - if not NULL, cut dendrogram tree at the height h;
#'   k - if not NULL, cut dendrogram tree to generate k clusters; if both k and h are not NULL, k will be used; If both k and h are NULL, the algorithm will find the largest k that contains no more than num.singleton singleton clusters; 
#'   num.singleton - the number of singleton clusters allowed;
#'   distance.method - the method to calculate cell distance, possible values include pearson, spearman, euclidean;
#'   linkage.method - linkage method, possible values include average, complete, ward.D2;
#'   do.shift - if TRUE, shit column mean to 0.
#'
#'      
#' When clustering.method is graph, the function utilizes the graph based method described in Shekhar et al., 2016, which first constructs a kNN graph of cells and then applies a community detection algorithm to partition cells.
#' Possible parameters include:
#'   num.nn - The number of nearest neighbors considered during the kNN graph construction; 
#'   do.jaccard - If TRUE, weigh kNN graph edges using Jaccard similarity;
#'   community.method - The community detection method, possible values include louvain, infomap, walktrap, spinglass, edge_betweenness, label_prop, optimal, fast_greedy.
#'    
#'    
#' When clustering.method is tight, the function utilizes tightClust::tight.clust() to perform tight clustering; please refer to ?tightClust::tight.clust for parameters and more informaiton
#' 
#' 
#' When clustering.method is consensus, the function utilizes ConsensusClusterPlus::ConsensusClusterPlus() to perform consensus clustering; please refer to ?ConsensusClusterPlus::ConsensusClusterPlus for parameters and more informaiton
#' 
#' 
#' When clustering.method is kmeans, the function utilizes stats::kmeans() to perform k-means clustering; please refer to ?kmeans for parameters and more informaiton
#' 
#' 
#' When clustering.method is pam, the function utilizes cluster::pam() to perform Partitioning Around Medoids clustering; please refer to ?cluster::pam for parameters and more information
#' 
setGeneric("cluster.assignment", function(object, feature.type="gene", rds.dims=1:3,
                                          update.cellgroup=T,
                                          clustering.method="hc",
                                          verbose=T, ...) standardGeneric("cluster.assignment"))
#' @export
setMethod("cluster.assignment","sincera",
          function(object, feature.type="gene", rds.dims=1:3,
                   update.cellgroup=T,
                   clustering.method="hc",
                   verbose=T, ...) {

            ft.valid <- c("gene","pca", "tsne")
            if (!(feature.type %in% ft.valid)) {
              stop("Please select a valid feature type: ", paste(ft.valid, collapse=", ", sep=""))
            }
            if (feature.type=="gene") {
              if (length(object@genes.forclustering)<2) {
                stop("The feature.type is set to \'gene\'. But too few genes were selected for clustering.\nPlease use cluster.geneSelection() to select more genes or use setGenesForClustering() to set the \'genes.forclustering\' slot")
              }
            } else if (feature.type=="pca") {
              if (dim(object@pca$rds)[1]==0) {
                stop("The feature.type is set to \'pca\'. But no PCA dimensions were found. \nPlease run doPCA() or use setPCAScores() to set the \'pca.scores\' slot\n")
              } else {
                rds.dims.notfound <- rds.dims[which(!(paste("PC", rds.dims, sep="") %in% colnames(object@pca$rds)))]
                if (length(rds.dims.notfound)>0) {
                  stop("The following PCA dimensions are not found: ", paste(rds.dims.notfound, collapse = ",", sep=""))
                }
              }
            } else if (feature.type=="tsne") {
              if (dim(object@tsne$rds)[1]==0) {
                stop("The feature.type is set to \'tsne\'. But no tSNE dimensions were found. \nPlease run doTSNE() or use setTSNE() to set the rds element of tsne slot\n")
              } else {
                rds.dims.notfound <- rds.dims[which(!(paste("tSNE", rds.dims, sep="") %in% colnames(object@tsne$rds)))]
                if (length(rds.dims.notfound)>0) {
                  stop("The following tSNE dimensions are not found: ", paste(rds.dims.notfound, collapse = ",", sep=""))
                }
              }
            }


            cmethods <- c("tight", "hc","consensus","ihc","graph","pam", "kmeans")
            cmethod.idx <- pmatch(clustering.method, cmethods)
            if (is.na(cmethod.idx)) {
              stop("invalid clustering method")
            }
            clustering.method <- cmethods[cmethod.idx]


            if (clustering.method=="tight") {
              if (!require(tightClust)) {
                stop("The package 'tightClust' is required for tight clustering.")
              }
            } else if (clustering.method=="consensus") {
              if (!require(ConsensusClusterPlus)) {
                stop("The package 'ConsensusClusterPlus' is required for consensus clustering.")
              }
            } else if (clustering.method=="pam") {
              if (!require(cluster)) {
                stop("The package 'cluster' is required for pam clustering.")
              }
            } else if (clustering.method=="graph") {
              if (!require(RANN)) {
                stop("The package 'RANN' is required for graph based clustering.")
              }
            }

            data.use <- NULL
            if (feature.type=="gene") {
              data.use <- getExpression(object, scaled=TRUE, genes=getGenesForClustering(object))
              if (all(is.na(data.use))) {
                stop("The scaled expression profiles of the genes for clustering were empty\n")
              }
            } else if (feature.type=="pca") {
              data.use <- t(object@pca$rds[, rds.dims])
            } else if (feature.type=="tsne") {
              data.use <- t(object@tsne$rds[, rds.dims])
            }


            clusters <- NULL
            hc.obj <- NULL
            hc.k <- 1

            if (clustering.method=="tight") {

              if (verbose) {
                cat("\tUsing tight clustering to find ", as.numeric(target), " cell clusters \n", sep="")
              }

              ret <- NULL
              clusters <- clustering.tight(data.use, ...)

              object <- setCellMeta(object, name="CLUSTER", value=clusters)

            } else if (clustering.method=="consensus") {

              if (verbose) {
                cat("Using consensus clustering", sep="")
              }

              clusters <- clustering.consensus(data.use, ...)

              object <- setCellMeta(object, name="CLUSTER", value=clusters)

            } else if (clustering.method=="hc") {

              if (verbose) {
                cat("Using hierarchical clustering", sep="")
              }

              ret <- clustering.hc(data.use, ...)

              object@hc.obj <- ret$hc.obj
              object@hc.k <- ret$kk

              object <- setCellMeta(object, name="CLUSTER", value=ret$clusters)

              plotHC(object)

            } else if (clustering.method=="graph") {

              clusters <- clustering.graph(data.use, ...)
              object <- setCellMeta(object, name="CLUSTER", value=clusters)

            } else if (clustering.method=="pam") {

              clusters <- clustering.pam(data.use, ...)
              object <- setCellMeta(object, name="CLUSTER", value=clusters)

            } else if (clustering.method=="kmeans") {

              clusters <- clustering.kmeans(data.use, ...)
              object <- setCellMeta(object, name="CLUSTER", value=clusters)

            }


            if (update.cellgroup==TRUE) {
              object <- copyCellMeta(object, from="CLUSTER", to="GROUP")
            }

            return(object)
          }
)


#' Compare and visualze multiple clustering assignments
#'
#'
#' @param m The data.frame encoding different clustering results. Rows are cells, Columns are clustering results
#' @param consistency.thresh For a cell, if less than consistency.thresh percent of its neighbors in clustering A remains in clustering B, the instability of the clustering assignment of the cell will increase
#' @param show.cell If TRUE, show cell names in the plot
#' @param bar.width The width of bar
#' @param font.size The font size of text in the plot
#' @return The data frame containing cell clustering assignments and instability scores
#' @export
#'
clusterings.compare <- function(m, consistency.thresh=0.5, do.plot=T, show.cell=FALSE, bar.width=1, font.size=12) {

  #description: plots cluster tracking plot
  #input: m - matrix where rows are k, columns are samples, and values are cluster assignments.
  mt <- apply(m, 2, function(x) factor(x))
  colnames(mt) <- colnames(m)
  rownames(mt) <- rownames(m)
  for (i in 1:dim(mt)[2]) {mt <- mt[order(mt[, i]), ]}



  mnames <- colnames(mt)
  mnames.pairs <- combn(mnames, 2)

  mi <- data.frame(Cell=1:dim(mt)[1], Instability=0)

  clusterings.sim <- matrix(NA, nrow=length(mnames), ncol=length(mnames))
  rownames(clusterings.sim) <- colnames(clusterings.sim) <- mnames

  clusterings.sim.helper <- function(x,y) {

    ret <- 0

    x.groups <- unique(x)
    y.groups <- unique(x)

    jsim <- matrix(0, ncol=length(y.groups), nrow=length(x.groups))
    colnames(jsim) <- paste("C", y.groups, sep="")
    rownames(jsim) <- paste("C", x.groups, sep="")


    for (i in 1:length(x.groups)) {
      for (j in 1:length(y.groups)) {
        i.set <- which(x==x.groups[i])
        j.set <- which(y==y.groups[j])
        jsim[i, j] <- length(intersect(i.set, j.set))/length(union(i.set, j.set))
      }
    }

    ret <- sum(jsim)/max(dim(jsim))

    return(ret)
  }

  for (j in 1:dim(mt)[1]) {

    for (i in 1:dim(mnames.pairs)[2]) {
      i.m1 <- mnames.pairs[1, i]
      i.m2 <- mnames.pairs[2, i]

      i.n1 <- length(unique(mt[, i.m1]))
      i.n2 <- length(unique(mt[, i.m2]))

      if (i.n1 < i.n2) {
        tmp <- i.m1
        i.m1 <- i.m2
        i.m2 <- tmp
      }

      j.m1 <- mt[j, i.m1]
      j.m2 <- mt[j, i.m2]

      j.set1 <- which(mt[, i.m1]==j.m1)
      j.set2 <- which(mt[, i.m2]==j.m2)
      j.set1in2 <- which(j.set1 %in% j.set2)

      j.set1in2.p <- length(j.set1in2)/length(j.set1)

      if (j.set1in2.p < consistency.thresh) {
        mi$Instability[j] <- mi$Instability[j]+1
      }

      i.m1.clustering <- mt[, i.m1]
      i.m2.clustering <- mt[, i.m2]


      clusterings.sim[i.m1, i.m2] <- clusterings.sim[i.m2, i.m1] <- clusterings.sim.helper(i.m1.clustering, i.m2.clustering)
    }

  }

  mi$Instability <- mi$Instability/dim(mnames.pairs)[2]


  if (do.plot) {

    viz <- as.data.frame(mt)
    viz$ID <- 1:dim(viz)[1]
    viz <- melt(viz, id.vars="ID")
    colnames(viz) <- c("ID","Clustering", "Cluster")
    g <- ggplot(viz, aes(x=ID, y=0.5)) + facet_wrap(~Clustering, ncol=1, strip.position = "left")
    #g <- g + ggtitle("Clustering comparison")
    #g <- g + geom_point(aes(x=Cell, y=1, col=Cluster))
    g <- g + geom_bar(stat="identity",aes(fill=Cluster, col=Cluster), width=bar.width)
    g <- g + ylab("Method") + xlab("Cell")
    g <- g + scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE)
    g <- g + theme_bw() + theme(panel.grid = element_blank())
    g <- g + theme(text=element_text(size=font.size))
    if (show.cell) {
      g <- g + theme(axis.text.x=element_text(angle = 90, vjust=0.5, hjust = 1))
    } else {
      g <- g + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
    }
    g <- g + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())
    g <- g + theme(panel.spacing = unit(0.1, "lines"))
    #g <- g + theme(panel.border = element_blank())
    g1 <- g


    g <- ggplot(mi, aes(x=Cell, y=Instability)) + geom_line(col="grey") + geom_area(fill="grey")
    g <- g + ggtitle("Clustering comparison") + ylab("Co-assignment Instability")
    g <- g + theme_bw() + theme(panel.grid = element_blank())
    g <- g + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    #g <- g + theme(panel.border = element_blank())
    g2 <- g


    instability.dist <- data.frame(table(factor(round(mi$Instability, 2))))
    colnames(instability.dist) <- c("Instability","Count")
    instability.dist$Percent <- round(instability.dist$Count/sum(instability.dist$Count), 3)

    g <- ggplot(instability.dist, aes(x=Instability, y=Percent, fill=Instability))
    g <- g + geom_bar( stat="identity")
    g <- g + ggtitle("Distribution of cell cluster assignment instability")
    g <- g + sincera_theme()
    g3 <- g

    library(pheatmap)
    library("RColorBrewer")
    colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
    g4 <- pheatmap(clusterings.sim, silent=T, col=colors)$gtable


    library(gridExtra)
    library(grid)
    grid.arrange(g2, g3, g1, g4, ncol=2)

  }

  mt <- as.data.frame(mt)
  mt$Instability <- mi$Instability

  return(mt)
}




#' Perform PCA based dimension reduction
#'
#' @param object A sincera object
#' @param genes The set of genes used for PCA
#' @param use.scaled If TRUE, use the scaled expression for PCA
#' @param use.fast If TRUE, use gmodels::fast.prcomp; otherwise, use stats::prcomp
#' @return The update sincera object with PCA results in the pca slot
#'
setGeneric("doPCA", function(object, genes=NULL, use.scaled=T, use.fast=T, ...) standardGeneric("doPCA"))
#' @export
setMethod("doPCA","sincera",
          function(object, genes=NULL, use.scaled=T, use.fast=T, ...) {

            cat("\nPerforming dimension reduction using PCA\n")

            expr=getExpression(object, scaled=use.scaled)

            pc.genes = getGenes(object)
            if (!is.null(genes)) {
              genes <- genes[which(genes %in% pc.genes)]
              if (length(genes)<3) stop("Too few genes selected for PCA")
              pc.genes <- genes
            }

            #Remove genes with zero variation
            pc.genes.var = apply(expr[pc.genes,],1,function(x) var(x))
            pc.genes = pc.genes[pc.genes.var>0]
            if (length(pc.genes)<3) stop("Too few non-zero-variance genes for PCA")
            expr = expr[pc.genes,]

            if (use.fast==T) {
              pca.obj = fast.prcomp(t(expr), ...)
            } else {
              pca.obj <- prcomp(t(expr), ...)
            }
            object <- setPCA(object, name="rds", value=pca.obj$x)
            object <- setPCA(object, name="loadings", value=pca.obj$rotation)
            object <- setPCA(object, name="sdev", value=pca.obj$sdev)

            return(object)
          }
)


#' Perform tSNE based dimension reduction
#'
#' Perform tSNE analysis on gene expression space or pca space
#'
#' @param object A sincera object
#' @param n The number of tSNE dimensions to be extracted
#' @param use.scaled If TRUE, use the scaled expression for tSNE
#' @param use.fast If TRUE, use Rtsne::Rtsne; otherwise, use tsne::tsne
#' @param max.iter The maximum number of iteration
#' @param seed The seed of randomness
#' @param feature.type The input feature type for tSNE reduction, possible values include gene and pca
#' @param dims If feature.type is pca, the pca components for tSNE
#' @param genes If feature.type is genes, the set of genes for tSNE; if NULL, set to all genes
#' @param use.scaled If TRUE, use the scaled expression values of genes for tSNE; only take effect when the feature.type is gene
#' @return The update sincera object with tsne results in the tsne slot
#'
setGeneric("doTSNE", function(object, n=2, use.fast=T, max.iter=2000, seed=0, feature.type="pca", dims=1:3, genes=NULL, use.scaled=F, ...) standardGeneric("doTSNE"))
#' @export
setMethod("doTSNE","sincera",
          function(object, n=2, use.fast=T, max.iter=2000, seed=0, feature.type="pca", dims=1:3, genes=NULL, use.scaled=F, ...) {

            cat("\nPerforming dimension reduction using tSNE\n")

            if (feature.type=="gene") {
              data.use=getExpression(object, scaled=use.scaled)
              pc.genes = getGenes(object)
              if (!is.null(genes)) {
                genes <- genes[which(genes %in% pc.genes)]
                if (length(genes)<3) stop("Too few genes selected for tSNE")
                pc.genes <- genes
              }
              #Remove genes with zero variation
              pc.genes.var = apply(data.use[pc.genes,],1,function(x) var(x))
              pc.genes = pc.genes[pc.genes.var>0]
              if (length(pc.genes)<3) stop("Too few non-zero-variance genes for tSNE")
              data.use = data.use[pc.genes,]
              data.use <- t(data.use)
            } else if (feature.type=="pca") {
              data.use <- getPCA(object, name="rds")
              data.use <- data.use[, dims]
            }

            set.seed(seed)
            if (use.fast==T) {
              tsne.obj=Rtsne(as.matrix(data.use),dims=n,...)
              tsne.obj <- tsne.obj$Y
            } else {
              tsne.obj=tsne(data.use,k=n, max_iter = max.iter,...)
            }
            tsne.obj <- data.frame(tsne.obj, check.names=FALSE)
            rownames(tsne.obj) <- rownames(data.use)
            colnames(tsne.obj) <- paste("tSNE", 1:n, sep="")
            object <- setTSNE(object, name="rds", value= tsne.obj)

            return(object)
          }
)



#' Validating clustering results using permutation analysis
#'
#' @param object A sincera object
#' @param n The number of random permutations
#' @param distance.method The method to calculate cell expression profile distance, possible values include pearson and euclidean
#'
setGeneric("cluster.permutation.analysis", function(object, n=20, distance.method="euclidean", ...) standardGeneric("cluster.permutation.analysis"))
#' @export
setMethod("cluster.permutation.analysis","sincera",
          function(object, n=20, distance.method="euclidean", ...) {

            es <- getES(object)
            es <- es[getGenesForClustering(object), ]
            cluster.permutation.analysis.old(es, group.by="GROUP", n=n, distance.method=distance.method, log.base=2, verbose=TRUE)

            return(object)
          }
)


#' Detect cluster specific differentially expressed genes
#'
#' @param object A sincera object
#' @param groups The cell groups included in the differential expression analysis; if NULL, set to all groups defined in the object
#' @param genes The set of genes included in the analysis; if NULL, set to all genes in the object
#' @param method The method for differential test, possible values include welch - one tailed welch's t-test, wilcoxon - one-tailed wilcoxon rank sum test
#' @param do.fdr If TRUE, do the FDR correction
#' @param thresh The threshold for significance
#' @return The update sincera object with diff test results in the difftests slot and the significant diff genes in the diffgenes slot
#'
setGeneric("cluster.diffgenes", function(object, groups=NULL, genes=NULL, method="welch", do.fdr=FALSE, thresh=0.05, ...) standardGeneric("cluster.diffgenes"))
#' @export
setMethod("cluster.diffgenes","sincera",
          function(object, groups=NULL, genes=NULL, method="welch", do.fdr=FALSE, thresh=0.05, ...) {

            cellgroup <- getCellMeta(object, name="GROUP")
            if (is.null(groups)) {
              groups <- sort(unique(cellgroup))
            }
            ng <- length(groups)

            if (is.null(genes)) {
                genes <- getGenes(object)
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


            expr <- getExpression(object)
            expr <- expr[which(rownames(expr) %in% genes), ]

            results <- data.frame(SYMBOL=rownames(expr))
            rownames(results) <- results$SYMBOL

            cat("\nPerforming differential expression using ", method, " test for groups: ", sep="")

            idx.all <- 1:getCellNum(object)
            if (method=="welch") {
                for(g in groups) {
                    idx.g <- which(cellgroup == g)
                    idx.ng <- setdiff(idx.all, idx.g)
                    if (length(idx.g)>=2 & length(idx.ng)>=2) {
                        g.name <- paste(method, ".", g, sep="")
                        results[, g.name] <- NA
                        results[, g.name] <- apply(expr, 1, FUN=welch.test.helper, idx.i=idx.g, idx.o=idx.ng)
                    } else {
                        warning(paste("Skip group ", g, " since it only contains one cell\n", sep=""))
                    }
                    cat(" ", g)
                }

            } else if (method=="wilcox") {
                for(g in groups) {
                  idx.g <- which(cellgroup == g)
                  idx.ng <- setdiff(idx.all, idx.g)
                  if (length(idx.g)>=2 & length(idx.ng)>=2) {
                    g.name <- paste(method, ".", g, sep="")
                    results[, g.name] <- NA
                    results[, g.name] <- apply(expr, 1, FUN=wilcoxon.test.helper, idx.i=idx.g, idx.o=idx.ng)
                  } else {
                    warning(paste("Skip group ", g, " since it only contains one cell\n", sep=""))
                  }
                  cat(" ", g)
                }
            }


            if (do.fdr==T) {
              test.cols <- grep(paste("^",method, sep=""), colnames(results), value=T)

              fdr.in <- NULL
              if (length(test.cols)==1) {
                fdr.in <- data.frame(v1=results[, test.cols], check.names=FALSE)
                colnames(fdr.in)[1] <- test.cols
              } else if (length(test.cols)>1){
                fdr.in <- results[, test.cols]
              }
              if (!is.null(fdr.in)) {
                fdr <- apply(fdr.in, 2, function(x) p.adjust(x, method="BH"))
                colnames(fdr) <- paste(test.cols, ".fdr", sep="")
                rownames(fdr) <- rownames(results)
                results <- cbind(results, fdr)
              }

            }

            object <- setDiffTest(object, value=results, method=method)

            diffgenes <- getDiffGenes(object, groups=groups, method=method, use.fdr=do.fdr, thresh=thresh, print.summary = T)
            object <- setDiffGenes(object, value=diffgenes)

            return(object)
          }
)

#' Perform cell type enrichment
#'
#' Perform cell type enrichment
#'
#' @param object A sincera oject
#' @param species The species of the data, possible values include MUSMU - mouse, HOMSA - human
#' @param id.type The type of gene identifier, possible values include SYMBOL - Entrez SYMBOL, EG - Entrez ID, ENSEMBL - Ensembl ID
#' @return The updated sincera object with enrichment results in the cte slot
#'
setGeneric("celltype.enrichment", function(object, species="MUSMU", id.type="SYMBOL",  ...) standardGeneric("celltype.enrichment"))
#' @export
setMethod("celltype.enrichment","sincera",
          function(object, species="MUSMU", id.type="SYMBOL",  ...) {

            ES <- getES(object)

            cat("\nUsing differentially expressed genes for cell type enrichment analysis\n")

            diffgenes <- getDiffGenes(object, print.summary = TRUE)

            # selecting cluster specific differentially expressed genes for cell type enrichment analysis
            groups <- sort(unique(getCellMeta(object, name="GROUP")))
            prefix="use_for_celltype_"
            for (i in groups) {
              i.de.name <- paste(prefix, i, sep="")
              i.diffgenes <- c()
              i.diffgenes <- as.character(diffgenes$SYMBOL[which(diffgenes$GROUP == i)])
              fData(ES)[,i.de.name] <- 0
              fData(ES)[which(rownames(fData(ES)) %in% i.diffgenes),i.de.name] <- 1
            }

            # initialize knowledge base for cell type enrichment analysis

            genome <- NULL
            associations <- NULL
            KB <- NULL

            genome <- rownames(object@rexprs)
            associations <- getAssociationTable(object)


            #' Contruct the knowledge base for cell type enrichment analysis
            #'
            #' associations (data.frame) a data frame containing the gene and cell type associations
            #' genome (character) the full set of genes in a scRNA-seq data
            #' species (character) the species of the genome: MUSMU - mouse, HOMSA - human
            #' id.type (character) the type of ids of the genes in the genome: ENSEMBL - ENSEMBL gene id, SYMBOL -  Entrez Gene Symbol, EG - Entrez Gene Id
            #' return a list of 3 items: celltype.gene.association - genome related associations, celltype.genome.count - genome-wide cell type associations, genome - symbol mapping for the genome
            KB <- celltype.enrichment.initKB(associations, genome, species=species, id.type=id.type)

            cat("\n")
            #' Contruct the knowledge base for cell type enrichment analysis
            #'
            #' ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
            #' KB (list) the knowledge base prepared by celltype.enrichment.initKB()
            #' group.by (character) the name of the column that contains the cluster information
            #' groups (character) the clusters for cell type enrichment
            #' species (character) the species of the genome: MUSMU - mouse, HOMSA - human
            #' id.type (character) the type of ids of the genes in the genome: ENSEMBL - ENSEMBL gene id, SYMBOL -  Entrez Gene Symbol, EG - Entrez Gene Id
            #' celltype.enrichment.prefix (character) the prefix of columns encoding the cluster-specific gene list for cell type enrichment analysis
            ret <- celltype.enrichment.old(ES, KB, group.by="GROUP", groups=NULL, species=species, id.type=id.type, celltype.enrichment.prefix="use_for_celltype_", verbose=TRUE)

            cat("\nTop 5 enriched cell type annotations per cluster:\n")
            for (i in 1:length(ret)) {
              cat("\nCluster", names(ret)[i],":\n")
              tmp <- ret[[i]]
              rownames(tmp) <- NULL
              print(tmp[1:5, c("TYPE","Fisher.PV")])
              cat("\n")
            }

            object <- setCellTypeEnrichment(object, groups=names(ret), ret)

            return(object)
          }
)

#' Performing rank-aggregation based cell type validation
#'
#' Performing rank-aggregation based cell type validation
#'
#' @param object A sincera object
#' @param use.scaled If TRUE, use zscore scaled data
#' @param thresh The expression threshold to include cells in individual rankings
#' @return The updated sincera object
#'
setGeneric("celltype.validation", function(object, use.scaled=T, thresh=0, ...) standardGeneric("celltype.validation"))
#' @export
setMethod("celltype.validation","sincera",
          function(object, use.scaled=T, thresh=0, ...) {

              ESz <- getES(object)

              if (use.scaled) {
                exprs(ESz) <- as.matrix(getExpression(object, scaled = T))
              }

              markers <- getCellTypeMarkers(object)

              types <- getCellType(object)
              types <- types[which(types %in% as.character(markers$TYPE))]
              groups <- as.character(names(types))

              prefix <- "use_as_marker_"

              for (i in 1:length(types)) {
                  i.type <- types[i]
                  i.group <- as.character(names(i.type))

                  i.markers <- as.character(markers$SYMBOL[which(markers$TYPE==as.character(i.type))])
                  i.name <- paste(prefix, i.group, sep="")
                  fData(ESz)[, i.name] <- 0
                  fData(ESz)[which(rownames(fData(ESz)) %in% i.markers), i.name] <- 1
              }

              #' Cell type validation using the expression patterns of known biomarkers
              #'
              #' ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
              #' group.by (character) the name of the column that contains the cluster information
              #' groups (character) the clusters for cell type validation
              #' threshold (numeric) the threshold of expression; 0 means that markers will not be used to rank cells where their expression < 0
              #' marker.prefix (character) the prefix of columns encoding the cell type markers
              celltype.validation.old(ESz, group.by="CLUSTER", groups=groups, marker.prefix=prefix, threshold=thresh)


            return(object)
          }
)


########## clustering algorithms #####################


clustering.tight <- function(x, ...) {
  ret <- NULL

  ret <- tight.clust(x, ...)

  clusters <- ret$cluster
  names(clusters) <- colnames(x)

  return(clusters)
}


clustering.consensus <- function(x, min.area.increase=0.2, ...) {


  ret <- NULL
  ret <- ConsensusClusterPlus(x, ..)

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
  names(clusters) <- colnames(data.use)

  return(clusters)
}


clustering.hc <- function(x, h=NULL, k=NULL, num.singleton=0, distance.method="pearson", linkage.method="average",  do.shift=TRUE) {

  if (do.shift) {
    x <- apply(x, 2, function(y) y-mean(y))
  }

  if (distance.method=="euclidean") {
    dd <- dist(t(x), method=distance.method)
  } else  {
    dd <- as.dist((1-cor(x, method=distance.method))/2)
  }

  hc <- hclust(dd, method=linkage.method)

  cell_order <- hc$labels[hc$order]


  if (is.null(k)) {
    if (is.null(h) || h<=0) {
      kk <- 1
      for (i in 2:dim(x)[2]) {
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

  names(clusters) <- colnames(x)

  return(list(clusters=clusters, kk=kk, hc.obj=hc))

}

# From Shekhar et al., Cell 2016
clustering.graph <- function(x, do.jaccard=T, num.nn=30, community.method="louvain") {

  weights=NULL;
  if (do.jaccard){
    weights=TRUE;
  }

  Adj = get_edges(t(x), nn=num.nn, do.jaccard=do.jaccard)

  g=igraph::graph.adjacency(Adj, mode = "undirected", weighted=weights)
  
  if (community.method=="louvain") graph.out = igraph::cluster_louvain(g)
  if (community.method=="infomap") graph.out = igraph::cluster_infomap(g)
  if (community.method=="fast_greedy") graph.out = igraph::cluster_fast_greedy(g)
  # if (community.method=="leading_eigen") graph.out = igraph::cluster_leading_eigen(g)
  if (community.method=="optimal") graph.out = igraph::cluster_optimal(g)
  if (community.method=="spinglass") graph.out = igraph::cluster_spinglass(g)
  if (community.method=="label_prop") graph.out = igraph::cluster_label_prop(g)
  if (community.method=="edge_betweenness") graph.out = igraph::cluster_edge_betweenness(g)
  if (community.method=="walktrap") graph.out = igraph::cluster_walktrap(g)

  clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
  names(clust.assign) = graph.out$names
  k=order(table(clust.assign), decreasing = TRUE)
  new.levels = rep(1,length(unique(graph.out$membership)))
  new.levels[k] = 1:length(unique(graph.out$membership))
  levels(clust.assign) = new.levels
  clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))

  return(clust.assign)
}

# From Shekhar et al., Cell 2016
# Build a nearest neighbor graph with or without edge weights, and return an adjacency matrix
get_edges=function(X,nn=30,do.jaccard=TRUE) {
  nearest=RANN::nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
  print("Found nearest neighbors")
  nearest$nn.idx = nearest$nn.idx[,-1]
  nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim = 1*(nearest$nn.dists >= 0 )


  edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
  edges$B = edges$C; edges$C=1

  #Remove repetitions
  edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))

  if (do.jaccard){

    NN = nearest$nn.idx
    jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )

    edges$C = jaccard_dist
    edges = subset(edges, C != 0)
    edges$C = edges$C/max(edges$C)
  }

  Adj = matrix(0, nrow=nrow(X), ncol=nrow(X))
  rownames(Adj) = rownames(X); colnames(Adj) = rownames(X)
  Adj[cbind(edges$A,edges$B)] = edges$C
  Adj[cbind(edges$B,edges$A)] = edges$C
  return(Adj)

}


clustering.pam <- function(x, k=4, ...) {
  require(cluster)
  pam.obj <- pam(t(x), k=k, ...)
  clusters <- pam.obj$clustering
  return(clusters)
}


clustering.kmeans <- function(x, k=4, ...) {
  km.obj <- kmeans(t(x), centers=k, ...)
  clusters <- km.obj$cluster
  return(clusters)
}

