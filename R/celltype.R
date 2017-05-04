

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
#' @param m The data.frame encoding different clustering results. Rows are cells, Columns are clustering results. Column names denote different clustering assignments.
#' @param consistency.thresh For a cell, if less than consistency.thresh percent of its neighbors in clustering A remains in clustering B, the instability of the clustering assignment of the cell will increase
#' @param show.cell If TRUE, show cell names in the plot
#' @param bar.width The width of bar
#' @param font.size The font size of text in the plot
#' @return The data frame containing cell clustering assignments and instability scores
#' @export
#' @details
#' The similarity between two clustering assignments is measured based on the metric proposed in Torres et al., International Journal of Electrical, Computer & Systems Engineer, 2009.
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
#' @param use.scaled If TRUE, use the scaled expression for PCA. 
#' @param use.fast If TRUE, use gmodels::fast.prcomp(); otherwise, use stats::prcomp()
#' @return The update sincera object with PCA results in the pca slot
#' @details 
#' When use.fast is TRUE, please refer to ?gmodels::fast.prcomp for more parameters to control the PCA
#' When use.fast is FALSE, please refer to ?stats:prcomp for more parameters to control the PCA
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
              pca.obj = gmodels::fast.prcomp(t(expr), ...)
            } else {
              pca.obj <- stats::prcomp(t(expr), ...)
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
#' @details
#' When feature.type is "gene", genes with zero variance will be excluded prior to calling tSNE functions.
#' When use.fast is TRUE, please refer to ?Rtsne::Rtsne for more parameters to control tSNE reduction.
#' When use.fast is FALSE, please refer to ?tsne::tsne for more parameters to control tSNE reduction.
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
              tsne.obj=Rtsne::Rtsne(as.matrix(data.use),dims=n,...)
              tsne.obj <- tsne.obj$Y
            } else {
              tsne.obj=tsne::tsne(data.use,k=n, max_iter = max.iter,...)
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
            cluster.permutation.analysis.1(es, group.by="GROUP", n=n, distance.method=distance.method, log.base=2, verbose=TRUE)

            return(object)
          }
)
#' Permutation Analysis for determining significance of cluster assignments
#'
#' @param ES (ExpressionSet) an ExpressionSet object containing the single cell RNA-seq data
#' @param group.by (character) the name of the column that contains the cluster information
#' @param n (numeric) the number of permutations
#' @param distance.method (character) the distance method: pearson or spearman - (1-correlation)/2; euclidean - euclidean distance
#' @param log.base (numeric) the base of logorithm; if log.base <=1 or log.base is NULL, no log transformations will be applied
#' @param verbose (logical)
#' @return NULL
cluster.permutation.analysis.1 <- function(ES, group.by="GROUP", n=20, distance.method="euclidean", log.base=2, verbose=TRUE) {
  
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
pa.helper <- function(ES, group.by="GROUP", log.base=2, cs, a, distance.method="euclidean") {
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

# to be updated
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




#' Perform cell type enrichment
#'
#' Perform cell type enrichment
#'
#' @param object A sincera oject
#' @param species The species of the data, possible values include MUSMU - mouse, HOMSA - human
#' @param id.type The type of gene identifier, possible values include SYMBOL - Entrez SYMBOL, EG - Entrez ID, ENSEMBL - Ensembl ID
#' @param do.plot If TRUE, plot top enriched cell types per group at the end of the analysis
#' @param top.k The number of top enriched cell types to be plotted at the end of the analysis
#' @param verbose If TRUE, print verbose messages
#' @return The updated sincera object with enrichment results in the cte slot
#'
setGeneric("celltype.enrichment", function(object, species="MUSMU", id.type="SYMBOL", do.plot=T, top.k=5, verbose=T, ...) standardGeneric("celltype.enrichment"))
#' @export
setMethod("celltype.enrichment","sincera",
          function(object, species="MUSMU", id.type="SYMBOL", do.plot=T, top.k=5, verbose=T, ...) {
            
            if (verbose) {
              cat("Sincera: cell type enrichment analysis ... \n")
            }
            
            cat("\nUsing differentially expressed genes for cell type enrichment analysis\n")
            diffgenes <- getDiffGenes(object, print.summary = TRUE)
            
            #wd <- paste(getwd(), dir.delim, "sincera.celltype.enrichment.", getTimestamp(), dir.delim, sep="")
            #dir.create(wd)
            

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
            # ret <- celltype.enrichment.old(ES, KB, group.by="GROUP", groups=NULL, species=species, id.type=id.type, celltype.enrichment.prefix="use_for_celltype_", verbose=TRUE)
            
            
            
            
            groups <- sort(unique(as.character(getCellMeta(object, name="GROUP"))))
            
            
            x <- KB$celltype.gene.association
            types.count <- KB$celltype.genome.count
            genome.type <- paste(species,".", id.type, sep="")
            
            ret <- list()
            
            for (i in groups) {
              if (verbose) {
                cat("Enriching cell types for group", i, "...")
              }
              #i.de.name <- paste(celltype.enrichment.prefix, i, sep="")
              #i.de <- rownames(fData(ES))[which(fData(ES)[, i.de.name]==1)]
              i.de <- as.character(diffgenes$SYMBOL[which(diffgenes$GROUP == i)])
              
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
                  
                  
                  i.types.count$Hits[j] <- paste(i.j.de, collapse=",")

                  if (FALSE) { #%%
                    if (genome.type =="MUSMU.ENSEMBL") {
                      i.j.de.originalid <- i.j.de
                    } else {
                      i.j.de.originalid <- unique(KB$genome[which(KB$genome[,"MUSMU.ENSEMBL"] %in% i.j.de),genome.type])
                    }
                    i.j.de.symbols <- unique(fData(ES)[which(rownames(fData(ES)) %in% i.j.de.originalid), GENE.SYMBOL.LABEL]) #%%
                    i.types.count$Hits_MUSMU_ENSEMBL[j] <- paste(i.j.de, collapse=",")
                    i.types.count$Hits_SYMBOL[j] <- paste(i.j.de.symbols, collapse=",")
                  }
                  
                  
                  
                  h <- h+1
                }
              }
              
              i.types.count <- i.types.count[order(i.types.count$Fisher.PV), ]
              ret[[as.character(i)]] <- i.types.count
              
              # write.table(i.types.count, file=paste(wd, i, "-celltype-enrichment.txt", sep=""), sep="\t", col.name=T, row.name=F)
              if (verbose) {
                cat("done\n")
              }
            }
           
            cat("\nTop 5 enriched cell type annotations per cluster:\n")
            
            for (i in 1:length(ret)) {
              cat("\nCluster", names(ret)[i],":\n")
              tmp <- ret[[i]]
              rownames(tmp) <- NULL
              
              i.viz <- tmp[1:5, c("TYPE","Fisher.PV")]
              
              print(i.viz)
              
              cat("\n")
            }
            
            
           
            object <- setCellTypeEnrichment(object, groups=names(ret), ret)
            
            plotCellTypeEnrichment(object, top.k=top.k)

            cat("Please use getCellTypeEnrichment() to assess full enrichment results.")
            
            return(object)
          }
)

#' Contruct the knowledge base for cell type enrichment analysis
#'
#' @param associations (data.frame) a data frame containing the gene and cell type associations
#' @param genome (character) the full set of genes in a scRNA-seq data
#' @param species (character) the species of the genome: MUSMU - mouse, HOMSA - human
#' @param id.type (character) the type of ids of the genes in the genome: ENSEMBL - ENSEMBL gene id, SYMBOL -  Entrez Gene Symbol, EG - Entrez Gene Id
#' @param verbose (logical)
#' @return a list of 3 items: celltype.gene.association - genome related associations, celltype.genome.count - genome-wide cell type associations, genome - symbol mapping for the genome
#'
celltype.enrichment.initKB <- function(associations, genome, species="MUSMU", id.type="ENSEMBL", verbose=T) {
  
  species.supported <- c("MUSMU", "HOMSA")
  s.id <- pmatch(species, species.supported)
  if (is.na(s.id)) {
    stop("invalid species")
  }
  species <- species.supported[s.id]
  
  id.supported <- c("ENSEMBL", "SYMBOL", "EG")
  id.id <- pmatch(id.type, id.supported)
  if (is.na(id.id)) {
    stop("invalid id type")
  }
  id.type <- id.supported[id.id]
  
  # convert the genome to MUSMU ENSEMBL
  genome.mapping=NULL
  
  if (!(species=="MUSMU" & id.type=="ENSEMBL")) {
    
    if (species=="MUSMU") {
      
      if (!require(AnnotationDbi)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("AnnotationDbi")
      }
      if (!require(org.Mm.eg.db)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("org.Mm.eg.db")
      }
      
      if (id.type=="SYMBOL") {
        y <- org.Mm.egSYMBOL2EG
        mapped_seqs <- mappedkeys(y)
        yy <- as.list(y[mapped_seqs])
        symbol2eg <- unlist(yy[genome])
        eg2ensembl <- idConverter(symbol2eg, srcSpecies="MUSMU", destSpecies="MUSMU", srcIDType="EG", destIDType="ENSEMBL")
        
        df1 <- data.frame(MUSMU.SYMBOL=names(symbol2eg), MUSMU.EG=as.character(symbol2eg))
        df2 <- data.frame(MUSMU.EG=names(eg2ensembl), MUSMU.ENSEMBL=as.character(eg2ensembl))
        
        genome.mapping <- merge(df1, df2, by.x="MUSMU.EG", by.y="MUSMU.EG")
      } else if (id.type=="EG") {
        eg2ensembl <- idConverter(genome, srcSpecies="MUSMU", destSpecies="MUSMU", srcIDType="EG", destIDType="ENSEMBL")
        genome.mapping <- data.frame(MUSMU.EG=names(eg2ensembl), MUSMU.ENSEMBL=as.character(eg2ensembl))
      }
    } else if (species=="HOMSA") {
      
      if (!require(AnnotationDbi)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("AnnotationDbi")
      }
      if (!require(org.Mm.eg.db)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("org.Mm.eg.db")
      }
      if (!require(hom.Hs.inp.db)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("hom.Hs.inp.db")
      }
      if (!require(hom.Mm.inp.db)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("hom.Mm.inp.db")
      }
      if (!require(org.Hs.eg.db)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("org.Hs.eg.db")
      }
      
      if (id.type=="SYMBOL") {
        y <- org.Hs.egSYMBOL2EG
        mapped_seqs <- mappedkeys(y)
        yy <- as.list(y[mapped_seqs])
        symbol2eg <- unlist(yy[genome])
        eg2ensembl <- idConverter(symbol2eg, srcSpecies="HOMSA", destSpecies="MUSMU", srcIDType="EG", destIDType="ENSEMBL")
        # MouseEGs = inpIDMapper(symbol2eg, "HOMSA", "MUSMU",  srcIDType="EG",
        #                         destIDType="EG", keepMultGeneMatches=FALSE, keepMultProtMatches=FALSE,
        #                         keepMultDestIDMatches = TRUE)
        
        
        df1 <- data.frame(HOMSA.SYMBOL=names(symbol2eg), HOMSA.EG=as.character(symbol2eg))
        df2 <- data.frame(HOMSA.EG=names(eg2ensembl), MUSMU.ENSEMBL=as.character(eg2ensembl))
        
        genome.mapping <- merge(df1, df2, by.x="HOMSA.EG", by.y="HOMSA.EG")
        
      } else if (id.type=="EG") {
        eg2ensembl <- idConverter(genome, srcSpecies="HOMSA", destSpecies="MUSMU", srcIDType="EG", destIDType="ENSEMBL")
        genome.mapping <- data.frame(HOMSA.EG=names(eg2ensembl), MUSMU.ENSEMBL=as.character(eg2ensembl))
      } else if (id.type=="ENSEMBL") {
        ensembl2ensembl <- idConverter(genome, srcSpecies="HOMSA", destSpecies="MUSMU", srcIDType="ENSEMBL", destIDType="ENSEMBL")
        genome.mapping <- data.frame(HOMSA.ENSEMBL=names(eg2ensembl), MUSMU.ENSEMBL=as.character(eg2ensembl))
      }
    }
  } else {
    genome.mapping <- data.frame(MUSMU.ENSEMBL=genome)
  }
  
  if (!is.null(genome.mapping) > 0 && !is.null(associations)) {
    associations <- associations[which(associations$GENE.ID %in% genome.mapping$MUSMU.ENSEMBL),]
    genome.mapping <- genome.mapping[which(genome.mapping$MUSMU.ENSEMBL %in% associations$GENE.ID), ]
    # genome association
    celltype.genome.count <- data.frame(table(associations$ANNOTATION.ID))
    colnames(celltype.genome.count) <- c("TYPE","COUNT")
  }
  
  return(list(celltype.gene.association=associations, celltype.genome.count=celltype.genome.count, genome=genome.mapping))
}


#' Performing rank-aggregation based cell type validation using marker expression
#'
#' Performing rank-aggregation based cell type validation using marker expression
#'
#' @param object A sincera object
#' @param use.scaled If TRUE, use zscore scaled data
#' @param thresh The expression threshold to include cells in individual rankings
#' @return The updated sincera object
#'
setGeneric("celltype.validation", function(object, use.scaled=T, thresh=0,  ...) standardGeneric("celltype.validation"))
#' @export
setMethod("celltype.validation","sincera",
          function(object, use.scaled=T, thresh=0, ...) {
            
              
            cat("Use \"TYPE\" metadata as the source of cell type assignments\n")
            
            markers <- getCellTypeMarkers(object)
            if ((!is.data.frame(markers)) | (dim(markers)[1]<2) | (dim(markers)[2]<2)) {
              stop("Invalid cell type marker definition. Please use getCellTypeMarkers() to check the current marker definition or use setCellTypeMarkers() to update it.")
            }
            
            types <- sort(unique(as.character(markers[, 1])))
            
            types <- types[which(types %in% getCellMeta(object, "TYPE"))]
            if (length(types)==0) {
              stop("The defined markers cannot find any corresponding cell types in the object.\nPlease use getCellTypeMarkers() and getCellType() to check the defined markers and cell types in the object.")
            }
            cat("The algorithm is validating the following cell types that are defined in the object and have defined markers: ", paste(types, collapse=", ", sep=""), "\n", sep="")
            
            cells <- getCells(object)
            
            type_validation <- data.frame(Name=NULL, Score=NULL, PREDICTION=NULL, GD=NULL, TYPE=NULL)
            
            n <- length(types)
            plot.dims <- getPlotDims(n)
            
            par.defaults <- par(no.readonly=TRUE)
            par(mfrow=c(plot.dims$nrow, plot.dims$ncol))
            
            for (i in 1:length(types)) {
              cat("\tValidating cell type", types[i], ":")
              i.markers <- as.character(markers$SYMBOL[which(markers$TYPE==types[i])])
              if (length(i.markers)<2) {
                cat("skipped due to only one marker was defined.\n")
              } else {
                i.x <- getExpression(object, scaled=use.scaled, genes=i.markers)
                i.list <- list()
                
                for (j in 1:dim(i.x)[1]) {
                  i.j.cells <- i.x[j,]
                  i.j.cells <- i.j.cells[i.j.cells > thresh]
                  if(length(i.j.cells) > 0) {
                    i.j.cells <- sort(i.j.cells, decreasing=T)
                    i.list[[as.character(rownames(i.x)[j])]] <- names(i.j.cells)
                  }
                }
                
                i.cells <- aggregateRanks(glist = i.list, N = length(cells))
                
                i.cells.not.idx <- which(!(cells %in% i.cells$Name))
                if (length(i.cells.not.idx) > 0) {
                  i.cells.not <- cells[i.cells.not.idx]
                  i.cells.not.df <- data.frame(Name=i.cells.not, Score=1)
                  rownames(i.cells.not.df) <- i.cells.not
                  i.cells <- rbind(i.cells, i.cells.not.df)
                }
                
                i.cells$PREDICTION <- -log(i.cells$Score)
                i.cells$PREDICTION <- i.cells$PREDICTION/max(i.cells$PREDICTION)
                i.cells$GD <- 0
                
                i.cells.gold_standard <- cells[which(getCellMeta(object, "TYPE") %in% types[i])]
                i.p.idx <- which(i.cells$Name %in% i.cells.gold_standard)
                
                if (length(i.p.idx) > 0) {
                  i.cells$GD[i.p.idx] <- 1
                }
                i.pred <- prediction(i.cells$PREDICTION, i.cells$GD)
                i.perf <- performance(i.pred,"auc")
                i.auc <- i.perf@y.values[[1]]
                i.perf <- performance(i.pred, "tpr", "fpr")
                
                #tiff(file=paste(wd, "celltype.validation.", types[i], ".tiff",sep=""), width = 1100, height =1200, units = "px", res = 300, compression="lzw", bg = "white")
                
                plot(i.perf, lwd=6, main=paste("ROC of ", as.character(types[i]), " (AUC:", round(i.auc,digit=6), ")", sep=""), cex.lab=1.2, colorize=T)
                axis(side=2, at=seq(0,1,0.1), font=2, lwd=3)
                axis(side=1, at=seq(0,1,0.1), font=2, lwd=3)
                box(lwd=3)
                
                #dev.off()
                
                i.cells$TYPE <- as.character(types[i])
                type_validation <- rbind(type_validation, i.cells)
                
                cat("Done\n")
              }
              
              
            }
            
            par(par.defaults)
              
            
            colnames(type_validation) <- c("CELL", "PVALUE", "SCORE", "TYPE.ASSIGNMENT", "TYPE")
            
            
            object <- setCellTypeValidation(object, value=type_validation)
              
            
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

