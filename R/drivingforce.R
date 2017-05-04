#' Select candidate TFs for network inference and TF ranking
#'
#' @param object A sincera object
#' @param groups A vector of cell groups; if NULL, use all defined cell groups
#' @param tflists A list containing user-provided candidate TFs for each group specified in groups parameter
#' @param diff.method The method used to perform the differential expression analysis
#' @param diff.thresh The pvalue for determining differentially expressed genes
#' @param min.expression The threshold expression value for determining expressed genes
#' @param pct.thresh The percentage threshold for determining genes commonly expressed in a cell group
#' @return An updated sincera object, use getDF with name="tfs" to access the selected TFs
#' @details 
#' Please use the tflists parameter to provide pre-selected TFs
#' Otherwise, the function will select candidate TFs for each group using the following criteria:
#' (1) TFs differentially expressed in the cell group (pvalue < diff.thresh); 
#' Or (2) TFs expressed greater than a certain percentage of the group cells (>min.expression in at least pct.thresh percentage of the group cells)
#
setGeneric("drivingforce.selectTFs", function(object, groups=NULL, tflists=NULL, diff.method="welch", diff.thresh=0.05, min.expression=5, pct.thresh=0.8, ...) standardGeneric("drivingforce.selectTFs"))
#' @export
setMethod("drivingforce.selectTFs","sincera",
          function(object, groups=NULL, tflists=NULL, diff.method="welch", diff.thresh=0.05, min.expression=5, pct.thresh=0.8, ...) {

            cat("\nSelect candidate TFs for cluster:")

            tfs <- getTFs(object)
            if (length(tfs)==0) {
                stop("Please use setTFs function to add TFs to sincera")
            }

            if (is.null(groups)) {
                groups <- sort(unique(getCellMeta(object, name="GROUP")))
            }

            if (is.null(tflists)) {

                tflists <- list()

                pcts <- GeneStats(dp=getExpression(object),
                                   ident=getCellMeta(object, name="GROUP"),
                                   groups=groups,
                                   min.expression=min.expression,
                                   stats=c("ident.pct"))

                diffgenes <- getDiffGenes(object, groups=groups, method=diff.method, use.fdr=FALSE, thresh=diff.thresh, print.summary = FALSE)

                for (i in 1:length(groups)) {

                    i.diff <- c()
                    i.common <- c()

                    i.diff <- as.character(diffgenes$SYMBOL[which(diffgenes$GROUP==groups[i])])
                    i.diff <- intersect(i.diff, tfs)

                    i.common <- rownames(pcts)[which(as.numeric(pcts[, paste("ident.pct.", groups[i], sep="")])>pct.thresh)]
                    i.common <- intersect(i.common, tfs)

                    tflists[[groups[i]]] <- union(i.diff, i.common)

                    tflists[[groups[i]]] <-  intersect(tflists[[groups[i]]], getGenesForClustering(object))
                }

            } else {
                if (length(groups) != length(tflists)) {
                    stop("Inconsistent number of cell groups and the size of TF lists.")
                }
            }

            for (i in 1:length(groups)) {
                object <- setDF(object, group=groups[i], name="tfs", value=tflists[[groups[i]]])
                cat(" ", groups[i], " (", length(tflists[[groups[i]]]), ") ", sep="")
            }
            cat("\n")

            return(object)
          }
)

#' Select candidate TGs for network inference and TF ranking
#'
#' @param object A sincera object
#' @param groups A vector of cell groups; if NULL, use to all defined cell groups
#' @param tglists A list containing user-provided candidate TGs for each group specified in the groups parameter
#' @param diff.method The method used to perform the differential expression analysis
#' @param diff.thresh The pvalue for determining differentially expressed genes
#' @return An updated sincera object, use getDF with name="tgs" to access the selected TGs
#' @details 
#' Please use the tglists parameter to provide pre-selected tgs
#' Otherwise, the function will select candidate TGs for each group using the following criteria: Genes differentially expressed in the cell group (pvalue < diff.thresh)
#'
setGeneric("drivingforce.selectTGs", function(object, groups=NULL, tglists=NULL, diff.method="welch", diff.thresh=0.01, ...) standardGeneric("drivingforce.selectTGs"))
#' @export
setMethod("drivingforce.selectTGs","sincera",
          function(object, groups=NULL, tglists=NULL, diff.method="welch", diff.thresh=0.01, ...) {

              cat("\nSelect candidate TGs for cell group:")

              if (is.null(groups)) {
                groups <- sort(unique(getCellMeta(object, name="GROUP")))
              }

              if (is.null(tglists)) {
                  tglists <- list()

                  diffgenes <- getDiffGenes(object, groups=groups, method=diff.method, use.fdr=FALSE, thresh=diff.thresh, print.summary = FALSE)

                  for (i in 1:length(groups)) {
                      i.diff <- as.character(diffgenes$SYMBOL[which(diffgenes$GROUP==groups[i])])
                      tglists[[groups[i]]] <- i.diff
                  }

              } else {
                  if (length(groups) != length(tflists)) {
                      stop("Inconsistent number of cell groups and the size of TF lists.")
                  }
              }

              for (i in 1:length(groups)) {
                  object <- setDF(object, group=groups[i], name="tgs", value=tglists[[groups[i]]])
                  cat(" ", groups[i], " (", length(tglists[[groups[i]]]), ") ", sep="")
              }
              cat("\n")

              return(object)

          }
)

#' Inferring significant TF->TF and TF->TG interactions based on gene expression patterns
#'
#' @param object A sincera object
#' @param groups The cell groups for network inference
#' @param use.scaled If TRUE, use the zscore scaled expression values for the inference
#' @param use.allcells If TRUE, the expression values of all cells will be used for inference; otherwise, only the expression values in the cell group will be used for inference.
#' @return An update sincera object, use getDF with name="edges" to access the significance of all interactions
#'
setGeneric("drivingforce.inferTRN", function(object, groups=NULL, use.scaled=T, use.allcells=T, ...) standardGeneric("drivingforce.inferTRN"))
#' @export
setMethod("drivingforce.inferTRN","sincera",
          function(object, groups=NULL, use.scaled=T, use.allcells=T, ...) {

              if (is.null(groups)) {
                groups <- sort(unique(getCellMeta(object, name="GROUP")))
              }

              for (i in 1:length(groups)) {
                  i.tfs <- getDF(object, group=groups[i], name="tfs")
                  i.tgs <- getDF(object, group=groups[i], name="tgs")
                  i.cells <- getCells(object, groups=groups[i])
                  if (use.allcells==T) {
                    i.expr <- getExpression(object, scaled=use.scaled, genes=union(i.tfs, i.tgs))
                  } else {
                    i.expr <- getExpression(object, scaled=use.scaled, genes=union(i.tfs, i.tgs), cells=i.cells)
                  }

                  i.tfs.idx <- which(rownames(i.expr) %in% i.tfs)
                  i.tgs.idx <- which(rownames(i.expr) %in% i.tgs)

                  # i.tfs <- rownames(i.expr)[i.tfs.idx]
                  # i.tgs <- rownames(i.expr)[i.tgs.idx]

                  i.expr <- t(i.expr)

                  cat("Inferring TRN for cell group ", groups[i], " using ", length(i.tgs.idx), " TGs and ", length(i.tfs.idx), " TFs: ", sep="")


                  edges <- data.frame(TF=NULL, TG=NULL, s1=NULL)

                  z <- 1
                  cpt <- 10
                  for (u in i.tgs.idx) {
                      # progress report
                      if ( ((z/length(i.tgs.idx)*100))>=cpt ) {
                          cat(cpt, "% ",sep="")
                          cpt=cpt+10
                      }
                      S1 <- network.inference(i.expr, method="ls",predPosition=i.tfs.idx, target=u,lag=0)
                      #write(paste(colnames(i.exprs)[u],paste(S1$S1ls[,],collapse="\t"),collapse="\t",sep="\t"),file=mat.file, append=T)

                      i.edges <- data.frame(TF=colnames(i.expr)[i.tfs.idx], TG=colnames(i.expr)[u], s1=as.numeric(S1$S1ls[1,]), stringsAsFactors = FALSE)
                      edges <- rbind(edges, i.edges)
                      rm(i.edges)
                      z <- z+1
                  }
                  cat("\n")

                  object <- setDF(object, group=groups[i], name="edges", value=edges)
                  rm(edges)
              }

            return(object)
          }
)

#' Extract largest connected component (LCC) of the cell group specific TRN
#'
#' @param object A sincera object
#' @param groups The cell groups for LCC extraction
#' @param thresh The threshold for significant interactions
#' @return An updated sincera object, use getDF with name="lcc" to access the lcc
#'
setGeneric("drivingforce.getLCC", function(object, groups=NULL, thresh=0.05, ...) standardGeneric("drivingforce.getLCC"))
#' @export
setMethod("drivingforce.getLCC","sincera",
          function(object, groups=NULL, thresh=0.05, ...) {

              if (is.null(groups)) {
                groups <- sort(unique(getCellMeta(object, name="GROUP")))
              }

              for (i in 1:length(groups)) {

                  edges <- getDF(object, group=groups[i], name="edges")
                  edges <- edges[which(edges$s1<thresh),]
                  g <- igraph::graph.data.frame(edges, directed=FALSE)
                  g <- igraph::simplify(g)

                  tfs <- getDF(object, group=groups[i], name="tfs")
                  tf.idx <- which(names(igraph::V(g)) %in% tfs)
                  tg.idx <- setdiff(1:length(igraph::V(g)), tf.idx)
                  g <- igraph::set_vertex_attr(g, name="Type", index=tf.idx, value="TF")
                  g <- igraph::set_vertex_attr(g, name="Type", index=tg.idx, value="TG")
                  igraph::V(g)$color <- "grey"
                  igraph::V(g)$color[which(V(g)$Type == "TF")] <- "red"

                #  plot(g, layout=layout.fruchterman.reingold, vertex.label=NA)

                  object <- setDF(object, group=groups[i], name="trn", value=g)


                  cl <- igraph::clusters(g, mode="weak")
                  lcc <- igraph::induced.subgraph(g, which(cl$membership == which.max(cl$csize)))

                  cat("The LCC of ", groups[i], ": ", length(igraph::V(lcc)), " nodes, ", length(igraph::E(lcc)), " edges\n", sep="")

                  #plot(lcc, layout=layout.fruchterman.reingold, vertex.label=NA)

                  if (length(igraph::V(lcc))/(length(igraph::V(g))) < 0.8) {
                      warning("\nCell group ", groups[i], ": less than 80% of tfs and tgs are in the LCC. Please consider increase threshold.", sep="")
                  }

                  object <- setDF(object, group=groups[i], name="lcc", value=lcc)
              }

            return(object)
          }
)

#' Rank TFs based on their importance to the LCC measured by six metrics
#'
#' @param object A sincera object
#' @param groups The cell groups for TF ranking
#' @param metrics The TF importance metrics that will be used for the ranking, possible values include DC, CC, BC, DFC, DCC, DDC.
#' @return An updated sincera object, use getDF with name="tfranks" to access the lcc
#'
setGeneric("drivingforce.rankTFs", function(object, groups=NULL, metrics=c("DC","CC","BC","DFC","DCC","DDC"), ...) standardGeneric("drivingforce.rankTFs"))
#' @export
setMethod("drivingforce.rankTFs","sincera",
          function(object, groups=NULL, metrics=c("DC","CC","BC","DFC","DCC","DDC"), ...) {

              if (is.null(groups)) {
                groups <- sort(unique(getCellMeta(object, name="GROUP")))
              }

              for (i in 1:length(groups)) {
                g <- getDF(object, group=groups[i], name="lcc")
                im <- getNodesImportance(g, measurements=metrics, componentType="weak", pathType="all", directed=FALSE, normalized=TRUE)
                object <- setDF(object, group=groups[i], name="im", value=im)

                tfs <- getDF(object, group=groups[i], name="tfs")

                im$TF <- 0
                im$TF[which(rownames(im) %in% tfs)] <- 1
                im <- im[which(im$TF == 1), metrics]

                im.r <- as.data.frame(apply(im, 2, function(y) rank(-y, ties.method="min")))
                avg.r <- apply(im.r, 1, mean)
                im.r$AVG.RANK <- as.numeric(rank(avg.r, ties.method="min"))

                im.r <- im.r[order(im.r$AVG.RANK),]

                im.r <- cbind(TF=rownames(im.r), im.r)

                print(im.r)

                object <- setDF(object, group=groups[i], name="tfranks", value=im.r)
              }

            return(object)
          }
)



#' Inferring first order conditional dependencies between TFs and TGs
#' Adapted from G1DBN package
#'
#' @param data (matrix) the expression data
#' @param method (character) regression method: ls, turkey, huber
#' @param prePosition (numeric) indices of regulators
#' @param targetPosition (numeric) indices of regulatory targets
#' @param lag (numeric) time lag; 0 for non-time series data
#' @return a list of score matrices
network.inference <-function(data, method='ls',predPosition=NULL, targetPosition=NULL, lag=0) {
  
  ## ===============================================
  ## INITIALIZING
  ## _______________________________________________
  
  data<-as.matrix(data)
  n=dim(data)[1] # nb of time points
  p=dim(data)[2] # nb of genes
  
  ## If predictor or target genes position is specified,
  ## we just work of them
  if(is.null(predPosition)){
    pred <- 1:p}
  else {
    pred <- predPosition}
  d <- length(pred)
  
  if(is.null(targetPosition)){
    tar <- 1:p}
  else {
    tar <- targetPosition}
  r <- length(tar)
  
  ## The score matrix
  S1ls=NULL
  S1tukey=NULL
  S1huber=NULL
  
  if('ls' %in% method){
    S1ls<-matrix(0,r,d)
  }
  if('tukey' %in% method){
    S1tukey<-matrix(0,r,d)
  }
  if('huber' %in% method){
    S1huber<-matrix(0,r,d)
  }
  
  ## ===============================================
  ## BUILDING THE SCORE MATRICES
  ## _______________________________________________
  
  ## Print the total number of vertices
  
  #cat("Treating", r ,"vertices:\n")
  cpt=10
  
  for (i in 1:r){
    ## Print percentage of accomplished vertices
    if ( ((i/r)*100)>=cpt ) {
      #cat(cpt, "% ",sep="")
      cpt=cpt+10
    }
    
    
    ## The regression model is
    ## Y = X
    
    ## Y is the vector containing the target genes
    ## on which the regression will be performed
    ## The time point 1 is removed
    if (1+lag >= n) {
      stop();
    }
    y<-data[(1+lag):n,tar[i]]
    
    for (j in  1:(d-1)){
      
      ## for all k > j
      for (k in c(1:d)[-c(1:j)]){
        
        ## X is a matrix with two columns containing the
        ## predicted gene on which the regression will
        ## be performed. The two columns contain
        ## respectively the data corresponding to the jth
        ## and the kth tested gene.
        ## The time n is removed
        x<-data[1:(n-lag),c(pred[j],pred[k])]
        
        ## =====================================
        ## ESTIMATION...
        ##
        ## Three estimators are available :
        ## - Least square
        ## - Huber
        ## - Tuckey
        
        ## =====================================
        ## LEAST SQUARE ESTIMATOR
        if('ls' %in% method){
          lm.3<-lm(y~x)
          prob<-abs(summary(lm.3)$coef[,"Pr(>|t|)"])
          ## coefficient aij(k) : aij given k
          S1ls[i,j]<-max(prob[2],S1ls[i,j],na.rm=TRUE)
          ## coefficient aik(j) : aik given j
          S1ls[i,k]<-max(prob[3],S1ls[i,k],na.rm=TRUE)
        }
        
        ## =====================================
        ## TUKEY'S ESTIMATOR
        if('tukey' %in% method){
          bisq.3<-rlm(y~x,method='MM')
          prob<-2*(1-pt(abs(summary(bisq.3)$coef[,"t value"]),n-4))
          ## coefficient aij(k) : aij given k
          S1tukey[i,j]<-max(prob[2],S1tukey[i,j],na.rm=TRUE)
          ## coefficient aik(j) : aik given j
          S1tukey[i,k]<-max(prob[3],S1tukey[i,k],na.rm=TRUE)
        }
        
        ## =====================================
        ## HUBER'S ESTIMATOR
        if('huber' %in% method ){
          hub.3<-rlm(y~x)
          prob<-2*(1-pt(abs(summary(hub.3)$coef[,"t value"]),n-4))
          ## coefficient aij(k) : aij given k
          S1huber[i,j]<-max(prob[2],S1huber[i,j],na.rm=TRUE)
          ## coefficient aik(j) : aik given j
          S1huber[i,k]<-max(prob[3],S1huber[i,k],na.rm=TRUE)
        }
        
        
      } # end k
    } # end j
  } # end i
  #cat("\n")
  
  ## The score matrices are return
  list(S1ls=S1ls,S1tukey=S1tukey,S1huber=S1huber)
}

#' Calculate the importance of a node in a graph
#
#' @param g an igraph object
#' @param measurement DFC, DCC, DDC, DC, CC, BC
#' @param componentType character string, either weak or strong; for directed graphs, weak implies weakly, strong strongly connected components to search.
#' @param pathType defined the types of the paths used for measuring the distance in directed graphs.
#' @param directed whether the graph is directed or undirected;
#' @param normalized whether to rescale the results to [0,1]
#' @return a importance matrix

getNodesImportance <- function(
  graph, nodes=NULL, measurements=c("DFC","DCC","DDC","DC","CC","BC"),
  componentType="weak",
  pathType="all",
  directed=FALSE,
  normalized=TRUE)
{
  cat("\tCalculating importance measurements: ")
  if (length(nodes)==0) {
    nodes=V(graph)
  }
  im <- as.data.frame(matrix(0, length(nodes), length(measurements)))
  rownames(im) <- nodes$name
  colnames(im) <- measurements
  
  for (m in measurements) {
    if (m == "DFC") {
      for (i in nodes) {
        rg <- delete.vertices(graph, i)
        n <- vcount(rg)
        cls <- clusters(rg, mode=componentType) # mode: strong, weak
        im[i,m]<- cls$no/n
      }
      cat("DFC ")
    } else if (m == "DCC") {
      for (i in nodes) {
        rg <- delete.vertices(graph, i)
        n <- vcount(rg)
        cls <- clusters(rg, mode=componentType) # mode: strong, weak
        s <- 0
        for (cs in cls$csize) {
          s <- s + cs*(cs-1)
        }
        im[i,m]<- 1 - s/(n*(n-1))
      }
      cat("DCC ")
    } else if (m == "DDC") {
      for (i in nodes) {
        rg <- delete.vertices(graph, i)
        n <- vcount(rg)
        spl <- shortest.paths(rg, mode=pathType) #mode: in, out, all
        s <- 0
        for ( si in 1:dim(spl)[1]) {
          for (sj in 1:dim(spl)[2]) {
            if (si == sj) {
            } else {
              if(spl[si,sj] == Inf) {
                s <- s + (1/n)
              } else {
                s <- s + (1/spl[si,sj])
              }
            }
          }
        }
        im[i,m]<- 1 - s/(n*(n-1))
      }
      cat("DDC ")
    } else if (m == "DC") {
      dgrs <- degree(graph,nodes,mode=pathType,loops=TRUE,normalized=normalized) # mode: in, out, all
      im[,m] <- dgrs
      cat("DC ")
    } else if (m == "CC") {
      clsns <- closeness(graph,nodes,mode=pathType,normalized=normalized) # mode: in, out, all
      im[,m] <- clsns
      cat("CC ")
    } else if (m == "BC") {
      btwns <- betweenness(graph, nodes, directed=directed, normalized=normalized)
      im[,m] <- btwns
      cat("BC ")
    } else {
      stop(cat("Unrecognized Importance Measurements:", m))
    }
  }
  cat("\n")
  return(im)
}


