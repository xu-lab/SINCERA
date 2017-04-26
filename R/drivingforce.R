

#' Select candidate TFs for network inference and TF ranking
#'
#' Please use the \"tflists\" parameter to provide pre-selected TFs
#' Otherwise, the function will select candidate TFs for each group using the following criteria:
#' * TFs differentially expressed in the cell group (pvalue < \"diff.thresh\")
#' * TFs expressed greater than a certain percentage of the group cells (>\"min.expression\" in at least \"pct.thresh\" percentage of the group cells)
#'
#'  @param object A sincera object
#'  @param groups A vector of cell groups; if NULL, use all defined cell groups
#'  @param tflists A list containing user-provided candidate TFs for each group specified in \"groups\"
#'  @param diff.method The method used to perform the differential expression analysis
#'  @param diff.thresh The pvalue for determining differentially expressed genes
#'  @param min.expression The threshold expression value for determining expressed genes
#'  @param pct.thresh The percentage threshold for determining genes commonly expressed in a cell group
#'  @return An updated sincera object, use getDF with name="tfs" to access the selected TFs
#'
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
#' Please use the \"tglists\" parameter to provide pre-selected tgs
#' Otherwise, the function will select candidate TGs for each group using the following criteria:
#' * Genes differentially expressed in the cell group (pvalue < \"diff.thresh\")
#'
#'  @param object A sincera object
#'  @param clusters A vector of cell groups; if NULL, use all defined cell groups
#'  @param tglists A list containing user-provided candidate TGs for each group specified in \"groups\"
#'  @param diff.method The method used to perform the differential expression analysis
#'  @param diff.thresh The pvalue for determining differentially expressed genes
#'  @return An updated sincera object, use getDF with name="tgs" to access the selected TGs
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
#' @param clusters The cell groups for LCC extraction
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
                  g <- graph.data.frame(edges, directed=FALSE)
                  g <- simplify(g)

                  tfs <- getDF(object, group=groups[i], name="tfs")
                  tf.idx <- which(names(V(g)) %in% tfs)
                  tg.idx <- setdiff(1:length(V(g)), tf.idx)
                  g <- set_vertex_attr(g, name="Type", index=tf.idx, value="TF")
                  g <- set_vertex_attr(g, name="Type", index=tg.idx, value="TG")
                  V(g)$color <- "grey"
                  V(g)$color[which(V(g)$Type == "TF")] <- "red"

                #  plot(g, layout=layout.fruchterman.reingold, vertex.label=NA)

                  object <- setDF(object, group=groups[i], name="trn", value=g)


                  cl <- clusters(g, mode="weak")
                  lcc <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))

                  cat("The LCC of ", groups[i], ": ", length(V(lcc)), " nodes, ", length(E(lcc)), " edges\n", sep="")

                  #plot(lcc, layout=layout.fruchterman.reingold, vertex.label=NA)

                  if (length(V(lcc))/(length(V(g))) < 0.8) {
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
#' @param clusters The cell groups for TF ranking
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

