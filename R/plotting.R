
#' 2D plot of cells in the PCA or tSNE reduced dimensional spaces
#'
#' 2D plot of cells using the PCA or tSNE components
#'
#' @param object A sincera object
#' @param feature.type The type of components, possible values include pca, tsne
#' @param dim1 The component for x axis
#' @param dim2 The component for y axis
#' @param color.by The cell metadata used to color cells
#' @param pt.size The size of points (cells) in the plot
#' @param do.label If TRUE, show group label in the plot
#' @param label.size The size of the group label
#'
setGeneric("plotRDS", function(object, feature.type="pca", dim1=1, dim2=2, color.by="GROUP", pt.size=5, do.label=F, label.size=10, ...) standardGeneric("plotRDS"))
#' @export
setMethod("plotRDS","sincera",
          function(object, feature.type="pca", dim1=1, dim2=2, color.by="GROUP", pt.size=5, do.label=F, label.size=10, ...) {

              if (feature.type=="pca") {

                  dim1 <- paste("PC", dim1, sep="")
                  dim2 <- paste("PC", dim2, sep="")

                  dims <- object@pca$rds

                  if (!(dim1 %in% colnames(dims))) stop("Dimension ", dim1, " cannot be found in the PCA space.")
                  if (!(dim2 %in% colnames(dims))) stop("Dimension ", dim2, " cannot be found in the PCA space.")

                  viz <- data.frame(getPCA(object, name="rds"), Group=factor(getCellMeta(object, name=color.by)), stringsAsFactors = FALSE, check.names=FALSE)

              } else if (feature.type=="tsne") {

                  dim1 <- paste("tSNE", dim1, sep="")
                  dim2 <- paste("tSNE", dim2, sep="")

                  dims <- object@tsne$rds

                  if (!(dim1 %in% colnames(dims))) stop("Dimension ", dim1, " cannot be found in the PCA space.")
                  if (!(dim2 %in% colnames(dims))) stop("Dimension ", dim2, " cannot be found in the PCA space.")

                  viz <- data.frame(getTSNE(object, name="rds"), Group=factor(getCellMeta(object, name=color.by)), stringsAsFactors = FALSE, check.names=FALSE)
              }
            
              colnames(viz)[1:2] <- c("x","y")
     
              g <- ggplot(viz, aes(x=x, y=y, col=Group))
              g <- g + geom_point(size=pt.size)
              if (do.label==T) {
                #viz %>% dplyr::group_by(Group) %>% summarize(x = median(x), y = median(y)) -> centers
                centers <- data.frame(x=tapply(viz$x, viz$Group, median), y=tapply(viz$y, viz$Group, median))
                centers$Group <- rownames(centers)
                g <- g + geom_point(data = centers, aes(x=x, y=y), size=0, alpha=0) + geom_text(data=centers, aes(label=Group), col="black", size = label.size)
              }
              g <- g + xlab(dim1) + ylab(dim2)
              g <- g + sincera_theme()
              print(g)
          }
)


#' Plot the standard deviation of principal components
#'
#' @param object A sincera object
#' @param num.pcs The number of principal components to be plotted
#'
setGeneric("plotPCASD", function(object, num.pcs=20, ...) standardGeneric("plotPCASD"))
#' @export
setMethod("plotPCASD","sincera",
          function(object, num.pcs=20, ...) {

            viz <- data.frame(Cell=1:length(object@pca$sdev), Sdev=object@pca$sdev)
            if (num.pcs>dim(viz)[1]) num.pcs=dim(viz)[1]
            viz <- viz[1:num.pcs, ]

            g <- ggplot(viz, aes(x=Cell, y=Sdev)) + geom_point() + geom_line(col="grey")
            g <- g + ggtitle("Standard deviation of principal components") + xlab("Principal components") + ylab("Standard Deviation")
            g <- g + theme(axis.text.x = element_blank())
            g <- g + theme_bw()
            g <- g + theme(panel.grid=element_blank())
            print(g)

          }
)


#' Barplot of the cell type enrichment results
#'
#' Barplot of the cell type enrichment results
#'
#' @param object A sincera object
#' @param groups The groups whose cell type enrichment results will be plotted
#' @param top.k The number of most enriched cell types to be plotted for each cell group 
#'
setGeneric("plotCellTypeEnrichment", function(object, groups=NULL, top.k=5, ...) standardGeneric("plotCellTypeEnrichment"))
#' @export
setMethod("plotCellTypeEnrichment","sincera",
          function(object, groups=NULL, top.k=5, ...) {
            
            if (length(object@cte)==0) stop("No cell type enrichment results are available. Please run celltype.enrichment() first.")
            
            if (is.null(groups)) {
              groups <- names(object@cte)
            } 
            
            groups.notfound <- groups[which(!(groups %in% names(object@cte)))]
            
            if (length(groups.notfound)>0) {
               if (length(groups.notfound)==length(groups)) {
                 stop("No cell groups have cell type enrichment results")
               } else {
                 warning("The cell type enrichment results for the following groups are not available: ", paste(groups.notfound, collapse=",", sep=""), "\n")
               }
            }
            groups <- setdiff(groups, groups.notfound)
            
            viz <- data.frame(Type=NULL, Pvalue=NULL, Group=NULL, Order=NULL)
            for (i in 1:length(groups)) {
              i.idx <- which(names(object@cte) == groups[i])
              if (length(i.idx)>0) {
                tmp <- object@cte[[i.idx]]
                rownames(tmp) <- NULL
                
                if (top.k>dim(tmp)[1]) top.k = dim(tmp)[1]
                
                i.viz <- tmp[1:top.k, c("TYPE","Fisher.PV")]
                
                
                i.viz$Group <- groups[i]
                i.viz$Order <- factor(top.k:1)
                
                i.viz$Pvalue <- -log(i.viz$Fisher.PV)

                viz <- rbind(viz, i.viz)
              }
            }
            
            # order bar per facet based pvalue
            viz$n <- as.numeric(factor(viz$Group))
            viz$TYPE <- as.character(viz$TYPE)
            viz <- ddply(viz, .(Group,TYPE), transform,
                         x=paste(c(rep(' ',n-1), TYPE), collapse=''))
            viz$x <- factor(viz$x, levels=viz[order(viz$Fisher.PV, decreasing=T),"x"])
            
            g <- ggplot(viz, aes(x=x, y=Pvalue, fill=Group)) + facet_wrap(~Group, scales="free") 
            g <- g + geom_bar(stat="identity") + coord_flip()
            g <- g + ggtitle("Cell type enrichment for cell groups") 
            g <- g + xlab("Cell type annotation") + ylab("-log(Pvalue)")
            g <- g + sincera_theme() 
            g <- g + theme(panel.border = element_blank()) + theme(axis.line = element_line())
            print(g)
            
          }
)


#' Plot heatmap
#'
#' Plot heatmap
#' Acknowledgment: Implementation was based on the DoHeatmap function in Seurat
#'
#' @param object A sincera object
#' @param genes The expression patterns of genes to be shown in the heatmap If NULL, set to all genes
#' @param cells The cells to be shown in the heatmap. If NULL, set to all cells
#' @param scaled If TRUE, use scaled expression
#' @param do.log2 IF TRUE, apply log2 transformatoin to the expression
#' @param order.by.group IF TRUE, order cells based on their GROUP information
#' @param minmax The minimum and maximum values for color scale
#' @param show.labcol If TRUE, show the column lables
#' @param show.labRow If TRUE, show the row labels
#'
setGeneric("plotHeatmap", function(object, genes=NULL, cells=NULL, scaled=T, do.log2=FALSE,
                                   order.by.group=T, minmax=c(-1,1),
                                   show.labCol=FALSE, show.labRow=FALSE, ...) standardGeneric("plotHeatmap"))
#' @export
setMethod("plotHeatmap","sincera",
          function(object, genes=NULL, cells=NULL, scaled=T, do.log2=FALSE,
                   order.by.group=T, minmax=c(-1,1),
                   show.labCol=FALSE, show.labRow=FALSE, ...) {

              cells.df <- data.frame(cell=getCells(object), cluster=factor(getCellMeta(object, name="GROUP")))
              rownames(cells.df) <- cells.df$cell

              expr <- getExpression(object, scaled=scaled)
              if (do.log2) {
                  expr <- log(expr+1, 2)
              }

              if (is.null(cells)) {
                cells <- getCells(object)
              }
              cells.notfound <- cells[which(!(cells %in% rownames(cells.df)))]
              if (length(cells.notfound)>0) {
                  stop(length(cells.notfound), " cells not found:", paste(cells.notfound, split=","))
              }
              cells.df <- cells.df[cells, ]

              col.lab <- colnames(expr)
              colsep.use=NULL
              if (order.by.group) {
                  cells.df <- cells.df[order(cells.df$cluster), ]
                  col.lab=rep("",length(cells))
                  col.lab[round(cumsum(table(cells.df$cluster))-table(cells.df$cluster)/2)+1]=levels(cells.df$cluster)
                  colsep.use=cumsum(table(cells.df$cluster))
              }
              expr <- expr[, rownames(cells.df)]


              if (is.null(genes)) {
                  genes <- getGenes(object)
              }
              genes.notfound <- genes[which(!(genes %in% rownames(expr)))]
              if (length(genes.notfound)>0) {
                  stop(length(genes.notfound), " genes not found: ", paste(genes.notfound, split=","))
              }
              expr <- expr[genes, ]

              expr[expr<minmax[1]] <- minmax[1]
              expr[expr>minmax[2]] <- minmax[2]

              row.lab <- rownames(expr)

              if (!show.labRow) row.lab=rep("",length(genes))
              if (!show.labCol) col.lab=rep("",length(cells))

              myPalette=
                  function (low = "white", high = c("green", "red"), mid = NULL,
                            k = 50)
                  {
                      low <- col2rgb(low)/255
                      high <- col2rgb(high)/255
                      if (is.null(mid)) {
                          r <- seq(low[1], high[1], len = k)
                          g <- seq(low[2], high[2], len = k)
                          b <- seq(low[3], high[3], len = k)
                      }
                      if (!is.null(mid)) {
                          k2 <- round(k/2)
                          mid <- col2rgb(mid)/255
                          r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1],
                                                                    len = k2))
                          g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2],
                                                                    len = k2))
                          b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3],
                                                                    len = k2))
                      }
                      rgb(r, g, b)
                  }

              pyCols=myPalette(low = "green",high = "red",mid = "black")

              heatmap.2(expr, Rowv=NA, Colv=NA, trace="none", col=pyCols, colsep = colsep.use, labRow = row.lab, labCol = col.lab, density.info = "none")

              #return(object)
          }
)


#' Plot the dendrogram tree of hierarchical clustering
#'
#' @param object A sincera object with the hc.obj slot stored the dendrogram
#' @param horiz If TRUE, plot the dendrogram tree with tips turned right
#' @param show.labels If TRUE, show the labels of leafs
#' @param do.radial IF TRUE, plot the dendrogram tree using radial layout
#' @details 
#' When do.radial is FALSE, the plotting uses plot(); one can use pdf(), tiff(), or other similar functions to save plots to files when needed
#' when do.radial is TRUE, the plotting uses ggplot2:ggplot(); one can use ggsave() to save plots to files.
#'
setGeneric("plotHC", function(object, horiz=FALSE, show.labels=FALSE, do.radial=FALSE, ...) standardGeneric("plotHC"))
#' @export
setMethod("plotHC","sincera",
          function(object, horiz=FALSE, show.labels=FALSE, do.radial=FALSE, ...) {

              library(dendextend)

              ret <- getHC(object)

              dend <- as.dendrogram(ret$hc.obj)

              clusters <- as.numeric(ret$hc.clusters[ret$hc.cellorder])

              dend <- color_branches(dend, k=length(unique(clusters)), groupLabels = factor(unique(clusters), levels=unique(clusters)) )#, clusters=clusters)

              if (do.radial) {
                ggd1 <- as.ggdend(dend)
                g <- ggplot(ggd1, horiz = horiz, labels=show.labels)
                g <- g + scale_y_reverse(expand = c(0.2, 0))
                g <- g + coord_polar(theta="x")
                print(g)
              } else {
                if (show.labels==FALSE) { labels(dend) <- rep("", length(clusters))}
                plot(dend, horiz=horiz)
              }

            #  detach("package:dendextend", unload=T, character.only = T)
          }
)




#' Plot gene expression patterns
#'
#' Plot gene expression patterns
#'
#' @param object A sincera object
#' @param genes The expression patterns of genes to be plotted
#' @param use.scaled If TRUE, use the scaled expression values for plotting
#' @param do.log2 If TRUE, apply log2 transformation to the values
#' @param do.order If TRUE, order cells based on their GROUP information
#' @param show.jitter If TRUE, in the violin and boxplot, plot individual cells
#' @param plots The type of plots to be shown, possible values include cell, violin, boxplot, pca, and tsne
#' @param dims.use A numeric vector containing the indices of the two dimensions for pca or tsne plot
#' @param pt.size The size of points in the plot
#' @param fontsize The size of font in the plot
#'
setGeneric("plotMarkers", function(object, genes, use.scaled=T, do.log2=FALSE, do.order=T, show.jitter=T, plots=c("cell", "boxplot", "pca","tsne"), dims.use=1:2, pt.size=1,  font.size=8, ...) standardGeneric("plotMarkers"))
#' @export
setMethod("plotMarkers","sincera",
          function(object, genes, use.scaled=T, do.log2=FALSE, do.order=T, show.jitter=T, plots=c("cell", "boxplot", "pca","tsne"), dims.use=1:2, pt.size=1, font.size=8,  ...) {

              genes <- unique(genes)
              genes.notfound <- genes[which(!(genes %in% getGenes(object)))]
              if (length(genes.notfound)>0) {
                if (length(genes.notfound) == length(genes)) {
                  stop("No expression profiles were found\n")
                } else {
                  warning(paste("The following expression profiles were not found: ", paste(genes.notfound, collapse = ", ", sep=""), sep=""))
                }
              }

              genes <- setdiff(genes, genes.notfound)

              n<-length(genes)
              if (n==0) {
                stop("No genes found")
              }
              pdim <- getPlotDims(n)
              nrow <- pdim$nrow
              ncol <- pdim$ncol

              expr <- getExpression(object, scaled=use.scaled, genes=genes)
              if (do.log2==TRUE) {
                expr <- log(expr+1,2)
              }

              if (n>1) {
                viz <- data.frame(Cell=getCells(object),
                                  Group=factor(getCellMeta(object, name="GROUP")),
                                  t(expr)
                )
              } else {
                viz <- data.frame(Cell=getCells(object),
                                  Group=factor(getCellMeta(object, name="GROUP"))
                )
                viz[, genes] <- as.numeric(expr)
              }
              viz <- melt(viz, id.vars=c("Cell","Group"))
              colnames(viz)[3:4] <- c("Gene","Expression")

              if ("cell" %in% plots) {

                  viz.1 <- viz
                  if (do.order) {
                      viz.1$CID <- paste(viz$Group, ".", viz$Cell, sep="")
                  }

                  g <- ggplot(viz.1, aes_string(x="CID", y="Expression")) + facet_wrap(~Gene, scale="free", nrow=nrow, ncol=ncol)
                  g <- g + geom_point(aes_string(col="Group"), size=pt.size)
                  g <- g + sincera_theme()
                  g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
                  g <- g + theme(panel.border = element_blank())
                  g <- g + theme(axis.line.y=element_line(), axis.line.x=element_line())
                  g <- g + theme(text=element_text(size=font.size))
                  g <- g + theme(axis.title.x=element_text(size=rel(1.2)), axis.title.y=element_text(size=rel(1.2)), strip.text.x=element_text(size=rel(1.2)))
                  print(g)
                  pause()
              }

              if ("violin" %in% plots) {
                  g <- ggplot(viz, aes_string(x="Group", y="Expression")) + facet_wrap(~Gene, scale="free", nrow=nrow, ncol=ncol)
                  if (show.jitter==TRUE) {
                    g <- g + geom_jitter(height=0, alpha=0.3, size=pt.size)
                  }
                  g <- g + geom_violin(aes_string(fill="Group"), alpha=0.7, scale="width",adjust=0.75)
                  g <- g + sincera_theme()
                  g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
                  g <- g + theme(panel.border = element_blank())
                  g <- g + theme(axis.line.y=element_line(), axis.line.x=element_line())
                  g <- g + theme(text=element_text(size=font.size))
                  g <- g + theme(axis.title.x=element_text(size=rel(1.2)), axis.title.y=element_text(size=rel(1.2)), strip.text.x=element_text(size=rel(1.2)))
                  print(g)
                  pause()
              }

              if ("boxplot" %in% plots) {
                g <- ggplot(viz, aes_string(x="Group", y="Expression")) + facet_wrap(~Gene, scale="free", nrow=nrow, ncol=ncol)
                if (show.jitter==TRUE) {
                  g <- g + geom_jitter(height=0, alpha=0.3, size=pt.size)
                }
                g <- g + geom_boxplot(aes_string(fill="Group"), alpha=0.7, outlier.shape = NA)
                g <- g + sincera_theme()
                g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
                g <- g + theme(panel.border = element_blank())
                g <- g + theme(axis.line.y=element_line(), axis.line.x=element_line())
                g <- g + theme(text=element_text(size=font.size))
                g <- g + theme(axis.title.x=element_text(size=rel(1.2)), axis.title.y=element_text(size=rel(1.2)), strip.text.x=element_text(size=rel(1.2)))
                print(g)
                pause()
              }
              
              if ("pca" %in% plots) {
                
                Dim1 <- Dim2 <- NULL
                
                Dim1=getPCA(object, name="rds")[,dims.use[1]]
                Dim2=getPCA(object, name="rds")[,dims.use[2]]
                
                
                if (!is.null(Dim1) & !is.null(Dim2)) {
                  
                  if (n>1) {
                    viz <- data.frame(Cell=getCells(object),
                                      Dim1=Dim1,
                                      Dim2=Dim2,
                                      Group=factor(getCellMeta(object, name="GROUP")),
                                      t(expr)
                    )
                  } else {
                    viz <- data.frame(Cell=getCells(object),
                                      Dim1=Dim1,
                                      Dim2=Dim2,
                                      Group=factor(getCellMeta(object, name="GROUP"))
                    )
                    viz[, genes] <- as.numeric(expr)
                  }
                  viz <- melt(viz, id.vars=c("Cell","Dim1","Dim2","Group"))
                  colnames(viz)[5:6] <- c("Gene","Expression")
                  
                  gs <- list()
                  for (i in 1:length(genes)) {
                    i.viz <- subset(viz, Gene==genes[i])
                    if (use.scaled==T) {
                      i.viz$Col <- i.viz$Expression
                    } else {
                      i.viz$Col <- factor(as.numeric(as.factor(cut(as.numeric(i.viz$Expression),breaks = 2))))
                    }
                    g <- ggplot(i.viz, aes(x=Dim1, y=Dim2, col=Expression)) #+ facet_wrap(~Gene+Time, ncol=ntime)
                    g <- g + geom_point(size=pt.size) 
                    g <- g + ggtitle(paste(genes[i]))
                    if(use.scaled==T) {
                      g <- g + scale_color_gradient2(low="grey80", mid="grey",  high="red", midpoint=0)
                    } else {
                      g <- g + scale_color_gradientn(colors=c("grey80", "red")) #, guide=FALSE)
                      #g <- g + guide_legend(guides=FALSE)
                    }
                    g <- g + xlab(paste("PC", dims.use[1], sep="")) + ylab(paste("PC", dims.use[2], sep=""))
                    g <- g + theme(legend.text = element_text(size=8), plot.title = element_text(size=14))
                    g <- g + sincera_theme()
                    #g <- g + guides(color=guide_legend(title="", ncol=2, override.aes=list(size=5)))
                    #g <- g + theme(legend.position = "bottom", legend.text=element_text(size=12))
                    gs[[i]] <- g
                    
                  }
                  grid.arrange(grobs=gs) 
                  pause()
                } else {
                  warning(paste("Invalid", rds.use, "dimensions for rds plot"))
                }
              }
              
              
              if ("tsne" %in% plots) {
                
                Dim1 <- Dim2 <- NULL
                
                Dim1=getTSNE(object, name="rds")[,dims.use[1]]
                Dim2=getTSNE(object, name="rds")[,dims.use[2]]
                
                
                if (!is.null(Dim1) & !is.null(Dim2)) {
                  
                  if (n>1) {
                    viz <- data.frame(Cell=getCells(object),
                                      Dim1=Dim1,
                                      Dim2=Dim2,
                                      Group=factor(getCellMeta(object, name="GROUP")),
                                      t(expr)
                    )
                  } else {
                    viz <- data.frame(Cell=getCells(object),
                                      Dim1=Dim1,
                                      Dim2=Dim2,
                                      Group=factor(getCellMeta(object, name="GROUP"))
                    )
                    viz[, genes] <- as.numeric(expr)
                  }
                  viz <- melt(viz, id.vars=c("Cell","Dim1","Dim2","Group"))
                  colnames(viz)[5:6] <- c("Gene","Expression")
                  
                  gs <- list()
                  for (i in 1:length(genes)) {
                    i.viz <- subset(viz, Gene==genes[i])
                    if (use.scaled==T) {
                      i.viz$Col <- i.viz$Expression
                    } else {
                      i.viz$Col <- factor(as.numeric(as.factor(cut(as.numeric(i.viz$Expression),breaks = 2))))
                    }
                    g <- ggplot(i.viz, aes(x=Dim1, y=Dim2, col=Expression)) #+ facet_wrap(~Gene+Time, ncol=ntime)
                    g <- g + geom_point(size=pt.size) 
                    g <- g + ggtitle(paste(genes[i]))
                    if(use.scaled==T) {
                      g <- g + scale_color_gradient2(low="grey80", mid="grey",  high="red", midpoint=0)
                    } else {
                      g <- g + scale_color_gradientn(colors=c("grey80", "red")) #, guide=FALSE)
                      #g <- g + guide_legend(guides=FALSE)
                    }
                    g <- g + xlab(paste("tSNE", dims.use[1], sep="")) + ylab(paste("tSNE", dims.use[2], sep=""))
                    g <- g + theme(legend.text = element_text(size=8), plot.title = element_text(size=14))
                    g <- g + sincera_theme()
                    #g <- g + guides(color=guide_legend(title="", ncol=2, override.aes=list(size=5)))
                    #g <- g + theme(legend.position = "bottom", legend.text=element_text(size=12))
                    gs[[i]] <- g
                    
                  }
                  grid.arrange(grobs=gs) 
                  pause()
                } else {
                  warning(paste("Invalid", rds.use, "dimensions for rds plot"))
                }
              }


          }
)








