#' The sincera Class
#'
#' The object storing all information associated with the sincera analysis, including data, annotations, analyes, etc.
#' It is recommended to use the defined get and set methods to access and manipulate the data in the sincera object.
#'
#' Information is stored in slots. Key slots include:
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{data}:}{\code{"ExpressionSet"}, the scRNA-seq data in ExpressionSet format, encoding the expression matrix, gene, and cell metadata, including sample and condition information. In most of the SINCERA functions, cell grouping is based on the "GROUP" metadata.  }
#'    \item{\code{genes.forclustering}:}{\code{"vector"}, the set of genes used for the cell cluster identification  }
#'    \item{\code{difftests}:}{\code{"list"}, results of differential expression test of all genes in the data using different methods}
#'    \item{\code{siggenes}:}{\code{"data.frame"}, a data frame containing the significant differentially expressed genes   }
#'    \item{\code{cte}:}{\code{"list"}, results of cell type enrichment  }
#'    \item{\code{ctmap}:}{\code{"data.frame"}, a data frame containing the mapping between cell clusters and cell types  }
#'    \item{\code{markers}:}{\code{"data.frame"}, a data frame containing known cell type markers }
#'    \item{\code{siggenes}:}{\code{"data.frame"}, a data frame containing the signature genes based on signature prediction   }
#'    \item{\code{sigvalidate}:}{\code{"data.frame"}, a data frame containing the results of signature validation   }
#'    \item{\code{df}:}{\code{"list"}, a list containing the results of driving force analysis   }
#'    \item{\code{associations}:}{\code{"data.frame"}, a data frame containing the gene and cell type association table   }
#'    \item{\code{tfs}:}{\code{"vector"}, a vector containing a set of all transcription factors/cofactors   }
#'    \item{\code{projname}:}{\code{"character"}, a string describing the data and analysis  }
#'
#'}
#' @name sincera
#' @rdname sincera
#' @aliases sincera-class
#' @exportClass sincera
sincera <- setClass(Class="sincera", slots =
                    c(rexprs="data.frame", rsampleinfo="vector", # raw data
                      data="ANY", sexprs="matrix",               # data used for downstream analysis

                      genes.forclustering="character",           # genes for clustering
                      hc.obj="ANY",                              # An object of class hclust
                      hc.k="numeric",                            # Number of clusters obtained from hierarchical clustering

                      pca="list",                                # PCA results
                      tsne="list",                               # tSNE results

                      # diff test
                      difftests="list",                          # Differential expression test results of all genes
                      diffgenes="data.frame",                    # Differentially expressed genes

                      cte = "list",                              # cell type enrichment results

                      ctmap = "data.frame",                      # mapping between cell types and clusters; column 1: cluster id; column 2: cell type
                      markers="data.frame",                      # marker information for cell types; column 1: cell type; column 2: gene

                      siggenes="data.frame",                     # signature prediction results

                      sigvalidate="data.frame",                  # signature validation results

                      df="list",                                 # driving force results

                      processing.log="list",                     # execution log

                      # knowledge base
                      associations="ANY",
                      tfs="ANY",

                      projname="character")
                    )



#' Construct a sincera object
#'
#' Construct a sincera object from data in memory or files on disk
#'
#' @param exprmatrix The expression matrix; rows are genes, columns are cells
#' @param exprfile   The path to the file containing the expression matrix; either exprmatrix or exprfile must be set; if both provided, exprmatrix will be used
#' @param samplevector A vector containing the batch information of cells
#' @param samplefile The path to the file containing the batch information of cells; if both samplevector and samplefile are set to NULL, all cells will be assigned as "sample1"
#' @param projname A string describing the data and analysis
#' @return The created sincera object
#' @export
#'
construct <- function(exprmatrix=NULL, exprfile=NULL,
                      samplevector=NULL, samplefile=NULL,
                      projname=NULL) {

  if (is.null(exprmatrix) & is.null(exprfile)) {
    stop("Please provide either the expression matrix or the path to the expression file")
  }

  if (!is.null(exprmatrix)) {
	  rexprs <- exprmatrix
  } else {
	  rexprs <- read.table(file=exprfile, sep="\t", header=T, row.names=1, check.names=FALSE)
  }

  if (!is.null(samplevector)) {
	  rsampleinfo <- samplevector
  } else {
	  if (!is.null(samplefile)) {
		  rsampleinfo <- as.character(read.table(file=samplefile, sep="\t", header=F, check.names=FALSE)$V1)
	  } else {
		  rsampleinfo <- rep("sample1", dim(rexpr)[2])
	  }
  }

  if (!(dim(rexprs)[2]) == length(rsampleinfo)) {
    stop("Inconsistent number of cells and sample info")
  }

  genenames <- rownames(rexprs)
  cellnames <- colnames(rexprs)

  if (is.null(projname)) {
    projname <- paste("sincera-", format(Sys.time(), "%b%d_%H_%M_%S"), sep="")
  }

  # create an sincera object and initialize using raw data
  object <- sincera(rexprs=rexprs, rsampleinfo=rsampleinfo, projname=projname)


  cells <- data.frame(CELL=cellnames, GROUP=rsampleinfo, SAMPLE=rsampleinfo, CLUSTER=rsampleinfo)
  rownames(cells) <- cellnames

  genes <- data.frame(SYMBOL=genenames)
  rownames(genes) <- genenames

  object@data <- new("ExpressionSet", exprs=as.matrix(rexprs), phenoData=new("AnnotatedDataFrame", data=cells), featureData=new("AnnotatedDataFrame", data=genes))

  object@sexprs <- matrix(NA, ncol=length(cellnames), nrow=length(genenames))
  rownames(object@sexprs) <- genenames
  colnames(object@sexprs) <- cellnames

  object@genes.forclustering <- unique(genenames)

  object@cte <- list()

  object@ctmap <- data.frame(GROUP=sort(unique(rsampleinfo)), TYPE=rep("undefined", length(unique(rsampleinfo))), stringsAsFactors = FALSE)
  object@markers <- data.frame(TYPE=NULL, SYMBOL=NULL, stringsAsFactors = FALSE)

  object@df <- list()

  object@associations <- data.frame()
  object@tfs <- c()

  return(object)
}


#' get expression matrix from a sincera object
#'
#' get expression values of a set of genes in a set of cells from a sincera object
#'
#' @param object A sincera object
#' @param scaled If TRUE, return the zscore scaled expression values
#' @param genes The set of genes to be returned; if NULL, set to all genes in the object
#' @param cells The set of cells to be returned; if NULL, set to all cells in the object
#' @return A matrix containing the expression values; rows are genes, columns are cells.
#' @export
#'
getExpression <- function(object, scaled=FALSE, genes=NULL, cells=NULL) {
  ret <- NULL
  if (scaled==TRUE) {
    ret <- object@sexprs
  } else {
    ret <- exprs(object@data)
  }

  if (is.null(genes)) genes <- rownames(ret)
  if (is.null(cells)) cells <- colnames(ret)

  # TODO: error handling: genes and cells

  ret <- ret[genes, cells]
  return(ret)
}


#' Update expression matrix in a sincera object
#'
#' update expression values of a set of genes in a set of cells in a sincera object
#'
#' @param object A sincera object
#' @param value An expression matrix, rows are genes, columns are cells, rownames are gene names, colnames are cell names
#' @param scaled If TRUE, the update will be performed on the scaled expression matrix, object@sexprs; otherwise, it will be performed on exprs(object@data)
#' @export
setExpression <- function(object, value, scaled=FALSE) {

  cells <- colnames(value)
  genes <- rownames(value)

  cells.notfound <- cells[which(!(cells %in% getCells(object)))]
  genes.notfound <- genes[which(!(genes %in% getGenes(object)))]
  if (length(cells.notfound)>0) {
    stop(paste(length(cells.notfound), "cells were not found:", paste(cells.notfound, collapse = ","), "\n"))
  }
  if (length(genes.notfound)>0) {
    stop(paste(length(genes.notfound), "genes were not found:", paste(genes.notfound, collapse = ","), "\n"))
  }

  if (scaled==TRUE) {
    if (is.null(object@sexprs)) {
      warning("The scaled expression matrix was undefined. An all-zero matrix has been created for the update.\n")
    }
    object@sexprs <- matrix(0, nrow=getGeneNum(object), ncol=getCellNum(object))
    rownames(object@sexprs) <- getGenes(object)
    colnames(object@sexprs) <- getCells(object)
    object@sexprs[genes, cells] <- value
  } else {
    exprs(object@data)[genes, cells] <- value
  }
  return(object)
}


#' get the names of cells in specified cell groups in a sincera object
#'
#' get the names of cells in specified cell groups in a sincera object
#'
#' @param object A sincera object
#' @param groups The cell groups; if NULL, set to all cell groups defined in the sincera object
#' @return A character vector containing cell names
#' @export
#'
getCells <- function(object, groups=NULL) {
    ret <- NULL
    allgroups <- sort(unique(getCellMeta(object, name="GROUP")))
    if (is.null(groups)) {
      groups <- allgroups
    } else {
      groups.notfound <- groups[which(!(groups %in% allgroups))]
      if (length(groups.notfound)>0) {
        stop(paste("The following groups were not found: ", paste(groups.notfound, collapse=", "), sep=""))
      }
    }

    ret <- as.character(rownames(pData(object@data))[which(getCellMeta(object, name="GROUP") %in% groups)])
    return(ret)
}


#' set the names of all cells in a sincera object
#'
#' set the names of all cells in a sincera object
#'
#' @param object A sincera object
#' @param value The names of all cells
#' @return An updated sincera object
#' @export
#'
setCellNames <- function(object, value) {

    if (length(value) != dim(exprs(object@data))[1]) {
      stop("The size of \'value\' is inconsistent with the number of cells\n")
    }
    if (any(is.na(value))) {
      warning("The \'value\' contains missing values\n")
    }
    colnames(exprs(object@data)) <- value
    rownames(pData(object@data)) <- value
    return(object)
}

#' Get the number of cells in a sincera object
#'
#' Get the number of cells in a sincera object
#'
#' @param object A sincera object
#' @return A numeric value representing the number of cells in the sincera object
#' @export
#'
getCellNum <- function(object) {
    ret <- 0
    ret <- dim(exprs(object@data))[2]
    return(ret)
}

#' Get the names of all genes in a sincera object
#'
#' @param object A sincera object
#' @return A character vector of the names of all genes in the sincera object
#' @export
#'
getGenes <- function(object) {
    ret <- NULL
    ret <- as.character(rownames(fData(object@data)))
    return(ret)
}

#' Update the gene names in a sincera object
#'
#' @param object A sincera object
#' @param value A character vector containing the new names of all genes
#' @return The updated sincera object
#' @export
#'
setGeneNames <- function(object, value) {
  if (length(value) != dim(exprs(object@data))[1]) {
    stop("The size of \'value\' is inconsistent with the number of genes in the object\n")
  }
  if (any(is.na(value))) {
    warning("The \'value\' contains missing values.")
  }
  rownames(exprs(object@data)) <- value
  rownames(fData(object@data)) <- value
  return(object)
}

#' Get the number of the genes in a sincera object
#'
#' @param object A sincera object
#' @return The number of genes
#' @export
#'
getGeneNum <- function(object) {
  ret <- 0
  ret <- dim(exprs(object@data))[1]
  return(ret)
}


#' Add or update meta data of cells in a sincera object
#'
#' Add or update meta data of cells in a sincera object
#'
#' @param object A sincera object
#' @param name The name of the meta data
#' @param value The vector containing the values of the meta data of each cell (should be of the same length as the number of cells in the object)
#' @return The sincera object with the updated cell meta data
#' @export
#'
setCellMeta <- function(object, name, value) {
  if (length(value) != getCellNum(object)) {
    stop("The size of \'value\' is different from the number of cells in the object\n")
  }
  if (!(name %in% colnames(pData(object@data)))) {
    cat("\nAdding cell meta data \'", name, "\'\n", sep="")
  } else {
    cat("\nUpdating cell meta data \'", name, "\'\n", sep="")
  }
  pData(object@data)[, name] <- value
  return(object)
}

#' Get the meta data of cells in a sincera object
#'
#' Get the meta data of cells in a sincera object
#'
#' @param object A sincera object
#' @param name The name of the meta data
#' @return A vector containing the value of the meta data of all cells
#' @export
#'
getCellMeta <- function(object, name) {
  ret <- NULL
  if (!(name %in% colnames(pData(object@data)))) {
    stop("Undefined cell metadata \'", name, "\'\nHere are the defined cell meta data: ", paste(colnames(pData(object@data)), collapse = ", "),"\n")
  } else {
    ret <- pData(object@data)[, name]
  }
  return(ret)
}


#' Add or update a cell meta data using a predefined cell meta data
#'
#' @param object A sincera object
#' @param from The name of the predefined cell meta data
#' @param to The name of the meta data to be added or updated
#' @return The update sincera object
#' @export
#'
copyCellMeta <- function(object, from, to) {
  if (!(from %in% colnames(pData(object@data)))) {
    stop("Undefined cell metadata \'", from, "\'\nHere are the defined cell meta data: ", paste(colnames(pData(object@data)), collapse = ", "),"\n")
  }
  object <- setCellMeta(object, name=to, value=getCellMeta(object, name=from))
  cat("\nCopied meta data \'", from, "\' to \'", to, "\'\n", sep="")
  return(object)
}


#' Add or update meta data of genes in a sincera object
#'
#' @param object A sincera object
#' @param name The name of the meta data
#' @param value A vector containing the value of the meta data of each gene (should be of the same length as the number of genes in the object)
#' @return The sincera object with the updated gene meta data
#' @export
#'
setGeneMeta <- function(object, name, value) {
  if (length(value) != getGeneNum(object)) {
    stop("The size of \'value\' is different from the number of cells in the object\n")
  }
  if (!(name %in% colnames(fData(object@data)))) {
    cat("\nAdding gene meta data \'", name, "\'\n", sep="")
  } else {
    cat("\nUpdating gene meta data \'", name, "\'\n", sep="")
  }
  fData(object@data)[, name] <- value
  return(object)
}

#' Get a gene meta data
#'
#' @param object A sincera object
#' @param name The name of the meta data
#' @return A vector encoding the value of the meta data of genes
#' @export
#'
getGeneMeta <- function(object, name) {
  ret <- NULL
  if (!(name %in% colnames(fData(object@data)))) {
    stop("Undefined gene metadata \'", name, "\'\nHere are the defined gene meta data: ", paste(colnames(fData(object@data)), collapse = ", "),"\n")
  } else {
    ret <- fData(object@data)[, name]
  }
  return(ret)
}


#' Add or update a gene meta data using a predefined gene meta data
#'
#' @param object A sincera object
#' @param from The name of the predefined gene meta data
#' @param to The name of the meta data to be added or updated
#' @return The update sincera object
#' @export
#'
copyGeneMeta <- function(object, from, to) {
  if (!(from %in% colnames(fData(object@data)))) {
    stop("Undefined gene metadata \'", from, "\'\nHere are the defined gene meta data: ", paste(colnames(fData(object@data)), collapse = ", "),"\n")
  }
  object <- setGeneMeta(object, name=to, value=getGeneMeta(object, name=from))
  return(object)
}



#' Map cell clusters to cell types
#'
#' Map cell clusters to cell types
#' @param object A sincera object
#' @param do.reset If TRUE, reset the cell types of all cells to "undefined"
#' @param groups A set of cell groups
#' @param types The cell types of cell groups
#' @return A sincera object with the updated cell type and cell group mapping
#' @export
#'
setCellType <- function(object, do.reset=FALSE, groups=NULL, types=NULL) {
    if (do.reset) {
        groups <- sort(unique(getCellMeta(object, name="GROUP")))
        object@ctmap <- data.frame(GROUP=groups, TYPE=rep("undefined", length(groups)), stringsAsFactors = FALSE)

    } else {
        if (length(groups) != length(types)) {
            stop("Inconsistent number of clusters and cell types.")
        }
        for (i in 1:length(groups)) {
            idx <- which(object@ctmap$GROUP == groups[i])
            if (length(idx)==0) {
                object@ctmap <- rbind(object@ctmap, data.frame(GROUP=groups[i], TYPE=types[i]))
            } else if (length(idx)==1) {
                object@ctmap$TYPE[idx] <- types[i]
            }
        }
    }
    return(object)
}


#' Get the cell types of given cell groups
#'
#' Get the cell types of given cell groups
#'
#' @param object A sincera object
#' @param groups The set of cell groups whose cell types will be returned; if NULL, set to all cell groups defined in the "GROUP" metadata
#' @return A character vector of cell types
#' @export
#'
getCellType <- function(object, groups=NULL) {
    if (is.null(groups)) {
      groups <- object@ctmap$GROUP
    }
    ret <- c()
    for (i in groups) {
        i.idx <- which(object@ctmap$GROUP == i)
        if (length(i.idx)>0) {
            ret <- c(ret, as.character(object@ctmap$TYPE[i.idx]))
        } else {
            ret <- c(ret, "NotFound")
        }
    }
    names(ret) <- groups
    return(ret)
}

#' Get cell type marker information defined in a sincera object
#'
#' Get cell type marker information defined in a sincera object
#'
#' @param object A sincera object
#' @return A data frame containing the markers and their corresponding cell types. The first column is cell type, and the second column is marker gene name.
#' @export
#'
getCellTypeMarkers <- function(object) {
    ret <- NULL
    ret <- object@markers
    return(ret)
}


#' Add cell type marker information into sincera
#'
#' Add cell type marker information into sincera
#'
#' @param object A sincera object
#' @param value A data frame containing the markers and their corresponding cell types. The first column is cell type, and the second column is marker gene name.
#' @param types
#'
#' @export
setCellTypeMarkers <- function(object, value) {

    value <- value[, 1:2]
    colnames(value) <- c("TYPE", "SYMBOL")
    object@markers <- value

    return(object)
}


#' Get the set of genes used for cell cluster identification
#'
#' Return the value of the genes.forclustering slot
#'
#' @param object A sincera object
#' @export
#'
getGenesForClustering <- function(object) {
  ret <- NULL
  ret <- object@genes.forclustering
  return(ret)
}


#' Update the set of genes used for cell cluster identification
#'
#' Update the value of the genes.forclustering slot
#'
#' @param object A sincera object
#' @param value A vector containing the names of genes for cell cluster identification
#' @export
#'
setGenesForClustering <- function(object, value) {
  object@genes.forclustering <- value
  return(object)
}

#' Update the PCA results in sincera
#'
#' update the pca slot
#'
#' @param object A sincera object
#' @param name The name of item to be updated, possible values include rds - principal components, loadings, sdev
#' @param value The new value of the item.
#' @return The updated sincera object
#' @export
#'
setPCA <- function(object, name, value) {
  if (!(name %in% c("rds","loadings","sdev"))) {
    stop("Invalid name. Please choose one from 'rds', 'loadings', 'sdev'.")
  }
  object@pca[[name]] <- value
  return(object)
}

#' Get the PCA results stored in a sincera object
#'
#' @param object A sicnera object
#' @param name The name of the item to be returned, possible values include rds, loadings, sdev
#' @return The value of the specified PCA results
#' @export
#'
getPCA <- function(object, name) {
  ret <- NULL
  if (name=="rds") {
  ret <- object@pca$rds
  } else if (name=="loadings") {
    ret <- object@pca$loadings
  } else if (name=="sdev") {
    ret <- object@pca$sdev
  } else {
    stop("Invalid name. Please choose one from 'rds', 'loadings', 'sdev'")
  }

  return(ret)
}


#' @export
getTSNE <- function(object, name) {
  ret <- NULL
  if (name=="rds") {
    ret <- object@tsne$rds
  } else {
    stop("Invalid name. Please choose one from 'rds', 'loadings', 'sdev'")
  }

  return(ret)
}

#' @export
setTSNE <- function(object, name, value) {
  if (!(name %in% c("rds"))) {
    stop("Invalid name. Please choose one from 'rds'.")
  }
  object@tsne[[name]] <- value
  return(object)
}


#' @export
setCellTypeEnrichment <- function(object, groups=NULL, value) {

    if (is.null(groups)) {
      groups <- sort(unique(getCellMeta(object, name="GROUP")))
    }

    if (length(groups) != length(value)) {
        stop("Inconsistent number of groups and vlaues!")
    }

    for (i in 1:length(groups)) {
        object@cte[[groups[i]]] <- value[[i]]
    }

    return(object)
}

#' @export
getCellTypeEnrichment <- function(object, groups=NULL) {
    ret <- NULL
    if (is.null(groups)) {
      groups <- sort(unique(getCellMeta(object, name="CLUSTER")))
    }
    if (length(object@cte)==0) {
        stop("No cell type enrichment data. Please run celltype.enrichment function first.")
    }
    cnames <- names(object@cte)
    ret <- object@cte[intersect(cnames, groups)]
    return(ret)
}


#' @export
getDiffTest <- function(object, method="welch") {
    ret <- NULL
    if (length(object@difftests) < 1 ) {
        stop("No differential expression results available. Please run cluster.diffgenes() first.")
    } else {
        if (!(method %in% names(object@difftests))) {
            stop("No differential expression results using method=", method, " available. Please run cluster.diffgenes with method=", method, " first", sep="")
        } else {
            ret <- object@difftests[[method]]
        }
    }
    return(ret)
}

#' @export
setDiffTest <- function(object, value, method) {
    object@difftests[[method]] <- value
    return(object)
}


#' @export
setDiffGenes <- function(object, value) {
    object@diffgenes <- value
    return(object)
}

#' @export
getDiffGenes <- function(object, groups=NULL, method="welch", use.fdr=FALSE, thresh=0.01, print.summary=T) {

    cellgroup <- getCellMeta(object, name="GROUP")
    if (is.null(groups)) {
        groups <- sort(unique(cellgroup))
    }
    ng <- length(groups)

    if (use.fdr) {
        diffgenes <- data.frame(GROUP=NULL, SYMBOL=NULL, PVALUE=NULL, FDR=NULL)
    } else {
        diffgenes <- data.frame(GROUP=NULL, SYMBOL=NULL, PVALUE=NULL)
    }

    results <- getDiffTest(object, method=method)

    for (g in groups) {
        g.name <- paste(method, ".", g, sep="")
        g.fdr.name <- paste(g.name, ".fdr", sep="")

        g.diff.idx <- c()
        if (use.fdr) {
            if (g.fdr.name %in% colnames(results)) {
                g.diff.idx <- which(results[, g.fdr.name] < thresh)
            }
        } else {
            if (g.name %in% colnames(results)) {
                g.diff.idx <- which(results[, g.name] < thresh)
            }
        }

        if (length(g.diff.idx)>0) {
            if (use.fdr) {
                g.diffgenes <- data.frame(GROUP=g, SYMBOL=results$SYMBOL[g.diff.idx], PVALUE=results[g.diff.idx, g.name], FDR=results[g.diff.idx, g.fdr.name])
            } else {
                g.diffgenes <- data.frame(GROUP=g, SYMBOL=results$SYMBOL[g.diff.idx], PVALUE=results[g.diff.idx, g.name])
            }

            diffgenes <- rbind(diffgenes, g.diffgenes)
        }
    }

    criterion.name <- "pvalue"
    if (use.fdr) criterion.name <- "fdr"

    if (print.summary) {
        cat("\nNumber of diff genes (", criterion.name, "<", thresh, ") per cluster\n")
        print(table(diffgenes$GROUP))
    }

    return(diffgenes)
}

#' @export
setSigGenes <- function(object, value) {
    object@siggenes <- value
    return(object)
}

#' @export
getSigGenes <- function(object, groups=NULL, print.summary=T) {
    ret <- NULL
    if (is.null(groups)) {
      groups <- sort(unique(getCellMeta(object, name="GROUP")))
    }
    ret <- object@siggenes
    ret <- ret[which(ret$GROUP %in% groups), ]

    if (print.summary) {
        cat("The number of signature genes per group:\n")
        print(table(ret$GROUP))
    }

    return(ret)
}

#' @export
setSigValidation <- function(object, value) {
  object@sigvalidate <- value
  return(object)
}

#' @export
getSigValidation <- function(object) {
  ret <- NULL
  ret <- object@sigvalidate
  return(NULL)
}

#' @export
getHC <- function(object) {
  ret <- NULL
  hc.obj <- object@hc.obj
  hc.k <- object@hc.k
  hc.cellorder <- hc.obj$labels[hc.obj$order]
  hc.clusters <- cutree(hc.obj, k=hc.k)
  names(hc.clusters) <- hc.obj$labels
  ret <- list(hc.obj=hc.obj, hc.k=hc.k, hc.cellorder=hc.cellorder, hc.clusters=hc.clusters)
  return(ret)
}

#' @export
setHC <- function(object, value) {
  object@hc.obj <- value$hc.obj
  object@hc.k <- value$hc.k
  return(object)
}

# convert S4 object to ExpressionSet object for SINCERA v1 functions
#' @export
getES <- function(object, genes=NULL, cells=NULL) {
  ret <- NULL

  if (is.null(genes)) genes <- rownames(fData(object@data))
  if (is.null(cells)) cells <- rownames(pData(object@data))

  # TODO: error handling

  ret <- object@data[genes, cells]
  return(ret)
}

#' @export
setES <- function(object, value) {
    object@data <- value
    return(object)
}

#' @export
getTFs <- function(object) {
  ret <- NULL
  ret <- object@tfs
  return(ret)
}

#' @export
setTFs <- function(object, value) {
  if (is.null(value) | length(value)==0) {
    stop("Invalid input\n")
  }
  if (any(is.na(value))) {
    stop("The input contains missing values\n")
  }
  object@tfs <- value
  cat("\n", length(value), "TFs were added.\n")
  return(object)
}

#' @export
getAssociationTable <- function(object) {
  ret <- NULL
  ret <- object@associations
  return(ret)
}

#' @export
setAssociationTable <- function(object, value) {
  object@associations <- value
  return(object)
}


# df name cluster tfs tgs edges metrics tfranking

# name=c("tfs", "tgs", "edges", "trn", "lcc", "metrics", "tfranking")
#' @export
setDF <- function(object, group, name, value) {

    if (length(object@df)==0) {
        object@df[[1]] <- list(group=group, tfs=NULL, tgs=NULL, edges=NULL, lcc=NULL, metrics=NULL, tfranking=NULL)
        object@df[[1]][[name]] <- value
        names(object@df)[1] <- group
    }
    cnames <- names(object@df)
    id <- which(cnames==group)

    if (length(id) >0) {

        object@df[[id]][[name]] <- value

    } else {
        object@df[[length(cnames)+1]] <- list(group=group, tfs=NULL, tgs=NULL, edges=NULL, lcc=NULL, metrics=NULL, tfranking=NULL)
        object@df[[length(cnames)+1]][[name]] <- value
        names(object@df)[length(cnames)+1] <- group
    }

    return(object)
}

#' @export
getDF <- function(object, group, name) {
    ret <- NULL

    if (length(object@df)==0) {
        stop("Please use selectTFs and selectTGs to set the tfs and tgs for cluster ", group, sep="")
    }

    cnames <- names(object@df)
    id <- which(cnames==group)

    if (length(id)>0) {
        ret <- (object@df[[id]])[[name]]
    } else {
        stop("Cluster", group, "not found")
    }

    return(ret)
}





#' Generate gene level global or per group statistics
#'
#' @param dp The expression matrix: rows are genes, columns are cells
#' @param ident The vector of cell identity
#' @param group The identity groups to be included in the calculation
#' @param genes The set of genes to be included in the calculation; if NULL, will be automatically set to all genes in the expression matrix
#' @param min.exp The threshold vaule of gene expression
#' @param min.gsize The minimum size of a group
#' @param min.avg The minimum of average expression values
#' @param stats A vector containing the statistics to be calculated for each gene, including
#' @return A data.frame containing the results
#' @export
GeneStats <- function(dp, ident, groups=NULL, genes=NULL, min.expression=1, min.gsize=2, min.avg=1,
                      stats=c("var", "specificity", "min", "max", "avg", "peak1", "peak2",
                              "welch", "ident.fc","ident.avg", "ident.min", "ident.max", "ident.cnt", "ident.pct")) {

    if (is.null(groups)) {
        groups <- sort(unique(as.character(ident)))
    }

    if (is.null(genes)) {
        genes <- as.character(rownames(dp))
    }

    ## statistics helpers
    t_helper <- function(x, idx1, idx2) {
        x <- as.numeric(x)
        p <- NA
        if (var(x)==0) {
            p<-1
        } else {
            p<-t.test(x[i.idx], x[j.idx], alternative="greater", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value
        }
        return(p)
    }
    avg_helper <- function(x, idx) {
        x <- as.numeric(x)
        avg <- mean(x[idx])
        return(avg)
    }
    specificity_helper <- function(x) {
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

    dp <- dp[genes, ]

    dd <- data.frame(gene=rownames(dp))
    rownames(dd) <- dd$gene

    if ("var" %in% stats) dd$var <- apply(dp, 1, var)
    if ("cv" %in% stats) dd$cv <- apply(dp, 1, function(x) ifelse(mean(x)==0, 0, sd(x)/mean(x)))
    if ("specificity" %in% stats) dd$specificity <- apply(dp, 1, specificity_helper)
    if ("avg" %in% stats) dd$avg <- apply(dp, 1, mean)
    if ("min" %in% stats) dd$min <- apply(dp, 1, min)
    if ("max" %in% stats) dd$max <- apply(dp, 1, max)
    if ("peak1" %in% stats) {
        dd$peak1.cluster <- NA
        dd$peak1.value <- NA
        for (i in 1:dim(dd)[1]) {
            i.order <- order(as.numeric(dp[i, ]), decreasing=T)
            i.peak1.idx <- i.order[1]
            dd$peak1.cluster[i] <- ident[i.peak1.idx]
            dd$peak1.value[i] <- dp[i, i.peak1.idx]
        }
    }
    if ("peak2" %in% stats) {
        dd$peak2.cluster <- NA
        dd$peak2.value <- NA
        for (i in 1:dim(dd)[1]) {
            i.order <- order(as.numeric(dp[i, ]), decreasing=T)
            i.peak2.idx <- i.order[2]
            dd$peak2.cluster[i] <- ident[i.peak2.idx]
            dd$peak2.value[i] <- dp[i, i.peak2.idx]
        }
    }

    for (i in 1:length(groups)) {

        i.idx <- which(ident == groups[i])
        j.idx <- which(ident != groups[i])

        if ("welch" %in% stats) {
            if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
                i.colname <- paste("welch.",groups[i],sep="")
                i.t <- apply(dp, 1, FUN=t_helper, idx1=i.idx, idx2=j.idx)
                dd[, i.colname] <- as.numeric(i.t)
            }
        }

        if ("ident.fc" %in% stats) {
            if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
                i.colname <- paste("ident.fc.",groups[i],sep="")
                i.avg <- apply(dp, 1, FUN=avg_helper, idx=i.idx)
                j.avg <- apply(dp, 1, FUN=avg_helper, idx=j.idx)
                i.avg[which(i.avg < min.avg)] <- min.avg
                j.avg[which(j.avg < min.avg)] <- min.avg
                i.fc <- i.avg/j.avg
                dd[, i.colname] <- as.numeric(i.fc)
            }
        }

        if ("ident.avg" %in% stats) {
            if (length(i.idx)>=min.gsize) {
                i.colname <- paste("ident.avg.",groups[i],sep="")
                i.avg <- apply(dp, 1, FUN=avg_helper, idx=i.idx)
                dd[, i.colname] <- as.numeric(i.avg)
            }
        }

        if ("ident.min" %in% stats) {
            if (length(i.idx)>=min.gsize) {
                i.colname <- paste("ident.min.",groups[i],sep="")
                i.avg <- apply(dp[, i.idx], 1, function(x) min(as.numeric(x)))
                dd[, i.colname] <- as.numeric(i.avg)
            }
        }

        if ("ident.max" %in% stats) {
            if (length(i.idx)>=min.gsize) {
                i.colname <- paste("ident.max.",groups[i],sep="")
                i.avg <- apply(dp[, i.idx], 1, function(x) max(as.numeric(x)))
                dd[, i.colname] <- as.numeric(i.avg)
            }
        }


        if ("ident.cnt" %in% stats) {
            if (length(i.idx) >= min.gsize) {
                i.colname <- paste("ident.cnt.",groups[i],sep="")
                i.cnt <- rowSums(dp[, i.idx] > min.expression)
                dd[, i.colname] <- as.numeric(i.cnt)
            }
        }

        if ("ident.pct" %in% stats) {
            if (length(i.idx) >= min.gsize) {
                i.colname <- paste("ident.pct.",groups[i],sep="")
                i.pct <- rowSums(dp[, i.idx] > min.expression)/length(i.idx)
                dd[, i.colname] <- as.numeric(i.pct)
            }
        }

    }

    ret <- data.frame(gene=dd$gene)
    rownames(ret) <- ret$gene
    for (s in stats) {
        s.cols <- grep(paste("^",s, sep=""), colnames(dd), value = TRUE)
        ret[, s.cols] <- dd[, s.cols]
    }
    return(ret)
}




###########################################################################################################
####################                                                                   ####################
####################                     auxilliary functions                          ####################
####################                                                                   ####################
###########################################################################################################




pause <- function() {
    cat("Press [enter] to continue")
    line <- readline()
}


getPlotDims <- function(n) {
    nrow <- 1
    if (n>4) nrow <- 2
    if (n>8) nrow <- 3
    if (n>12) nrow <- 4
    if (n>16) nrow <- n
    ncol <- ceiling(n/nrow)
    return(list(nrow=nrow, ncol=ncol))
}


sincera_theme <- function(font_size=8) {

  theme(text=element_text(size=font_size)) +
    theme(plot.title = element_text(vjust=2)) +
    theme(axis.title.x = element_text(colour = "black", size=font_size), axis.title.y = element_text(colour = "black", vjust=2, size=font_size)) +
    theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) +
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_rect(colour="black", fill=NA)) + #axis.line = element_line(),
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(plot.margin = unit(c(10,10,10,10), "pt"))
}
