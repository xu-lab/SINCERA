#' The sincera Class
#'
#' The object storing all information associated with the sincera analysis, including data, annotations, analyes, etc.
#' It is recommended to use the defined get and set methods to access and manipulate the data in the sincera object.
#'
#' Information is stored in slots. Key slots include:
#'
#' @section Slots:
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
#' }
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
                      ctv = "data.frame",                        # cell type validation results

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
  
  object@ctv <- data.frame(CELL=NULL, PVALUE=NULL, SCORE=NULL, CLUSTER.ASSIGNMENT=NULL, TYPE=NULL)

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
    ret <- Biobase::exprs(object@data)
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
#' @param scaled If TRUE, the update will be performed on the scaled expression matrix, object@sexprs; otherwise, it will be performed on Biobase::exprs(object@data)
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
    Biobase::exprs(object@data)[genes, cells] <- value
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

    if (length(value) != dim(Biobase::exprs(object@data))[1]) {
      stop("The size of \'value\' is inconsistent with the number of cells\n")
    }
    if (any(is.na(value))) {
      warning("The \'value\' contains missing values\n")
    }
    colnames(Biobase::exprs(object@data)) <- value
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
    ret <- dim(Biobase::exprs(object@data))[2]
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
  if (length(value) != dim(Biobase::exprs(object@data))[1]) {
    stop("The size of \'value\' is inconsistent with the number of genes in the object\n")
  }
  if (any(is.na(value))) {
    warning("The \'value\' contains missing values.")
  }
  rownames(Biobase::exprs(object@data)) <- value
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
  ret <- dim(Biobase::exprs(object@data))[1]
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
setCellMeta <- function(object, name, value, verbose=T) {
  if (length(value) != getCellNum(object)) {
    stop("The size of \'value\' is different from the number of cells in the object\n")
  }
  if (verbose) {
    if (!(name %in% colnames(pData(object@data)))) {
      cat("\nAdding cell meta data \"", name, "\"\n", sep="")
    } else {
      cat("\nUpdating cell meta data \"", name, "\"\n", sep="")
    }
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
copyCellMeta <- function(object, from, to, verbose=T) {
  if (!(from %in% colnames(pData(object@data)))) {
    stop("Undefined cell metadata \'", from, "\'\nHere are the defined cell meta data: ", paste(colnames(pData(object@data)), collapse = ", "),"\n")
  }
  object <- setCellMeta(object, name=to, value=getCellMeta(object, name=from))
  if (verbose) cat("\nCopied meta data \'", from, "\' to \'", to, "\'\n", sep="")
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
setGeneMeta <- function(object, name, value, verbose=T) {
  if (length(value) != getGeneNum(object)) {
    stop("The size of \'value\' is different from the number of cells in the object\n")
  }
  if (verbose) {
    if (!(name %in% colnames(fData(object@data)))) {
      cat("\nAdding gene meta data \"", name, "\"\n", sep="")
    } else {
      cat("\nUpdating gene meta data \"", name, "\"\n", sep="")
    }
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
#' @param do.reset If TRUE, reset the cell types of all cells to "undefined". All the other parameters will be ignored.
#' @param groups A set of cell groups
#' @param types The cell types of cell groups
#' @return A sincera object with the updated cell type and cell group mapping
#' @export
#'
setCellType <- function(object, do.reset=FALSE, groups=NULL, types=NULL) {
  
    if (do.reset) {
      
      if (FALSE) {
        groups <- sort(unique(getCellMeta(object, name="GROUP")))
        object@ctmap <- data.frame(GROUP=groups, TYPE=rep("undefined", length(groups)), stringsAsFactors = FALSE)
      }
      
      if (!("TYPE" %in% colnames(pData(object@data)))) {
        object <- setCellMeta(object, name="TYPE", value=rep("undefined", getCellNum(object)))
        warning("Cell type was not defined previously. All cells were assigned to \"undefined\" cell type")
      } else {
        object <- setCellMeta(object, name="TYPE", value=rep("undefined", getCellNum(object)))
        cat("All cells were reset to \"undefined\" cell type")
      }
        
    } else {
      if (length(groups) != length(types)) {
          stop("The numbers of items in the \"groups\" and \"types\" parameters are inconsistent")
      }
      if (!("TYPE" %in% colnames(pData(object@data)))) {
        object <- setCellMeta(object, name="TYPE", value=rep("undefined", getCellNum(object)))
        warning("Cell type was not defined previously. All cells were assigned to \"undefined\" cell type")
      }
      for (i in 1:length(groups)) {
        
          if (FALSE) {
            idx <- which(object@ctmap$GROUP == groups[i])
            if (length(idx)==0) {
                object@ctmap <- rbind(object@ctmap, data.frame(GROUP=groups[i], TYPE=types[i]))
            } else if (length(idx)==1) {
                object@ctmap$TYPE[idx] <- types[i]
            }
          }
          
          i.idx <- which(getCellMeta(object, name="GROUP") == groups[i])
          if (length(i.idx)==0) {
            stop("Cell group \"", groups[i], "\" was not defined.")
          } else {
            celltype <- as.character(getCellMeta(object, name="TYPE"))
            celltype[i.idx] <- types[i]
            object <- setCellMeta(object, name="TYPE", value=factor(celltype))
          }
      }
      cat("\nThe updated cell type definitions are as follows:\n")
      print(getCellType(object))
      cat("\n")
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
#' @details
#' This function takes "GROUP" metadata as cell group information, and returns mapping between the "GROUP" and "TYPE" metadata.
#' An error message will be returned if the "TYPE" metadata does not exist.
#' If any cell groups specified in the groups parameter are not defined in the "GROUP" metadata, an warning will be raised.
#' If any cell groups mapped to multiple cell types, an warning will be raised.
#' 
#'
getCellType <- function(object, groups=NULL) {
  
    if (FALSE) {
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
    }
    
    
    if (!("TYPE" %in% colnames(pData(object@data)))) {
      stop("Cell type information was not defined in the object. Please run setCellType() first.")
    }
    
    if (is.null(groups)) {
      groups2 <- sort(unique(as.character(getCellMeta(object, name="GROUP"))))
    } else {
      groups2 <- groups
    }
    
    j <- 1
    ret2 <- c()
    groups.multimapped <- c()
    groups.uniquemapped <- c()
    groups.notmapped <- c()
    celltype <- as.character(getCellMeta(object, name="TYPE"))
    for (i in 1:length(groups2)) {
      i.types <- unique(celltype[which(getCellMeta(object, "GROUP")==groups2[i])])
      if (length(i.types)<1) {
        groups.notmapped <- c(groups.notmapped, groups2[i])
      } else if (length(i.types)==1) {
        groups.uniquemapped <- c(groups.uniquemapped, groups2[i])
        ret2[j] <- i.types
        names(ret2)[j] <- groups2[i]
        j <- j+1
      } else {
        groups.multimapped <- c(groups.multimapped, groups2[i])
        ret2[j] <- paste(i.types, collapse=", ", sep="")
        names(ret2)[j] <- groups2[i]
        j <- j+1
      }
    }
    if (length(groups.notmapped)>0) {
      warning("The following groups are not found in the object", paste(groups.notmapped, collapse=", "), "\n") 
    }
    if (length(groups.multimapped)>0) {
      warning("The following groups mapped to multiple celltypes", paste(groups.multimapped, collpase=", "), "\n")
    }
    
    
    return(ret2)
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
#' @export
#' @details 
#' The value is a data frame with the first column containing the cell types and the second column containing the markers, e.g.,
#'  TYPE   SYMBOL
#'   AT2   SFTPB
#'   AT2   ABCA3
#'   AT2   SLC34A2
#'   AT2   LPCAT1
#'  Basal  KRT5
#'  Basal  KRT14
#'  
setCellTypeMarkers <- function(object, value) {

    if (is.data.frame(value) & dim(value)[2]>=2 & dim(value)[1]>0) { 
      value <- value[, 1:2]
      colnames(value) <- c("TYPE", "SYMBOL")
      object@markers <- value
    } else {
      stop("Invalid input format. Please make sure \"value\" is a data frame with the first column containing the cell types and the second column containing the markers.")
    }

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
setCellTypeValidation <- function(object, value) {
  
  if (!is.data.frame(value)) {
    stop("Invalid input value. Please refer to ?setCellTypeValidation for input format.")
  }
  object@ctv <- value
  return(object)
}

#' @export
getCellTypeValidation <- function(object) {
  ret <- NULL
  if (dim(object@ctv)[1]==0) {
    stop("No cell type validation data. Please run celltype.validation function first.")
  }
  ret <- object@ctv
  return(ret)
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
#' @param dp The expression matrix: rows are genes, columns are cells.
#' @param ident A character/factor vector encoding the identities (e.g., cluster, type, sample, condition, etc.) of cells in the expression matrix. The length of the vector should be the same as the number of cells in the expression matrix and the order of cell identities should be the same as the order of cells in the expression matrix.
#' @param groups The identity groups to be included in the calculation.
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



#' consensus maximization
#' @param x - normalized model prediction
consensus_maximization <- function(x, prior=NULL, alpha=4, epslon=0.1, max.iter=100) {
  
  dir.delim <- "/"
  
  cat("sincera: consensus maximization...\n")
  
  wd <- paste(getwd(), dir.delim, "sincera.consensus.maximization.", getTimestamp(), dir.delim, sep="")
  dir.create(wd)
  
  fs <- x
  
  # fs.type <- data.frame(NAME=colnames(fs), TYPE=model.type)
  
  A <- matrix(0, nrow=dim(fs)[1], ncol=dim(fs)[2]*2)
  rownames(A) = rownames(fs)
  colnames(A) = paste("q", seq(1,dim(A)[2],1), sep="")
  
  for (j in 0:(dim(fs)[2]-1)) {
    for (i in 1:dim(fs)[1]) {
      A[i,(2*j+1)] <- fs[i, (j+1)]
      A[i,(2*j+2)] <- 1-fs[i, (j+1)]
    }
  }
  
  if (is.null(prior)) {
    prior= rep(1, length=dim(A)[2])
    idx <- seq(2, dim(A)[2], by=2)
    prior[idx] <- 0
  }
  
  Y <- matrix(0, nrow=dim(A)[2], ncol=2)
  rownames(Y) <- colnames(A)
  colnames(Y) <- c("P","N")
  
  Y[,1] <- prior
  Y[,2] <- 1-Y[,1]
  
  
  if (FALSE) {
    type2.p <- 1
    type2.n <- 0.5
    
    for (j in 0:(dim(fs)[2]-1) ) {
      if (fs.type$TYPE[j+1] == 1) {
        Y[(2*j+1),] <- c(1,0)
        Y[(2*j+2),] <- c(0, 1)
      } else if (fs.type$TYPE[j+1] == 2) {
        Y[(2*j+1),] <- c(type2.p,1-type2.p)
        Y[(2*j+2),] <- c(type2.n,1-type2.n)
      }
    }
  }
  
  # start integration
  Q <- matrix(0, nrow=dim(Y)[1], ncol=dim(Y)[2])
  rownames(Q) <- rownames(Y)
  colnames(Q) <- colnames(Y)
  
  U <- matrix(0, nrow=dim(A)[1], ncol=dim(Y)[2])
  rownames(U) <- rownames(A)
  colnames(U) <- colnames(Y)
  
  # init U
  U[,1] <- fs[,2]
  U[,2] <- 1-U[,1]
  
  iter <- 1
  
  data4viz <- data.frame(Iteration=1:(max.iter-1), Delta=NA)
  
  delta <- Inf #
  while ( (delta>epslon) & (iter < max.iter)) {
    U.1 <- U
    # update Q
    dv <- matrix(0, nrow=dim(Q)[1], ncol=dim(Q)[1])
    diag(dv) <- apply(A, 2, sum)
    kv <- matrix(0, nrow=dim(Q)[1], ncol=dim(Q)[1])
    diag(kv) <- 1
    Q <- solve(dv + alpha*kv) %*% (t(A) %*% U + (alpha*kv) %*% Y)
    # update U
    dn <- matrix(0, nrow=dim(A)[1], ncol=dim(A)[1])
    diag(dn) <- apply(A, 1, sum)
    U <- solve(dn) %*% A %*% Q
    delta <- U-U.1
    delta <- norm(as.matrix(delta), "f")
    
    data4viz$Delta[iter] <- delta
    
    iter <- iter + 1
  }
  
  # visualization
  ggplot(data4viz[1:(iter-1),], aes(x=Iteration, y=Delta)) + geom_line(colour="blue") + geom_hline(yintercept=epslon, colour="red") + xlab("Iteration") + ggtitle("Consensus Maximization") + sincera_theme(font_size=12)
  ggsave(filename=paste(wd, "CM.pdf", sep=""))
  
  
  dU <- cbind(RID=rownames(fs), U)
  dU <- cbind(dU, x)
  dU <- dU[order(dU$P, decreasing=T),]
  
  write.table(dU, file=paste(wd, "CM.U-alpha", alpha,"-epslon", epslon, "-max.iter", max.iter, ".txt", sep=""), sep="\t", col.names=T, row.names=F)
  
  cat("sincera: consensus maximization completed\n\n")
}


###########################################################################################################
####################                                                                   ####################
####################                     auxilliary functions                          ####################
####################                                                                   ####################
###########################################################################################################

# return Sys.time in specified format
getTimestamp <- function(format.str='%m%d%Y_%Hh%Mm%Ss'){
  return(format(Sys.time(),format.str))
}


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
