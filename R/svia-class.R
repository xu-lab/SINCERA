svia <- setClass(
  "svia",
  slots = c(
    raw.data = "ANY",
    data = "ANY",
    scale.data = "ANY",
    cellmeta = "data.frame",
    ident="vector",
    tsne = "data.frame"
  )
)


setMethod(
  f = "show",
  signature = "svia",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      #"in project",
      #object@project.name,
      "\n",
      nrow(x = object@data),
      "genes across",
      ncol(x = object@data),
      "samples.\n"
    )
    invisible(x = NULL)
  }
)


getExpression <- function(object, type="norm") {
  if (type == "raw") {
    return(object@raw.data)
  } else if (type == "scaled") {
    return(object@scale.data)
  } else {
    return(object@data)
  }
}

getCellMeta <-  function(object, name) {
  metas <- colnames(object@cellmeta)
  if (name %in% metas) {
    return(as.character(object@cellmeta[, name]))
  } else {
    stop("Cell meta not found")
  }
}

getCellMetas <-  function(object, name) {
  metas <- colnames(object@cellmeta)
  return(metas)
}

getCellName <- function(object) {
  return(colnames(getExpression(object)))
}

getCellIdent <- function(object) {
  return(object@ident)
}

getGeneMeta <- function(object, name="symbol") {
  if (name=="symbol") {
    return(rownames(getExpression(object)))
  }
}

getTSNE <- function(object) {
  object@tsne
}

convert2Svia <- function(object, from="seurat1") {
  if (from=="seurat1") {
    obj.new <- svia(raw.data=as(object=as.matrix(object@raw.data), Class = "dgCMatrix"),
                    data=as(object=as.matrix(object@data), Class = "dgCMatrix"),
                    scale.data=object@scale.data,
                    cellmeta=object@data.info,
                    ident=as.character(object@ident),
                    tsne=data.frame(object@tsne.rot)
                    )
  } else if (from=="seurat2") {
    obj.new <- svia(raw.data=object@raw.data,
                    data=object@data,
                    scale.data=object@scale.data,
                    cellmeta=object@meta.data,
                    ident=as.character(object@ident),
                    tsne=data.frame(object@dr$tsne@cell.embeddings)
                    )
  }
  
  return(obj.new)
}

if (FALSE) {
  sdps.old <- sdps
  sdps.new <- list()
  for (i in 1:length(sdps)) {
    sdps.new[[i]] <- list()
    for (j in 1:length(sdps[[i]])) {
      sdps.new[[i]][[j]] <- convert2Svia(sdps[[i]][[j]])
    }
    names(sdps.new[[i]]) <- names(sdps[[i]])
  }
  names(sdps.new) <- names(sdps)
  saveRDS(sdps.new, file="sdps.new.rds")
}
