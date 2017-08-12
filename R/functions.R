# from seurat
extract_field=function(string,field=1,delim="_") {
  fields=as.numeric(unlist(strsplit(as.character(field),",")))
  if (length(fields)==1)  return(strsplit(string,delim, fixed=TRUE)[[1]][field])
  return(paste(strsplit(string,delim, fixed=TRUE)[[1]][fields],collapse = delim))
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

prefiltering.cells <- function(x, y, mt.genes, min.exp=0, 
               min.ngene=500, max.ngene=5000, 
               min.libsize=500, max.libsize=10000,
               min.mt=0, max.mt=0.1,
               plot.individual=T,
               draw.line=T,
               pt.size=0.5
               ) {
  
  ret <- data.frame(Cell=colnames(x), Group=y, nGene=0, nUMI=0, MT=0)
  rownames(ret) <- ret$Cell
  ret$nGene <- colSums(x>min.exp)
  ret$nUMI <- colSums(x)
  mt.genes <- mt.genes[which(mt.genes %in% rownames(x))]
  if (length(mt.genes)>0) {
    ret$MT <- colSums(x[mt.genes, ])/ret$nUMI
  }
  
  cat("\nNumber of expressed genes per cell\n")
  print(tapply(ret$nGene, ret$Group, summary))
  
  cat("\nLibrary size per cell\n")
  print(tapply(ret$nUMI, ret$Group, summary))
  
  ngroups <- length(unique(y))
  
  gs <- list() 
  i <- 1
  
  g <- ggplot(data=ret, aes(x=nGene, y=nUMI))
  g <- g + geom_point(size=pt.size, col="grey60")
  if (draw.line==T) {
    g <- g + geom_vline(xintercept=min.ngene, size=1, linetype="dashed", col="red")
    g <- g + geom_vline(xintercept=max.ngene, size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=min.libsize, size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=max.libsize, size=1, linetype="dashed", col="red")
  }
  g <- g + ggtitle("nGenes vs. nCounts")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  if (ngroups>1) {
    g <- ggplot(data=ret, aes(x=nGene, y=nUMI))
    g <- g + facet_wrap(~Group, scales="free")
    if (draw.line==T) {
      g <- g + geom_vline(xintercept=min.ngene, size=1, linetype="dashed", col="red")
      g <- g + geom_vline(xintercept=max.ngene, size=1, linetype="dashed", col="red")
      g <- g + geom_hline(yintercept=min.libsize, size=1, linetype="dashed", col="red")
      g <- g + geom_hline(yintercept=max.libsize, size=1, linetype="dashed", col="red")
    }
    g <- g + geom_point(size=pt.size, col="grey60")
    g <- g + ggtitle("nGenes vs. nCounts")
    g <- g + sincera_theme()
    gs[[i]] <- g
    i <- i+1
  }
  
  g <- ggplot(data=ret, aes(x=nGene, y=MT))
  g <- g + geom_point(size=pt.size, col="grey60")
  if (draw.line==T) {
    g <- g + geom_vline(xintercept=min.ngene, size=1, linetype="dashed", col="red")
    g <- g + geom_vline(xintercept=max.ngene, size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=min.mt, size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=max.mt, size=1, linetype="dashed", col="red")
  }
  g <- g + ggtitle("nGenes vs. pMT")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  if (ngroups>1) {
    g <- ggplot(data=ret, aes(x=nGene, y=MT))
    g <- g + facet_wrap(~Group, scales="free")
    if (draw.line==T) {
      g <- g + geom_vline(xintercept=min.ngene, size=1, linetype="dashed", col="red")
      g <- g + geom_vline(xintercept=max.ngene, size=1, linetype="dashed", col="red")
      g <- g + geom_hline(yintercept=min.mt, size=1, linetype="dashed", col="red")
      g <- g + geom_hline(yintercept=max.mt, size=1, linetype="dashed", col="red")
    }
    g <- g + geom_point(size=pt.size, col="grey60")
    g <- g + ggtitle("nGenes vs. pMT")
    g <- g + sincera_theme()
    gs[[i]] <- g
    i <- i+1
  }
  
  
  
  g <- ggplot(data=ret, aes(x=nGene))
  g <- g + geom_histogram(binwidth = 100, fill="grey",  col="black")
  if (draw.line==T) {
    g <- g + geom_vline(xintercept=min.ngene, size=1,  linetype="dashed", col="red")
    g <- g + geom_vline(xintercept=max.ngene, size=1,  linetype="dashed", col="red")
  }
  g <- g + ggtitle("nGenes")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  if (ngroups>1) {
    g <- ggplot(data=ret, aes(x=nGene))
    g <- g + facet_wrap(~Group, scales="free")
    g <- g + geom_histogram(binwidth = 100, fill="grey",  col="black")
    if (draw.line==T) {
      g <- g + geom_vline(xintercept=min.ngene, size=1,  linetype="dashed", col="red")
      g <- g + geom_vline(xintercept=max.ngene, size=1,  linetype="dashed", col="red")
    }
    g <- g + ggtitle("nGenes")
    g <- g + sincera_theme()
    gs[[i]] <- g
    i <- i+1
  }
  
  g <- ggplot(data=ret, aes(x=Group, y=nGene))
  if (draw.line==T) {
    g <- g + geom_hline(yintercept=min.ngene, size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=max.ngene, size=1, linetype="dashed", col="red")
  }
  g <- g + geom_jitter(height=0, size=pt.size, col="grey60") + geom_boxplot(outlier.shape = NA)
  g <- g + ggtitle("nGenes")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  g <- ggplot(data=ret, aes(x=Group, y=log(nGene+1,2)))
  if (draw.line==T) {
    g <- g + geom_hline(yintercept=log(min.ngene+1,2), size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=log(max.ngene+1,2), size=1, linetype="dashed", col="red")
  }
  g <- g + geom_jitter(height=0, size=pt.size, col="grey60") + geom_boxplot(outlier.shape = NA)
  g <- g + ggtitle("nGenes")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1

  g <- ggplot(data=ret, aes(x=nUMI))
  g <- g + geom_histogram(binwidth = 500, fill="grey",  col="black")
  if (draw.line==T) {
    g <- g + geom_vline(xintercept=min.libsize, size=1, linetype="dashed", col="red")
    g <- g + geom_vline(xintercept=max.libsize, size=1, linetype="dashed", col="red")
  }
  g <- g + ggtitle("nCounts")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  if (ngroups>1) {
    g <- ggplot(data=ret, aes(x=nUMI))
    g <- g + facet_wrap(~Group, scales="free")
    g <- g + geom_histogram(binwidth = 500, fill="grey",  col="black")
    if (draw.line==T) {
      g <- g + geom_vline(xintercept=min.libsize, size=1, linetype="dashed", col="red")
      g <- g + geom_vline(xintercept=max.libsize, size=1, linetype="dashed", col="red")
    }
    g <- g + ggtitle("nCounts")
    g <- g + sincera_theme()
    gs[[i]] <- g
    i <- i+1
  }
  
  g <- ggplot(data=ret, aes(x=Group, y=nUMI))
  if (draw.line==T) {
    g <- g + geom_hline(yintercept=min.libsize, size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=max.libsize, size=1, linetype="dashed", col="red")
  }
  g <- g + geom_jitter(height=0, size=pt.size, col="grey60") + geom_boxplot(outlier.shape = NA)
  g <- g + ggtitle("nCounts")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  g <- ggplot(data=ret, aes(x=Group, y=log(nUMI+1, 2)))
  if (draw.line==T) {
    g <- g + geom_hline(yintercept=log(min.libsize+1,2), size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=log(max.libsize+1,2), size=1, linetype="dashed", col="red")
  }
  g <- g + geom_jitter(height=0, size=pt.size, col="grey60") + geom_boxplot(outlier.shape = NA)
  g <- g + ggtitle("nCounts")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  g <- ggplot(data=ret, aes(x=MT))
  g <- g + geom_histogram(binwidth = 0.01, fill="grey",  col="black")
  if (draw.line==T) {
    g <- g + geom_vline(xintercept=min.mt, size=1, linetype="dashed", col="red")
    g <- g + geom_vline(xintercept=max.mt, size=1, linetype="dashed", col="red")
  }
  g <- g + ggtitle("pMT")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
  
  if (ngroups>1) {
    g <- ggplot(data=ret, aes(x=MT))
    g <- g + facet_wrap(~Group, scales="free")
    g <- g + geom_histogram(binwidth = 0.01, fill="grey",  col="black")
    if (draw.line==T) {
      g <- g + geom_vline(xintercept=min.mt, size=1, linetype="dashed", col="red")
      g <- g + geom_vline(xintercept=max.mt, size=1, linetype="dashed", col="red")
    }
    g <- g + ggtitle("pMT")
    g <- g + sincera_theme()
    gs[[i]] <- g
    i <- i+1
  }
  
  g <- ggplot(data=ret, aes(x=Group, y=MT))
  if (draw.line==T) {
    g <- g + geom_hline(yintercept=min.mt, size=1, linetype="dashed", col="red")
    g <- g + geom_hline(yintercept=max.mt, size=1, linetype="dashed", col="red")
  }
  g <- g + geom_jitter(height=0, size=pt.size, col="grey60") + geom_boxplot(outlier.shape = NA)
  g <- g + ggtitle("pMT")
  g <- g + sincera_theme()
  gs[[i]] <- g
  i <- i+1
 
  ret.selected <- subset(ret, nGene>=min.ngene & nGene<=max.ngene & nUMI>=min.libsize & nUMI <= max.libsize & MT >= min.mt & MT <=max.mt)
  ret$Prefiltering <- 0
  ret$Prefiltering[which(ret$Cell %in% ret.selected$Cell)] <- 1
  
  print("Cell prefiltering:")
  print(table(ret$Prefiltering))
  print(table(ret$Group, ret$Prefiltering))
  
  if (plot.individual==T) {
    for (g in gs) print(g)
  } else {
    grid.arrange(grobs=gs, ncol=4)
  }
  return(ret)
}


prefiltering.genes <- function(x, min.exp=0, min.cells=3, min.counts=1) {
  cat("\nSelecting genes with acceptable abundancy for downstream analysis\n")
  sgenes <- c()
  num.cells= rowSums(x > min.exp)            
  sgenes=names(num.cells[which(num.cells>=min.cells)])
  cat("\t", length(sgenes), " genes expressed (>", min.exp, " UMI) in at least ", min.cells, "/", dim(x)[2], " cells\n", sep="")
  genes.use = names(which(rowSums(x) >= min.counts))
  cat("\t", length(genes.use), " genes have at least ", min.counts, " UMIs across ", dim(x)[2], " cells\n", sep="")
  sgenes = intersect(genes.use, sgenes)
  cat("\tFinally, ", length(sgenes), "/", dim(x)[1], " genes passed the abundancy criteria\n", sep="")
  return(sgenes)
}

NormalizeCounts <- function(x, genes=NULL, sf=0, log.base=2, count.pseudo=1) {
  cat("\nNormalizing counts\n")
  dp <- x
  if (is.null(genes)) {
    genes <- as.character(rownames(dp))
  }  
  
  col.sums=colSums(dp[genes, ])
  if (sf<=0) {
    sf = median(col.sums)
  }
  norm_counts = as.data.frame(sf* scale(dp, center=FALSE, scale=col.sums), check.names=FALSE)
  cat("\tnormalized to library size of each cell\n")
  cat("\trescaled with sf=", sf, "\n")
  lognorm_counts <- as.data.frame(log(norm_counts + count.pseudo, log.base), check.names=FALSE)
  cat("\tadded ", count.pseudo, " pseudo count, and log", log.base, " transformed\n", sep="")
  return(list(norm=norm_counts, lognorm=lognorm_counts, libsize=col.sums, sf=sf))
}

enrichment <- function(params, use.scaled=T) {
  
  gene <- params$gene
  if ((length(gene)>0) & (length(gene)<max.genes+1)) {
    
    if (length(gene)==1) {
      gene <- c(gene, "Actb")
    }
    
    if (use.scaled) {
      gs <- GeneStats(dp=params$sdp@scale.data[gene, ], ident=params$sdp@ident, groups=NULL, 
                    min.exp=params$min.exp, stats=c("ident.avg","ident.cnt","ident.pct", "ident.recall"))
    } else {
      gs <- GeneStats(dp=params$sdp@data[gene, ], ident=params$sdp@ident, groups=NULL, 
                      min.exp=params$min.exp, stats=c("ident.avg","ident.cnt","ident.pct", "ident.recall"))
      
    }
    gene.stats <- GetSigs(gs, groups=as.character(unique(params$sdp@ident)), criteria=c("ident.avg","ident.cnt","ident.pct", "ident.recall"))
    
     gene.stats$ident.avg <- round(gene.stats$ident.avg, 2)
     gene.stats$ident.pct <- round(gene.stats$ident.pct, 2)
     gene.stats$ident.recall <- round(gene.stats$ident.recall, 2)
    
    colnames(gene.stats) <- c("Gene","Group","Cluster.AverageExp", "Cluster.Expressed", "Cluster.Frequency", "Cluster.Recall")
    
    gene.stats <- gene.stats[which(gene.stats$Gene %in% params$gene), ]
    
    gene.stats <- unique(gene.stats)
    
    gene.stats$Total.Expressed <- rowSums(params$sdp@data[as.character(gene.stats$Gene), ]>params$min.exp)
    gene.stats$Total.Frequency <- gene.stats$Total.Expressed/length(which(substr(as.character(params$sdp@ident), 1, 6) != "Contam"))
    
    
    tmp <- table(params$sdp@ident)
    typedist <- data.frame(Type=as.character(names(tmp)), Count=as.numeric(tmp), stringsAsFactors = FALSE)
    typedist$Percent <- round(100*typedist$Count/sum(typedist$Count),2)
    typedist[dim(typedist)[1]+1, ] <- c("Total", sum(typedist$Count), 100)
    
    ub <- 1000
    
    #score1 = recall/percent
    gene.stats$Score1 <- 0
    for (i in 1:dim(gene.stats)[1]) {
      if (gene.stats$Total.Expressed[i] == 0) {
        gene.stats$Score1[i] <- ifelse(gene.stats$Cluster.Recall[i]>0, ub, 0)
      } else {
        gene.stats$Score1[i] <- gene.stats$Cluster.Recall[i]/(as.numeric(typedist$Percent[which(typedist$Type == as.character(gene.stats$Group[i]))])/100)
      }
    }
    
    #score2 = frequency*prior
    gene.stats$Score2 <- gene.stats$Cluster.Frequency*gene.stats$Total.Frequency
    
    gene.stats$Score3 <- 1
    for (i in 1:dim(gene.stats)[1]) {
      if (gene.stats$Total.Expressed[i] > 0) {
        i.a <- gene.stats$Cluster.Expressed[i]
        i.b <- (as.numeric(typedist$Count[which(typedist$Type == as.character(gene.stats$Group[i]))])) - i.a
        i.c <- gene.stats$Total.Expressed[i]
        i.d <- length(which(substr(as.character(params$sdp@ident), 1, 6) != "Contam")) - i.c
        i.data <- matrix(c(i.a, i.b, i.c, i.d), nrow = 2) 
        gene.stats$Score3[i] <- fisher.test(i.data, alternative="greater")$p.value
      }
    }
    
     gene.stats$Total.Frequency <- round(gene.stats$Total.Frequency, 2)
     gene.stats$Score1 <- round(gene.stats$Score1, 2)
     gene.stats$Score2 <- round(gene.stats$Score2, 2)
     gene.stats$Score3 <- round(gene.stats$Score3, 3)
  }
  
  return(gene.stats)
  
}



