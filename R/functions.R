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



GS.all <- function(dp, ident, groups=NULL, genes=NULL, percluster=T, min.exp=1, min.gsize=2, min.avg=1,
                   stats=c("t", "bimod", "binom", "t.fdr", "bimod.fdr","binom.fdr",  "efsize", "fc","ident.avg", "ident.q0","ident.q25","ident.q50", "ident.q75","ident.q100", "ident.cnt", "ident.pct", "ident.recall")) {
  if (is.null(groups)) {
    groups <- sort(unique(ident))
  }
  
  if (is.null(genes)) {
    genes <- rownames(dp)
  }
  
  ret <- data.frame(gene=genes)
  rownames(ret) <- ret$gene
  for (i in 1:length(groups)) {
    ident.stats <- stats[which(substring(stats, 1, 6)=="ident.")]
    if (length(ident.stats)>0) {
      i.ret <- GS.one(dp=dp, ident=ident, cluster.1=groups[i], cluster.2=NULL, 
                    min.exp=min.exp, min.gsize=min.gsize, min.avg=min.avg, stats=ident.stats)
      ret <- data.frame(ret, i.ret[rownames(ret), -1], check.names=FALSE) 
    }
    pair.stats <- stats[which(substring(stats, 1, 6) != "ident.")]
    if (length(pair.stats)>0) {
      for (s in pair.stats) {
        if (percluster==T) {
          i.gs <- groups[-i]
          i.ret <- data.frame(gene=rownames(dp))
          rownames(i.ret) <- i.ret$gene
          for (j in 1:length(i.gs)) {
            i.j.ret <- GS.one(dp=dp, ident=ident, cluster.1=groups[i], cluster.2=i.gs[j], 
                            min.exp=min.exp, min.gsize=min.gsize, min.avg=min.avg, stats=s)
            i.ret[, colnames(i.j.ret)[2]] <- as.numeric(i.j.ret[, 2]) #, check.names=FALSE)
          }
          if (s %in% c("t","binom", "bimod")) {
            ret[, paste(s, ".", groups[i], sep="")] <- apply(i.ret[rownames(ret), -1], 1, max)
            ret <- data.frame(ret, i.ret[rownames(ret), -1], check.names=FALSE)
          } else if (s %in% c("efsize","fc")) {
            ret[, paste(s, ".", groups[i], sep="")] <- apply(i.ret[rownames(ret), -1], 1, min)
            ret <- data.frame(ret, i.ret[rownames(ret), -1], check.names=FALSE)
          }
        } else {
          i.ret <- GS.one(dp=dp, ident=ident, cluster.1=groups[i], cluster.2=NULL, 
                          min.exp=min.exp, min.gsize=min.gsize, min.avg=min.avg, stats=s)
          ret[, colnames(i.ret)[2]] <- as.numeric(i.ret[, 2])
        }
      }
    }
  }
  
  dd <- data.frame(gene=ret$gene)
  for (s in stats) {
    s.cols <- grep(paste("^", s, ".", sep=""), colnames(ret), value = TRUE)
    dd <- cbind(dd, ret[, s.cols])
  }
  return(dd)
  return(ret)
}


GS.one <- function(dp, ident, cluster.1, cluster.2=NULL, genes=NULL, min.exp=1, min.gsize=2, min.avg=1,
                      stats=c("t","bimod", "binom", "t.fdr", "bimod.fdr","binom.fdr", "sSeq",  "efsize", "fc","ident.avg", "ident.q0","ident.q25","ident.q50", "ident.q75","ident.q100", "ident.cnt", "ident.pct", "ident.recall")) {

	
  if (is.null(genes)) genes <- as.character(rownames(dp))    
  
  t_helper <- function(x, idx1, idx2) {
    x <- as.numeric(x)
    p <- NA
    if (var(x[c(idx1, idx2)])==0) {
      p<-1
    } else {
      tryCatch({
          p<-t.test(x[idx1], x[idx2], alternative="greater", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value
        }, warning=function(w) {
          p <- 1
        }, error=function(e) {
          p <- 1
        }, finally={
          
        }
      )
      
    }
    return(p)
  }
  
  bimod_helper <- function(x, idx1, idx2) {
    p<-NA
    x <- as.numeric(x)
    p <- diffLRT(x[idx1], x[idx2])
    return(p)
  }
  
  # the probability of seeing n1 or more cells in cluster cells for a gene given its probability in the comparing cells
  binom_helper <- function(x, idx1, idx2, min.exp=0) {
    x <- as.numeric(x)
    n1 <- sum(x[idx1]>min.exp)
    n2 <- max(sum(x[idx2]>min.exp), 1)
    p <- pbinom(n1, length(idx1), n2/length(idx2), lower.tail = FALSE) + dbinom(n1, length(idx1), n2/length(idx2))
    return(p)
  }
  
  efsize_helper <- function(x, idx1, idx2, min.exp=0) {
    x <- as.numeric(x)
    n1 <- sum(x[idx1]>min.exp)
    n2 <- max(sum(x[idx2]>min.exp), 1)
    efs <- n1*length(idx2)/(n2*length(idx1))
    return(efs)
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
  
  dp <- dp[genes, ]

  dd <- data.frame(gene=rownames(dp))
  rownames(dd) <- dd$gene
  
  if (FALSE) {
	  dd$nexpressed <- rowSums(dp>min.exp)
	  dd$var <- apply(dp, 1, var)
	  dd$specificity <- apply(dp, 1, specificity_helper)
	  dd$avg <- apply(dp, 1, mean)
	  dd$min <- apply(dp, 1, min)
	  dd$max <- apply(dp, 1, max)
	  	  
	  dd$peak1.cluster <- NA
	  dd$peak1.value <- NA
	  dd$peak2.cluster <- NA
	  dd$peak2.value <- NA
	  	  
	  for (i in 1:dim(dd)[1]) {
		i.order <- order(as.numeric(dp[i, ]), decreasing=T)
		i.peak1.idx <- i.order[1]
		i.peak2.idx <- i.order[2]
		dd$peak1.cluster[i] <- as.character(ident[i.peak1.idx])
		dd$peak2.cluster[i] <- as.character(ident[i.peak2.idx])
		dd$peak1.value[i] <- dp[i, i.peak1.idx]
		dd$peak2.value[i] <- dp[i, i.peak2.idx]
	  }
  }
  
  #for (i in 1:length(groups)) {
    
  i.idx <- which(ident == cluster.1)
	if (is.null(cluster.2)) {
		#cluster.2 <- "rest"
		j.idx <- which(ident != cluster.1)
	} else {
		j.idx <- which(ident == cluster.2)
	}
        
  if ("t" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("t.", cluster.1, sep="")
      } else {
        i.colname <- paste("t.", cluster.1, ".", cluster.2, sep="")
      }
      i.t <- apply(dp, 1, FUN=t_helper, idx1=i.idx, idx2=j.idx)
      dd[, i.colname] <- as.numeric(i.t)
    }
  }
  if ("t.fdr" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("t.fdr.", cluster.1, sep="")
      } else {
        i.colname <- paste("t.fdr.", cluster.1, ".", cluster.2, sep="")
      }
      i.t <- apply(dp, 1, FUN=t_helper, idx1=i.idx, idx2=j.idx)
      i.t <- as.numeric(i.t)
      i.t <- p.adjust(i.t, method="fdr")
      dd[, i.colname] <- i.t
    }
  }
  
  if ("bimod" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("bimod.", cluster.1, sep="")
      } else {
        i.colname <- paste("bimod.", cluster.1, ".", cluster.2, sep="")
      }
      i.t <- apply(dp, 1, FUN=bimod_helper, idx1=i.idx, idx2=j.idx)
      dd[, i.colname] <- as.numeric(i.t)
    }
  }
  if ("bimod.fdr" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("bimod.fdr.", cluster.1, sep="")
      } else {
        i.colname <- paste("bimod.fdr.", cluster.1, ".", cluster.2, sep="")
      }
      i.t <- apply(dp, 1, FUN=bimod_helper, idx1=i.idx, idx2=j.idx)
      i.t <- as.numeric(i.t)
      i.t <- p.adjust(i.t, method="fdr")
      dd[, i.colname] <- i.t
    }
  }
  
  if ("binom" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("binom.", cluster.1, sep="")
      } else {
        i.colname <- paste("binom.", cluster.1, ".", cluster.2, sep="")
      }
      i.t <- apply(dp, 1, FUN=binom_helper, idx1=i.idx, idx2=j.idx, min.exp=min.exp)
      dd[, i.colname] <- as.numeric(i.t)
    }
  }
  
  if ("binom.fdr" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("binom.fdr.", cluster.1, sep="")
      } else {
        i.colname <- paste("binom.fdr.", cluster.1, ".", cluster.2, sep="")
      }
      i.t <- apply(dp, 1, FUN=binom_helper, idx1=i.idx, idx2=j.idx, min.exp=min.exp)
      i.t <- as.numeric(i.t)
      i.t <- p.adjust(i.t, method="fdr")
      dd[, i.colname] <- i.t
    }
  }
	
 if ("sSeq" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("sSeq.", cluster.1, sep="")
	i.colname.2 <- paste("log2FC.", cluster.1, sep="")
      } else {
        i.colname <- paste("sSeq.", cluster.1, ".", cluster.2, sep="")
	i.colname.2 <- paste("log2FC.", cluster.1, ".", cluster.2, sep="")
      }
      #i.t <- apply(dp, 1, FUN=binom_helper, idx1=i.idx, idx2=j.idx, min.exp=min.exp)
      #i.t <- as.numeric(i.t)
      #i.t <- p.adjust(i.t, method="fdr")
      conds <- rep("A", dim(dp)[2])
      conds[j.idx] <- "B"
      rs.2 <- nbTestSH.2(dp, conds=conds)
      dd[, i.colname] <- rs.2$pval   
      dd[, i.colname.2] <- rs.2$rawLog2FoldChange
    }
 }
  
  if ("efsize" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("efsize.", cluster.1, sep="")
      } else {
        i.colname <- paste("efsize.", cluster.1, ".", cluster.2, sep="")
      }
      i.t <- apply(dp, 1, FUN=efsize_helper, idx1=i.idx, idx2=j.idx, min.exp=min.exp)
      dd[, i.colname] <- as.numeric(i.t)
    }
  }
  
  if ("fc" %in% stats) {
    if (length(i.idx)>=min.gsize & length(j.idx)>=min.gsize) {
      if (is.null(cluster.2)) {
        i.colname <- paste("fc.", cluster.1, sep="")
      } else {
        i.colname <- paste("fc.", cluster.1, ".", cluster.2, sep="")
      }
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
      i.colname <- paste("ident.avg.", cluster.1, sep="")
      i.avg <- apply(dp, 1, FUN=avg_helper, idx=i.idx)
      dd[, i.colname] <- as.numeric(i.avg)
    }
  }
  
  if ("ident.min" %in% stats) {
    if (length(i.idx)>=min.gsize) {
      i.colname <- paste("ident.min.", cluster.1, sep="")
      i.avg <- apply(dp[, i.idx], 1, function(x) min(as.numeric(x)))
      dd[, i.colname] <- as.numeric(i.avg)
    }
  }
  
  if ("ident.max" %in% stats) {
    if (length(i.idx)>=min.gsize) {
      i.colname <- paste("ident.max.", cluster.1, sep="")
      i.avg <- apply(dp[, i.idx], 1, function(x) max(as.numeric(x)))
      dd[, i.colname] <- as.numeric(i.avg)
    }
  }
  
  if ("ident.q0" %in% stats) {
    if (length(i.idx)>=min.gsize) {
      i.colname <- paste("ident.q0.", cluster.1, sep="")
      i.q <- apply(dp[, i.idx], 1, function(x) quantile(as.numeric(x), probs=c(0)))
      dd[, i.colname] <- as.numeric(i.q)
    }
  }
  
  if ("ident.q25" %in% stats) {
    if (length(i.idx)>=min.gsize) {
      i.colname <- paste("ident.q25.", cluster.1, sep="")
      i.q <- apply(dp[, i.idx], 1, function(x) quantile(as.numeric(x), probs=c(0.25)))
      dd[, i.colname] <- as.numeric(i.q)
    }
  }
  
  if ("ident.q50" %in% stats) {
    if (length(i.idx)>=min.gsize) {
      i.colname <- paste("ident.q50.", cluster.1, sep="")
      i.q <- apply(dp[, i.idx], 1, function(x) quantile(as.numeric(x), probs=c(0.50)))
      dd[, i.colname] <- as.numeric(i.q)
    }
  }
  
  if ("ident.q75" %in% stats) {
    if (length(i.idx)>=min.gsize) {
      i.colname <- paste("ident.q75.", cluster.1, sep="")
      i.q <- apply(dp[, i.idx], 1, function(x) quantile(as.numeric(x), probs=c(0.75)))
      dd[, i.colname] <- as.numeric(i.q)
    }
  }
  
  if ("ident.q100" %in% stats) {
    if (length(i.idx)>=min.gsize) {
      i.colname <- paste("ident.q100.", cluster.1, sep="")
      i.q <- apply(dp[, i.idx], 1, function(x) quantile(as.numeric(x), probs=c(1)))
      dd[, i.colname] <- as.numeric(i.q)
    }
  }

  if ("ident.cnt" %in% stats) {
    if (length(i.idx) >= min.gsize) {
      i.colname <- paste("ident.cnt.", cluster.1, sep="")
      i.cnt <- rowSums(dp[, i.idx] > min.exp)
      dd[, i.colname] <- as.numeric(i.cnt)
    }
  }
  
  if ("ident.pct" %in% stats) {
    if (length(i.idx) >= min.gsize) {
      i.colname <- paste("ident.pct.", cluster.1, sep="")
      i.pct <- rowSums(dp[, i.idx] > min.exp)/length(i.idx)
      dd[, i.colname] <- as.numeric(i.pct)
    }
  }
  
  if ("ident.recall" %in% stats) {
    if (length(i.idx) >= min.gsize) {
      i.colname <- paste("ident.recall.", cluster.1, sep="")
      i.recall <- rowSums(dp[, i.idx] > min.exp)/rowSums(dp>min.exp)
      dd[, i.colname] <- as.numeric(i.recall)
    }
  }
  #}
  
  if (FALSE) {
	  ret <- data.frame(gene=dd$gene)
	  for (s in stats) {
		s.cols <- grep(paste("^",s,".",sep=""), colnames(dd), value = TRUE)
		ret <- cbind(ret, dd[, s.cols])
	  }
  }
  return(dd)				  
					  
}


# op: 2 >, 1>=, 0==, -1 <=, -2 <
GetSigs <- function(gs, groups, criteria, thresh=NULL, op=NULL) {
  
  if(is.null(thresh)) thresh=rep(-Inf, length(criteria))
  if(is.null(op)) op=rep(1, length(criteria))
  
  getValid <- function(x, thresh, op) {
    ret <- c()
    if (op==0) {
      ret <- which(x==thrsh)
    } else if (op==1) {
      ret <- which(x>=thresh)
    } else if (op==2) {
      ret <- which(x>thresh)
    } else if (op==-1) {
      ret <- which(x<=thresh)
    } else if (op==-2) {
      ret <- which(x<thresh)
    }
    return(ret)
  }
  
  dd <- data.frame(gene=NULL, group=NULL)
  for (i in criteria) {
    dd[, i] <- NULL
  }
  
  for (i in 1:length(groups)) {
    i.g <- groups[i]
    i.criteria <- paste(criteria, ".", i.g, sep="")
    i.sig <- 1:dim(gs)[1]
    for (j in 1:length(i.criteria)) {
      #cat(i, ".", j,"\n")
      i.sig <- intersect(i.sig, getValid(as.numeric(gs[, i.criteria[j]]), thresh[j], op[j]))
    }
    if (length(i.sig)>0) {
      i.dd <- data.frame(gene=rownames(gs)[i.sig], group=i.g)
      i.dd[, criteria] <- gs[i.sig, i.criteria]
      dd <- rbind(dd, i.dd)
    }
  }
  
  return(dd)
  
}
                   
                   
# from Seurat for bimod test ###############
diffLRT = function(x,y,xmin=1) {
  lrtX=bimodLikData(x)
  lrtY=bimodLikData(y)
  lrtZ=bimodLikData(c(x,y))
  lrt_diff=2*(lrtX+lrtY-lrtZ)
  return(pchisq(lrt_diff,3,lower.tail = F))
}

bimodLikData=function(x,xmin=0) {
  x1=x[x<=xmin]
  x2=x[x>xmin]
  xal=minmax(length(x2)/length(x),min=1e-5,max=(1-1e-5))
  likA=length(x1)*log(1-xal)
  mysd=sd(x2)
  if(length(x2)<2) {
    mysd=1
  }
  likB=length(x2)*log(xal)+sum(dnorm(x2,mean(x2),mysd,log=TRUE))
  return(likA+likB)
}

minmax=function(data,min,max) {
  data2=data
  data2[data2>max]=max
  data2[data2<min]=min
  return(data2)
}

# A wrapper of MeanVarPlot function in Seurat to perform batch-aware highly variable gene selection
pc.genes.selection <- function(lognorm, sgenes, scells.run1601, scells.run1888, x.low.cutoff = 0.2, x.high.cutoff=5, y.cutoff=0.25) {
  
  sdp <- new("seurat", raw.data = lognorm[sgenes, scells.run1601])
  sdp <- Setup(sdp, min.cells = 0, min.genes = 0, is.expr=0, do.logNormalize = F, total.expr = 1e4, do.scale=TRUE, do.center=TRUE, project = db)
  sdp@scale.data <- t(scale(t(sdp@data), center=T, scale=T))
  sdp <- MeanVarPlot(sdp ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff, do.contour = F)
  vargenes.1601 <- sdp@var.genes
  
  sdp <- new("seurat", raw.data = lognorm[sgenes, scells.run1888])
  sdp <- Setup(sdp, min.cells = 0, min.genes = 0, is.expr=0, do.logNormalize = F, total.expr = 1e4, do.scale=TRUE, do.center=TRUE, project = db)
  sdp@scale.data <- t(scale(t(sdp@data), center=T, scale=T))
  sdp <- MeanVarPlot(sdp ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff, do.contour = F)
  vargenes.1888 <- sdp@var.genes
  
  
  sdp <- new("seurat", raw.data = lognorm[sgenes, c(scells.run1601, scells.run1888)])
  sdp <- Setup(sdp, min.cells = 0, min.genes = 0, is.expr=0, do.logNormalize = F, total.expr = 1e4, do.scale=TRUE, do.center=TRUE, project = db)
  sdp@scale.data <- t(scale(t(sdp@data), center=T, scale=T))
  #sdp@raw.data <- dp[sgenes, scells]
  #sdp@ident <- factor(sruns)
  #names(sdp@ident) <- scells
  #sdp@data.info <- data.frame(nGene=sdp@data.info$nGene, nUMI=sdp@data.info$nUMI, orig.ident=sdp@ident)
  sdp <- MeanVarPlot(sdp ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff, do.contour = F)
  vargenes.scells <- sdp@var.genes
  
  vargenes <- data.frame(gene=sort(unique(c(vargenes.1601,vargenes.1888,vargenes.scells))))
  vargenes$run1601 <-  vargenes$run1888 <- vargenes$allruns <- 0
  vargenes$run1601[which(vargenes$gene %in% vargenes.1601)] <- 1
  #vargenes$run1887[which(vargenes$gene %in% vargenes.1887)] <- 1
  vargenes$run1888[which(vargenes$gene %in% vargenes.1888)] <- 1
  vargenes$allruns[which(vargenes$gene %in% vargenes.scells)] <- 1
  rownames(vargenes) <- vargenes$gene
  table(vargenes$allruns, vargenes$run1601)
  nc <- rowSums(vargenes[, c("run1601","run1888")])
  pc.genes <- rownames(vargenes)[which(nc==2)]
  
  return(list(vargenes=vargenes, pc.genes=pc.genes))
}

# x - row is genes, column is cells, and gene expression has been per gene zscore normalized
pc.test <- function(x, B=1000, seed=NULL, min.ev=sqrt(3), max.r.ev=2, df=NULL, do.plot=T,...) {
	if (!is.null(seed)) set.seed(seed)
  n <- ncol(x)
  m <- nrow(x)
  
  if (is.null(df)) df <- n
  ev <- (prcomp(x)$sdev)
  evs <- as.data.frame(matrix(0, nrow=B, ncol=df))
  j <- 1
  for (i in 1:B) {
    xi <- t(apply(x, 1, sample, replace=FALSE))
    evi <- (prcomp(xi)$sdev)
    if (!any(is.na(evi[1:df]))) {
      evs[j, ]<-evi[1:df]
      j <- j + 1
    } else {
      cat("\n", i, " contains NAs\n", sep="")
    }
    if ((i %% 50)==0) cat("\r",i,"/",B, sep="")
  }
  cat("\n")
  if (j < (B+1)) {
    evs <- evs[1:(j-1),]
  }
  colnames(evs) <- paste("P", 1:dim(evs)[2], sep="")
  evsm <- melt(evs)
  colnames(evsm) <- c("PC","EV")
  
  r <- 3
  max.ev <- max(max(evsm[which(evsm$EV<=max.r.ev), "EV"]), min.ev)
  #if (do.ceiling) max.ev <- ceiling(max.ev)
  ev.p <- length(which(evsm$EV>max.ev))/length(evsm$EV)
  r <- length(which(ev>max.ev))
  
  cat(round(100*ev.p, digits = 2), "% random EVs > ", max.ev, "\n", sep="")
  cat(r, " PCs with EV greater than ", max.ev, "\n", sep="")
  
  if (do.plot==T) {
    plot(ev)
    g <- ggplot(data=evsm) + geom_density(aes(x=EV), fill="grey", col="black")
    g <- g + geom_vline(xintercept=max.ev, col="blue")
    g <- g + geom_text(aes( max.ev, 0.2, label = paste("threshold=", round(max.ev,3), sep=""), hjust = -0.1), col="blue", size = 5)
    ev <- sort(ev, decreasing=T)
    g <- g + geom_point(data=data.frame(x=ev[1:r],y=0.5), aes(x=x, y=y), col=2, size=2)
    #for (i in 1:r) {
    # g <- g + geom_vline(xintercept=ev[i], col="red")
    #}
     
    g <- g + theme_bw() 
    g <- g + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    print(g)
  }
  
  if (FALSE) {
    ret <- permutationPA(x, B, verbose=FALSE)
    cat(ret$r, " sig PCs based on permutationPA\n")
  }
  
  return(list(r=r, evsm=evsm))
}

# a wrapper of ComBat for batch correction
pc.batch.combat <- function(x, sample, do.plot=T, par.prior=TRUE, db) {
  library(caret)
  x1 <- t(sva::ComBat(t(x), sample, par.prior = par.prior))
  xm <- melt(data.frame(x, cell=rownames(x), sample=sample), id.vars=c("cell", "sample"))
  x1m <- melt(data.frame(x1, cell=rownames(x1), sample=sample), id.vars=c("cell", "sample"))
  xm$category <- "original"
  x1m$category <- "combat"
  viz <- rbind(xm, x1m)
  #g <- ggplot(viz, aes(x=sample, y=value, fill=sample)) + facet_wrap(~variable + category, scales="free") + geom_point(pch=21)
  #ggsave(file=paste(db, ".pc.combat.pdf", sep=""))
  return(x1)
}

		   
		   
getNormFactor.2 <- function(countsTable) {
  
  libsize <- colSums(countsTable)
  sf <- libsize/median(libsize)
  names(sf) <- colnames(countsTable)
  return(sf)
  
}

nbTestSH.2 <- function (countsTable, conds, condA = "A", condB = "B", numPart = 1, 
    SHonly = FALSE, propForSigma = c(0, 1), shrinkTarget = NULL, 
    shrinkQuantile = NULL, plotASD = FALSE, coLevels = NULL, 
    contrast = NULL, keepLevelsConsistant = FALSE, useMMdisp = FALSE, 
    addRawData = FALSE, shrinkVariance = FALSE, pairedDesign = FALSE, 
    pairedDesign.dispMethod = "per-pair", useFisher = FALSE, 
    Dispersions = NULL, eSlope = 0.05, lwd_ASD = 4.5, cex_ASD = 1.2) 
{
    if (is.null(coLevels)) {
        coLevels = data.frame(exp = rep(1, ncol(countsTable)))
    }
    for (j in 1:ncol(coLevels)) {
        cL = apply(coLevels, 1, paste, collapse = "_")
    }
    cL.tb = table(cL)
    if (sum(cL.tb < 2) > 0) {
        stop("Errors in 'coLevels'. All the levels must have paired comparisons.")
    }
    if (!is.null(contrast)) {
        if (length(contrast) != ncol(countsTable)) {
            stop(paste("Error: the length of contrast vector must equal", 
                "to the number of columns in countsTable."))
        }
        if (length(unique(conds[contrast > 0])) != 1) {
            stop(paste("Error: this package is currently only available", 
                "for the contrast between conditions, not the contrast", 
                "within conditions. Please revise the contrast vector,", 
                "such as c(1,1,-1,-1) for cond=c('A','A', 'B, 'B'),", 
                "instead of c(1,-1,1,-1)."))
        }
        countsTable = countsTable[, contrast != 0]
        conds = conds[contrast != 0]
        cL = cL[contrast != 0]
        contrast = contrast[contrast != 0]
        contrast.mat = t(matrix(contrast, nrow = length(contrast), 
            ncol = nrow(countsTable)))
        if (sum(abs(contrast) != 1) > 0) {
            countsTable = round(countsTable * abs(contrast.mat))
        }
    }
    # sf = getNormFactor(countsTable)
	sf = getNormFactor.2(countsTable)
    colA = conds == condA
    colB = conds == condB
    if (sum(colB[1:sum(colA)]) > 0) {
        stop(paste("re-order the columns of the countsTable so that all", 
            "samples in the same condition are in the adjacent columns.", 
            "For example, 'A A B B' is good, but 'A B A B' is bad."))
    }
    counts = as.matrix(countsTable)
    ng = nrow(counts)
    cntA = matrix(counts[, colA], ncol = sum(colA))
    cntB = matrix(counts[, colB], ncol = sum(colB))
    if (SHonly) {
        disp = nbinomTestForMatricesSH(countsA = cntA, countsB = cntB, 
            sizeFactorsA = sf[colA], sizeFactorsB = sf[colB], 
            numPart = numPart, SHonly = TRUE, propForSigma = propForSigma, 
            shrinkTarget = shrinkTarget, shrinkQuantile = shrinkQuantile, 
            cLA = cL[colA], cLB = cL[colB], keepLevelsConsistant, 
            useMMdisp, shrinkVariance = shrinkVariance, eSlope = eSlope, 
            plotASD = plotASD, lwd_ASD = lwd_ASD, cex_ASD = cex_ASD)
        return(disp)
    }
    else {
        t1 = Sys.time()
        pval0 = nbinomTestForMatricesSH(countsA = cntA, countsB = cntB, 
            sizeFactorsA = sf[colA], sizeFactorsB = sf[colB], 
            numPart = numPart, SHonly = FALSE, propForSigma = propForSigma, 
            shrinkTarget = shrinkTarget, shrinkQuantile = shrinkQuantile, 
            cLA = cL[colA], cLB = cL[colB], contrast = contrast, 
            keepLevelsConsistant = keepLevelsConsistant, useMMdisp, 
            shrinkVariance = shrinkVariance, pairedDesign = pairedDesign, 
            pairedDesign.dispMethod = pairedDesign.dispMethod, 
            useFisher = useFisher, Dispersions = Dispersions, 
            eSlope = eSlope, plotASD = plotASD, lwd_ASD = lwd_ASD, 
            cex_ASD = cex_ASD)
        print(Sys.time() - t1)
        dispMM = pval0$dispMM
        dispSH = pval0$dispSH
        Mean = pval0$mu
        pval = pval0$pval
    }
    cl.nm = sort(unique(c(cL[colA], cL[colB])))
    rM.A = rowMeans(as.matrix(counts[, colA]))
    rM.B = rowMeans(as.matrix(counts[, colB]))
    l2f = log2(rM.A/rM.B)
    rs = data.frame(Mean = Mean, rawMeanA = rM.A, rawMeanB = rM.B, 
        rawLog2FoldChange = l2f, dispMM = dispMM, dispSH = dispSH, 
        pval = pval, stringsAsFactors = FALSE)
    rownames(rs) = rownames(countsTable)
    if (addRawData) {
        rs = data.frame(rs, countsTable)
    }
    return(rs)
}
