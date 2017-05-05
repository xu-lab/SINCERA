# This is demonstration of using SINCERA to analysis E16.5 single cell RNA-seq data (Guo et al., 2015 & Du et al., 2015)
# Author: Minzhe Guo (minzhe.guo@cchmc.org)


library(SINCERA)

# Loading E16.5 data
data("E16.5")

ls()

# expression matrix
dim(expressions)

# cell information
head(cells)

# gene information
head(genes)

# The analysis starts with running the construct function to create an R S4 object,
# which will hold all the data and analysis results.
sc <- construct(exprmatrix=expressions,
                samplevector=paste("sample", cells$SAMPLE, sep=""))

# Identify and remove low quality cells.
# The key parameters of running this function include: “min.expression”,
# which specifies the minimum expression value for a gene to be considered as expressed,
# and “min.genes”, which specifies the lower bound of the number of expressed genes in a cell.
sc <- filterLowQualityCells(sc, min.expression=1, min.genes=1000, do.plot=T)

# Use filterContaminatedCells function to remove potential contaminated cells
# based on the co-expression of known marker genes of two distinct cell types,
# such as the co-expression of mouse lung epithelial marker Epcam and mouse lung endothelial cell marker Pecam1.
# Epithelial ENSMUSG00000045394    Epcam
# Endothelial ENSMUSG00000020717   Pecam1
sc <- filterContaminatedCells(sc, markers.1="ENSMUSG00000045394", markers.2="ENSMUSG00000020717")

# Filter out non- or low-expressed genes, as well as genes that are expressed in less than a certain number of cells per sample preparation.
# By default, genes expressed (>5 FPKM/TPM) in less than 2 cells will be filtered out by this function.
sc <- prefilterGenes(sc, pergroup=TRUE, min.expression=5, min.cells=2, min.samples=2)

# Use expr.minimum function to set a minimum expression value
sc <- expr.minimum(sc, value=0.01)

# Run batch.analysis function to identify batch differences.
sc <- batch.analysis(sc, analysis=c("q", "qq", "ma", "isccd", "distribution"), min.expression=1)

# Perform per-sample z-score transformation
sc <- normalization.zscore(sc, pergroup=TRUE)

# Determine the threshold value for the expression specificity filter, which will be used to select genes for cell cluster identification
rgenes <- rownames(genes)[genes$SYMBOL %in% mouse.ribosomal.genes]
specifity.thresh <- specificity.thresholdSelection(sc, rgenes)

# Select expression profiles that are potentially informative about cell types/states
# and remove genes that may increase noise in the cell type identification step
sc <- cluster.geneSelection(sc, method="specificity", pergroup=TRUE, min.samples=2, specifity.thresh=specifity.thresh)

sc <- doPCA(sc, genes=getGenesForClustering(sc), use.fast = T, center=T, scale.=T)

plotRDS(sc)
plotPCASD(sc)

sc <- doTSNE(sc, dims=1:9)
plotRDS(sc, feature.type="tsne")

# Run cluster.assignment function to assign cells to initial clusters.
# The default algorithm uses hierarchical clustering with average linkage,
# Pearson’s correlation based distance measurement, and z-score transformed expression values of the selected genes.
sc <- cluster.assignment(sc, h=0.5)

# save dendrogram with cell names to pdf file for 
pdf(file="HC.pdf", width=20, height=10) #adjust width or height based on the number of cells in the dendrogram
plotHC(sc, show.labels = T)
dev.off()


# visualize cell clusters in tSNE plot
plotRDS(sc, feature.type="tsne")


# The number of cells in each cluster
print(table(getCellMeta(sc, name="CLUSTER")))

# For reproducing the results in the paper, use the cluster ids in the paper
sc <- setCellMeta(sc, name="CLUSTER", value=paste("C", cells$CLUSTER, sep=""))

sc <- copyCellMeta(sc, "CLUSTER", "GROUP")

# Run plotMarkers function to check the quality of the obtained clustering scheme and
# inspect the expression patterns of a number of known markers across cell clusters.
plotMarkers(sc, genes=as.character(markers$SYMBOL[which(markers$TYPE=="PMP")]))
plotMarkers(sc, genes=as.character(markers$SYMBOL[which(markers$TYPE=="MyoFB")]))
plotMarkers(sc, genes=as.character(markers$SYMBOL[which(markers$TYPE=="Pericyte")]))
plotMarkers(sc, genes=as.character(markers$SYMBOL[which(markers$TYPE=="MatrixFB")]))
plotMarkers(sc, genes=as.character(markers$SYMBOL[which(markers$TYPE=="Endothelial")]))
plotMarkers(sc, genes=as.character(markers$SYMBOL[which(markers$TYPE=="Immune")]))
plotMarkers(sc, genes=as.character(markers$SYMBOL[which(markers$TYPE=="Epithelial")]))

# Heatmap of markers can also be used to inspect the quality of clustering scheme
plotHeatmap(sc, genes=as.character(markers$SYMBOL))


# Run the cluster.permutation.analysis function to perform a cluster membership permutation analysis to determine cluster significance.
sc <- cluster.permutation.analysis(sc)

# Identify differentially expressed genes in each cluster.
sc <- cluster.diffgenes(sc, genes=getGenesForClustering(sc), thresh=0.01)

diffgenes <- getDiffGenes(sc)

sc <- setAssociationTable(sc, value=associations.01112014)

# Run celltype.enrichment function to predict cell type for each cluster
# add gene and cell type association table
sc <- celltype.enrichment(sc, id.type="ENSEMBL")

# Map cell clusters to cell types
sc <- setCellType(sc, do.reset=T)
sc <- setCellType(sc, groups=paste("C",1:9, sep=""), types=c("PMP", "MyoFB","Pericyte","Undefined", "MatrixFB", "Undefined", "Endothelial","Immune","Epithelial"))

# Print the cell cluster and cell type mapping
print(getCellType(sc))

# Pre-collected cell type markers
print(markers)

# Use setCellTypeMarkers function to add the collected markers into SINCERA.
sc <- setCellTypeMarkers(sc, markers)

#	Run celltype.validation function to perform a rank-aggregation-based quantitative assessment of
# the consistency between mapped cell type and the expression pattern of known cell type marker genes.
sc <- celltype.validation(sc)

# Run the signature.prediction function to predict cell type signature genes.
sc <- signature.prediction(sc, use.logireg = TRUE)

# Get the predicted signature genes
siggenes <- getSigGenes(sc)

# Visualize the expression of the predicted signature genes across cell types (clusters)
plotHeatmap(sc, genes=as.character(siggenes$SYMBOL))

# Validate the signature prediction using a repeated random subsampling approach
sc <- signature.validation(sc, siggenes=siggenes, repeats=5)

# add TF knowledge
sc <- setTFs(sc, value=rownames(genes)[which(genes$TF==1)])

# Select candidate transcription factors for driving force analysis
sc <- drivingforce.selectTFs(sc, diff.thresh=0.001, pct.thresh=1)
# Select cell type specific differentially expressed genes or signature genes as candidate target genes
sc <- drivingforce.selectTGs(sc, diff.thresh=0.01)

# Infer a TRN using the cell type specific expression patterns of the selected candidate TFs and TGs.
sc <- drivingforce.inferTRN(sc, groups=c("C9"))
sc <- drivingforce.getLCC(sc, groups=c("C9"), thresh = 0.001)

# Rank TFs based on their importance to the inferred TRN
sc <- drivingforce.rankTFs(sc, groups=c("C9"))





