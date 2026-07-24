### Load library
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(DoubletFinder)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(sctransform)))

## ---------------------------
##
## Script name:  
##
## Purpose of script: filter low quality cell and generate report graph and statistics
##
## Author: Yuanzhen Zhu, Jiahao Xiang 
##
## Date Created: 2020-11-11, Data Modified: 2021-06-01, 2021-10-08
##
## Version: V2.1
##
## ---------------------------
##
## Notes: This script is universal and general used!
##
## ---------------------------

### Get the parameters

parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input raw matrix or 10X-like dir')
parser$add_argument('-F','--nf',help='HVG number set')
parser$add_argument('-D','--dim',help='PCA dim usage')
parser$add_argument('-P','--percentage',help='doublets percentage')
parser$add_argument('-R','--res',help='Map resolution usage')
parser$add_argument('-O','--out',help='out directory')
parser$add_argument('-S','--sample',help='sample sample name')
args = parser$parse_args()


#EC.data <- readRDS(args$input)
EC.data <- Read10X(data.dir = args$input,gene.column = 1)
EC.data <- EC.data[,-1]
nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 3000)
dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 25)
#doublets.percentage <- if(!is.null(args$percentage)) args$percentage else 0.05
#doublets.percentage <- as.numeric(doublets.percentage)
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 0.6)


### Creat Seurat object, basic filtering and statistics
EC <- CreateSeuratObject(EC.data, project = args$sample, min.cells = 3, min.features = 200)
EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "AOB78-")

meanT <- sum(EC@meta.data$nCount_RNA)/nrow(EC@meta.data)
meanG <- sum(EC@meta.data$nFeature_RNA)/nrow(EC@meta.data)
meanMT <- sum(EC[["percent.mt"]])/nrow(EC@meta.data)
info <- c(meanT, meanG, meanMT)

pdf(paste0(args$out, "/count_mt_sctransform.cor.pdf"))
FeatureScatter(EC, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf(paste0(args$out, "/count_gene_sctransform.cor.pdf"))
FeatureScatter(EC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("Gene") + labs(title="")+ theme(axis.text.x = element_blank())
p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("Transcript") + labs(title="")+ theme(axis.text.x = element_blank()) + ylim(0,15000)
p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
p <-p1|p2|p3
pdf(paste0(args$out,"/QC_before.dis.SCTransform.pdf"))
print(p)
dev.off()

ECmeta <- EC@meta.data[order(-EC@meta.data$nFeature_RNA),]
n95 <- as.numeric(as.integer(nrow(ECmeta) * 0.05))
n95features <- as.numeric(ECmeta[n95, "nFeature_RNA"])
#EC <- subset(EC, subset = nFeature_RNA > 400 & nFeature_RNA < n95features & percent.mt < 5)
EC <- subset(EC, subset = nFeature_RNA > 400 & nFeature_RNA < 2500 & percent.mt < 5)
p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("Gene") + labs(title="")+ theme(axis.text.x = element_blank())
p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("Transcript") + labs(title="")+ theme(axis.text.x = element_blank()) + ylim(0,15000)
p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
p <-p1|p2|p3
pdf(paste0(args$out,"/QC_after.dis.SCTransform.pdf"))
print(p)
dev.off()


meanT <- sum(EC@meta.data$nCount_RNA)/nrow(EC@meta.data)
meanG <- sum(EC@meta.data$nFeature_RNA)/nrow(EC@meta.data)
meanMT <- sum(EC[["percent.mt"]])/nrow(EC@meta.data)
info <- rbind(info, c(meanT, meanG, meanMT))
rownames(info) <- c("before", "after")
colnames(info) <- c("Transcripts", "Genes", "pct.mt")
write.table(info, paste0(args$out, "/QC.mean.SCTransform.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

EC <- SCTransform(EC, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = nf.usage)

pdf(paste0(args$out, "/VariableFeaturePlot.SCTransform.pdf"))
VariableFeaturePlot(EC, selection.method = "sct")
dev.off()

### PCA, statistics
DefaultAssay(EC) <- "SCT"
EC <- RunPCA(EC, npcs = dim.usage)
pdf(paste0(args$out,"/PCA_DimHeatmap.SCTransform.pdf"))
DimHeatmap(EC, dims = 1:dim.usage, cells = 100, balanced = TRUE)
dev.off()

pdf(paste0(args$out, "/PCA_ElbowPlot.SCTransform.pdf"))
ElbowPlot(EC)
dev.off()

### Choose suitable dim for UMAP & T-SNE: dim reduce unfavored signal, determine percent of variation associated with each PC; calculate cumulative percents for each PCs; Determine the difference between variation of PC and subsequent PC;
# pct <- EC[["pca"]]@stdev / sum(EC[["pca"]]@stdev) * 100
# cumu <- cumsum(pct)
# co1 <- which(cumu > 80 & pct < 5)[1]
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# pcs <- min(co1, co2)
# pdf(paste0(args$out, "/PCA_choose.SCTransform.pdf"))
# plot(EC[["pca"]]@stdev, xlab="PC", ylab = "SD explaineded (%)")
# abline(v=pcs, col="red")
# dev.off()
# print( paste0("PCA dim: ", as.character(dim.usage), "; Favored choose dim 1:", as.character(pcs)))

### UMAP and T-SNE
EC <- FindNeighbors(EC, dims = 1:(dim.usage-5))
EC <- FindClusters(EC, resolution = res.usage)
EC <- RunUMAP(EC, dims = 1:(dim.usage-5))
EC <- RunTSNE(EC, dims = 1:(dim.usage-5))

p1 <- DimPlot(EC, reduction = "umap")
p2 <- DimPlot(EC, reduction = "tsne")
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE.SCTransform.pdf"), plot = p, device = "pdf", width = 18, height = 9)

saveRDS(EC,paste(args$out,"/",args$sample,"_sctransformed_QCed.RDS",sep=""))


