require(Seurat)
require(ggplot2)
library(dplyr)
require(cowplot)
require(plyr)
require(glmGamPoi)

parser = argparse::ArgumentParser(description="Script to integrate scRNA data")
parser$add_argument('-I','--input', help='input raw matrix or 10X-like dir')
parser$add_argument('-F','--nf',help='per sample HVG number set')
parser$add_argument('-A','--acf',help='integrate anchor features')
parser$add_argument('-AD','--adim',help='Which dimensions to use from the CCA to specify the neighbor search space')
parser$add_argument('-IND','--indim',help='Number of dimensions to use in the anchor weighting procedure')
parser$add_argument('-TN','--tn',help='top markers plot heatmap')
parser$add_argument('-SP','--sp',help='per sample RunPCA npcs')
parser$add_argument('-D','--dim',help='PCA dim usage')
parser$add_argument('-R','--res',help='Map resolution usage')
parser$add_argument('-O','--out',help='out directory')
args = parser$parse_args()

nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 3000)
dim.usage <- as.numeric(if(!is.null(args$dim)) args$dim else 25)
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 0.6)
acf.usage <- as.numeric(if(!is.null(args$acf)) args$acf else 3000)
adim.usage <- as.numeric(if(!is.null(args$adim)) args$adim else 30)
tn.usage <- as.numeric(if(!is.null(args$tn)) args$tn else 20)
sp.usage <- as.numeric(if(!is.null(args$sp)) args$sp else 30)
indim.usage <- as.numeric(if(!is.null(args$indim)) args$indim else 30)

#循环建立path向量
path.list <- read.table(args$input)
for (line in path.list){line}
#循环读取rds文件,seurat对象名必须是 S1_lib2_QCed这样的格式
for (i in line){
abc=strsplit(i,split="/")
id <- abc[[1]][lengths(abc)]
del_id <- as.vector((strsplit(id,split = "\\.")[[1]]))[1]
assign(del_id ,readRDS(i))
}

#建立seurat list对象
id.list <- c()
for (i in grep(ls(),pattern = 'QCed',value=TRUE)) { id.list <- append(id.list,i)}
rds.list <- lapply(id.list, get)


rds.list <- lapply( X = rds.list, FUN = function(x){ 
       x <- SCTransform(x ,vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n = nf.usage)
})

features  <- SelectIntegrationFeatures(object.list = rds.list, nfeatures = acf.usage)
rds.list <- PrepSCTIntegration(object.list = rds.list, anchor.features = features)

rds.list <- lapply(X = rds.list, FUN = function(x) {
       x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = rds.list, normalization.method = "SCT", anchor.features = features, dims = 1:adim.usage, reduction = "rpca", k.filter = 200)
Bmor.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:indim.usage)


write.table(Bmor.integrated@assays[["integrated"]]@var.features, paste0(args$out,'/sct_variable.txt'),sep = '\t',row.names =FALSE, col.names =FALSE,quote =FALSE)
#看top30的高变基因有哪些
save.image(paste0(args$out,"/only_integrete.RData"))
#做pca，且画出前20PCA热图
Bmor.integrated <- RunPCA(Bmor.integrated, npcs = dim.usage)
pdf(paste0(args$out,"/sct_PCA_DimHeatmap.pdf"))
DimHeatmap(Bmor.integrated, dims = 1:20, cells = 100, balanced = TRUE)
dev.off()

pdf(paste0(args$out,"/sct_PCA_ElbowPlot.pdf"))
ElbowPlot(Bmor.integrated)
dev.off()

Bmor.integrated <- RunTSNE(Bmor.integrated, dims = 1:(dim.usage-5), reduction = "pca",check_duplicates = FALSE)
Bmor.integrated <- RunUMAP(Bmor.integrated, dims = 1:(dim.usage-5), reduction = "pca")
Bmor.integrated <- FindNeighbors(Bmor.integrated, reduction = "pca", dims = 1:(dim.usage-5))
Bmor.integrated <- FindClusters(Bmor.integrated, resolution = res.usage, n.start=10)

ncluster = as.numeric(nlevels(Bmor.integrated@active.ident))
best_color<- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFE4E1","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#1B4400","#4FC601","#3B5DFF","#FF2F80","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#DDEFFF","#7B4F4B","#A1C299","#0AA6D8","#00A087FF","#4DBBD5FF","#E64B35FF","#3C5488FF","#F38400","#A1CAF1", "#C2B280","#848482","#E68FAC", "#0067A5","#F99379", "#604E97","#F6A600", "#B3446C","#DCD300","#882D17", "#8DB600","#654522", "#E25822", "#2B3D26","#191970","#000080","#6495ED","#1E90FF","#00BFFF","#00FFFF","#FF1493","#FF00FF","#A020F0","#63B8FF","#008B8B","#54FF9F","#00FF00","#76EE00","#FFF68F","Yellow1","Gold1","DarkGoldenrod4","#FF6A6A","#FF8247","#FFA54F","#FF7F24","#FF3030","#FFA500","#FF7F00","#FF7256","#FF6347","#FF4500","#FF1493","#FF6EB4","#EE30A7","#8B008B")

saveRDS(Bmor.integrated, paste0(args$out,"/sct_integrate.RDS"))

p1 <- DimPlot(Bmor.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Bmor.integrated, reduction = "umap", label = T, pt.size=0.5,label.size = 5,cols = best_color[1:ncluster])
p3 <- DimPlot(Bmor.integrated, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(Bmor.integrated, reduction = "tsne", label = T, pt.size=0.5,label.size = 5,cols = best_color[1:ncluster])

pdf(paste0(args$out,"/sct_cell_cluster.umap.pdf"),height=10,width=20)
plot_grid(p1, p2)
dev.off()
pdf(paste0(args$out,"/sct_cell_cluster.tsne.pdf"),height=10,width=20)
plot_grid(p3, p4)
dev.off()
pdf(paste0(args$out, "/sct_cell_cluster_umi.pdf"),height=10,width=20)
p5 <- ggplot(as.data.frame(Bmor.integrated@reductions$umap@cell.embeddings), aes(x=UMAP_1, y=UMAP_2, color=Bmor.integrated@meta.data$nCount_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="UMIs") + scale_color_viridis_c(direction = -1)
p6 <- ggplot(as.data.frame(Bmor.integrated@reductions$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2, color=Bmor.integrated@meta.data$nCount_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="UMIs") + scale_color_viridis_c(direction = -1)
plot_grid(p5, p6)
dev.off()
pdf(paste0(args$out,"/sct_cell_cluster_gene.pdf"),height=10,width=20)
p7 <- ggplot(as.data.frame(Bmor.integrated@reductions$umap@cell.embeddings), aes(x=UMAP_1, y=UMAP_2, color=Bmor.integrated@meta.data$nFeature_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="genes") + scale_color_viridis_c(direction = -1)
p8 <- ggplot(as.data.frame(Bmor.integrated@reductions$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2, color=Bmor.integrated@meta.data$nFeature_RNA)) + geom_point(size=0.5) + theme(legend.position = "none") + theme_bw() + labs(color="genes") + scale_color_viridis_c(direction = -1)
plot_grid(p7, p8)
dev.off()

DefaultAssay(Bmor.integrated) <- "RNA"
       markers <- FindAllMarkers(Bmor.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay = "RNA")
write.table(markers,paste0(args$out,'/sct_markers.txt'),sep = '\t',row.names =FALSE, col.names =TRUE,quote =FALSE) 
top <- markers %>% group_by(cluster) %>% top_n(n = tn.usage, wt = avg_log2FC)
pdf(paste0(args$out,"/top_marker_heatmap.pdf"), width = 25, height = 20)
DoHeatmap( 
  Bmor.integrated,
  features = top$gene,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "integrated",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)

dev.off()

pdf(paste0(args$out,"/sct_neuron_glia.markers_tsne.pdf"), width = 15, height = 15)
FeaturePlot(Bmor.integrated, reduction = "tsne", features = c("repo", "GLaz","Syt1","nSyb"), pt.size = 0.1)
dev.off()
pdf(paste0(args$out,"/sct_neuron_glia.markers_umap.pdf"), width = 15, height = 15)
FeaturePlot(Bmor.integrated, reduction = "umap", features = c("repo", "GLaz","Syt1","nSyb"), pt.size = 0.1)
dev.off()

p1<- VlnPlot(Bmor.integrated, pt.size = 0, features = c("BMSK0001774", "BMSK0000580","BMSK0008066","BMSK0008068","trio","Eip75B","Eip93F","mub","Pka-C1","ort","Ih","qvr","bt"))
pdf(paste0(args$out,'/kc.pdf'),height=20,width=70)
plot_grid(p1)
dev.off()

p1<- VlnPlot(Bmor.integrated,pt.size = 0 ,features = c("repo", "GLaz","bdl","Syt1","nSyb","Mlc-c"))
pdf(paste0(args$out,'/neuron_gial.pdf'),height=20,width=60)
plot_grid(p1)
dev.off()

p3 <- VlnPlot(Bmor.integrated,pt.size = 0, features = c("BMSK0001207","ple","SerT","Tdc2","DAT","Ddc","Trh","BMSK0004580"))
pdf(paste0(args$out,'/ol_ple_sert_.pdf'),height=20,width=60)
plot_grid(p3)
dev.off()

p3 <- FeaturePlot(Bmor.integrated,feature = c("DAT","ple"))
pdf(paste0(args$out,'/ple_featureplot.pdf'),height=9,width=18)
plot_grid(p3)
dev.off()


p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("BMSK0008113","trp","ninaB","ninaC","Sulf1"))
pdf(paste0(args$out,'/trp.pdf'),height=15,width=50)
plot_grid(p4)
dev.off()

p5 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("Gs2", "BMSK0011682","BMSK0012704","Eaat1","Rh50","Gat","BMSK0011307"))
pdf(paste0(args$out,'/Astrocytes.pdf'),height=20,width=60)
plot_grid(p5)
dev.off()

p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("hml","Fer2LCH"))
pdf(paste0(args$out,'/hml.pdf'),height=15,width=30)
plot_grid(p4)
dev.off()
p8 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("amon","BMSK0006239", "BMSK0006241","BMSK0015259","BMSK0005225","CCAP"))
pdf(paste0(args$out,"/pep.pdf"),height=20,width=50)
plot_grid(p8)
dev.off()

p8 <- VlnPlot(Bmor.integrated,pt.size = 0, features = c("dpn","Imp","pros","ab","Syp"))
pdf(paste0(args$out,"/neuroblast.pdf"),height=20,width=40)
plot_grid(p8)
dev.off()

p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("hoe1","zyd"))
pdf(paste0(args$out,'/Cotex_glia.pdf'),height=15,width=20)
plot_grid(p4)
dev.off()
p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("hth","Oaz","acj6","kn"))
pdf(paste0(args$out,'/Olfactory_projectional_Neuron.pdf'),height=15,width=50)
plot_grid(p4)
dev.off()

p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("Mdr49","Indy"))
pdf(paste0(args$out,'/Surface_glia.pdf'),height=15,width=50)
plot_grid(p4)
dev.off()

p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("BMSK0010690","BMSK0010691","BMSK0000716"))
pdf(paste0(args$out,'subperineurial_glia.pdf'),height=15,width=50)
plot_grid(p4)
dev.off()

p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("SPARC","vkg"))
pdf(paste0(args$out,'/perineurial_glia.pdf'),height=15,width=30)
plot_grid(p4)
dev.off()
p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("egr","Tsf1","e","Idgf4"))
pdf(paste0(args$out,'/Ensheathing_glia.pdf'),height=15,width=30)
plot_grid(p4)
dev.off()
p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("BMSK0001420"))
pdf(paste0(args$out,'/lozenge.pdf'),height=15,width=20)
plot_grid(p4)
dev.off()
p4 <- VlnPlot(Bmor.integrated, pt.size = 0,features = c("bsh","Lim3","svp","BMSK0012283","Mip"))
pdf(paste0(args$out,'/ol.pdf'),height=15,width=50)
plot_grid(p4)
dev.off()

Bmor.integrated$sample <- vapply( colnames(Bmor.integrated), function(x) strsplit(x, split = "_C")[[1]][1], FUN.VALUE = character(1) )

Bmor.integrated$sample <- factor(Bmor.integrated$sample)

Bmor.integrated$single <- vapply( Bmor.integrated$single, function(x) gsub("S1_lib5","s10",x)[[1]][1], FUN.VALUE = as.character(1))
Bmor.integrated$single <- vapply( Bmor.integrated$single, function(x) gsub("S1_lib6","s10",x)[[1]][1], FUN.VALUE = as.character(1))


pdf(paste0(args$out,'/split_tsne.pdf'),height=20,width=180)
DimPlot(Bmor.integrated,reduction="tsne",split.by="single",pt.size= 1,label.size = 8,label=T,cols = best_color[1:ncluster])
dev.off()

pdf(paste0(args$out,'/split_umap.pdf'),height=20,width=180)
DimPlot(Bmor.integrated,reduction="umap",split.by="single",pt.size= 1,label.size = 8,label=T,cols = best_color[1:ncluster])
dev.off()




