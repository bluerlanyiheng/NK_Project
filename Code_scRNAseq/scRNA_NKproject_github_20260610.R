library(Seurat)
library(patchwork)
library(clustree)
library(tidyverse)
library(plyr)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(harmony)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
set.seed(101)
library(future)
library(scDblFinder)
library(glmGamPoi)
#BiocManager::install('glmGamPoi')

plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #50G
options(future.seed = TRUE)
set.resolutions <- seq(0.2, 1.0, by = 0.1)
setwd("/data/activate_data/lanyiheng/zhaoweixi/20260128_singlecellZNF205/Results/Final_figures_20260307")

# Doublet detection function
FunctionScDblFinder.RNA <- function(Seurat.object, clusters = NULL, samples = NULL, nfeatures = 2000, dims = 20, dbr = NULL, n.core = 10, seed = 12345){
  set.seed(seed)
  require(scDblFinder)
  require(BiocParallel)
  sce <- as.SingleCellExperiment(Seurat.object)
  sce <- scDblFinder(sce, clusters = clusters, samples = samples, nfeatures = nfeatures, dims = dims, BPPARAM = MulticoreParam(n.core))
  Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.class, "scDblFinder.class")
  Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.score, "scDblFinder.score")
  return(Seurat.object)
}

# Load data
samples <- c("Vector", "ZNF205")
base_dir <- "/data/activate_data/lanyiheng/zhaoweixi/20260128_singlecellZNF205/2260128_SC_ZNF205/data/matrix"
scRNA.list <- lapply(samples, function(x){
  scRNA.data <- Read10X(data.dir = file.path(base_dir, paste0("26011103_", x), "filtered_feature_bc_matrix"))
  scRNA.data <- CreateSeuratObject(counts = scRNA.data, project = x, min.cells = 3, min.features = 200)
  scRNA.data[["mt_ratio_RNA"]] <- PercentageFeatureSet(scRNA.data, pattern = "^MT-")
  scRNA.data[["rp_ratio_RNA"]] <- PercentageFeatureSet(scRNA.data, pattern = "^RPL|^RPS")
  scRNA.data <- subset(scRNA.data, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & mt_ratio_RNA < 5 & rp_ratio_RNA < 25)
  scRNA.data <- FunctionScDblFinder.RNA(Seurat.object = scRNA.data, nfeatures = 2000, dims = 30, n.core = 16)
  scRNA.data <- subset(scRNA.data, subset = scDblFinder.class == "singlet")
  return(scRNA.data)
})

# Merge objects
scRNA <- merge(scRNA.list[[1]], y = scRNA.list[2:length(scRNA.list)],
               add.cell.ids = c("Vector", "ZNF205"),
               project = "scRNA")
orig.ident <- mapvalues(scRNA$orig.ident, 
                        from = samples, 
                        to = c("Vector", "ZNF205"))
scRNA <- AddMetaData(object = scRNA, metadata = orig.ident, col.name = "orig.ident")
scRNA.pro <- subset(scRNA, subset = nFeature_RNA < 6000 & nCount_RNA < 15000)
saveRDS(scRNA.pro, file = "scRNA.pro.rds")

scRNA.pro <- readRDS("/data/activate_data/lanyiheng/zhaoweixi/20260128_singlecellZNF205/Results/Final_figures_20260305/scRNA.pro.NKcell.v4.20260305.rds")

# Filtering by PTPRC (CD45)
genes <- c("PTPRC")
cnt_Vector <- GetAssayData(scRNA.pro, assay = "RNA", layer = "counts.Vector")[genes, , drop = FALSE]
cnt_ZNF205 <- GetAssayData(scRNA.pro, assay = "RNA", layer = "counts.ZNF205")[genes, , drop = FALSE]
keep_Vector <- (cnt_Vector["PTPRC", ] > 0 )
keep_ZNF205 <- (cnt_ZNF205["PTPRC", ] > 0 )
cells_use <- c(colnames(cnt_Vector)[keep_Vector],
               colnames(cnt_ZNF205)[keep_ZNF205])
scRNA.pro <- subset(scRNA.pro, cells = cells_use)

# QC violin plots
pdf("1.FinalFiltering.QC.pdf")
VlnPlot(
  object = scRNA.pro,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 2,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
VlnPlot(
  object = scRNA.pro,
  features = c("mt_ratio_RNA", "rp_ratio_RNA"),
  ncol = 2,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
VlnPlot(
  object = scRNA.pro,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 2,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & yscale("log10", .format = TRUE) & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

# Batch correction, PCA, clustering, UMAP/tSNE
scRNA.pro <- NormalizeData(scRNA.pro)
scRNA.pro <- SCTransform(
  scRNA.pro,
  vst.flavor     = "v2",
  vars.to.regress= c("nCount_RNA","mt_ratio_RNA", "rp_ratio_RNA"),
  verbose        = FALSE
)

length(VariableFeatures(scRNA.pro))
scRNA.pro <- RunPCA(scRNA.pro, npcs = 50, features = VariableFeatures(scRNA.pro), verbose = FALSE)
pdf("2.pca_before_batch_effect.pdf")
DimPlot(object = scRNA.pro, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA.pro, ndims = 50)
dev.off()

options(future.seed = TRUE)
scRNA.pro <- FindNeighbors(scRNA.pro, dims = 1:30, reduction = "pca", verbose = FALSE)
scRNA.pro <- FindClusters(scRNA.pro, reduction = "pca",resolution = set.resolutions,dims = 1:30, verbose = FALSE)
scRNA.pro <- RunUMAP(scRNA.pro, dims = 1:20, n.neighbors = 50,reduction = "pca",verbose = FALSE)
scRNA.pro <- RunTSNE(scRNA.pro, dims = 1:20,reduction = "pca",verbose = FALSE)

pdf("3.umap&tsne_res.1.0.pdf")
DimPlot(scRNA.pro, reduction = "umap", group.by = "seurat_clusters")
DimPlot(scRNA.pro, reduction = "tsne", group.by = "seurat_clusters")
dev.off()

pdf("clutree.pdf")
clustree(scRNA.pro)
dev.off()
scRNA.pro$seurat_clusters <- scRNA.pro$SCT_snn_res.0.3

pdf("4.umap&tsne_res.0.3.pdf", height = 5, width = 6)
DimPlot(scRNA.pro, reduction = "umap",label = T,repel = TRUE, group.by = "SCT_snn_res.0.3", pt.size = 0.5,label.size = 4)+ coord_fixed()
DimPlot(scRNA.pro, reduction = "umap",label = T,repel = TRUE, group.by = "SCT_snn_res.0.3", split.by = "orig.ident", pt.size = 0.5,label.size = 4)+ coord_fixed()
DimPlot(scRNA.pro, reduction = "umap",label = T,repel = TRUE, group.by = "orig.ident", pt.size = 0.5,label.size = 4)+ coord_fixed()
dev.off()

saveRDS(scRNA.pro, file = "scRNA.pro.NKcell.v2.rds")

# Cluster marker plotting
DefaultAssay(scRNA.pro) <- "RNA"
pdf("5.cluster.FeaturePlot.pdf", height = 30, width = 14)
FeaturePlot(scRNA.pro, reduction = "umap", cols = c("lightgrey", "red"), split.by = "orig.ident", features = c("CD3D", "PTPRC", "CD4", "CD8A", "GZMB", "MKI67"), ncol = 2, pt.size = 0.6,label.size = 4)+ coord_fixed()
FeaturePlot(scRNA.pro, reduction = "tsne", cols = c("lightgrey", "red"), split.by = "orig.ident", features = c("CD3D", "PTPRC", "CD4", "CD8A", "GZMB", "MKI67"), ncol = 2, pt.size = 0.6,label.size = 4)+ coord_fixed()
dev.off()
pdf("5.qc.FeaturePlot.pdf", height = 26, width =8)
FeaturePlot(
  scRNA.pro,
  reduction  = "umap",
  features   = c('EGFR','EPCAM','KRT18','KRT19','LAG3', 'TIGIT',"IFNG", "PTPRC","NCAM1","CD300A","MKI67","nCount_RNA","mt_ratio_RNA", "rp_ratio_RNA","scDblFinder.score","nFeature_RNA"),
  cols       = c("lightgrey", "red"),
  pt.size    = 0.6,
  label.size = 4,
  ncol       = 2
)+ coord_fixed()
dev.off()

# Differential expression
scRNA.pro <- PrepSCTFindMarkers(scRNA.pro)
Idents(scRNA.pro) <- scRNA.pro$seurat_clusters
DEG.cluster <- FindAllMarkers(scRNA.pro, only.pos = TRUE,
                              group.by = "seurat_clusters",
                              assay = "SCT", min.pct = 0.25)
DEG.cluster.sig <- DEG.cluster[which(DEG.cluster$p_val_adj < 0.05 & DEG.cluster$avg_log2FC > 0.25),]
saveRDS(DEG.cluster.sig, file = "DEG.cluster.sig.rds")
write.table(DEG.cluster.sig, file = "DEG.cluster.sig.txt", quote = F, sep = "\t", row.names = F)

# Dotplots
pdf("8.cluster.Dotplot.topmakers.pdf", height = 8, width =24)
featureMarker <- c("PTPRC",
                   "KLRC1", "SRGAP3", "PLCB1", "GZMK", "KRT86", "GAS7", "FCER1G", "KLRB1", "LAG3", "CD3E",
                   "RGS9", "KLRC2", "CD2", "GNLY", "CD3G", "CD3D", "CD8B", "TRAC", "CD8A", "CAMK4", "MARCKS",
                   "BAIAP2L1", "ZNF90", "H3F3C", "MEF2C", "AIF1", "IFNG", "CCL4", "CCL3", "CD69", "CCL4L2",
                   "NFKBIA", "MKI67", "TOP2A", "PCLAF", "TYMS", "CDK1", "ASPM", "CLSPN", "UTRN", "PTPRJ",
                   "ARID1B", "NEAT1", "FTX", "CD300A", "FGFBP2", "SELL", "TGFBR3", "MYOM2", "CX3CR1"
)
DotPlot(scRNA.pro, features = featureMarker, assay = "SCT", group.by = "seurat_clusters", cols = c("#1e90ff", "#ff5a36"), dot.scale = 12) + theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5))
dev.off()

# Cluster renaming
scRNA.pro$seurat_clusters
clusters <- plyr::mapvalues(scRNA.pro$seurat_clusters,
                            from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"),
                            to = c("1", "2", "3", "1", "4", "5", "6", "7", "2", "8"))
scRNA.pro$cellType <- paste0("C", clusters)
scRNA.pro$cellType <- factor(scRNA.pro$cellType, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"))
scRNA.pro$seurat_clusters <- scRNA.pro$cellType
Idents(scRNA.pro) <- scRNA.pro$cellType

# UMAP/TSNE by cluster
pdf("7.umap_Final.splitraw2.pdf", height = 5, width = 10)
DimPlot(scRNA.pro, reduction = "umap",label = T,repel = TRUE, split.by = "orig.ident", group.by = "seurat_clusters", pt.size = 1,label.size = 8)+ NoLegend() + coord_fixed()
dev.off()
pdf("7.umap_Final.mergeraw2.pdf", height = 5, width = 5)
DimPlot(scRNA.pro, reduction = "umap",label = T,repel = TRUE, group.by = "seurat_clusters", pt.size = 1,label.size = 8)+ NoLegend()+ coord_fixed()
dev.off()

# Cluster proportions
cell_percentage <- scRNA.pro@active.ident
write.table(cell_percentage,file = 'merged_cell_percentage.csv',sep = ",")


cell_ratio <- prop.table(table(Idents(scRNA.pro)))
write.table(cell_ratio,file = 'cell_ratio_merged.csv',sep = ",")
Cellratio <- prop.table(table(Idents(scRNA.pro), scRNA.pro$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)
colnames(Cellratio)<-c("cellType","group","cellratio")

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


p1<-ggplot(Cellratio) + geom_bar(aes(x =group, y= cellratio, fill = cellType),stat = "identity",width = 0.7,linewidth = 0.5,colour = '#222222')+
  theme_classic() + labs(y = 'Ratio',x='Groups') + scale_fill_manual(values = allcolour)+ theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))
ggsave(p1,file="9.cell_ratio.pdf",width=3,height=8)

# Proportion change comparison
cellRatio <- as.data.frame(table(scRNA.pro@meta.data[,c("orig.ident", "cellType")]))
cellNumber <- as.data.frame(table(scRNA.pro@meta.data[,c("orig.ident")]))
cellNumber <- rep(cellNumber$Freq, times = length(levels(scRNA.pro$seurat_clusters)))
cellRatio$cellNumber <- cellNumber
cellRatio$ratio <- round(cellRatio$Freq / cellRatio$cellNumber * 100, 2)
ZNF205.change <- cellRatio$ratio[seq(2, nrow(cellRatio), by = 2)] / cellRatio$ratio[seq(1, nrow(cellRatio), by = 2)]
ZNF205.change <- data.frame(cellType = cellRatio$cellType[seq(2, nrow(cellRatio), by = 2)], change = ZNF205.change)
ZNF205.change$type <- "Increase"
ZNF205.change$type[which(ZNF205.change$change < 1)] <- "Decrease"
ZNF205.change$type <- factor(ZNF205.change$type, levels = c("Increase", "Decrease"))
ZNF205.change$change <- log2(ZNF205.change$change)
pdf("10.cell_ratio_change.pdf", height = 5, width = 4)
ggbarplot(ZNF205.change, x = "cellType", y = "change", fill = "type", color = "type", sort.val = "desc", palette = c("#E64B35", "#214DA9"), position = position_dodge(1), ylab = "Log2(fold change (ZNF205 vs Vector))", xlab = "") + rotate_x_text(90)
dev.off()

saveRDS(scRNA.pro, file = "scRNA.pro.NKcell.v3.rds")


# Cell number stats
stats <- table(scRNA.pro$celltype.group)
stats_df <- as.data.frame(stats)
colnames(stats_df) <- c("celltype.group", "cell_num")
write.table(
  stats_df,
  file = "celltype_group_cellcounts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Add module scores (gene signatures)
DefaultAssay(scRNA.pro) <- "SCT"
# (Gene sets are defined here for downstream AddModuleScore use)
# Cytokines, chemotaxis, activating/inhibitory receptors, proliferation, glycolysis, etc.
# ... (keep all gene set assignments as is; omitted here for clarity) ...

scRNA.pro <- AddModuleScore(scRNA.pro, list(Cytokines_Genes), name = "Cytokines", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Chemotaxis_receptors_Genes), name = "Chemotaxis_receptors", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Activating_receptors_Genes), name = "Activating_receptors", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Inhibitory_receptor_Genes), name = "Inhibitory_receptor", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Cytokine_receptor_Genes), name = "Cytokine_receptor", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(cytotoxicity_Genes), name = "Cytotoxicity", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(AP1_signaling_Genes), name = "AP1_signaling", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Cell_proliferation_Genes), name = "Proliferation", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Cytotoxicityall_Genes), name = "Cytotoxicityall", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Inflammatory_Genes), name = "Inflammatory", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Stress_Genes), name = "Stress", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Glycolysis_REACTOME_Genes), name = "Glycolysis_REACTOME", seed = 2,nbin = 20)
scRNA.pro <- AddModuleScore(scRNA.pro, list(Glycolysis_HALLMARK_Genes), name = "Glycolysis_HALLMARK", seed = 2,nbin = 20)

saveRDS(scRNA.pro, file = "scRNA.pro.NKcell.v4.rds")


# DE analysis: ZNF205 vs Vector within clusters
scRNA.pro@meta.data$group <- scRNA.pro$orig.ident
scRNA.pro$celltype.group <- paste(scRNA.pro@active.ident, scRNA.pro$group, sep = "_")
scRNA.pro$celltype <- Idents(scRNA.pro)
DefaultAssay(scRNA.pro) <- "RNA"
for(i in c("C1","C2", "C3", "C4", "C5","C6", "C7","C8")){
  cluster <- subset(scRNA.pro, idents=i)
  cluster <- JoinLayers(cluster)
  Idents(cluster) <- "group"
  diff_expr <- FindMarkers(cluster, ident.1="ZNF205", ident.2="Vector", assay="RNA", test.use = "wilcox", logfc.threshold=0, verbose=FALSE)
  write.table(diff_expr, file=paste0(i,"_diff_expr_ZNF205-vs-Vector.xls"), sep="\t", quote=FALSE)
}

# Merge DE results
all_df <- list()
for(i in c("C1","C2", "C3", "C4", "C5","C6", "C7","C8")){
  fname <- paste0(i,"_diff_expr_ZNF205-vs-Vector.xls")
  tmp <- read.table(fname, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  tmp$Cluster <- i
  all_df[[i]] <- tmp
}
final_df <- bind_rows(all_df)
write.table(final_df, file="AllCluster_diff_expr_ZNF205-vs-Vector.xls",
            sep="\t", quote=FALSE, row.names=FALSE)


saveRDS(scRNA.pro, file = "scRNA.pro.NKcell.v5.rds")


# Plotting gene expression and signature scores
DefaultAssay(scRNA.pro) <- "SCT"
Idents(scRNA.pro) <- scRNA.pro$cellType
plist <- VlnPlot(
  scRNA.pro,
  layer = "data",
  features = c("IFNG", "NKG7", "CST7", "TNFSF10", "GZMA", "GZMB", "GZMH", "GZMM", "GZMK", "GNLY", "PRF1", "CTSW"),
  pt.size = 0,
  ncol = 2,
  log = T,
  combine = T
)
plist_box <- lapply(plist, function(x){
  x + geom_boxplot(width=0.14, color="black", fill="white", outlier.shape=NA)
})

pdf("11.cytotoxicity_marker_vln.pdf", width=16, height=24)
wrap_plots(plist_box, ncol=2)
dev.off()


DefaultAssay(scRNA.pro) <- "SCT"
Idents(scRNA.pro) <- scRNA.pro$cellType
plist <- VlnPlot(
  scRNA.pro,
  layer = "data",
  features = c("KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2", "LILRB1", "LAG3", "PDCD1", "SIGLEC7", "CD300A", "CD96", "IL1RAPL1", "TIGIT", "HAVCR2"),
  pt.size = 0,
  ncol = 2,
  log = T,
  combine = T
)
plist_box <- lapply(plist, function(x){
  x + geom_boxplot(width=0.14, color="black", fill="white", outlier.shape=NA)
})

pdf("12.inhibitory_marker_vln.pdf", width=16, height=28)
wrap_plots(plist_box, ncol=2)
dev.off()

# Signature scores - per cluster violin
colnames(scRNA.pro@meta.data)
Cytotoxicity2_df <- scRNA.pro@meta.data[,c("celltype","Cytotoxicity1")]
colnames(Cytotoxicity2_df) <- c("CellType","Signature score")
Inflammatory2_df <- scRNA.pro@meta.data[,c("celltype","Inflammatory1")]
colnames(Inflammatory2_df) <- c("CellType","Signature score")
Stress2_df <- scRNA.pro@meta.data[,c("celltype","Stress1")]
colnames(Stress2_df) <- c("CellType","Signature score")
Cytokines2_df <- scRNA.pro@meta.data[,c("celltype","Cytokines1")]
colnames(Cytokines2_df) <- c("CellType","Signature score")
Inhibitory2_df <- scRNA.pro@meta.data[,c("celltype","Inhibitory_receptor1")]
colnames(Inhibitory2_df) <- c("CellType","Signature score")
Proliferation2_df <- scRNA.pro@meta.data[,c("celltype","Proliferation1")]
colnames(Proliferation2_df) <- c("CellType","Signature score")
Glycolysis_REACTOME2_df <- scRNA.pro@meta.data[,c("celltype","Glycolysis_REACTOME1")]
colnames(Glycolysis_REACTOME2_df) <- c("CellType","Signature score")
Glycolysis_HALLMARK2_df <- scRNA.pro@meta.data[,c("celltype","Glycolysis_HALLMARK1")]
colnames(Glycolysis_HALLMARK2_df) <- c("CellType","Signature score")

Cytotoxicity2_df$`Signature score` <- as.numeric(as.character(Cytotoxicity2_df$`Signature score`))
celltypes <- as.character(unique(Cytotoxicity2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Cytotoxicity_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Cytotoxicity2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("cytotoxicity Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-1.5, 2))
dev.off()

Inflammatory2_df$`Signature score` <- as.numeric(as.character(Inflammatory2_df$`Signature score`))
celltypes <- as.character(unique(Inflammatory2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Inflammatory_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Inflammatory2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Inflammatory Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-1, 1.8))
dev.off()

Stress2_df$`Signature score` <- as.numeric(as.character(Stress2_df$`Signature score`))
celltypes <- as.character(unique(Stress2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Stress_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Stress2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Stress Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-0.5, 1.5))
dev.off()

Cytokines2_df$`Signature score` <- as.numeric(as.character(Cytokines2_df$`Signature score`))
celltypes <- as.character(unique(Cytokines2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Cytokines_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Cytokines2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Cytokines Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-1, 2))
dev.off()

Inhibitory2_df$`Signature score` <- as.numeric(as.character(Inhibitory2_df$`Signature score`))
celltypes <- as.character(unique(Inhibitory2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Inhibitory_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Inhibitory2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Inhibitory Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-1, 1))
dev.off()

Proliferation2_df$`Signature score` <- as.numeric(as.character(Proliferation2_df$`Signature score`))
celltypes <- as.character(unique(Proliferation2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Proliferation_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Proliferation2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Proliferation Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-0.5, 1.5))
dev.off()

Glycolysis_REACTOME2_df$`Signature score` <- as.numeric(as.character(Glycolysis_REACTOME2_df$`Signature score`))
celltypes <- as.character(unique(Glycolysis_REACTOME2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Glycolysis_REACTOME_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Glycolysis_REACTOME2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Glycolysis_REACTOME Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-0.3, 0.3))
dev.off()

Glycolysis_HALLMARK2_df$`Signature score` <- as.numeric(as.character(Glycolysis_HALLMARK2_df$`Signature score`))
celltypes <- as.character(unique(Glycolysis_HALLMARK2_df$CellType))
my_comparisons <- combn(celltypes, 2, simplify = FALSE)
pdf("Glycolysis_HALLMARK_signature_score_celltype.pdf", width=5, height=5)
ggviolin(Glycolysis_HALLMARK2_df, x = "CellType", y = "Signature score", fill = "CellType",
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Glycolysis_HALLMARK Signature Score") + xlab("") +
  theme(title=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=30,size=10,color="black", vjust=1, hjust=1),
        axis.text.y=element_text(size=10,color="black"),
        axis.line=element_line(linewidth=0.5, color="black")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(colour="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5)) +
  ylim(c(-0.2, 0.2))
dev.off()


# By group within cell clusters and all merged
Glycolysis_REACTOME3df <- scRNA.pro@meta.data[,c("celltype","group","Glycolysis_REACTOME1")]
colnames(Glycolysis_REACTOME3df) <- c("CellType","Group","Signature score")

pdf("Glycolysis_REACTOME_signature_score_ZNF205vsVector.pdf",width=5,height=5)
ggviolin(Glycolysis_REACTOME3df, x = "CellType", y = "Signature score", color = "Group", palette=c("#4575B4","#D73027"),
         add = "boxplot", add.params = list(fill = "white"))+ylab("Glycolysis_REACTOME Signature Score")+xlab("")+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=30,size=10,color="black", vjust = 1, hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5),plot.title = element_text(hjust = 0.5))+stat_compare_means(aes(group=Group), label = "p.signif",label.y=0.25)+ylim(c(-0.2,0.3))
dev.off()

Glycolysis_HALLMARK3df <- scRNA.pro@meta.data[,c("celltype","group","Glycolysis_HALLMARK1")]
colnames(Glycolysis_HALLMARK3df) <- c("CellType","Group","Signature score")

pdf("Glycolysis_HALLMARK_signature_score_ZNF205vsVector.pdf",width=5,height=5)
ggviolin(Glycolysis_HALLMARK3df, x = "CellType", y = "Signature score", color = "Group", palette=c("#4575B4","#D73027"),
         add = "boxplot", add.params = list(fill = "white"))+ylab("Glycolysis_HALLMARK Signature Score")+xlab("")+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=30,size=10,color="black", vjust = 1, hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5),plot.title = element_text(hjust = 0.5))+stat_compare_means(aes(group=Group), label = "p.signif",label.y=0.2)+ylim(c(-0.1,0.2))
dev.off()

# By group merged
Glycolysis_REACTOME4df <- scRNA.pro@meta.data[,c("group","Glycolysis_REACTOME1")]
colnames(Glycolysis_REACTOME4df) <- c("Group","Signature score")

pdf("Glycolysis_REACTOME_signature_score_merged_ZNF205vsVector.pdf",width=2.2,height=3)
ggviolin(Glycolysis_REACTOME4df, x = "Group", y = "Signature score", color = "Group", palette = c("#4575B4","#D73027"),
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Glycolysis_REACTOME Signature Score") + xlab("") +
  theme(
    title = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 30, size = 10, color = "black", vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black")
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  stat_compare_means(
    label = "p.signif",
    method = "wilcox.test",
    label.y = 0.25
  ) +
  ylim(c(-0.2, 0.3))
dev.off()

write.csv(
  aggregate(`Signature score` ~ Group, Glycolysis_REACTOME4df, mean),
  "Glycolysis_REACTOME_score_group_mean.csv",
  row.names = FALSE
)

Glycolysis_HALLMARK4df <- scRNA.pro@meta.data[,c("group","Glycolysis_HALLMARK1")]
colnames(Glycolysis_HALLMARK4df) <- c("Group","Signature score")


pdf("Glycolysis_HALLMARK_signature_score_merged_ZNF205vsVector.pdf",width=2.2,height=3)
ggviolin(Glycolysis_HALLMARK4df, x = "Group", y = "Signature score", color = "Group", palette = c("#4575B4","#D73027"),
         add = "boxplot", add.params = list(fill = "white")) +
  ylab("Glycolysis_HALLMARK Signature Score") + xlab("") +
  theme(
    title = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 30, size = 10, color = "black", vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black")
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  stat_compare_means(
    label = "p.signif",
    method = "wilcox.test",
    label.y = 0.18
  ) +
  ylim(c(-0.1, 0.2))
dev.off()
write.csv(
  aggregate(`Signature score` ~ Group, Glycolysis_HALLMARK4df, mean),
  "Glycolysis_HALLMARK_score_group_mean.csv",
  row.names = FALSE
)

# Within C1 comparison
Glycolysis_REACTOME3df <- scRNA.pro@meta.data[,c("celltype","group","Glycolysis_REACTOME1")]
colnames(Glycolysis_REACTOME3df) <- c("CellType","Group","Signature score")
C1_df <- Glycolysis_REACTOME3df[Glycolysis_REACTOME3df$CellType == "C1", ]
pdf("Glycolysis_REACTOME_C1_signature_score_ZNF205vsVector.pdf",width=2.2,height=3)
ggviolin(C1_df, x = "CellType", y = "Signature score", color = "Group", palette=c("#4575B4","#D73027"),
         add = "boxplot",
         width = 0.8, add.params = list(fill = "white"))+ylab("Glycolysis_REACTOME Signature Score")+xlab("")+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=30,size=10,color="black", vjust = 1, hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5),plot.title = element_text(hjust = 0.5))+stat_compare_means(aes(group=Group), label = "p.signif",label.y=0.25)+ylim(c(-0.2,0.3))
dev.off()

Glycolysis_HALLMARK3df <- scRNA.pro@meta.data[,c("celltype","group","Glycolysis_HALLMARK1")]
colnames(Glycolysis_HALLMARK3df) <- c("CellType","Group","Signature score")
C1_df <- Glycolysis_HALLMARK3df[Glycolysis_HALLMARK3df$CellType == "C1", ]
pdf("Glycolysis_HALLMARK_C1_signature_score_ZNF205vsVector.pdf",width=2.2,height=3)
ggviolin(C1_df, x = "CellType", y = "Signature score", color = "Group", palette=c("#4575B4","#D73027"),
         add = "boxplot",
         width = 0.8, add.params = list(fill = "white"))+ylab("Glycolysis_HALLMARK Signature Score")+xlab("")+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=30,size=10,color="black", vjust = 1, hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5),plot.title = element_text(hjust = 0.5))+stat_compare_means(aes(group=Group), label = "p.signif",label.y=0.15)+ylim(c(-0.1,0.18))
dev.off()