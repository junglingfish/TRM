setwd('D:/ZJU-FISH/doctor/TRM/01scRNA/')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'

dir.create('results')
options(stringsAsFactors = F,check.bounds = F)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(scRNAtoolVis)
#
dir_name=list.dirs('GSE196756_RAW/',full.names = F,recursive = F)
dir_name
datalist=list()
#读取10x的数据创建CreateSeuratObject对象
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE196756_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  #细胞增加标签
  colnames(my.data)=paste0(dir_name[i],colnames(my.data))
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], min.cells = 3, min.features = 250)
  datalist[[i]]$Samples=dir_name[i]
  datalist[[i]]$type=substr(dir_name[i],1,1)
}
names(datalist)=dir_name
#批量计算线粒体和rRNA的含量
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}
#合并所有的数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#细胞数的统计
raw_cell=sce@meta.data
raw_count <- table(raw_cell$Samples)
raw_count
sum(raw_count)#12554

pearplot_befor<-VlnPlot(sce,group.by ='Samples', 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                        pt.size = 0, 
                        ncol = 4)
pearplot_befor
ggsave('results/pearplot_befor.pdf',pearplot_befor,height = 5,width = 8)
ggsave('results/pearplot_befor.jpg',pearplot_befor,height = 5,width = 10,dpi = 300)

#样本的颜色
sample_color<-pal_nejm(alpha = 0.5)(8)[1:4]
sample_color
Feature_ber1<-FeatureScatter(sce,feature1 = 'nFeature_RNA',
                             feature2 = 'nCount_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber2<-FeatureScatter(sce,feature1 = 'percent.mt',
                             feature2 = 'nCount_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber3<-FeatureScatter(sce,feature1 = 'percent.mt',
                             feature2 = 'nFeature_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber1=Feature_ber1+theme(legend.position = 'none')
Feature_ber2=Feature_ber2+theme(legend.position = 'none')

Feature_ber<-ggarrange(Feature_ber1,Feature_ber2,Feature_ber3,ncol = 3,nrow = 1,widths = c(1,1,1.2))
ggsave('results/Feature_cor.pdf',Feature_ber,height = 3,width = 10)
ggsave('results/Feature_cor.jpg',Feature_ber,height = 5,width = 17,dpi = 300)

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset =  
              nFeature_RNA < 5000 & 
              percent.mt < 15)
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_cell=sce@meta.data

clean_count <- table(clean_cell$Samples)
clean_count
sum(clean_count)
pearplot_after <- VlnPlot(sce,group.by ='Samples', 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                          pt.size = 0, 
                          ncol = 4)
pearplot_after
ggsave('results/pearplot_after.pdf',Feature_ber,height = 4,width = 10)
ggsave('results/pearplot_after.jpg',Feature_ber,height = 5,width = 15,dpi = 300)

#保存数据
save(datalist,file = 'datalist.RData')

load("datalist.RData")
#标准化
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
#筛选高变基因
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,#筛选前2000个高变，可修改的
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
#对全部的数据进行scale
sce <- ScaleData(sce, features =  rownames(sce))
#PCA降维，将高维降到低维
sce <- RunPCA(sce, features = VariableFeatures(sce)) 

elbowplot <- ElbowPlot(sce, ndims=50, reduction="pca") 
elbowplot
ggsave('results/elbowplot.pdf',height = 5,width = 5)

#可修改，选择合适的PC进行后续的聚类和降维
Dims <- 26
# sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
sce <- RunTSNE(sce,dims=1:Dims, reduction="pca")
# raw.umap<-DimPlot(sce,group.by='Samples',
#                   reduction="umap",
#                   label = "T", 
#                   pt.size = 0.2,
#                   label.size = 0)+
#   ggtitle('')
# raw.umap
# ggsave('results/raw.umap.pdf',raw.umap,height = 7,width = 7)

raw.tsne<-DimPlot(sce,group.by='Samples',
                  reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 0)+
  ggtitle('')
raw.tsne
ggsave('results/raw.tsne.pdf',raw.tsne,height = 6,width = 6)
#聚类
library(clustree)
sce <- FindNeighbors(sce, dims = 1:Dims)
sce <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1,.1))
)
colnames(sce@meta.data)
clustree(sce@meta.data, prefix = "RNA_snn_res.")

pdf('results/clust.snn_res.pdf',he=10,wi=10)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
dev.off()
#注释TRM细胞
#聚类分析
Resolution <- 0.8
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#亚群注释
DefaultAssay(sce) <- "RNA"

VlnPlot(sce,features = c('CD8A', 'CD69', 'ITGAE', 'ITGA1'),pt.size = 0,group.by = 'seurat_clusters',ncol = 2)
library(randomcoloR)
allcolour <- c(pal_npg(alpha = 0.8)(9),
               pal_igv(alpha = 0.8)(9),
               pal_jama(alpha = 0.8)(7),
               pal_jco(alpha = 0.8)(9),
               pal_nejm(alpha = 0.8)(8))
length(table(sce@active.ident))
#32
mycolor1 = allcolour[1:length(table(sce$seurat_clusters))]

figs2b<-FeaturePlot(sce,
                    features =  c('CD8A', 'CD69', 'ITGAE', 'ITGA1'),
                    pt.size = 0.15,reduction = 'tsne',ncol = 2)
figs2a<-DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
                reduction="tsne",
                label = "T", 
                pt.size = 0.2,
                label.size = 4) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

figs2ab<-ggarrange(figs2a,figs2b,nrow = 1,ncol = 2,widths = c(1,1))
figs2ab
table(sce$seurat_clusters)
save(sce,file = 'sce1.RData')

load('sce1.RData')
DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
        reduction="tsne",
        label = "T", 
        pt.size = 0.4,
        label.size = 3) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
#提取9,11,21,重新聚类
load('sce1.RData')
Idents(sce)='seurat_clusters'
sce<-subset(sce,idents =c(1, 12, 22))
#二次聚类
Resolution <- 0.6
DefaultAssay(sce) <- "RNA"

sce <- FindNeighbors(object = sce, dims = 1:10)
sce <- FindClusters(object = sce, resolution = Resolution)
DefaultAssay(sce) <- "RNA"

p <- VlnPlot(sce,features = c('CD8A', 'CD69', 'ITGAE', 'ITGA1'),pt.size = 0,group.by = 'seurat_clusters',ncol = 4)
ggsave('results/violinplot_4genes_TRM.pdf',p,height = 4,width = 10)


#重新降维
sce <- RunUMAP(sce, 
               dims=1:10, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)
sce <- RunTSNE(sce,
               dims=1:10,
               reduction="pca",
               perplexity=10,
               max_iter=1000)
figs2c<-DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
                reduction="umap",
                label = "T", 
                pt.size = 0.8,
                label.size = 3) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

figs2c <- DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
                  reduction="tsne",
                  label = "T",
                  pt.size = 0.8,
                  label.size = 6) +
  theme(axis.line = element_line(size=0.1, colour = "black"),
        #axis.text = element_blank(),
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

figs2d<-FeaturePlot(sce,
                    features =  c('CD8A', 'CD69', 'ITGAE', 'ITGA1'),
                    pt.size = 0.3,reduction = 'tsne',ncol = 2)
fig2cd<-ggarrange(figs2c,figs2d,nrow = 1,ncol = 2,widths = c(1,1))
fig2cd

ggsave('results/all_cluster.pdf', figs2ab, height = 4.5, width = 11)
ggsave('results/TRM_cluster.pdf', fig2cd, height = 4.5, width = 11)

figs2<-ggarrange(figs2ab,fig2cd,nrow = 2,ncol = 1)
ggsave('results/FigS3.pdf',figs2,height =15,width = 15)
FeaturePlot(sce,
            features =  c('CD8A', 'CD69', 'ITGAE', 'ITGA1'),
            pt.size = 0.1,reduction = 'tsne',ncol = 2)
#寻找差异基因时的差异倍数
Logfc = 0.5
#差异基因时最小的表达比例
Minpct = 0.35
DefaultAssay(sce) <- "RNA"
Idents(sce)<-'seurat_clusters'
sce<-JoinLayers(sce)
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct, only.pos = T)
sce.markers.all <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct, only.pos = F)
# make a polar plot
vol <- jjVolcano(diffData = sce.markers.all,
                 tile.col = corrplot::COL2('RdBu', 15)[4:12],
                 size  = 3.5,
                 fontface = 'italic',
                 polar = T)
ggsave("results/round_volocano.pdf", plot = vol, width = 8, height = 8)

sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)
write.table(sce.markers,'results/scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')
# write.table(sce.markers,'results/scRNA_marker_gene_TF.txt',quote = F,row.names = F,sep='\t')


### 选择前5个marker基因
sce$celltype <- Idents(sce) 
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)  

Top5 <- intersect(unique(Top5$gene),rownames(sce))
# 确保 Top5 的基因顺序，KRT13 放在第十和第十五个位置
Top5 <- c(setdiff(Top5, "KRT13")[1:9], "KRT13", setdiff(Top5, "KRT13")[-(1:9)])

# 创建 DotPlot
sc_marker_dotplot <- DotPlot(object = sce, features = Top5, 
                             group.by = "celltype", cols = c("#ffffff", "#9b0900"), scale = TRUE) +
  RotatedAxis() + 
  ggtitle("Top 5 Marker Genes") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +  # 黑色边框
  xlab('Genes') +
  ylab('Clusters')
sc_marker_dotplot
ggsave('results/sc_marker_dotplot.pdf',sc_marker_dotplot,height = 6,width = 10)

dot_plot <- jjDotPlot(object = sce,
                          gene = Top5,
                          anno = T,
                          plot.margin = c(3,1,1,1),
                          tree.pos = 'left',
                          same.pos.label = T,
                          yPosition = 15, 
                          cluster.order = c('0', '1', '2', '3', '4', '5'),
                          ytree = F)
ggsave("results/Top5_genes_dotplot_vision.pdf", plot = dot_plot, width = 14, height = 9)

library(Seurat)
library(ggplot2) # 确保加载 ggplot2 来调整热图的颜色
# 绘制热图并调整颜色及其他参数
heatmap_plot <- DoHeatmap(object = sce,
                          features = Top5,
                          label = TRUE,       # 显示标签
                          slot = "scale.data",
                          size = 5,          # 标签字体大小
                          angle = 45,        # 标签旋转角度
                          draw.lines = TRUE)  # 绘制分隔线

# 使用 scale_fill_gradientn 自定义颜色
heatmap_plot <- heatmap_plot + scale_fill_gradientn(colors = c("#4a79b6", "#e0f0f4", "#db4b43"))
ggsave("results/Top5_genes_heatmap_plot.pdf", plot = heatmap_plot, width = 7, height = 5)


#第二种
bubble.df=as.matrix(GetAssayData(sce, assay = "RNA", slot = "data")[Top5, ])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
bubble.df$CB=rownames(bubble.df)
bubble.df=merge(bubble.df,
                data.frame(CB=rownames(sce@meta.data),
                           celltype=sce@meta.data$seurat_clusters),
                by = "CB")
bubble.df$CB=NULL

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$celltype)) {
  bubble.df_small=bubble.df%>%filter(celltype==i)
  for (j in Top5) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    celltype_v=append(celltype_v,i)
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}
plotdf=data.frame(
  celltype=celltype_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v)
plotdf$celltype=factor(plotdf$celltype,levels = unique(as.character(sce.markers$cluster)))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(Top5)))
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
sc_marker_dotplot1<-plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
sc_marker_dotplot1
ggsave('results/sc_marker_dotplot1.pdf',sc_marker_dotplot1,height = 7,width = 9)

#绘图
save(sce,file = 'sce.rdata')
load('sce.rdata')
mycolor =pal_npg('nrc')(9)

fig1a = DimPlot(sce,group.by = 'Samples',
                reduction="tsne",
                label = "F", 
                pt.size = 1.5,
                label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')+guides(colour = guide_legend(ncol = 1))

fig1a

fig1b<-DimPlot(sce,cols=mycolor,group.by = 'seurat_clusters',
               reduction="tsne",split.by = 'type',
               label = "F", 
               pt.size = 1.5,
               label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
fig1b
ggsave('results/cluster_to_tumornormal.pdf',fig1b,height = 5,width = 7)
#
Idents(sce)='seurat_clusters'
library("ggplot2")
sample_clust<-as.matrix(table(sce$Samples,sce$seurat_clusters))
sample_clust=apply(sample_clust,1,function(x){return(x/sum(x))})
sample_clust=reshape2::melt(sample_clust)
colnames(sample_clust)<-c("cluster","Samples","proportion")
sample_clust$cluster=paste0('TRM_',sample_clust$cluster)
write.table(sample_clust,'results/sample_clust1.txt',quote = F,row.names = T,sep='\t')

clust_freq<-as.data.frame(table(sce$Samples))
colnames(clust_freq)=c('Samples','cell_num')
clust_freq=clust_freq[order(clust_freq$cell_num,decreasing = T),]
clust_freq$Samples=factor(clust_freq$Samples,levels = clust_freq$Samples)
sample_clust$Samples=factor(sample_clust$Samples,levels =clust_freq$Samples)

fig1e1<-ggplot(sample_clust,aes(x = Samples,y = proportion,fill=cluster))+
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +scale_fill_manual(values = c('#Fa7F72', '#FBB463', '#FBF8B4', "#80B1D3", '#8DD1C6', "#BDBADB"))+
  theme_bw() + 
  theme(axis.ticks.length = unit(0.1, 'cm'),
        legend.position = "left") +xlab('')+
  coord_flip()+scale_y_continuous(expand = expand_scale(mult = c(0, 0)))
fig1e1

fig1e2 <- ggplot(clust_freq, aes(x = Samples, y = cell_num, fill = Samples)) +
  geom_bar(stat = "identity") +
  ggtitle("") +
  theme_bw() + 
  scale_fill_manual(values = c('#DE582B', "#018A67", '#1868B2', '#F3A332')) +
  theme(axis.ticks.length = unit(0, 'cm'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip() +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0)), 
                     breaks = c(0, 500, 1000),
                     limits = c(0, max(clust_freq$cell_num) + 80))  # Ensure space for 1000 tick

fig1e2


fig1e3<-ggpubr::ggarrange(fig1e1,fig1e2,nrow = 1,ncol = 2,widths = c(2,1))
fig1e3
ggsave('results/cluster_to_cellfrom.pdf',fig1e3,height = 4,width = 5)

library(clusterProfiler)
library(org.Hs.eg.db)

# marker基因进行注释
ids = bitr(sce.markers$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db') # 将SYMBOL转成ENTREZID
sce.markers2 = merge(sce.markers, ids, by.x = 'gene', by.y = 'SYMBOL')

# 按照分组因子分组
gcSample = split(sce.markers2$ENTREZID, sce.markers2$cluster)

# KEGG富集分析
sce.markers2.enrich.res <- compareCluster(gcSample,
                                          fun = "enrichKEGG",
                                          organism = "hsa", pvalueCutoff = 0.05)

# 获取富集分析的摘要
summary_result <- summary(sce.markers2.enrich.res)

# 绘制dotplot并设置颜色渐变
fig1f <- dotplot(sce.markers2.enrich.res, showCategory = 5, # 设置显示的类别数量
                 color = "p.adjust") + 
  scale_color_gradient(low = "blue", high = "red") +  # 修改颜色渐变
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
# 显示图形
print(fig1f)
ggsave('results/KEGG_cluster.pdf',fig1f,height = 6,width = 5)

save(sce,file = 'sce2.RData')
#恶性和非恶性的区别
load('sce2.RData')
library(copykat)

copykat <- function (rawmat = rawdata, id.type = "S", cell.line = "no", 
                     ngene.chr = 0, LOW.DR = 0.05, UP.DR = 0.1, win.size = 25, 
                     norm.cell.names = "", KS.cut = 0.1, sam.name = "", distance = "euclidean", 
                     n.cores = 1) {
  start_time <- Sys.time()
  set.seed(1)
  sample.name <- paste(sam.name, "_copykat_", sep = "")
  print("running copykat v1.0.4")
  print("step1: read and filter data ...")
  rawmat <- GetAssayData(sce, assay = "RNA", slot = "counts")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", 
              sep = ""))
  # genes.raw <- apply(rawmat, 2, function(x) (sum(x > 0)))
  # if (sum(genes.raw > 200) == 0) 
  #   stop("none cells have more than 200 genes")
  # if (sum(genes.raw < 100) > 1) {
  #   rawmat <- rawmat[, -which(genes.raw < 200)]
  #   print(paste("filtered out ", sum(genes.raw <= 200), 
  #               " cells with less than 200 genes; remaining ", ncol(rawmat), 
  #               " cells", sep = ""))
  # }
  der <- apply(rawmat, 1, function(x) (sum(x > 0)))/ncol(rawmat)
  if (sum(der > LOW.DR) >= 1) {
    rawmat <- rawmat[which(der > LOW.DR), ]
    print(paste(nrow(rawmat), " genes past LOW.DR filtering", 
                sep = ""))
  }
  WNS1 <- "data quality is ok"
  if (nrow(rawmat) < 7000) {
    WNS1 <- "low data quality"
    UP.DR <- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  print("step 2: annotations gene coordinates ...")
  anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type)
  anno.mat <- anno.mat[order(anno.mat$abspos, decreasing = FALSE), 
  ]
  HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
  toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
  # if (length(toRev) > 0) {
  #   anno.mat <- anno.mat[-toRev, ]
  # }
  # ToRemov2 <- NULL
  # for (i in 8:ncol(anno.mat)) {
  #   cell <- cbind(anno.mat$chromosome_name, anno.mat[, i])
  #   cell <- cell[cell[, 2] != 0, ]
  #   if (length(as.numeric(cell)) < 5) {
  #     rm <- colnames(anno.mat)[i]
  #     ToRemov2 <- c(ToRemov2, rm)
  #   }
  #   else if (length(rle(cell[, 1])$length) < 23 | min(rle(cell[, 
  #                                                              1])$length) < ngene.chr) {
  #     rm <- colnames(anno.mat)[i]
  #     ToRemov2 <- c(ToRemov2, rm)
  #   }
  #   i <- i + 1
  # }
  # if (length(ToRemov2) == (ncol(anno.mat) - 7)) 
  #   stop("all cells are filtered")
  # if (length(ToRemov2) > 0) {
  #   anno.mat <- anno.mat[, -which(colnames(anno.mat) %in% 
  #                                   ToRemov2)]
  # }
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat <- log(sqrt(rawmat3) + sqrt(rawmat3 + 1))
  norm.mat <- apply(norm.mat, 2, function(x) (x <- x - mean(x)))
  colnames(norm.mat) <- colnames(rawmat3)
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c) {
    model <- dlm::dlmModPoly(order = 1, dV = 0.16, dW = 0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x <- x[2:length(x)]
    x <- x - mean(x)
  }
  test.mc <- parallel::mclapply(1:ncol(norm.mat), dlm.sm, 
                                mc.cores = 1)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), 
                            byrow = FALSE)
  colnames(norm.mat.smooth) <- colnames(norm.mat)
  print("step 4: measuring baselines ...")
  if (cell.line == "yes") {
    print("running pure cell line mode")
    relt <- baseline.synthetic(norm.mat = norm.mat.smooth, 
                               min.cells = 10, n.cores = 1)
    norm.mat.relat <- relt$expr.relat
    CL <- relt$cl
    WNS <- "run with cell line mode"
    preN <- NULL
  }
  else if (length(norm.cell.names) > 1) {
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% 
                                                    norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset", 
                sep = ""))
    if (NNN == 0) 
      stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")
    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% 
                                             norm.cell.names)], 1, median)
    print("baseline is from known input")
    d <- parallelDist::parDist(t(norm.mat.smooth), threads = 1, 
                               method = "euclidean")
    km <- 6
    fit <- hclust(d, method = "ward.D2")
    CL <- cutree(fit, km)
    while (!all(table(CL) > 5)) {
      km <- km - 1
      CL <- cutree(fit, k = km)
      if (km == 2) {
        break
      }
    }
    WNS <- "run with known normal"
    preN <- norm.cell.names
    norm.mat.relat <- norm.mat.smooth - basel
  }
  else {
    basa <- baseline.norm.cl(norm.mat.smooth = norm.mat.smooth, 
                             min.cells = 5, n.cores = 1)
    basel <- basa$basel
    WNS <- basa$WNS
    preN <- basa$preN
    CL <- basa$cl
    if (WNS == "unclassified.prediction") {
      Tc <- colnames(rawmat)[which(as.numeric(apply(rawmat[which(rownames(rawmat) %in% c("PTPRC", "LYZ", "PECAM1")), ], 2, mean)) > 1)]
      length(Tc)
      preN <- intersect(Tc, colnames(norm.mat.smooth))
      if (length(preN) > 5) {
        print("start manual mode")
        WNS <- paste("copykat failed in locating normal cells; manual adjust performed with ", 
                     length(preN), " immune cells", sep = "")
        print(WNS)
        basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% 
                                                 preN)], 1, mean)
      }else {
        basa <- baseline.GMM(CNA.mat = norm.mat.smooth, 
                             max.normal = 5, mu.cut = 0.05, Nfraq.cut = 0.99, 
                             RE.before = basa, n.cores = 1)
        basel <- basa$basel
        WNS <- basa$WNS
        preN <- basa$preN
      }
    }
    norm.mat.relat <- norm.mat.smooth - basel
  }
  DR2 <- apply(rawmat3, 1, function(x) (sum(x > 0)))/ncol(rawmat3)
  norm.mat.relat <- norm.mat.relat[which(DR2 >= UP.DR), ]
  anno.mat2 <- anno.mat[which(DR2 >= UP.DR), ]
  # ToRemov3 <- NULL
  # for (i in 8:ncol(anno.mat2)) {
  #   cell <- cbind(anno.mat2$chromosome_name, anno.mat2[, 
  #                                                      i])
  #   cell <- cell[cell[, 2] != 0, ]
  #   if (length(as.numeric(cell)) < 5) {
  #     rm <- colnames(anno.mat2)[i]
  #     ToRemov3 <- c(ToRemov3, rm)
  #   }
  #   else if (length(rle(cell[, 1])$length) < 23 | min(rle(cell[, 
  #                                                              1])$length) < ngene.chr) {
  #     rm <- colnames(anno.mat2)[i]
  #     ToRemov3 <- c(ToRemov3, rm)
  #   }
  #   i <- i + 1
  # }
  # if (length(ToRemov3) == ncol(norm.mat.relat)) 
  #   stop("all cells are filtered")
  # if (length(ToRemov3) > 0) {
  #   norm.mat.relat <- norm.mat.relat[, -which(colnames(norm.mat.relat) %in% 
  #                                               ToRemov3)]
  # }
  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]
  print("step 5: segmentation...")
  results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, bins = win.size, 
                      cut.cor = KS.cut, n.cores = 1)
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, 
                        bins = win.size, cut.cor = 0.5 * KS.cut, n.cores = 1)
  }
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, 
                        bins = win.size, cut.cor = 0.5 * 0.5 * KS.cut, n.cores = 1)
  }
  if (length(results$breaks) < 25) 
    stop("too few segments; try to decrease KS.cut; or improve data")
  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA, 2, function(x) (x <- x - 
                                                         mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)
  write.table(RNA.copycat, paste(sample.name, "CNA_raw_results_gene_by_cell.txt", 
                                 sep = ""), sep = "\t", row.names = FALSE, quote = F)
  print("step 6: convert to genomic bins...")
  Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat = RNA.copycat, 
                              n.cores = n.cores)
  uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
  print("step 7: adjust baseline ...")
  if (cell.line == "yes") {
    mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
    write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep = ""), 
                sep = "\t", row.names = FALSE, quote = F)
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(mat.adj),threads = 1, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), method = "ward.D")
    }
    saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", sep = ""))
    print("step 8: ploting heatmap ...")
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
    chr <- as.numeric(Aj$DNA.adj$chrom)%%2 + 1
    rbPal1 <- colorRampPalette(c("black", "grey"))
    CHR <- rbPal1(2)[as.numeric(chr)]
    chr1 <- cbind(CHR, CHR)
    if (ncol(mat.adj) < 3000) { h <- 10}else {h <- 15}
    col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                   seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
    if (distance == "euclidean") {
      jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
           height = h * 250, width = 4000, res = 100)
      heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x,threads = 1, method = distance), hclustfun = function(x) hclust(x,method = "ward.D"), ColSideColors = chr1, Colv = NA, 
                Rowv = TRUE, notecol = "black", col = my_palette, 
                breaks = col_breaks, key = TRUE, keysize = 1, 
                density.info = "none", trace = "none", cexRow = 0.1, 
                cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, 
                symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 10))
      dev.off()
    }
    else {
      jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
           height = h * 250, width = 4000, res = 100)
      heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 -cor(t(x), method = distance)), hclustfun = function(x) hclust(x,  method = "ward.D"), ColSideColors = chr1, Colv = NA, 
                Rowv = TRUE, notecol = "black", col = my_palette, 
                breaks = col_breaks, key = TRUE, keysize = 1, 
                density.info = "none", trace = "none", cexRow = 0.1, 
                cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, 
                symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), cex.main = 4, margins = c(10,10))
      dev.off()
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
    names(reslts) <- c("CNAmat", "hclustering")
    return(reslts)
  }
  else {
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(uber.mat.adj), threads = 1, method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(uber.mat.adj, method = distance)), method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)
    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.pred[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "nondiploid"
    names(com.pred) <- names(hc.umap)
    results.com.rat <- uber.mat.adj - apply(uber.mat.adj[, which(com.pred == "diploid")], 1, mean)
    results.com.rat <- apply(results.com.rat, 2, function(x) (x <- x-mean(x)))
    results.com.rat.norm <- results.com.rat[, which(com.pred == "diploid")]
    dim(results.com.rat.norm)
    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)
    adjN <- function(j) {a <- results.com.rat[, j]
    a[abs(a - base) <= 0.25 * cf.h] <- mean(a)
    a
    }
    mc.adjN <- parallel::mclapply(1:ncol(results.com.rat),  adjN, mc.cores = n.cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
    colnames(adj.results) <- colnames(results.com.rat)
    rang <- 0.5 * (max(adj.results) - min(adj.results))
    mat.adj <- adj.results/rang
    print("step 8: final prediction ...")
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(mat.adj), threads = n.cores, method = distance),
                    method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)
    saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", 
                              sep = ""))
    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.preN[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
    names(com.preN) <- names(hc.umap)
    if (WNS == "unclassified.prediction") {
      com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
      com.preN[which(com.preN == "nondiploid")] <- "c2:aneuploid:low.conf"
    }
    print("step 9: saving results...")
    res <- cbind(names(com.preN), com.preN)
    colnames(res) <- c("cell.names", "copykat.pred")
    write.table(res, paste(sample.name, "prediction.txt", 
                           sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name,"CNA_results.txt", sep = ""), sep = "\t", row.names = FALSE, quote = F)
    print("step 10: ploting heatmap ...")
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
    chr <- as.numeric(Aj$DNA.adj$chrom)%%2 + 1
    rbPal1 <- colorRampPalette(c("black", "grey"))
    CHR <- rbPal1(2)[as.numeric(chr)]
    chr1 <- cbind(CHR, CHR)
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]
    cells <- rbind(compreN_pred, compreN_pred)
    if (ncol(mat.adj) < 3000) {
      h <- 10
    }
    else {
      h <- 15
    }
    col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, 
                                                   -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                   seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
    if (distance == "euclidean") {
      jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
           height = h * 250, width = 4000, res = 100)
      heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x, threads = n.cores, method = distance), hclustfun = function(x) hclust(x,method = "ward.D"), ColSideColors = chr1, RowSideColors = cells,  Colv = NA, Rowv = TRUE, notecol = "black", col = my_palette,  breaks = col_breaks, key = TRUE, keysize = 1, ensity.info = "none", trace = "none", cexRow = 0.1, cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), cex.main = 4, margins = c(10,10))
      legend("topright", paste("pred.", names(table(com.preN)), sep = ""), pch = 15, col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex = 1)
      dev.off()
    }
    else {
      jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
           height = h * 250, width = 4000, res = 100)
      heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 - 
                                                                              cor(t(x), method = distance)), hclustfun = function(x) hclust(x,method = "ward.D"), ColSideColors = chr1, RowSideColors = cells, 
                Colv = NA, Rowv = TRUE, notecol = "black", col = my_palette, 
                breaks = col_breaks, key = TRUE, keysize = 1, 
                density.info = "none", trace = "none", cexRow = 0.1, 
                cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, 
                symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                cex.main = 4, margins = c(10, 10))
      legend("topright", paste("pred.", names(table(com.preN)),sep = ""), 
             pch = 15, col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex = 1)
      dev.off()
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), 
                   hcc)
    names(reslts) <- c("prediction", "CNAmat", "hclustering")
    return(reslts)
  }
}
copykat.test <- copykat(rawmat=sce@assays$RNA@counts,
                        id.type="S",
                        cell.line="no",
                        ngene.chr=5,
                        #每个染色体中至少有 5 个基因来计算 DNA 拷贝数
                        win.size=25,
                        #每个片段至少取 25 个基因
                        KS.cut=0.15,
                        #0-1,值越大灵敏度越低
                        sam.name="LUAD",
                        #随意固定一个名称
                        distance="euclidean",
                        n.cores=1
                        #并行计算
)
save(copykat.test,file = 'TRM_copykat.test.RData')

#读取CNV的结果
copykat.test<-read.delim('LUAD_copykat_prediction.txt',sep='\t',header = T)
head(copykat.test)
table(copykat.test$copykat.pred)
rownames(copykat.test)=copykat.test$cell.names
copykat.test=copykat.test[rownames(sce@meta.data),]
#添加分组
sce <- AddMetaData(sce, copykat.test$copykat.pred,col.name = "copykat.pred")
sce$copykat.pred[is.na(sce$copykat.pred)]<-'Unknown'
table(sce$copykat.pred)
# aneuploid   diploid 
#347       241
sce$copykat.pred=ifelse(sce$copykat.pred=='aneuploid','malignant','no_malignant')

save(sce,file = 'sce2.RData')
fig1h<-DimPlot(sce,cols=c('red','blue'),group.by = 'copykat.pred',
               reduction="tsne",
               label = "F", 
               pt.size = 0.5,
               label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
fig1h
ggsave(filename = 'results/malignant.pdf',plot = fig1h,he=5,wi=6)


fig1ef<-ggarrange(fig1e3,fig1f,fig1h,labels = c('D','E','F'),nrow = 1,ncol = 3,widths = c(1.2,1.3,1))

fig1ab<-ggarrange(fig1a,fig1b,nrow = 1,ncol=2,labels = c('A','B'),widths = c(1,1.5))
fig1=ggarrange(fig1ab,sc_marker_dotplot,fig1ef,labels = c('','C',''),nrow = 3,ncol = 1,heights = c(1.5,1,1.5))

ggsave(filename = 'results/FigS1.pdf',plot = fig1,he=15,wi=18)
ggsave(filename = 'results/FigS1.jpg',plot = fig1,he=15,wi=18)