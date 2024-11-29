setwd('D:/ZJU-FISH/doctor/TRM/06GSVA')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'

dir.create('results')

tcga.dat<-read.delim(paste0(routine_data, 'results/tumor_exp_cleaned.txt'),sep='\t',header = T,row.names = 1,check.names = F)
tcga.module<-read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)
colnames(tcga.module)
tcga.module.gene.exp=tcga.module[,3:10]
rownames(tcga.module.gene.exp)=tcga.module$Samples
library(GSVA)
library(GSEABase)
gmtFile='c2.cp.kegg.v7.5.1.symbols.gmt'
c2KEGG <- getGmt(gmtFile,
                 collectionType=BroadCollection(category="h"),
                 geneIdType=SymbolIdentifier())
gsvapar <- gsvaParam(as.matrix(tcga.dat), c2KEGG, maxDiff=TRUE)
tcga.kegg <- gsva(gsvapar)
# tcga.kegg <- gsva(as.matrix(tcga.dat),
#                   c2KEGG,
#                   method = 'ssgsea',
#                   min.sz = 10,
#                   max.sz = 500,
#                   verbose = TRUE)
save(tcga.kegg,file = 'tcga.kegg.Rdata')
load('tcga.kegg.Rdata')
tcga.kegg[1:4,1:4]
write.table(tcga.kegg,'results/tcga.kegg.enrichscore.txt',quote = F,row.names = T,sep='\t')
gene.kegg.dat<-cbind.data.frame(tcga.module.gene.exp,
                                t(tcga.kegg[,rownames(tcga.module.gene.exp)]))
gene.kegg.dat[1:4,1:7]
gene.kegg.cor <- Hmisc::rcorr(as.matrix(gene.kegg.dat))
gene.kegg.cor.r=reshape2::melt(gene.kegg.cor$r)
gene.kegg.cor.p=reshape2::melt(gene.kegg.cor$P)
colnames(gene.kegg.cor.r)=c('gene','pathway','cor')
colnames(gene.kegg.cor.p)=c('gene','pathway','p')
gene.kegg.cor.res<-merge(gene.kegg.cor.r,gene.kegg.cor.p,by=c('gene','pathway'))
head(gene.kegg.cor.res)
gene.kegg.cor.res=gene.kegg.cor.res[which(gene.kegg.cor.res$gene %in% colnames(tcga.module.gene.exp) & gene.kegg.cor.res$pathway %in% rownames(tcga.kegg)),]
head(gene.kegg.cor.res)
write.table(gene.kegg.cor.res,'results/gene.kegg.cor.res.txt',quote = F,row.names = F,sep = '\t')
gene.kegg.cor.res.fit=gene.kegg.cor.res[which(abs(gene.kegg.cor.res$cor)>0.5 & gene.kegg.cor.res$p<0.001),]
dim(gene.kegg.cor.res.fit)
sig_pathway<-as.character(unique(gene.kegg.cor.res.fit$pathway))
length(sig_pathway)#37
#相关性分析的热图
gene.kegg.cor1=gene.kegg.cor
rownames(gene.kegg.cor1$r)=gsub('KEGG_','',rownames(gene.kegg.cor1$r))
colnames(gene.kegg.cor1$r)=gsub('KEGG_','',colnames(gene.kegg.cor1$r))

rownames(gene.kegg.cor1$P)=gsub('KEGG_','',rownames(gene.kegg.cor1$P))
colnames(gene.kegg.cor1$P)=gsub('KEGG_','',colnames(gene.kegg.cor1$P))
sig_pathway1=gsub('KEGG_','',sig_pathway)


pdf('results/Fig1A.pdf',height =12,width = 9)
corrplot::corrplot(as.matrix(t(gene.kegg.cor1$r[colnames(tcga.module.gene.exp),sig_pathway1])), 
                   p.mat = as.matrix(t(gene.kegg.cor1$P[colnames(tcga.module.gene.exp),sig_pathway1])),
                   mar = c(0,0,1,1),
                   col=colorRampPalette(c("#00468b", "white","#ed0000"))(100),
                   tl.srt = 90,
                   tl.cex = 1,
                   tl.col = 'black',
                   tl.offset = 0.5,
                   cl.pos = c("b","r","n")[1], 
                   cl.align.text = 'l',
                   cl.length = 5,
                   cl.ratio = 0.1,
                   cl.cex = 0.8,
                   addgrid.col = 'white',
                   method = 'color',
                   insig = 'label_sig',
                   sig.level=c(0.001,0.01,0.05),
                   pch.cex=1,
                   is.corr=T,
                   xpd=T)

dev.off()
#热图
anno_col=as.data.frame(tcga.module.gene.exp)
colnames(anno_col)
tcga.kegg1=tcga.kegg
tcga.kegg1[1:4,1:4]
rownames(tcga.kegg1)=gsub('KEGG_','',rownames(tcga.kegg1))
cols=ggsci::pal_jco()(9)
library(pheatmap)
pdf('results/Fig1B.pdf',height = 12,width = 12)
pheatmap(mat = as.matrix(tcga.kegg1[sig_pathway1,rownames(anno_col)]),
         scale = 'row',show_colnames = F,show_rownames = T,
         annotation_col = anno_col,
         cluster_cols = T,cluster_rows =T,
         color = colorRampPalette(c("#3C5488FF", "white","#E64B35FF"))(100)
         ,breaks = unique(c(seq(-2, 2, length=100))),
         annotation_colors =data.frame(BTK=colorRampPalette(c("white",cols[1]))(20),
                                       CLEC3B=colorRampPalette(c("white",cols[2]))(20),
                                       ANLN=colorRampPalette(c("white",cols[3]))(20),
                                       CD302=colorRampPalette(c("white",cols[4]))(20),
                                       CYP4B1=colorRampPalette(c("white",cols[5]))(20),
                                       ECT2=colorRampPalette(c("white",cols[6]))(20),
                                       GRIA1=colorRampPalette(c("white",cols[7]))(20),
                                       HMMR=colorRampPalette(c("white",cols[8]))(20)))
dev.off()
