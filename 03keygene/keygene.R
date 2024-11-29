setwd('D:/ZJU-FISH/doctor/TRM/03keygene/')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_2 = 'D:/ZJU-FISH/doctor/TRM/02TRMgenecluster/'

dir.create('results')
tcga.dat<-read.delim(paste0(routine_data, 'GEO/GSE53625/GSE53625_count_loged.txt'),sep='\t',header = T,row.names = 1,check.names = F)
# tcga.group<-read.delim(paste0(routine_data, 'GEO/GSE53625/GSE53625_count_group.txt'), sep='\t',header = T)
tcga.TRM.score<-read.delim(paste0(routine_2, 'results/tcga.TRM.score.txt'),sep='\t',header = T)
tcga.cli<-read.csv(paste0(routine_data, 'clinic/GEO_clinic_data.csv'), check.names = F)
#差异分析
limma_DEG=function(exp,group,ulab,dlab){
  library(limma)
  ind1=which(group==ulab)
  ind2=which(group==dlab)
  sml <- c(rep('A',length(ind1)),rep('B',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  
  design <- model.matrix(~fl+0)
  colnames(design) <- levels(fl)
  cont.matrix<-makeContrasts(contrasts='A-B',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
  return(tT)
}
#火山图
volcano_plot=function(logfc,pvalue,cutFC=1,cutPvalue=0.05
                      ,colors=c("#BF1D2D",'grey',"#293890")
                      ,ylab='-log10(FDR)',
                      leg='Group',
                      xlab='log2(FoldChange)'){
  library(ggplot2)
  diff.name=rep('None',length(logfc))
  diff.name[which(logfc>cutFC&pvalue<cutPvalue)]='Up'
  diff.name[which(logfc< -cutFC&pvalue<cutPvalue)]='Down'
  dat.input=data.frame(logFC=logfc,FDR=pvalue,change=diff.name)
  p1 <- ggplot(data = dat.input, 
               aes(x = logFC, 
                   y = -log10(FDR)))+
    geom_point(alpha=0.6, size=1, aes(color=change))+
    scale_color_manual(values=colors,limits = c("Down",'None', "Up"),name=leg)+
    geom_vline(xintercept=c(-cutFC,cutFC),lty=2,col="black",lwd=0.8)+
    geom_hline(yintercept = -log10(cutPvalue),lty=2,col="black",lwd=0.8)+
    ylab(ylab)+xlab(xlab)+
    theme_get()+
    theme(legend.background = element_rect(fill = NA, colour = NA))
  return(p1)
}
#差异分析
tcga.diff<-limma_DEG(exp = tcga.dat[,tcga.TRM.score$Samples],
                     group = tcga.TRM.score$group,
                     ulab = 'tumor',dlab = 'normal')
fc.fit=1;p.fit=0.05
fig1a<-volcano_plot(logfc = tcga.diff$logFC,
                    pvalue = tcga.diff$adj.P.Val,
                    cutFC = fc.fit,
                    cutPvalue = p.fit,
                    colors = c('#293890','grey','#BF1D2D'),
                    ylab = '-log10(FDR)',
                    leg='',
                    xlab='log2(FoldChange)')
fig1a
ggsave('results/volcano.plot.pdf',fig1a,height = 4,width = 5)
#
write.table(tcga.diff,'results/tcga.diff.txt',quote = F,row.names = T,sep='\t')
tcga.diff.fit<-tcga.diff[which(abs(tcga.diff$logFC)>fc.fit & tcga.diff$adj.P.Val<p.fit),]
table(tcga.diff.fit$logFC>0)
#FALSE  TRUE 
#1006   725
sam_T<-tcga.TRM.score[which(tcga.TRM.score$group=='tumor'),"Samples"]
tcga.gene.exp<-t(tcga.dat[rownames(tcga.diff.fit),sam_T])
tcga.TRM.exp<-merge(data.frame(Samples=rownames(tcga.gene.exp),tcga.gene.exp),
                    tcga.TRM.score[,c("Samples","TRM_2","TRM_5")],
                    by='Samples')
write.table(tcga.TRM.exp,'results/tcga.TRM.exp.txt',quote = F,row.names = F,sep='\t')
#相关性分析
tcga.TRM.exp=tcga.TRM.exp[,-1]
tcga.TRM.gene.cor<-Hmisc::rcorr(as.matrix(tcga.TRM.exp))
tcga.TRM.gene.cor.r<-reshape2::melt(tcga.TRM.gene.cor$r)
tcga.TRM.gene.cor.p<-reshape2::melt(tcga.TRM.gene.cor$P)
colnames(tcga.TRM.gene.cor.r)=c('TRM','gene','cor')
colnames(tcga.TRM.gene.cor.p)=c('TRM','gene','p')
tcga.TRM.gene.cor.r=tcga.TRM.gene.cor.r[tcga.TRM.gene.cor.r$TRM %in% c("TRM_2","TRM_5") & tcga.TRM.gene.cor.r$gene %in% rownames(tcga.diff.fit),]
tcga.TRM.gene.cor.p=tcga.TRM.gene.cor.p[tcga.TRM.gene.cor.p$TRM %in% c("TRM_2","TRM_5") & tcga.TRM.gene.cor.p$gene %in% rownames(tcga.diff.fit),]
tcga.TRM.gene.cor.res<-merge(tcga.TRM.gene.cor.r,tcga.TRM.gene.cor.p,by=c('TRM','gene'))
head(tcga.TRM.gene.cor.res)
tcga.TRM.gene.cor.res.fit<-tcga.TRM.gene.cor.res[abs(tcga.TRM.gene.cor.res$cor)>0.4 & tcga.TRM.gene.cor.res$p<0.05,]
dim(tcga.TRM.gene.cor.res.fit)
#867
tcga.TRM.gene.cor.res.fit$TRM=as.character(tcga.TRM.gene.cor.res.fit$TRM)
tcga.TRM.gene.cor.res.fit$gene=as.character(tcga.TRM.gene.cor.res.fit$gene)
table(tcga.TRM.gene.cor.res.fit$TRM,tcga.TRM.gene.cor.res.fit$cor>0)
#         neg pos
# TRM_1     4  196
# TRM_3   167  500
write.table(tcga.TRM.gene.cor.res.fit,'results/tcga.TRM.gene.cor.res.fit.txt',quote = F,row.names = F,sep='\t')
#富集分析
library("clusterProfiler")
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)

module_genes=as.character(unique(tcga.TRM.gene.cor.res.fit$gene))
gene = bitr(module_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(gene,2)

#GO富集分析
ego <- enrichGO(gene = gene$ENTREZID,
                OrgDb = 'org.Hs.eg.db', 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
write.table(ego@result,file="results/GO.txt",sep="\t",quote=F,row.names = F)
#KEGG富集
R.utils::setOption("clusterProfiler.download.method",'auto') 
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05)
write.table(kk@result,file="results/kegg.txt",sep="\t",quote=F,row.names = F)

library(dplyr)
library("ggpubr")
keggdata <- data.frame(ONTOLOGY = 'KEGG', kk@result[ , -c(1, 2)])

enrich<-rbind.data.frame(ego@result,
                         keggdata)
table(enrich$ONTOLOGY)
write.table(enrich,'results/enrich.txt',quote = F,row.names = T,sep='\t')
data=enrich %>% group_by(ONTOLOGY) %>% slice_head(n=10)
data$ONTOLOGY=factor(data$ONTOLOGY,levels = c('BP','CC','MF','KEGG'))
pdf('results/enrich.pdf',height = 15,width = 12)
ggbarplot(data, x="Description", y="Count", fill = "ONTOLOGY", color = "white",
          xlab="",
          orientation = "horiz",   #图形横向显示
          palette = c("#DE582B", "#1868B2", "#018A67", "#F3A332"),    #设置颜色方案
          legend = "right",    #图例位置
          sort.val = "asc",    #GO排序方式
          sort.by.groups=TRUE)+    #根据GO分组进行排序
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))+
  xlab('')
dev.off()