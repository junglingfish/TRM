setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

#TCGA表达谱
tcga.dat<-read.delim(paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt'),sep='\t',header = T,row.names = 1,check.names = F)

tcga.dat[1:4,1:4]
tcga.module<-read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)
colnames(tcga.module)
tcga.module.gene.exp=tcga.module[,3:10]
rownames(tcga.module.gene.exp)=tcga.module$Samples
#Estimate
immu_estimate=function(exp,platform='illumina',isTCGA=T){#affymetrix
  library(estimate)
  inf=tempfile()
  ouf=tempfile()
  ouf2=tempfile()
  writeMatrix=function(dat,outpath,row=T,header=T){
    if(row){
      write.table(cbind(Tag=row.names(dat),dat),file=outpath,sep="\t",quote = F,row.names=F,col.names = header)
    }else{
      write.table(dat,file=outpath,sep="\t",quote = F,row.names=F,col.names = header)
    }
  }
  writeMatrix(exp,outpath = inf)
  outputGCT(inf, ouf)
  estimateScore(ouf,ouf2,platform=platform)
  est.score=t(read.csv(ouf2,sep = '\t',row.names = 1,check.names = F,skip = 2))
  est.score=est.score[-1,]
  rnames=row.names(est.score)
  if(isTCGA){
    rnames=gsub('\\.','-',row.names(est.score))
  }
  est.score=apply(est.score, 2, as.numeric)
  row.names(est.score)=rnames
  
  return(est.score)
}
tcga.esti=immu_estimate(exp = tcga.dat,platform = 'illumina',isTCGA = F)
save(tcga.esti,file = 'results/tcga.esti.RData')
cor_point <- function(dat, xlab, ylab, method, leg, cols, title_size = 20, axis_text_size = 18, legend_text_size = 18) {
  library(ggpubr)
  colnames(dat) <- c('name1', 'name2', 'name3')
  
  p <- ggplot(dat, aes(x = name1, y = name2, color = name3)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor(method = method, aes(x = name1, y = name2, color = name3)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size = title_size),
          axis.text = element_text(size = axis_text_size),
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_text_size)) +
    xlab(xlab) + ylab(ylab) +
    scale_color_manual(values = cols) +
    labs(colour = leg)
  
  return(p)
}

load('results/tcga.esti.RData')
head(tcga.esti)

# tcga.esti <- read.csv('immunescores.csv', row.names = 1)
tcga.esti1=cbind.data.frame(tcga.esti,tcga.module.gene.exp[rownames(tcga.esti),])
tcga.esti1=reshape2::melt(tcga.esti1,id=colnames(tcga.esti))
head(tcga.esti1)
fig1a<-cor_point(dat = tcga.esti1[,c("ImmuneScore","value","variable")],
                 xlab = 'ImmuneScore',ylab = 'Gene Expression',
                 method = 'spearman',leg = 'Gene',
                 cols = ggsci::pal_jco(alpha = 0.8)(9))
ggsave('results/Fig1a.pdf',height = 9,width = 9)
#
esti.gene=cbind.data.frame(tcga.esti,tcga.module.gene.exp[rownames(tcga.esti),])
esti.gene.cor <- Hmisc::rcorr(as.matrix(esti.gene))
pdf('results/Fig1a1.pdf', height = 5, width = 5)
corrplot::corrplot(t(as.matrix(esti.gene.cor$r[colnames(tcga.esti), colnames(tcga.module.gene.exp)])), 
                   p.mat = t(as.matrix(esti.gene.cor$P[colnames(tcga.esti), colnames(tcga.module.gene.exp)])),
                   mar = c(0,0,1,1),
                   col = colorRampPalette(c('#34499d', 'white', '#e92428'))(100),
                   tl.srt = 60,
                   tl.cex = 1,
                   tl.col = 'black',
                   tl.offset = 0.5,
                   cl.pos = c("b", "r", "n")[1], 
                   cl.align.text = 'l',
                   cl.length = 5,
                   cl.ratio = 0.1,
                   cl.cex = 0.8,
                   addgrid.col = 'white',
                   method = 'color',
                   insig = 'label_sig',
                   sig.level = c(0.001, 0.01, 0.05),
                   pch.cex = 1,
                   is.corr = TRUE,
                   xpd = TRUE)

dev.off()

#高低表达的比较
head(esti.gene)
esti.gene1=esti.gene
colnames(esti.gene1)
esti.gene1$ABLIM1=ifelse(esti.gene1$ABLIM1>median(esti.gene1$ABLIM1),'High','Low')
esti.gene1$CSTB=ifelse(esti.gene1$CSTB>median(esti.gene1$CSTB),'High','Low')
esti.gene1$MAL2=ifelse(esti.gene1$MAL2>median(esti.gene1$MAL2),'High','Low')
esti.gene1$SLPI=ifelse(esti.gene1$SLPI>median(esti.gene1$SLPI),'High','Low')
esti.gene1$FANCB=ifelse(esti.gene1$FANCB>median(esti.gene1$FANCB),'High','Low')
esti.gene1$IFT57=ifelse(esti.gene1$IFT57>median(esti.gene1$IFT57),'High','Low')
esti.gene1$KIF11=ifelse(esti.gene1$KIF11>median(esti.gene1$KIF11),'High','Low')
esti.gene1$PPEF1=ifelse(esti.gene1$PPEF1>median(esti.gene1$PPEF1),'High','Low')

head(esti.gene1)
sig_point<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggplot(data = dat,aes(x = group, 
                           y = gene, 
                           fill = group))+ 
    scale_fill_manual(values = c("#e92428", "#34499d")) + #用自定义颜色填充
    geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                size = 0.8, color="black") +
    geom_boxplot(notch = FALSE, outlier.size = -1, 
                 color="black", lwd=0.8, alpha = 0.7) +
    geom_point(shape = 21, size=2, # 点的性状和大小
               position = position_jitterdodge(), # 让点散开
               color="black", alpha = 1) +
    theme_classic() + 
    ylab("ImmuneScore") +
    xlab("") +
    theme(axis.text.x = element_text(hjust = 1, size = 16,face = "bold.italic"),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none",
          axis.title = element_text(size = 18,face = "bold.italic"),
          axis.text = element_text(size = 16)) + 
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
  return(pp)
}
fig1b<-list()
fig1b[[1]]<-sig_point(dat = esti.gene1[,c("ABLIM1","ImmuneScore")],
                      leg = 'ABLIM1',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[1]]
fig1b[[2]]<-sig_point(dat = esti.gene1[,c("CSTB","ImmuneScore")],
                      leg = 'CSTB',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[2]]
fig1b[[3]]<-sig_point(dat = esti.gene1[,c("MAL2","ImmuneScore")],
                      leg = 'MAL2',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[3]]
fig1b[[4]]<-sig_point(dat = esti.gene1[,c("SLPI","ImmuneScore")],
                      leg = 'SLPI',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[4]]
fig1b[[5]]<-sig_point(dat = esti.gene1[,c("FANCB","ImmuneScore")],
                      leg = 'FANCB',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[5]]
fig1b[[6]]<-sig_point(dat = esti.gene1[,c("IFT57","ImmuneScore")],
                      leg = 'IFT57',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[6]]
fig1b[[7]]<-sig_point(dat = esti.gene1[,c("KIF11","ImmuneScore")],
                      leg = 'KIF11',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[7]]
fig1b[[8]]<-sig_point(dat = esti.gene1[,c("PPEF1","ImmuneScore")],
                      leg = 'PPEF1',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[8]]
# fig1b[[9]]<-sig_point(dat = esti.gene1[,c("SLC4A9","ImmuneScore")],
#                       leg = 'SLC4A9',ylab = 'ImmuneScore',
#                       palette = ggsci::pal_lancet()(9)[3:4])
# fig1b[[9]]
fig1b<-ggpubr::ggarrange(plotlist = fig1b,nrow = 2,ncol = 4)
ggsave('results/Fig1b.pdf',fig1b,height = 9,width = 18)
#CIBERSORT
source('CIBERSORT.R')
# tcga.ciber=CIBERSORT("ref.txt", paste0(routine_result, 'tumor_exp_cleaned.txt'), perm=100, QN=TRUE)
# write.table(tcga.ciber,'results/tcga.ciber.txt',quote = F,sep = '\t',row.names = T)
# save(tcga.ciber,file = 'results/tcga.ciber.RData')

load('results/tcga.ciber.RData')
tcga.ciber[1:4,1:4]
#
ciber.gene=cbind.data.frame(tcga.ciber[,1:22],tcga.module.gene.exp[rownames(tcga.ciber),])
ciber.gene.cor <- Hmisc::rcorr(as.matrix(ciber.gene))
pdf('results/Fig1c.pdf',height = 6,width = 12)
corrplot::corrplot(t(as.matrix(ciber.gene.cor$r[colnames(tcga.ciber)[1:22],colnames(tcga.module.gene.exp)])), 
                   p.mat = t(as.matrix(ciber.gene.cor$P[colnames(tcga.ciber)[1:22],colnames(tcga.module.gene.exp)])),
                   mar = c(0,0,1,1),
                   col=colorRampPalette(c('#34499d', 'white','#e92428'))(100),
                   tl.srt = 60,
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
#mcpcount
immu_MCPcounter=function(exp,isTCGA=T){
  mcpEstimates=MCPcounter::MCPcounter.estimate(exp,featuresType="HUGO_symbols",
                                               probesets=read.table('immu_mcp_probes.txt',sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                                               genes=read.table('immu_mcp_genes.txt',sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
  )
  mcpEstimates=t(mcpEstimates)
  if(isTCGA){
    rnames=gsub('\\.','-',row.names(mcpEstimates)) 
    row.names(mcpEstimates)=rnames
  }
  return(mcpEstimates)
}
tcga.mcp=immu_MCPcounter(exp = tcga.dat,isTCGA = T)
head(tcga.mcp)
library(ggcorrplot) 
mcp.gene=cbind.data.frame(tcga.mcp,tcga.module.gene.exp[rownames(tcga.mcp),])
write.table(mcp.gene,'results/mcp.gene.txt',quote = F,row.names = T,sep='\t')
pmtcars <- cor_pmat(mcp.gene)
cormtcars <- round(cor(mcp.gene), 3)
pdf('results/fig1d.pdf',height = 7,width = 7)
ggcorrplot(cormtcars[colnames(tcga.mcp),colnames(tcga.module.gene.exp)],
           hc.order = F,  #分等级聚类重排矩阵
           ggtheme = ggplot2::theme_void(base_size = 15), #主题修改
           colors = c('#34499d','white','#e92428'), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
           lab = T,lab_size = 3,    #相关系数文本字体大小
           tl.cex = 12,             #坐标轴字体大小
           p.mat = pmtcars[colnames(tcga.mcp),colnames(tcga.module.gene.exp)],         #添加显著性信息
           sig.level = 0.01,        #显著性水平
           pch = 4,                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
           pch.cex = 10)  
dev.off()