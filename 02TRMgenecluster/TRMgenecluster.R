setwd('D:/ZJU-FISH/doctor/TRM/02TRMgenecluster/')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
dir.create('results')

library(GSVA)
library(ggrain)

dir.create('results')
bioSurvival=function(OS,OS.time,riskscore,labs,palette,leg){
  dat=data.frame(OS=OS,OS.time=OS.time,riskscore=riskscore)
  library(survival)
  library(survminer)
  res.cut <- surv_cutpoint(dat,
                           time = "OS.time", 
                           event = "OS", 
                           variables = c("riskscore"))
  cut_va=as.numeric(res.cut$cutpoint[1])
  print(cut_va)
  dat$risk=ifelse(dat$riskscore>=cut_va,'High','Low')
  ggplotKM<-function(time,status,group,labs,palette,leg){
    library(ggplot2)
    library(survival)
    dat1=data.frame(time=time,status=status,group=group)
    colnames(dat1)=c('time','status','groups')
    sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
    surv=survminer::ggsurvplot(sf, data = dat1, 
                               palette = palette, 
                               pval = TRUE,
                               surv.median.line='hv'
                               #,conf.int = T
                               ,conf.int.style ='step'
                               , pval.coord=c(0, 0.2), #Add p-value 
                               risk.table = TRUE, 
                               legend.title = leg
                               ,legend.labs =labs
                               ,conf.int=T
    )
    p1=surv$plot+theme_pubr()+
      theme(axis.text.y=element_text(family="Times",face="plain"),
            #axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
            legend.position=c(1,1),
            legend.justification=c(1,1),
            legend.background = element_rect(fill = NA, colour = NA),
            legend.title = element_text(family="Times",face="plain"),
            legend.text = element_text(family="Times",face="plain"))
    p2=surv$table+theme_pubr()+
      theme(axis.text.y=element_text(family="Times",face="plain"),
            plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
            plot.title=element_blank(),
            legend.position=c(1,1), 
            legend.justification=c(1,1),
            legend.title = element_text(family="Times",face="plain"),
            legend.text = element_text(family="Times",face="plain"))+
      xlab('Year')
    
    g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
    return(p1)
  }
  p=ggplotKM(time = dat$OS.time,status=dat$OS,group=dat$risk,labs=labs,palette=palette,leg=leg)
  return(p)
}
#
tcga.dat<-read.delim(paste0(routine_data, 'GEO/GSE53625/GSE53625_count_loged.txt'),sep='\t',header = T,row.names = 1,check.names = F)
tcga.group<-read.delim(paste0(routine_data, 'GEO/GSE53625/GSE53625_count_group.txt'),sep='\t',header = T)
#TRM的marker基因
cell.module.gene=read.delim('D:/ZJU-FISH/doctor/TRM/01scRNA/results/scRNA_marker_gene.txt',sep='\t',header = T)
cell.module.gene=cell.module.gene[,c("gene","cluster")]
head(cell.module.gene)
cell.module.gene$cluster=paste0('TRM_',cell.module.gene$cluster)
#ssGSEA计算TCGA样本的TRM评分
pathway.score<-function(exp,gene){
  ssGSEAScore_by_genes<-function(gene.exp,genes){
    #library('GSVA')
    #library(GSEABase)
    #all.list=list()
    gs=GSEABase::GeneSet(setName='GeneSet', 
                         setIdentifier=paste0("101"),
                         geneIds=unique(genes),
                         GSEABase::SymbolIdentifier()) 
    gsc <- GSEABase::GeneSetCollection(list(gs))
    fl <- tempfile()
    GSEABase::toGmt(gsc, fl)
    cgeneset=GSEABase::getGmt(fl)
    gsvaP <- ssgseaParam(
      exprData = as.matrix(gene.exp),
      geneSets = cgeneset,
      # assay = NA_character_,
      # annotation = NA_character_,
      minSize = 1,
      maxSize = 5000,
      # alpha = 0.25,
      # normalize = TRUE
    )
    ssGSEA.geneset <- gsva(gsvaP, verbose=TRUE)
    # ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp),
    #                              cgeneset,method='ssgsea',
    #                              min.sz=1, max.sz=5000, 
    #                              verbose=TRUE)
    #detach('package:GSVA')
    #detach('package:GSEABase')
    #row.names(ssGSEA.geneset)
    return(ssGSEA.geneset)
  }
  pathway_score<-data.frame()
  for (i in unique(gene[,2])){
    gene_set=gene[gene[,2]==i,1]
    score=ssGSEAScore_by_genes(exp,gene_set)
    rownames(score)=i
    pathway_score=rbind.data.frame(pathway_score,score)
  }
  return(t(pathway_score))
}
tcga.TRM.score<-pathway.score(exp = tcga.dat,gene = cell.module.gene)
head(tcga.TRM.score)
library(ggplot2)
library(ggsignif)
sig_violin<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  p<-ggplot(dat, aes(group, gene, fill = group, color = group)) +
    geom_rain(alpha = .5, rain.side = 'l',
              boxplot.args = list(color = "black", outlier.shape = NA),
              boxplot.args.pos = list(
                position = ggpp::position_dodgenudge(x = .1), width = 0.1
              )) +
    theme_classic() +
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2') +
    guides(fill = 'none', color = 'none') +
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(color=leg,fill=leg)+scale_color_manual(values = palette)+theme_classic()+
    scale_fill_manual(values = palette)
  return(p)
}
tcga.TRM.score.group<-merge(tcga.group,
                            data.frame(Samples=rownames(tcga.TRM.score),tcga.TRM.score),
                            by='Samples')
write.table(tcga.TRM.score.group,'results/tcga.TRM.score.txt',quote = F,row.names = T,sep='\t')

fig1a<-sig_violin(dat = tcga.TRM.score.group[,c("group","TRM_0")],
                  leg = 'Groups',ylab = 'TRM 0 Score',palette = c('#6888F5','#D77071'))
fig1a
fig1b<-sig_violin(dat = tcga.TRM.score.group[,c("group","TRM_1")],
                  leg = 'Groups',ylab = 'TRM 1 Score',palette = c('#6888F5','#D77071'))
fig1b
fig1c<-sig_violin(dat = tcga.TRM.score.group[,c("group","TRM_2")],
                  leg = 'Groups',ylab = 'TRM 2 Score',palette = c('#6888F5','#D77071'))
fig1c
fig1d<-sig_violin(dat = tcga.TRM.score.group[,c("group","TRM_3")],
                  leg = 'Groups',ylab = 'TRM 3 Score',palette = c('#6888F5','#D77071'))
fig1d
fig1e<-sig_violin(dat = tcga.TRM.score.group[,c("group","TRM_4")],
                  leg = 'Groups',ylab = 'TRM 4 Score',palette = c('#6888F5','#D77071'))
fig1e
fig1f<-sig_violin(dat = tcga.TRM.score.group[,c("group","TRM_5")],
                  leg = 'Groups',ylab = 'TRM 5 Score',palette = c('#6888F5','#D77071'))
fig1f
# fig1g<-sig_violin(dat = tcga.TRM.score.group[,c("Type","TRM_6")],
#                   leg = 'Groups',ylab = 'TRM_6 Score',palette = c('#66b131','#fd2615'))
# fig1g
# fig1h<-sig_violin(dat = tcga.TRM.score.group[,c("Type","TRM_7")],
#                   leg = 'Groups',ylab = 'TRM_7 Score',palette = c('#66b131','#fd2615'))
# fig1h
# fig1i<-sig_violin(dat = tcga.TRM.score.group[,c("Type","TRM_8")],
#                   leg = 'Groups',ylab = 'TRM_8 Score',palette = c('#66b131','#fd2615'))
# fig1i
fig1<-ggarrange(fig1a,fig1b,fig1c,fig1d,fig1e,fig1f,nrow = 2,ncol = 3,common.legend = T)
fig1
ggsave('results/Fig1.pdf',fig1,height = 6,width = 10)
library(ggrain)
fig1a+fig1b+fig1c+fig1d+fig1e+fig1f+plot_layout(guides = 'collect') &
  theme_bw()

#KM曲线
# 筛选出 group 列属于 "tumor" 的行
tcga.TRM.score.tumor <- tcga.TRM.score.group[tcga.TRM.score.group$group == "tumor", ]
# 修改 Sample 列，保留 "_" 之前的内容
tcga.TRM.score.tumor$Samples <- sub("_.*", "", tcga.TRM.score.tumor$Samples)

tcga.cli<-read.csv(paste0(routine_data, 'clinic/GEO_clinic_data.csv'))
tcga.TRM.cli<-merge(tcga.cli,
                    data.frame(Samples=tcga.TRM.score.tumor$Samples,tcga.TRM.score.tumor),
                    by='Samples')
write.table(tcga.TRM.cli,'results/tcga.TRM.cli.txt',quote = F,row.names = F,sep='\t')
fig2a<-bioSurvival(OS = tcga.TRM.cli$status,
                   OS.time = tcga.TRM.cli$time,
                   riskscore = tcga.TRM.cli$TRM_0,labs = c('High','Low'),
                   palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_0')
fig2a
fig2b<-bioSurvival(OS = tcga.TRM.cli$status,
                   OS.time = tcga.TRM.cli$time,
                   riskscore = tcga.TRM.cli$TRM_1,labs = c('High','Low'),
                   palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_1')
fig2b
fig2c<-bioSurvival(OS = tcga.TRM.cli$status,
                   OS.time = tcga.TRM.cli$time,
                   riskscore = tcga.TRM.cli$TRM_2,labs = c('High','Low'),
                   palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_2')
fig2c
fig2d<-bioSurvival(OS = tcga.TRM.cli$status,
                   OS.time = tcga.TRM.cli$time,
                   riskscore = tcga.TRM.cli$TRM_3,labs = c('High','Low'),
                   palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_3')
fig2d
fig2e<-bioSurvival(OS = tcga.TRM.cli$status,
                   OS.time = tcga.TRM.cli$time,
                   riskscore = tcga.TRM.cli$TRM_4,labs = c('High','Low'),
                   palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_4')
fig2e
fig2f<-bioSurvival(OS = tcga.TRM.cli$status,
                   OS.time = tcga.TRM.cli$time,
                   riskscore = tcga.TRM.cli$TRM_5,labs = c('High','Low'),
                   palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_5')
fig2f
# fig2g<-bioSurvival(OS = tcga.TRM.cli$status,
#                    OS.time = tcga.TRM.cli$time,
#                    riskscore = tcga.TRM.cli$TRM_6,labs = c('High','Low'),
#                    palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_6')
# fig2g
# fig2h<-bioSurvival(OS = tcga.TRM.cli$status,
#                    OS.time = tcga.TRM.cli$time,
#                    riskscore = tcga.TRM.cli$TRM_7,labs = c('High','Low'),
#                    palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_7')
# fig2h
# fig2i<-bioSurvival(OS = tcga.TRM.cli$status,
#                    OS.time = tcga.TRM.cli$time,
#                    riskscore = tcga.TRM.cli$TRM_8,labs = c('High','Low'),
#                    palette = c("#ed0000", "#00468b"),leg = 'GSE53625 TRM_8')
# fig2i
fig2<-ggarrange(fig2a,fig2b,fig2c,fig2d,fig2e,fig2f, nrow = 2,ncol = 3)
fig2
ggsave('results/tcga.TRM.KM.pdf',fig2,height = 6,width = 10)
#临床特征表达的比较
fig3a=list()
fig3a[[1]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicT","TRM_0")],
                       leg = 'T.Stage',ylab = 'TRM 0 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3a[[1]]
fig3a[[2]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicN","TRM_0")],
                       leg = 'N.Stage',ylab = 'TRM 0 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3a[[2]]

fig3a[[3]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicStage","TRM_0")],
                       leg = 'Stage',ylab = 'TRM 0 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3a[[3]]
fig3a=ggarrange(plotlist = fig3a,ncol = 3,nrow = 1,widths = 1)
fig3a
#
fig3b=list()
fig3b[[1]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicT","TRM_1")],
                       leg = 'T.Stage',ylab = 'TRM 1 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3b[[1]]
fig3b[[2]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicN","TRM_1")],
                       leg = 'N.Stage',ylab = 'TRM 1 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3b[[2]]

fig3b[[3]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicStage","TRM_1")],
                       leg = 'Stage',ylab = 'TRM 1 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3b[[3]]
fig3b=ggarrange(plotlist = fig3b,ncol = 3,nrow = 1,widths = c(1,1,1))
fig3b
#
fig3c=list()
fig3c[[1]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicT","TRM_2")],
                       leg = 'T.Stage',ylab = 'TRM 2 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3c[[1]]
fig3c[[2]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicN","TRM_2")],
                       leg = 'N.Stage',ylab = 'TRM 2 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3c[[2]]

fig3c[[3]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicStage","TRM_2")],
                       leg = 'Stage',ylab = 'TRM 2 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3c[[3]]
fig3c=ggarrange(plotlist = fig3c,ncol = 3,nrow = 1,widths = 1)
fig3c
#
fig3d=list()
fig3d[[1]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicT","TRM_3")],
                       leg = 'T.Stage',ylab = 'TRM 3 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3d[[1]]
fig3d[[2]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicN","TRM_3")],
                       leg = 'N.Stage',ylab = 'TRM 3 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3d[[2]]

fig3d[[3]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicStage","TRM_3")],
                       leg = 'Stage',ylab = 'TRM 3 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3d[[3]]
fig3d=ggarrange(plotlist = fig3d,ncol = 3,nrow = 1,widths = 1)
fig3d
#
fig3e=list()
fig3e[[1]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicT","TRM_4")],
                       leg = 'T.Stage',ylab = 'TRM 4 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3e[[1]]
fig3e[[2]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicN","TRM_4")],
                       leg = 'N.Stage',ylab = 'TRM 4 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3e[[2]]

fig3e[[3]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicStage","TRM_4")],
                       leg = 'Stage',ylab = 'TRM 4 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3e[[3]]
fig3e=ggarrange(plotlist = fig3e,ncol = 3,nrow = 1,widths = 1)
fig3e
fig3f=list()
fig3f[[1]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicT","TRM_5")],
                       leg = 'T.Stage',ylab = 'TRM 5 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3f[[1]]
fig3f[[2]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicN","TRM_5")],
                       leg = 'N.Stage',ylab = 'TRM 5 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3f[[2]]

fig3f[[3]]<-sig_violin(dat = tcga.TRM.cli[,c("pathologicStage","TRM_5")],
                       leg = 'Stage',ylab = 'TRM 5 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3f[[3]]
fig3f=ggarrange(plotlist = fig3f,ncol = 3,nrow = 1,widths = 1)
fig3f
fig3<-ggarrange(fig3a,fig3b,fig3c,fig3d,fig3e,fig3f,nrow = 6,ncol = 1,labels = '')
ggsave('results/tcga.TRM.cli.pdf',fig3,height = 18,width = 10)
