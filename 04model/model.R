setwd('D:/ZJU-FISH/doctor/TRM/04model/')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_2 = 'D:/ZJU-FISH/doctor/TRM/03keygene/'

library(ggpubr)
library(ggsci)

dir.create('results')
tcga.cli<-read.csv(paste0(routine_data, 'clinic/GEO_clinic_data.csv'), check.names = F)
tcga.dat<-read.delim(paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt'),sep='\t',header = T,row.names = 1,check.names = F)
tcga.caf.gene.cor.res.fit<-read.delim(paste0(routine_2, 'results/tcga.TRM.gene.cor.res.fit.txt'),sep='\t',header = T)
mycols <- pal_npg('nrc',alpha = 0.6)(9)
mycols <- rep(mycols,14)
cox_batch<-function(dat,time,event){
  coxRun<-function(dat){
    library(survival)
    colnames(dat)=c('time','status','AS')  
    dat=dat[which(!is.na(dat[,1])&!is.na(dat[,3])&!is.na(dat[,2])),]
    #print(nrow(dat))
    if(nrow(dat)<10){
      print(paste0('Sample Num is small:',nrow(dat)))
      return(c(NA,NA,NA,NA))
    }
    #if(quantile(dat[,3])['25%']==quantile(dat[,3])['50%']) return(c(NA,NA,NA,NA))
    fmla <- as.formula("Surv(time, status) ~AS")
    if(table(dat[,2])[1]>1&table(dat[,2])[2]>1){
      cox <- survival::coxph(fmla, data = dat)
      re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
      return(re)
    }else{
      return(c(NA,NA,NA,NA))
    }
  }
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('p.value','HR','Low 95%CI','High 95%CI')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}
tcga.cox<-cox_batch(dat = tcga.dat[as.character(unique(tcga.caf.gene.cor.res.fit$gene)),tcga.cli$Samples],
                    time = tcga.cli$time,event = tcga.cli$status)
cox.p<-0.05#可修改
table(tcga.cox$p.value<cox.p)

tcga.cox.fit<-tcga.cox[which(tcga.cox$p.value<cox.p),]
dim(tcga.cox.fit)#219
dim(tcga.cox)
tcga.cox$coef=log(tcga.cox$HR)
tcga.cox$Gene=rownames(tcga.cox)
tcga.cox$type=rep('None',nrow(tcga.cox))
tcga.cox$type[which(tcga.cox$p.value<cox.p & tcga.cox$coef>0)]='Risk'
tcga.cox$type[which(tcga.cox$p.value<cox.p & tcga.cox$coef<0)]='Protective'
table(tcga.cox$type)

write.table(tcga.cox,'results/tcga.cox.txt',quote = F,row.names = T,sep='\t')
tcga.cox$label <- c(names(tcga_risk$module.gene),rep(NA,nrow(tcga.cox)-8))
library(ggVolcano)
library(latex2exp)
library(ggrepel)
library(ggplot2)
fig1a <- ggplot(data = tcga.cox,
                aes(x = coef,
                    y = -log10(p.value)))+
  geom_point(alpha=0.4, size=3.5, aes(color=type))+
  scale_color_manual(values=c("#00468b",'grey',"#ed0000"),
                     limits = c("Protective",'None', "Risk"),name='State')+
  geom_hline(yintercept = -log10(cox.p),lty=4,col="black",lwd=0.8)+
  ylab('-log10(pvalue)')+xlab('Cox coefficient')+
  theme_bw()+
  theme(axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
        legend.text=element_text(face="plain", family="Times", colour="black"),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black" ),#设置图例的总标题的字体属性
        legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill = NA, colour = NA)
  )
data_selected <- tcga.cox[names(tcga_risk$module.gene),]
fig <- fig1a + geom_label_repel(data=data_selected,
                                aes(label=rownames(data_selected)))

fig
ggsave('results/tcga.sig.cox.pdf',fig,height = 4.5,width = 5)
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  crbind2DataFrame=function(dat){
    print(class(dat))
    if(class(dat)=='table'){
      if(!is.na(ncol(dat))){
        dat=apply(dat,2,function(x){
          return(x)
        })
      }
    }
    if(class(dat)!='data.frame'){
      dat1=as.data.frame(as.matrix(dat))
    }else{
      dat1=dat
    }
    #print(head(dat1))
    for(i in 1:ncol(dat1)){
      dat1[,i]=as.character(dat1[,i])
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      if(sum(is.na(as.numeric(dt)))==0){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
    return(dat1)  
  }
  library(glmnet)
  set.seed(3)
  options(ggrepel.max.overlaps=Inf)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=sig.coef,lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  crbind2DataFrame=function(dat){
    print(class(dat))
    if(class(dat)=='table'){
      if(!is.na(ncol(dat))){
        dat=apply(dat,2,function(x){
          return(x)
        })
      }
    }
    if(class(dat)!='data.frame'){
      dat1=as.data.frame(as.matrix(dat))
    }else{
      dat1=dat
    }
    #print(head(dat1))
    for(i in 1:ncol(dat1)){
      dat1[,i]=as.character(dat1[,i])
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      if(sum(is.na(as.numeric(dt)))==0){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
    return(dat1)  
  }
  
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=cox1$coefficients,model=mult_results))
}
tcga.module.dat<- tcga.dat[rownames(tcga.cox.fit),tcga.cli$Samples]
rownames(tcga.module.dat)=gsub('-','__',rownames(tcga.module.dat))
tcga.lasso<-get_riskscore.lasso(dat = t(tcga.module.dat[,tcga.cli$Samples]),
                                os=tcga.cli$status,
                                os.time = tcga.cli$time,
                                labels='')

length(tcga.lasso$lasso.gene)#9
tcga.lasso$plot
ggsave('results/lasso.pdf',tcga.lasso$plot,height = 5,width = 10)
tcga.lasso$lambda.min
#0.06533129
tcga_risk<-get_riskscore(dat=as.data.frame(t(tcga.module.dat[names(tcga.lasso$lasso.gene),tcga.cli$Samples])),
                         os=tcga.cli$status,
                         os.time=tcga.cli$time,
                         step=F,#不进行逐步回归
                         direction=c("both", "backward", "forward")[1])

length(tcga_risk$module.gene)#9
names(tcga_risk$module.gene)
#"ANGPTL7" "C6"      "CSRP1"   "EXPH5"   "F2RL2"   "KCNMA1" 
# "MAGEC3"  "MAMDC2"  "SLC4A9" 
tcga_risk$model
#"0.093*ANGPTL7+0.15*C6+0.121*CSRP1+-0.08*EXPH5+0.12*F2RL2+0.014*KCNMA1+-0.373*MAGEC3+0.143*MAMDC2+-0.188*SLC4A9"

#柱状图
gene.coef=as.data.frame(tcga_risk$module.gene)
gene.coef=data.frame(Gene=rownames(gene.coef),Coef=gene.coef[,1])
gene.coef$Type=ifelse(gene.coef$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
library(dplyr)
fig1d=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  coord_flip() +
  scale_fill_manual(values=ggsci::pal_lancet('lanonc',alpha =0.6)(9)[c(7,1)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Cox coefficient") +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")

fig1d
ggsave('results/fig1d.pdf',fig1d,height = 5,width = 5)
#AUC和ROC
risk_plot<-function(time,event,riskscore,
                    group,mk,labs=c('High','Low'),
                    palette=c("#00468b", "#ed0000", "#42b540")){
  dat=data.frame(time=time,status=event,riskscore=riskscore,group=group)
  #ROC曲线
  ROC_rt=timeROC::timeROC(T=dat$time, 
                          delta=dat$status,
                          marker=dat$riskscore, cause=1,
                          weighting='marginal',
                          times=mk, 
                          ROC=TRUE,iid = T)
  p.dat=data.frame()
  for(i in which(ROC_rt$times>0)){
    lbs=paste0('AUC at ',mk[i],' years: ',round(ROC_rt$AUC[i],2))
    p.dat=rbind.data.frame(p.dat,data.frame(V1=ROC_rt$FP[,i],V2=ROC_rt$TP[,i],Type=lbs))
  }
  #colnames(p.dat)=c('V1','V2','Type')
  p.dat=as.data.frame(p.dat)
  p.dat$Type=as.factor(p.dat$Type)
  roc_plot=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))+
    stat_smooth(aes(colour=Type),se = FALSE, size = 1)+
    theme_bw()+xlab('False positive fraction')+
    ylab('True positive fraction') +
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_text(family="Times",face="plain"),
          axis.title.x=element_text(family="Times",face="plain"),
          axis.title.y=element_text(family="Times",face="plain"),
          plot.title=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
          legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  #KM曲线
  library(ggplot2)
  library(survival)
  dat1=dat[,c("time","status","group")]
  colnames(dat1)=c('time','status','groups')
  sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
  surv=survminer::ggsurvplot(sf, data = dat1, 
                             palette = c("#ed0000", "#00468b"), 
                             pval = TRUE,
                             surv.median.line='hv'
                             #,conf.int = T
                             ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             linetype = 'solid',
                             legend.title = 'Group'
                             ,legend.labs =labs
                             ,conf.int=T
  )
  p1=surv$plot+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
          legend.position=c(1,1),
          legend.justification=c(1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  p2=surv$table+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
          plot.title=element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  gg<-ggpubr::ggarrange(g2,roc_plot,ncol = 2,nrow = 1,labels = '')
  return(gg)
}
#开始分析
tcga_risk$result$Risk=ifelse(tcga_risk$result$riskscorez>0,'High','Low')
tcga_roc_km<-risk_plot(time=tcga_risk$result$time,
                       event=tcga_risk$result$status,
                       riskscore=tcga_risk$result$riskscorez,
                       group=tcga_risk$result$Risk,
                       mk=c(1,3,5),labs=c('High','Low'),
                       palette=c("#00468b", "#ed0000", "#42b540"))
tcga_roc_km
ggsave('results/gse53625.auc.km.pdf',tcga_roc_km,height = 5,width = 10)
#TCGA验证
gse31210.dat<-read.delim(paste0(routine_data, 'results/TCGA_tumor_cleaned.txt'),sep='\t',header = T,row.names = 1, check.names = F)
gse31210.cli<-read.csv(paste0(routine_data, 'clinic/TCGA_clinic_data.csv'), row.names = 1)
gse31210.dat <- as.data.frame(gse31210.dat)
colnames(gse31210.dat) <- rownames(gse31210.cli)
gse31210_dat_m<-cbind.data.frame(OS.time=gse31210.cli$time,
                                 OS=gse31210.cli$status,
                                 t(gse31210.dat[intersect(names(tcga_risk$module.gene),rownames(gse31210.dat)),]))
gse31210_dat_m[1:4,1:4]
gse31210_risk<-get_riskscore(dat=gse31210_dat_m[,-c(1,2)],
                             os=gse31210_dat_m$OS,
                             os.time=gse31210_dat_m$OS.time,
                             step=F,#不进行逐步回归
                             direction=c("both", "backward", "forward")[1])
gse31210_risk$result$Risk=ifelse(gse31210_risk$result$riskscorez>0,'High','Low')
gse31210_roc_km<-risk_plot(time=gse31210_risk$result$time,
                           event=gse31210_risk$result$status,
                           riskscore=gse31210_risk$result$riskscorez,
                           group=gse31210_risk$result$Risk,
                           mk=c(1, 2, 3),labs=c('High','Low'),
                           palette=ggsci::pal_lancet()(10)[c(2,1)])
gse31210_roc_km
ggsave('results/tcga.auc.km.pdf',gse31210_roc_km,height = 5,width = 10)
# #
# gse3141.dat<-read.delim('../03.data.pre/GSE3141/results/gse3141.dat.txt',sep='\t',header = T,row.names = 1)
# gse3141.cli<-read.delim('../03.data.pre/GSE3141/results/gse3141.cli.txt',sep='\t',header = T)
# gse3141_dat_m<-cbind.data.frame(OS.time=gse3141.cli$OS.time,
#                                 OS=gse3141.cli$OS,
#                                 t(gse3141.dat[intersect(names(tcga_risk$module.gene),rownames(gse3141.dat)),gse3141.cli$Accession]))
# gse3141_dat_m[1:4,1:4]
# 
# gse3141_risk<-get_riskscore(dat=gse3141_dat_m[,-c(1,2)],
#                             os=gse3141_dat_m$OS,
#                             os.time=gse3141_dat_m$OS.time,
#                             step=F,#不进行逐步回归
#                             direction=c("both", "backward", "forward")[1])
# gse3141_risk$result$Risk=ifelse(gse3141_risk$result$riskscorez>0,'High','Low')
# gse3141_roc_km<-risk_plot(time=gse3141_risk$result$time/365,
#                           event=gse3141_risk$result$status,
#                           riskscore=gse3141_risk$result$riskscorez,
#                           group=gse3141_risk$result$Risk,
#                           mk=c(1,3,5),labs=c('High','Low'),
#                           palette=ggsci::pal_lancet()(10)[c(2,1)])
# gse3141_roc_km
# ggsave('results/gse3141.auc.km.pdf',gse3141_roc_km,height = 5,width = 10)
# #GSE37745
# gse37745.dat<-read.delim('../03.data.pre/GSE37745/results/GSE37745.exp.txt',sep='\t',header = T,row.names = 1)
# gse37745.cli<-read.delim('../03.data.pre/GSE37745/results/GSE37745.cli.txt',sep='\t',header = T)
# gse37745_dat_m<-cbind.data.frame(OS.time=gse37745.cli$OS.time,
#                                  OS=gse37745.cli$OS,
#                                  t(gse37745.dat[intersect(names(tcga_risk$module.gene),rownames(gse37745.dat)),gse37745.cli$Samples]))
# gse37745_dat_m[1:4,1:4]
# 
# gse37745_risk<-get_riskscore(dat=gse37745_dat_m[,-c(1,2)],
#                              os=gse37745_dat_m$OS,
#                              os.time=gse37745_dat_m$OS.time,
#                              step=F,#不进行逐步回归
#                              direction=c("both", "backward", "forward")[1])
# gse37745_risk$result$Risk=ifelse(gse37745_risk$result$riskscorez>0,'High','Low')
# gse37745_roc_km<-risk_plot(time=gse37745_risk$result$time/365,
#                            event=gse37745_risk$result$status,
#                            riskscore=gse37745_risk$result$riskscorez,
#                            group=gse37745_risk$result$Risk,
#                            mk=c(1,3,5),labs=c('High','Low'),
#                            palette=ggsci::pal_lancet()(10)[c(2,1)])
# gse37745_roc_km
# ggsave('results/gse37745.auc.km.pdf',gse37745_roc_km,height = 5,width = 10)
# #GSE50081
# gse50081.dat<-read.delim('../03.data.pre/GSE50081//results/GSE50081.exp.txt',sep='\t',header = T,row.names = 1)
# gse50081.cli<-read.delim('../03.data.pre/GSE50081/results/GSE50081.cli.txt',sep='\t',header = T)
# gse50081_dat_m<-cbind.data.frame(OS.time=gse50081.cli$OS.time,
#                                  OS=gse50081.cli$OS,
#                                  t(gse50081.dat[intersect(names(tcga_risk$module.gene),rownames(gse50081.dat)),gse50081.cli$Samples]))
# gse50081_dat_m[1:4,1:4]
# 
# gse50081_risk<-get_riskscore(dat=gse50081_dat_m[,-c(1,2)],
#                              os=gse50081_dat_m$OS,
#                              os.time=gse50081_dat_m$OS.time,
#                              step=F,#不进行逐步回归
#                              direction=c("both", "backward", "forward")[1])
# gse50081_risk$result$Risk=ifelse(gse50081_risk$result$riskscorez>0,'High','Low')
# gse50081_roc_km<-risk_plot(time=gse50081_risk$result$time,
#                            event=gse50081_risk$result$status,
#                            riskscore=gse50081_risk$result$riskscorez,
#                            group=gse50081_risk$result$Risk,
#                            mk=c(1,3,5),labs=c('High','Low'),
#                            palette=ggsci::pal_lancet()(10)[c(2,1)])
# gse50081_roc_km
# ggsave('results/gse50081.auc.km.pdf',gse50081_roc_km,height = 5,width = 10)
# #GSE68465
# gse68465.dat<-read.delim('../03.data.pre/GSE68465/results/GSE68465.exp.txt',sep='\t',header = T,row.names = 1)
# gse68465.cli<-read.delim('../03.data.pre/GSE68465/results/GSE68465.cli.txt',sep='\t',header = T)
# gse68465_dat_m<-cbind.data.frame(OS.time=gse68465.cli$OS.time,
#                                  OS=gse68465.cli$OS,
#                                  t(gse68465.dat[intersect(names(tcga_risk$module.gene),rownames(gse68465.dat)),gse68465.cli$Samples]))
# gse68465_dat_m[1:4,1:4]
# 
# gse68465_risk<-get_riskscore(dat=gse68465_dat_m[,-c(1,2)],
#                              os=gse68465_dat_m$OS,
#                              os.time=gse68465_dat_m$OS.time,
#                              step=F,#不进行逐步回归
#                              direction=c("both", "backward", "forward")[1])
# gse68465_risk$result$Risk=ifelse(gse68465_risk$result$riskscorez>0,'High','Low')
# gse68465_roc_km<-risk_plot(time=gse68465_risk$result$time/365,
#                            event=gse68465_risk$result$status,
#                            riskscore=gse68465_risk$result$riskscorez,
#                            group=gse68465_risk$result$Risk,
#                            mk=c(1,3,5),labs=c('High','Low'),
#                            palette=ggsci::pal_lancet()(10)[c(2,1)])
# gse68465_roc_km
# ggsave('results/gse68465.auc.km.pdf',gse68465_roc_km,height = 5,width = 10)
tcga_risk1=tcga_risk$result
gse31210_risk1=gse31210_risk$result
gse3141_risk1=gse3141_risk$result
gse37745_risk1=gse37745_risk$result
gse50081_risk1=gse50081_risk$result
gse68465_risk1=gse68465_risk$result

write.table(tcga_risk1,'results/gse53625.risk.txt',quote = F,row.names = F,sep='\t')
write.table(gse31210_risk1,'results/tcga_risk.txt',quote = F,row.names = F,sep='\t')
# write.table(gse3141_risk1,'results/gse3141_risk.txt',quote = F,row.names = F,sep='\t')
# write.table(gse37745_risk1,'results/gse37745_risk.txt',quote = F,row.names = F,sep='\t')
# write.table(gse50081_risk1,'results/gse50081_risk.txt',quote = F,row.names = F,sep='\t')
# write.table(gse68465_risk1,'results/gse68465_risk.txt',quote = F,row.names = F,sep='\t')