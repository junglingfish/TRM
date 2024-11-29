setwd('D:/ZJU-FISH/doctor/TRM/08PDL')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')


tcga.module<-read.delim(paste0(routine_4, 'results/gse53625.risk.txt'),sep='\t',header = T)
colnames(tcga.module)
module.gene=colnames(tcga.module)[3:10]
module.gene
library(ggpubr)
library(colorspace)
sig_boxplot<-function(dat,leg,ylab,palette=c("#80b8df", "#EF767A")){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggplot(dat,aes(x=group,y=gene,fill=group))+
    geom_violin(width =0.8,color='black',size=1)+
    theme_classic() + 
    theme(text = element_text(size=10, colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 10),
          axis.title.y = element_text(color = 'black', size = 12),
          axis.line = element_line(size = 1))+ 
    theme(legend.position="none") +  
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
                 geom = "pointrange", color = "black", size=1)+
    scale_fill_manual(values = palette) +  # 使用传入的调色板
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(color=leg)
  return(pp)
}
#计算风险得分
library(survival)
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
ggplotKM<-function(time,status,group,labs,palette){
  library(ggplot2)
  library(survival)
  dat1=data.frame(time=time,status=status,group=group)
  colnames(dat1)=c('time','status','groups')
  sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
  surv=survminer::ggsurvplot(sf, data = dat1, 
                             palette = c("#ed0000", "#00468b"), 
                             pval = TRUE,
                             surv.median.line='hv'
                             # ,conf.int = T
                             ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
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
  return(g2)
}

#1、IMV210数据集####
imv210.exp<-read.delim('IMvigor210/IMvigor210_Counts2TPM.txt',sep='\t',header = T)
IMvigor210_anno <- read.delim('IMvigor210/IMvigor210_entrez2gene.txt', header = T)
imv210.exp=merge(IMvigor210_anno[,c("entrez_id","symbol")],
                 imv210.exp,by.x='entrez_id',by.y='genes')
imv210.exp[1:4,1:4]
imv210.exp=imv210.exp[,-1]
imv210.exp=aggregate(.~symbol,imv210.exp,mean)
imv210.exp[1:4,1:4]
rownames(imv210.exp)=imv210.exp$symbol
imv210.exp=imv210.exp[-1,-1]
range(imv210.exp)
imv210.exp <- log2(imv210.exp + 1)
IMvigor210_cli <- read.delim('IMvigor210/IMvigor210_cli.txt', header = T)
IMvigor210_cli <- IMvigor210_cli[, c("os", "censOS",
                                     "Best.Confirmed.Overall.Response",
                                     "IC.Level", "TC.Level", "Immune.phenotype",
                                     "FMOne.mutation.burden.per.MB",
                                     "Neoantigen.burden.per.MB","TCGA.Subtype")]
colnames(IMvigor210_cli) <- c('OS.time', 'OS', 'Response', 
                              "IC.Level", "TC.Level", "Immune.phenotype",
                              'TMB', 'NEO','Stage')
table(IMvigor210_cli$Response)
IMvigor210_cli=IMvigor210_cli[which(IMvigor210_cli$Response != 'NE'),]
IMvigor210_cli$Response[IMvigor210_cli$Response=='CR'|IMvigor210_cli$Response=='PR']<-'CR/PR'
IMvigor210_cli$Response[IMvigor210_cli$Response=='PD'|IMvigor210_cli$Response=='SD']<-'PD/SD'
imv210.exp=imv210.exp[,rownames(IMvigor210_cli)]
imv.risk<-get_riskscore(dat = as.data.frame(t(imv210.exp[intersect(rownames(imv210.exp),module.gene),])),
                        os = IMvigor210_cli$OS,
                        os.time = IMvigor210_cli$OS.time,
                        step = F,direction = 'both')
imv.risk$result$Risk=ifelse(imv.risk$result$riskscorez>0,'High','Low')
fig1a<-ggplotKM(time = imv.risk$result$time/30,
                status =imv.risk$result$status,
                group = imv.risk$result$Risk,
                labs = c('High','Low'),
                palette = c("#00468b", "#ed0000") )
fig1a
imv.risk.cli<-merge(imv.risk$result[,c("Samples","riskscorez","Risk")],
                    data.frame(Samples=rownames(IMvigor210_cli),IMvigor210_cli),
                    by='Samples')
write.table(x = imv.risk.cli,'results/imv.risk.cli.txt',quote = F,row.names = F,sep='\t')
fig1b=sig_boxplot(dat = imv.risk.cli[,c("Response","riskscorez")],
                  leg = 'Response',
                  ylab = 'RiskScore',
                  palette = c("#80b8df", "#EF767A"))
fig1b

response.risk=prop.table(table(imv.risk.cli$Response,imv.risk.cli$Risk),margin=2)

response.risk=reshape2::melt(response.risk)
colnames(response.risk)<-c("type","Risk","Percentage")

response.risk$Percentage<-round(response.risk$Percentage,digits=2)
write.table(response.risk,'results/response.risk.txt',quote = F,row.names = F,sep='\t')

fig1c=ggplot(response.risk,aes(x=Risk,y=Percentage,fill=type))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  scale_fill_manual(values = c("CR/PR" = "#80b8df", "PD/SD" = "#EF767A"))  # 替换为你的类型和颜色
fig1c
#I+II
imv.risk.cli1=imv.risk.cli[which(imv.risk.cli$Stage=='I'|imv.risk.cli$Stage =='II'),]
fig1d<-ggplotKM(time = imv.risk.cli1$OS.time/30,
                status =imv.risk.cli1$OS,
                group = imv.risk.cli1$Risk,
                labs = c('High','Low'),
                palette = c("#00468b", "#ed0000") )
fig1d
imv.risk.cli2=imv.risk.cli[which(imv.risk.cli$Stage=='III'|imv.risk.cli$Stage=='IV'),]

fig1e<-ggplotKM(time = imv.risk.cli2$OS.time/30,
                status =imv.risk.cli2$OS,
                group = imv.risk.cli2$Risk,
                labs = c('High','Low'),
                palette = c("#00468b", "#ed0000") )
fig1e
#GSE78220
GSE78220_cli <- read.delim('GSE78220/table_2.xls',sep=',',header = T)

GSE78220_cli1=data.frame(Samples=GSE78220_cli$Accession,
                         Title=GSE78220_cli$Title,
                         Rresponse=GSE78220_cli$anti.pd.1.response,
                         OS.time=GSE78220_cli$overall.survival..days.,
                         OS=GSE78220_cli$vital.status)
GSE78220_cli1 <- na.omit(GSE78220_cli1)
table(GSE78220_cli1$OS)
GSE78220_cli1$OS <- ifelse(GSE78220_cli1$OS == 'Alive', 0, 1)
rownames(GSE78220_cli1) <- GSE78220_cli1$Title

GSE78220_exp <- openxlsx::read.xlsx('GSE78220/GSE78220_PatientFPKM.xlsx',
                                    sheet = 1)
rownames(GSE78220_exp) <- GSE78220_exp$Gene
GSE78220_exp <- GSE78220_exp[, -1]
colnames(GSE78220_exp) <- stringr::str_split_fixed(colnames(GSE78220_exp), '\\.', 3)[, 1]
boxplot(GSE78220_exp[, 1:5])
GSE78220_exp=log2(GSE78220_exp+1)

gse.risk<-get_riskscore(dat = as.data.frame(t(GSE78220_exp[intersect(rownames(GSE78220_exp),module.gene),GSE78220_cli1$Title])),
                        os = GSE78220_cli1$OS,
                        os.time = GSE78220_cli1$OS.time,
                        step = F,direction = 'both')
gse.risk$result$Risk=ifelse(gse.risk$result$riskscorez>0,'High','Low')
fig1f<-ggplotKM(time = gse.risk$result$time/30,
                status =gse.risk$result$status,
                group = gse.risk$result$Risk,
                labs = c('High','Low'),
                palette = c("#00468b", "#ed0000") )
fig1f
gse.risk.cli<-merge(gse.risk$result[,c("Samples","riskscorez","Risk")],
                    data.frame(Samples=rownames(GSE78220_cli1),Response=GSE78220_cli1$Rresponse),
                    by='Samples')
table(gse.risk.cli$Response)
gse.risk.cli$Response[gse.risk.cli$Response=='Complete Response'|gse.risk.cli$Response=='Partial Response']<-'CR/PR'
gse.risk.cli$Response[gse.risk.cli$Response=='Progressive Disease']<-'PD/SD'

write.table(x = gse.risk.cli,'results/gse.risk.cli.txt',quote = F,row.names = F,sep='\t')
gse.risk.cli$Response <- factor(gse.risk.cli$Response,levels = c('CR/PR','PD/SD'))
fig1g=sig_boxplot(dat = gse.risk.cli[,c("Response","riskscorez")],
                  leg = 'Response',
                  ylab = 'RiskScore',
                  palette = c("#80b8df", "#EF767A"))
fig1g
response.risk=prop.table(table(gse.risk.cli$Response,gse.risk.cli$Risk),margin=2)

response.risk=reshape2::melt(response.risk)
colnames(response.risk)<-c("type","Risk","Percentage")

response.risk$Percentage<-round(response.risk$Percentage,digits=2)
write.table(response.risk,'results/gse.response.risk.txt',quote = F,row.names = F,sep='\t')

fig1h=ggplot(response.risk,aes(x=Risk,y=Percentage,fill=type))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  scale_fill_manual(values = c("CR/PR" = "#80b8df", "PD/SD" = "#EF767A"))   # 替换为你的类型和颜色
fig1h
fig1<-ggarrange(fig1a,fig1d,fig1e,fig1b,fig1c,fig1f,fig1g,fig1h,
                nrow = 2,ncol = 4,labels = '')
fig1
ggsave('results/Fig1.pdf',fig1,height = 12,width = 20)