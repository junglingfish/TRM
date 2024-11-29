#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("GSEABase")

# install.packages("ggpubr")
# install.packages("reshape2")

setwd('D:/ZJU-FISH/doctor/TRM/09immunofunction')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

#引用包
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)

expFile=paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt')           #表达数据文件
gmtFile="immune.gmt"           #免疫数据集文件
riskFile = paste0(routine_4, 'results/gse53625.risk.txt')       #风险文件

#读取表达输入文件，并对输入文件处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
# rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

#读取免疫基因集文件
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssgsea分析
# ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
gsvapar <- gsvaParam(mat, geneSet, maxDiff=TRUE)
ssgseaScore <- gsva(gsvapar)


#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)

# #去除正常样品
# group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# data=t(data[,group==0])
# rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
# data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F)
# risk=read.csv(riskFile, row.names = 1, check.names = F)
rownames(risk) <- risk$Samples
#合并数据
data <- t(data)
sameSample=intersect(row.names(data), risk$Samples)
data=data[sameSample,,drop=F]
risk=risk[sameSample,"Risk",drop=F]
rt1=cbind(data, risk)

#对免疫相关功能绘制箱线图
data=melt(rt1, id.vars=c("Risk"))
colnames(data)=c("Risk","Type","Score")
data$risk=factor(data$Risk, levels=c("Low","High"))
p=ggboxplot(data, x="Type", y="Score", color = "risk",
            xlab="",ylab="Score",add = "none",palette = c("#34499d","#e92428") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

#输出图片文件
pdf(file="immFunction.pdf", width=11, height=6)
print(p)
dev.off()