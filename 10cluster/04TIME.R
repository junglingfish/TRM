#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results/TIME/')

#???ð?
library(limma)
library(ggpubr)
library(ggplot2)

cluFile="results/cluster_risk.csv"          #???͵Ľ????ļ?
scoreFile="immunescores+riskscore_TRM.csv"      #΢?????????ļ?

score_df=read.csv(scoreFile, row.names=1)
score_df=as.matrix(score_df)
score <- matrix(as.numeric(score_df), nrow = nrow(score_df), ncol = ncol(score_df))  # 转换后重建矩阵
rownames(score) <- rownames(score_df)
colnames(score) <- colnames(score_df)
# score=avereps(score)

#??ȡ?????ļ????????????ļ?????
cluster=read.csv(cluFile)
rownames(cluster) = cluster$Samples

#??Ʒȡ????
sameSample=intersect(row.names(score), cluster$Samples)
score1=score[sameSample,,drop=F]
cluster1=cluster[sameSample,"Cluster",drop=F]
data=cbind(score1, cluster1)
# data$Cluster <- as.factor(data$Cluster)
library(ggrain)
#???ñȽ???
type=levels(factor(data[,"Cluster"]))
data$Cluster=factor(data$Cluster, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#??????ɫ
bioCol=c("#ff5104","#0395cf")
bioCol=bioCol[1:length(unique(data$Cluster))]

# boxplot=ggplot(data,aes(x = Cluster, 
#                         y = ImmuneScore, 
#                         fill = Cluster))


for(i in colnames(data)[1:(ncol(data)-1)]){
 # i = "StromalScore"
  boxplot=ggplot(data,aes(x = Cluster, 
                          y = !!sym(i), 
                          fill = Cluster)) +
  scale_fill_manual(values = c("#ff5104","#0395cf")) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color="black") +
  geom_boxplot(notch = FALSE, outlier.size = -1,
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() +
  ylab(")") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 1, size = 12,face = "bold.italic"),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15,face = "bold.italic"),
        axis.text = element_text(size = 10)) +
  stat_compare_means(comparisons=my_comparisons)
  
  #????ͼ??
  pdf(file=paste0(i, ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}

# library(ggplot2)
# 
# # 假设 data 是你的数据框，并包含 Cluster 和 ImmuneScore 列
# data$Cluster <- as.factor(data$Cluster)  # 转换 Cluster 为因子
# 
# # 创建箱线图
# boxplot <- ggplot(data, aes(x = Cluster, y = ImmuneScore, fill = Cluster)) +
#   geom_boxplot() +  # 添加箱线图
#   theme_bw() +  # 使用白色背景主题
#   labs(title = "Boxplot of Immune Score by Cluster", 
#        x = "Cluster", 
#        y = "Immune Score")
# 
# # 打印箱线图
# print(boxplot)
