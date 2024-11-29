setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results/TIME/')

###box plot
library(tidyverse)
library(pheatmap)

b <- read.csv('results/cluster_clinic_data.csv', row.names = 1)
a <- read.csv("geo_cibersort_result.csv", row.names = 1)

# 假设 a 和 b 是你的矩阵
# 将 b 的 group 列提取为向量
group_vector <- b[match(rownames(a), rownames(b)), "Cluster"]
# 将 group 列添加到 a
result <- cbind(a, group = group_vector)

write.csv(result, 'results/TIME/cibersort_group_result.csv')

result <- read.csv('results/TIME/cibersort_group_result.csv', row.names = 1)
result <- result %>% rownames_to_column("sample") #把行名变成一列，并加上“sample”
library(ggsci)
library(tidyr)
library(ggpubr)
result <- result[, -seq(ncol(result) - 3, ncol(result) - 1)]
result$group <- as.factor(result$group)
b <- gather(result,key=CIBERSORT,value = Proportion,-c(group,sample))

# 绘制箱线图
p <- ggboxplot(b, x = "CIBERSORT", y = "Proportion",
               fill = NULL,  # 填充颜色为空
               color = "group",  # 边框颜色根据 group
               palette = c("2" = "#34499d", "1" = "#e92428"),  # 定义边框颜色
               size = 0.6,  # 设置边框的粗细
               width = 0.5) +  # 调整箱线图的宽度
  stat_boxplot(geom = "errorbar",
               aes(ymin = after_stat(ymax), color = group), 
               position = position_dodge(0.8),  # 确保误差线对齐
               width = 0.5, size = 0.6) +  # 调整 width 与箱线图一致
  stat_boxplot(geom = "errorbar",
               aes(ymax = after_stat(ymin), color = group), 
               position = position_dodge(0.8),  # 确保误差线对齐
               width = 0.5, size = 0.6) +  # 调整 width 与箱线图一致
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 8,  # 调整显著性标注的大小
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", ""))) +
  geom_boxplot(aes(color = group), 
               position = position_dodge(0.8), 
               size = 0.6,
               width = 0.5, 
               outlier.colour = NA) +  # 设置箱线图的间隔
  # geom_jitter(aes(color = group), size = 1.5, width = 0.2) +  # 添加点，颜色根据 group
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),  # 调整 x 轴标注的字体大小
        axis.text.y = element_text(size = 15),  # 调整 y 轴标注的字体大小
        axis.title.y = element_text(size = 20))  # 调整 y 轴标题的字体大小

# 显示图形
print(p)

# 使用 ggsave() 函数保存图形
ggsave(filename = "results/TIME/cibersort_boxplot.pdf",  plot = p, width = 12, height = 6, units = "in")



cluster <- b
immune <- a
sameSample=intersect(row.names(cluster), row.names(immune))
cluster=cluster[sameSample, "Cluster", drop=F]
immune=immune[sameSample, , drop=F]
data=cbind(cluster, immune)
data <- data[,-c(24:26)]

outTab=data.frame()
sigCell=c("Cluster")
for(i in colnames(data)[2:ncol(data)]){
  if(sd(data[,i])<0.001){next}
  if(length(levels(factor(data[,"Cluster"])))>2){
    test=kruskal.test(data[,i] ~ data[,"Cluster"])
  }else{
    test=wilcox.test(data[,i] ~ data[,"Cluster"])
  }
  pvalue=test$p.value
  if(pvalue<1){
    outTab=rbind(outTab,cbind(immune=i, pvalue))
    sigCell=c(sigCell, i)
  }
}
write.table(file="results/TIME/immuneCor.txt", outTab, sep="\t", quote=F, row.names=F)

#??ͼ????
data=data[,sigCell]
#???͵?ע??
data$Cluster <- factor(data$Cluster,levels = c('1','2'))
data=data[order(data[,"Cluster"]),]

annCol=data[,1,drop=F]
annCol[,"Cluster"]=factor(annCol[,"Cluster"], unique(annCol[,"Cluster"]))
data=t(data[,(2:ncol(data))])
#???????͵?ע??
annRow=sapply(strsplit(rownames(data),"_"), '[', 2)
annRow=as.data.frame(annRow)
row.names(annRow)=row.names(data)
colnames(annRow)=c("Methods")
annRow[,"Methods"]=factor(annRow[,"Methods"], unique(annRow[,"Methods"]))
gapCol=as.vector(cumsum(table(annCol[,"Cluster"])))
gapRow=as.vector(cumsum(table(annRow[,"Methods"])))

#??????ͼע?͵???ɫ
bioCol=c("#0066FF","#FF9900","#00DB00","#FF0000","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(annCol[,"Cluster"]))]
bioCol <- c("#D20A13", '#11AA4D')
Cluster=bioCol
names(Cluster)=levels(factor(annCol[,"Cluster"]))
ann_colors=list(Cluster=Cluster)

#??ͼ???ӻ?
pdf("results/TIME/immHeatmap.pdf", width=7, height=5)
pheatmap(data,
         annotation=annCol,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#4DBBD5FF",5), "white", rep("#E64B35FF",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         gaps_col=gapCol,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=8,
         fontsize_row=8,
         fontsize_col=8)
dev.off()