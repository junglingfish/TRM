setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'

TPM <- read.table(paste0(routine_data, 'GEO/GSE53625/GSE53625_tumor_count_loged.txt'), check.names = F)
#计算行平均值，按降序排列
index=order(rowMeans(TPM[,-1]),decreasing = T)
#调整表达谱的基因顺序
TPM_ordered=TPM[index,]

#在线读取probesets、genes文件的运行方式
library(MCPcounter)
estimate<- MCPcounter.estimate( 
  TPM_ordered, 
  featuresType='HUGO_symbols',  
  probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep='\t',stringsAsFactors=FALSE,colClasses='character'),
  genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep='\t',stringsAsFactors=FALSE,header=TRUE,colClasses='character',check.names=FALSE))

a <- data.frame(t(estimate))
# 假设 scores 是你的数据矩阵
# 获取所有行名
rownames_a <- rownames(a)
# 创建一个逻辑向量，标记符合条件的行
keep_rows <- logical(length(rownames_a))
for (i in 1:length(rownames_a)) {
  rowname <- rownames_a[i]
  
  if (grepl("^TCGA", rowname) && substr(rowname, 15, 15) != "6") {
    # 保留前12位字符并将'.'替换为'-'
    rownames_a[i] <- gsub("\\.", "-", substr(rowname, 1, 12))
    keep_rows[i] <- TRUE
  } else if (grepl("^ec", rowname)) {
    # 保留'_'之前的内容
    rownames_a[i] <- sub("_.*", "", rowname)
    keep_rows[i] <- TRUE
  } else {
    keep_rows[i] <- FALSE
  }
}
# 更新行名
rownames(a) <- rownames_a
# 筛选符合条件的行
a <- a[keep_rows, ]
write.csv(a, paste0(routine_fig7, 'MPcounter_result.csv'))


###box plot

library(tidyverse)
b <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)
a <- read.csv("results/MPcounter_result.csv", row.names = 1)

# 假设 a 和 b 是你的矩阵
# 将 b 的 group 列提取为向量
group_vector <- b[match(rownames(a), b$Samples), "Risk"]
# 将 group 列添加到 a
result <- cbind(a, group = group_vector)

write.csv(result, 'results/MPcounter_group_result.csv')

result <- read.csv('results/MPcounter_group_result.csv', row.names = 1)
result$group <- ifelse(result$group == "High", "Low", "High")
result <- result %>% rownames_to_column("sample") #把行名变成一列，并加上“sample”
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(result,key=MPcounter,value = Abundance,-c(group,sample))

# 绘制箱线图
p <- ggboxplot(b, x = "MPcounter", y = "Abundance",
               fill = NULL,  # 填充颜色为空
               color = "group",  # 边框颜色根据 group
               palette = c("Low" = "#34499d", "High" = "#e92428"),  # 定义边框颜色
               size = 0.4,  # 设置边框的粗细
               width = 0.5) +  # 调整箱线图的宽度
  stat_boxplot(geom = "errorbar",
               aes(ymin = after_stat(ymax), color = group), 
               position = position_dodge(0.8),  # 确保误差线对齐
               width = 0.5, size = 0.4) +  # 调整 width 与箱线图一致
  stat_boxplot(geom = "errorbar",
               aes(ymax = after_stat(ymin), color = group), 
               position = position_dodge(0.8),  # 确保误差线对齐
               width = 0.5, size = 0.4) +  # 调整 width 与箱线图一致
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 8,  # 调整显著性标注的大小
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", ""))) +
  geom_boxplot(aes(color = group), 
               position = position_dodge(0.8), 
               size = 0.4,
               width = 0.5, 
               outlier.colour = NA) +  # 设置箱线图的间隔
  # geom_jitter(aes(color = group), size = 1.5, width = 0.2) +  # 添加点，颜色根据 group
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),  # 调整 x 轴标注的字体大小
        axis.text.y = element_text(size = 15),  # 调整 y 轴标注的字体大小
        axis.title.y = element_text(size = 20))  # 调整 y 轴标题的字体大小

# 将图形保存为 PDF 文件
pdf("results/MCPcounter.pdf", width = 8, height = 8)  # 调整 width 和 height 根据需要更改图形尺寸

print(p)  # 显示图形
dev.off()  # 关闭 PDF 设备