setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

library(estimate)
library(readxl)


ESCCExpr <- paste0(routine_data, 'GEO/GSE53625/GSE53625_tumor_count_loged.txt')
# read.table(ESCCExpr)
filterCommonGenes(input.f = ESCCExpr, 
                  output.f = 'ESCC_exp.gct', 
                  id = "GeneSymbol")
estimateScore(input.ds = 'ESCC_exp.gct',
              output.ds = "ESCC_estimate_score.gct", 
              platform = "affymetrix")
scores <- read.table("ESCC_estimate_score.gct",skip = 2,header = T)
rownames(scores) <- scores[,1]
scores <- t(scores[,3:ncol(scores)])

# 假设 scores 是你的数据矩阵
# 获取所有行名
rownames_scores <- rownames(scores)
# 创建一个逻辑向量，标记符合条件的行
keep_rows <- logical(length(rownames_scores))
for (i in 1:length(rownames_scores)) {
  rowname <- rownames_scores[i]
  
  if (grepl("^TCGA", rowname) && substr(rowname, 15, 15) != "6") {
    # 保留前12位字符并将'.'替换为'-'
    rownames_scores[i] <- gsub("\\.", "-", substr(rowname, 1, 12))
    keep_rows[i] <- TRUE
  } else if (grepl("^ec", rowname)) {
    # 保留'_'之前的内容
    rownames_scores[i] <- sub("_.*", "", rowname)
    keep_rows[i] <- TRUE
  } else {
    keep_rows[i] <- FALSE
  }
}
# 更新行名
rownames(scores) <- rownames_scores
# 筛选符合条件的行
scores <- scores[keep_rows, ]
# write.csv(scores, paste0(routine_fig7, 'scores.csv'))
write.csv(scores, 'results/GEO_scores.csv')

# riskscore <- read.csv(paste0(routine_fig7, 'riskscore.csv'), check.names = F)
riskscore <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)

# 假设 scores 和 riskscore 是你的数据矩阵
# 提取 group 列，使用 match 函数确保顺序一致
group_column <- riskscore[match(rownames(scores), riskscore$Samples), "Risk"]
# 将 group 列添加到 scores 矩阵
scores <- cbind(scores, group = group_column)

write.csv(scores, 'results/immunescores+riskscore_TRM.csv')


##boxplot
scores <- read.csv('results/immunescores+riskscore_TRM.csv', row.names = 1)
# scores <- read_xlsx('results/immunescores+riskscore_TRM.xlsx', sheet = 1)
# scores <- read.csv(paste0(routine_fig7, 'GEO_immunescores+riskscore_TRM.csv'), row.names = 1)

# namescore <- scores[[1]]
# scores <- as.data.frame(scores[-1])
# rownames(scores) <- namescore

library(ggplot2)
library(dplyr)
library(ggpubr)

# 假设 scores 矩阵包含 ImmuneScore 和 group 列
# 将 scores 转换为数据框，以便于 ggplot 使用
scores_df <- as.data.frame(scores)

# 绘制箱线图
p <- ggplot(scores_df, aes(x = group, y = ImmuneScore)) +
  geom_jitter(aes(color = group), size = 2, alpha = 0.6) +
  geom_boxplot(aes(color = group), fill = NA, outlier.colour = NA, size = 1.5) +  # 边框颜色
  stat_boxplot(geom = "errorbar",
               aes(ymin = after_stat(ymax), color = group), 
               width = 0.3, size = 1.5) +
  stat_boxplot(geom = "errorbar",
               aes(ymax = after_stat(ymin), color = group), width = 0.3, size = 1.5) +
  scale_color_manual(values = c("Low" = "#34499d", "High" = "#e92428")) +  # 边框颜色
  labs(x = NULL, y = "Immune Score") +  # 取消 x 轴标题
  scale_x_discrete(labels = c("Low" = "Low risk group", "High" = "High risk group")) +  # 改变 x 轴标注
  theme_minimal() +
  theme(
    text = element_text(size = 14),  # 设置字体大小
    panel.grid.major = element_blank(),  # 取消网格线
    panel.grid.minor = element_blank(),  # 取消小网格线
    axis.line.x = element_line(color = "black", size = 1.5),  # 添加 x 轴黑色实线
    axis.line.y = element_line(color = "black", size = 1.5),  # 添加 y 轴黑色实线
    axis.ticks = element_line(color = "black", size = 1),  # 显示 x 和 y 轴的刻度线
    axis.title.y = element_text(color = "black", size = 20),  # 设置 y 轴标题颜色和大小
    axis.text.x = element_text(color = "black", size = 15, angle = 45, hjust = 1),  # 设置 x 轴标尺颜色和大小
    axis.text.y = element_text(color = "black", size = 15),  # 设置 y 轴标尺颜色和大小
    axis.ticks.length = unit(5, "mm"),  # 设置刻度线长度为 5 毫米
    plot.title = element_text(hjust = 0.5, face = "bold")  # 居中加粗标题
  ) +
  ggtitle("Immune Score by Risk Group")
  # stat_compare_means(aes(group = list(c("Low", "High"))),
  #                    method = "wilcox.test",
  #                    label = "p.signif",
  #                    size = 5,  # 调整显著性标注的大小
  #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.1, 1),
  #                                       symbols = c("***", "**", "*", "")))

# 添加显著性标注
p <- p + geom_signif(
  comparisons = list(c("Low", "High")),  # 根据实际分组的名称调整
  map_signif_level = function(p) {
    if (p <= 0.001) {
      return("***")
    } else if (p <= 0.01) {
      return("**")
    } else if (p <= 0.05) {
      return("*")
    } else {
      return("NS.")
    }
  },
  y_position = c(2600),  # 根据数据的实际范围调整
  textsize = 7,  # 设置显著性标注的字体大小
  size = 1.5  # 设置标注的线条宽度
)

# 将图形保存为 PDF 文件
pdf("results/immunescore_TRM.pdf", width = 6, height = 8)  # 调整 width 和 height 根据需要更改图形尺寸

print(p)  # 显示图形
dev.off()  # 关闭 PDF 设备


#####################################################################################################################
# 绘制箱线图
p <- ggplot(scores_df, aes(x = group, y = TumorPurity)) +
  geom_jitter(aes(color = group), size = 2, alpha = 0.6) +
  geom_boxplot(aes(color = group), fill = NA, outlier.colour = NA, size = 1.5) +  # 边框颜色
  stat_boxplot(geom = "errorbar",
               aes(ymin = after_stat(ymax), color = group), 
               width = 0.3, size = 1.5) +
  stat_boxplot(geom = "errorbar",
               aes(ymax = after_stat(ymin), color = group), width = 0.3, size = 1.5) +
  scale_color_manual(values = c("Low" = "#34499d", "High" = "#e92428")) +  # 边框颜色
  labs(x = NULL, y = "Immune Score") +  # 取消 x 轴标题
  scale_x_discrete(labels = c("Low" = "Low risk group", "High" = "High risk group")) +  # 改变 x 轴标注
  theme_minimal() +
  theme(
    text = element_text(size = 14),  # 设置字体大小
    panel.grid.major = element_blank(),  # 取消网格线
    panel.grid.minor = element_blank(),  # 取消小网格线
    axis.line.x = element_line(color = "black", size = 1.5),  # 添加 x 轴黑色实线
    axis.line.y = element_line(color = "black", size = 1.5),  # 添加 y 轴黑色实线
    axis.ticks = element_line(color = "black", size = 1),  # 显示 x 和 y 轴的刻度线
    axis.title.y = element_text(color = "black", size = 20),  # 设置 y 轴标题颜色和大小
    axis.text.x = element_text(color = "black", size = 15, angle = 45, hjust = 1),  # 设置 x 轴标尺颜色和大小
    axis.text.y = element_text(color = "black", size = 15),  # 设置 y 轴标尺颜色和大小
    axis.ticks.length = unit(5, "mm"),  # 设置刻度线长度为 5 毫米
    plot.title = element_text(hjust = 0.5, face = "bold")  # 居中加粗标题
  ) +
  ggtitle("TumorPurity by Risk Group")
# stat_compare_means(aes(group = list(c("Low", "High"))),
#                    method = "wilcox.test",
#                    label = "p.signif",
#                    size = 5,  # 调整显著性标注的大小
#                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.1, 1),
#                                       symbols = c("***", "**", "*", "")))

# 添加显著性标注
p <- p + geom_signif(
  comparisons = list(c("Low", "High")),  # 根据实际分组的名称调整
  map_signif_level = function(p) {
    if (p <= 0.001) {
      return("***")
    } else if (p <= 0.01) {
      return("**")
    } else if (p <= 0.05) {
      return("*")
    } else {
      return("NS.")
    }
  },
  y_position = c(1),  # 根据数据的实际范围调整
  textsize = 7,  # 设置显著性标注的字体大小
  size = 1.5  # 设置标注的线条宽度
)

# 将图形保存为 PDF 文件
pdf("results/TumorPurity_TRM.pdf", width = 6, height = 8)  # 调整 width 和 height 根据需要更改图形尺寸

print(p)  # 显示图形
dev.off()  # 关闭 PDF 设备