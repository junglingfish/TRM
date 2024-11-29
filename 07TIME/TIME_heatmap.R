library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
library(TCGAbiolinks)
library(homologene)
library(readxl)
library(org.Hs.eg.db)
library(LEA)
library('DESeq2')#加载包；

setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

###基因表达热图
library(pheatmap)
exp <- read.csv('results/ssgsea_result.csv', row.names = 1, check.names=FALSE)

immune <- read.csv('results/immunescores+riskscore_TRM.csv',row.names = 1)
# gene <- read.csv(paste0(routine_sfig1, '1_cox_6_genes.csv'),row.names = 1)

# 首先确保 `exp` 和 `group` 矩阵的行名内容一致
exp <- exp[rownames(immune), ]

# 将 `exp_selected` 合并到 `group` 中
group <- cbind(exp, immune)

exp_scaled <- t(scale(exp))

# 创建注释数据框，确保 Risk 列在第一列
annotation_col <- data.frame(
  Risk = factor(group$group, levels = c("Low", "High")),  # Risk 列，确保是因子类型
  TumorPurity = group$TumorPurity,
  ESTIMATEScore = group$ESTIMATEScore,
  ImmuneScore = group$ImmuneScore,
  StromalScore = group$StromalScore
)

# 设置行名
rownames(annotation_col) <- rownames(group)

# 定义注释颜色
annotation_colors <- list(
  Risk = c("Low" = "#1777bc", "High" = "#bf3a2a"),  # 为 Risk 列设置颜色
  TumorPurity = colorRampPalette(c("blue", "yellow"))(100),  # 数值型变量渐变色
  ESTIMATEScore = colorRampPalette(c("#6ac5ec", "#070f12", "#d5d100"))(100),  # 数值型变量渐变色
  ImmuneScore = colorRampPalette(c("darkblue", "lightblue"))(100),  # 数值型变量渐变色
  StromalScore = colorRampPalette(c("#54feff", "#040000", "#f20202"))(100)  # 数值型变量渐变色
)

# 调整注释列顺序，将 Risk 列移到最前面
annotation_col <- annotation_col[, c("ESTIMATEScore", "StromalScore", "TumorPurity", "ImmuneScore", "Risk")]


# 绘制热图
pheatmap(exp_scaled,  # 转置数据矩阵
         annotation_col = annotation_col,   # 样本分组注释
         annotation_colors = annotation_colors, # 注释颜色
         cluster_rows = F,              # 行（基因）聚类
         cluster_cols = FALSE,             # 列（样本）不聚类
         color = colorRampPalette(c("#1777bc", "white", "#bf3a2a"))(100), # 蓝-白-红渐变
         show_colnames = FALSE,            # 不显示列名
         show_rownames = TRUE,             # 显示行名
         annotation_legend = TRUE,         # 显示注释图例
         legend = TRUE,                    # 显示颜色标尺
         fontsize_row = 10,                # 调整行名字体大小
         fontsize_col = 10,                # 调整列名字体大小
         main = "") # 主标题
