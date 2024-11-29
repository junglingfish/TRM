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

setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results/Consensus/')

###基因表达热图
library(pheatmap)
exp <- read.table(paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt'), header = TRUE, row.names = 1, check.names=FALSE)

group <- read.csv('results/cluster_clinic_data.csv',row.names = 1)
risk <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'),sep='\t',header = T)
gene <- colnames(risk)[3:10]

##数据清洗
# 对 pathologicT 列中等于 '--' 的值替换为 'T0'
group$pathologicT[group$pathologicT == "'--"] <- "T1"

# 对 pathologicN 列中等于 '--' 的值替换为 'N0'
group$pathologicN[group$pathologicN == "'--"] <- "N0"
group$pathologicN[group$pathologicN == "NX"] <- "N0"

# 替换 pathologicStage 中的值
group$pathologicStage <- gsub("Stage IA|Stage IB", "I", group$pathologicStage)
group$pathologicStage <- gsub("Stage II|Stage IIA|Stage IIB", "II", group$pathologicStage)
group$pathologicStage <- gsub("Stage III|Stage IIIA|Stage IIIB|IIIA|IIIB|IIIC", "III", group$pathologicStage)
group$pathologicStage <- gsub("Stage IV|Stage IVA", "IV", group$pathologicStage)

# 对等于 '--' 的值替换为 'I'
group$pathologicStage[group$pathologicStage == "'--"] <- "I"

# 根据 group 中的 Samples 列来重新排列 risk 的行顺序
risk_ordered <- risk[match(rownames(group), risk$Samples), ]
# 确保两者的 Samples 列顺序一致后，将 Risk 列添加到 group 中
group$Risk <- risk_ordered$Risk

write.csv(group, 'groupdata_for_heatmap.csv')


group <- read.csv('groupdata_for_heatmap.csv', row.names = 1)
group$Cluster <- as.factor(group$Cluster)
# 1. 筛选关键基因
# 只保留基因表达矩阵 exp 中的关键基因（即 gene 的行名）
exp_filtered <- exp[rownames(exp) %in% gene, ]

# 2. 对表达矩阵进行标准化（z-score）
# 行标准化，确保每个基因的表达值具有可比性
exp_scaled <- t(scale(t(exp_filtered)))

# 3. 裁剪数据，限制在 -2 到 2 之间
exp_clipped <- pmin(pmax(exp_scaled, -2), 2)

# 确保样本按照分组排序
group_order <- order(group$Cluster)
exp_clipped <- exp_clipped[, group_order]

# 创建注释数据框，顺序为 group, TStage, NStage, TNM, gender
annotation_col <- data.frame(
  gender = group$gender[group_order],
  grade = group$Tumor.grade[group_order],
  NStage = group$pathologicN[group_order],
  TStage = group$pathologicT[group_order],
  TNM = group$pathologicStage[group_order],
  Risk = group$Risk[group_order],
  group = group$Cluster[group_order]
)
rownames(annotation_col) <- colnames(exp_clipped)

# 定义注释颜色
annotation_colors <- list(
  group = c("1" = "#bf3a2a", "2" = "#1777bc"),
  Risk = c('High' = '#e46870', 'Low' = '#76bd77'),
  TStage = c("T1" = "blue", "T2" = "cyan", "T3" = "#7bcc00", "T4" = "#f4aa02"),
  NStage = c("N0" = "#fd9e4c", "N1" = "#0bd492", 'N2' = '#fc85b7', 'N3' = '#abc601'),
  TNM = c("I" = "lightgreen", "II" = "yellow", "III" = "orange", "IV" = "red"),
  grade = c('well' = '#fa81f8', 'moderately' = '#09d6c4', 'poorly' = '#bdbe0c'),
  gender = c('male' = '#00dbbb', 'female' = '#f3948e')
)

# 绘制热图
# pdf("results/Consensus/cluster_heatmap.pdf", width = 8, height = 6.5) # 可根据需要调整宽度和高度
pheatmap(exp_clipped,
         annotation_col = annotation_col,   # 样本分组注释
         annotation_colors = annotation_colors, # 注释颜色
         cluster_rows = TRUE,              # 行（基因）聚类
         cluster_cols = FALSE,             # 列（样本）不聚类
         color = colorRampPalette(c("#1777bc", "white", "#bf3a2a"))(100), # 蓝-白-红渐变
         show_colnames = FALSE,            # 不显示列名
         show_rownames = TRUE,             # 显示行名
         annotation_legend = TRUE,         # 显示注释图例
         legend = TRUE,                    # 显示颜色标尺
         fontsize_row = 12,                # 调整行名字体大小
         fontsize_col = 12,                # 调整列名字体大小
         main = "") # 主标题
# dev.off()