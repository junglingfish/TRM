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
library(ggrepel)

setwd('D:/ZJU-FISH/doctor/TRM/06GSVA')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
routine_3 = 'D:/ZJU-FISH/doctor/TRM/03keygene/'
dir.create('results/difgene')

###########################################################################################
##GSE53625
exp <- read.table(paste0(routine_data, 'all_expr_batch.txt'), header = TRUE, row.names = 1)
group <- read.table(paste0(routine_data, 'all_group.txt'), header = TRUE, row.names = 1)
exp <- round(2^exp)

# 确保表达矩阵的列名与分组数据的行名一致
all(rownames(group) == colnames(exp)) # 检查是否匹配

# 确保group的分组信息为因子类型
group$group <- as.factor(group$group)

# 构建DESeq2的数据集
dds <- DESeqDataSetFromMatrix(countData = exp, colData = group, design = ~group)

# 过滤一些低表达量的基因
dds <- dds[rowSums(counts(dds)) > 1, ]

# 运行DESeq2的标准化和差异表达分析
dds1 <- DESeq(dds)

# 差异分析，注意：contrast中的group应根据你的实际分组信息进行修改，如 "Tumor" 和 "Normal"
res <- results(dds1, contrast = c('group', 'tumor', 'normal'))

# 查看标准化后的基本信息
summary(res)

# 输出结果为数据框
res <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

# 将结果保存到本地
write.table(res, file = 'results/difgene/deseq_diff.csv', col.names = NA, sep = ',', quote = FALSE)

# 筛选差异表达基因
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05), 'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05), 'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) < 1 | res1$padj >= 0.05), 'sig'] <- 'none'

write.table(res1, file = 'results/difgene/deseq_diff_upanddown.csv', col.names = NA, sep = ',', quote = FALSE)

# 输出筛选的差异基因总表
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'results/difgene/normal_tumor.DESeq2.select.csv', sep = ',', col.names = NA, quote = FALSE)

# 根据上调和下调分别输出
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'results/difgene/normal_tumor.DESeq2.up.csv', sep = ',', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'results/difgene/normal_tumor.DESeq2.down.csv', sep = ',', col.names = NA, quote = FALSE)

tcga.TRM.gene.cor.res.fit <- read.delim(paste0(routine_3, 'results/tcga.TRM.gene.cor.res.fit.txt'),sep='\t',header = T)
keygene <- tcga.TRM.gene.cor.res.fit$gene
res1_filtered <- res1_select[rownames(res1_select) %in% keygene, ]
res1_filtered$log2FoldChange = abs(res1_filtered$log2FoldChange)
write.table(res1_filtered, file = 'results/DESeq2.gene.for.circle.csv', sep = ',', col.names = NA, quote = FALSE)

# 使用ggplot2绘制火山图
library(ggplot2)

p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 2, alpha = 0.5) +  # 绘制散点图
  scale_color_manual(values = c('#b93a29', 'gray', '#0172b6'), limits = c('up', 'none', 'down')) +  # 自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'Normal vs Tumor', color = '') +  # 坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), # 背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  # 添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-6, 6) + ylim(0, 300)  # 定义刻度边界

# 保存为 PDF 文件
ggsave(filename = "results/difgene/volcano_plot.pdf", plot = p, device = "pdf", width = 6, height = 5)

###################################################################################################################################
# ##缺氧通路差异基因火山图绘制
# #genesets <- msigdbr(species = "Homo sapiens")
# genesets <- msigdbr(species = "Homo sapiens", category = "H")
# genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
# genesets <- split(genesets$gene_symbol, genesets$gs_name)

gene <- read.delim(paste0(routine_3, 'results/tcga.TRM.gene.cor.res.fit.txt'),sep='\t',header = T)
genemarker <- gene$gene

res <- read.csv('results/difgene/deseq_diff_upanddown.csv', row.names = 1)

res$log2FoldChange <- abs(res$log2FoldChange)

# 筛选res矩阵的行名，仅保留在genesets$HALLMARK_HYPOXIA中的基因
res_filtered <- res[rownames(res) %in% genemarker, ]

p <- ggplot(data = res_filtered, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 2.5) +  # 绘制散点图
  scale_color_manual(values = c('red', 'gray', 'red'), limits = c('up', 'none', 'down')) +  # 自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = '', color = '') +  # 坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), # 背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  # 添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(0, 3) + ylim(0, 250)  # 定义刻度边界

# 显示火山图
print(p)

# # 输出筛选的差异基因总表
# res1_select <- subset(res_filtered, sig %in% c('up', 'down'))
# write.table(res1_select, file = paste0(routine_fig2, 'DESeq2_select_TRM.csv'), sep = ',', col.names = NA, quote = FALSE)
