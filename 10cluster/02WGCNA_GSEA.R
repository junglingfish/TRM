# 加载所需的R包
library(clusterProfiler)  # 功能丰富的基因富集分析工具
library(data.table)       # 数据处理
library(DESeq2)           # 差异表达分析
library(dplyr)            # 数据处理
library(enrichplot)       # 富集分析可视化
library(fgsea)            # GSEA分析主程序
library(GSEABase)         # 提供GSEA基础结构和函数
library(GseaVis)          # GSEA可视化工具
library(ggplot2)         # 绘图处理
library(GSVA)            # 基因集变异分析
library(homologene)      # 基因同源信息
library(limma)           # 差异分析和线性模型
library(msigdbr)         # 包含基因集合，通常与GSEA分析共同使用
library(org.Hs.eg.db)    # 人类基因注释数据库
library(patchwork)       # 组合多个ggplot2图形
library(readxl)          # 读取Excel文件
library(TCGAbiolinks)    # TCGA数据下载和处理
library(LEA)             # 进行群体遗传学分析

setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results/WGCNA/')

###########################################################################################
##DEG_DESeq2差异分析
# exp <- read.table(paste0(routine_data_GEO, 'GSE53625_tumor_count_loged.txt'), header = TRUE, row.names = 1, check.names=FALSE)
exp <- read.table(paste0(routine_data, 'GEO/GSE53625/GSE53625_tumor_count_loged.txt'), header = TRUE, row.names = 1, check.names=FALSE)
group <- read.csv('results/consensus_class_results.csv', row.names = 1)
exp <- round(2^exp)

# 确保表达矩阵的列名与分组数据的行名一致
all(rownames(group) == colnames(exp)) # 检查是否匹配

# 确保group的分组信息为因子类型
group$Cluster <- as.factor(group$Cluster)

# 构建DESeq2的数据集
dds <- DESeqDataSetFromMatrix(countData = exp, colData = group, design = ~Cluster)

# 过滤一些低表达量的基因
dds <- dds[rowSums(counts(dds)) > 1, ]

# 运行DESeq2的标准化和差异表达分析
dds1 <- DESeq(dds)

# 差异分析，注意：contrast中的group应根据你的实际分组信息进行修改，如 "Tumor" 和 "Normal"
res <- results(dds1, contrast = c('Cluster', '1', '2'))

# 查看标准化后的基本信息
summary(res)

# 输出结果为数据框
res <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

# 将结果保存到本地
write.table(res, file = 'results/WGCNA/deseq_diff.csv', col.names = NA, sep = ',', quote = FALSE)

# 筛选差异表达基因
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01), 'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01), 'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) < 1 | res1$padj >= 0.01), 'sig'] <- 'none'

write.table(res1, file = 'results/WGCNA/deseq_diff_upanddown.csv', col.names = NA, sep = ',', quote = FALSE)

# 输出筛选的差异基因总表
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'results/WGCNA/cluster.DESeq2.select.csv', sep = ',', col.names = NA, quote = FALSE)

# 根据上调和下调分别输出
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'results/WGCNA/cluster.DESeq2.up.csv', sep = ',', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'results/WGCNA/cluster.DESeq2.down.csv', sep = ',', col.names = NA, quote = FALSE)

# 使用ggplot2绘制火山图
library(ggplot2)

p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 3, alpha = 0.5) +  # 绘制散点图
  scale_color_manual(values = c('#ba3d25', 'gray', '#127cba'), limits = c('up', 'none', 'down')) +  # 自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = '', color = '') +  # 坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), # 背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  # 添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-3, 3) + ylim(0, 40)  # 定义刻度边界

# 显示火山图
print(p)



DEG_DESeq2 <- read.csv('results/WGCNA/deseq_diff.csv',row.names = 1)
DEG_DESeq2 <- DEG_DESeq2[order(DEG_DESeq2$log2FoldChange, decreasing = TRUE), ]


library(org.Hs.eg.db) # human的OrgDB
library(clusterProfiler)
# ID转化
gene_entrezid <- bitr(geneID = rownames(DEG_DESeq2), 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", # 转成ENTREZID
                      OrgDb = "org.Hs.eg.db"
)

head(gene_entrezid)

gene_entrezid$logFC <- DEG_DESeq2$log2FoldChange[match(gene_entrezid$SYMBOL, rownames(DEG_DESeq2))]
genelist = gene_entrezid$logFC
names(genelist) = gene_entrezid$ENTREZID 

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

# # 使用 filter 函数筛选 gs_name 列中包含 'HYPOXIA' 的行
# m_t2g_hypoxia <- m_t2g %>%
#   filter(grepl("HYPOXIA", gs_name))

gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH"
)


#提取显著富集的基因集
g1 <- as.data.frame(gsea_res)
g1 <- subset(g1, pvalue < 0.05, abs(NES) > 1)
g1 <- g1[order(g1$NES,decreasing = T),]

write.csv(gsea_res, 'results/WGCNA/gsea_function_cluster.csv')

# clsaasic with pvalue
p1 <- gseaNb(object = gsea_res,
             geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
             lineSize = 1,
             # addGene = mygene,
             addPval = T,
             pvalX = 0.7,pvalY = 0.75,
             pvalSize = 4,
             pCol = 'black',
             pHjust = 0)
p1

# 保存图形为 PDF 文件
ggsave("results/WGCNA/GSEA_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf", plot = p1, width = 6, height = 5, units = "in")


# clsaasic with pvalue
p2 <- gseaNb(object = gsea_res,
             geneSetID = 'HALLMARK_P53_PATHWAY',
             lineSize = 1,
             # addGene = mygene,
             addPval = T,
             pvalX = 0.7,pvalY = 0.75,
             pvalSize = 4,
             pCol = 'black',
             pHjust = 0)
p2

# 保存图形为 PDF 文件
ggsave("results/WGCNA/GSEA_HALLMARK_P53_PATHWAY.pdf", plot = p2, width = 6, height = 5, units = "in")

# clsaasic with pvalue
p3 <- gseaNb(object = gsea_res,
             geneSetID = 'HALLMARK_INTERFERON_ALPHA_RESPONSE',
             lineSize = 1,
             # addGene = mygene,
             addPval = T,
             pvalX = 0.7,pvalY = 0.75,
             pvalSize = 4,
             pCol = 'black',
             pHjust = 0)
p3

# 保存图形为 PDF 文件
ggsave("results/WGCNA/GSEA_HALLMARK_INTERFERON_ALPHA_RESPONSE.pdf", plot = p3, width = 6, height = 5, units = "in")

