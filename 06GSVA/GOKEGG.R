setwd('D:/ZJU-FISH/doctor/TRM/06GSVA')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'

library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例

#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.csv('DESeq2_for_GOKEGG.csv', row.names = NULL)

#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
gene <- bitr(info$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
write.csv(gene, 'gene_to_ENTREZID.csv', row.names = F)

#GO
GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pvalueCutoff = 0.05,#设定p值阈值
             qvalueCutoff = 0.05,#设定q值阈值
             readable = T)
write.csv(GO, 'results/GO.csv', row.names = F)

##气泡图
GO_dataset <- read.csv('results/GO.csv')
# 假设 GeneRatio 列的格式是 "a/b"
GO_dataset$GeneRatio <- sapply(GO_dataset$GeneRatio, function(x) {
  # 将分数转换为小数
  if (grepl("/", x)) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    return(parts[1] / parts[2])  # 计算小数
  } else {
    return(as.numeric(x))  # 如果不是分数，直接转换为数值
  }
})

# GO_BP
GOBP_dataset <- GO_dataset %>% filter(ONTOLOGY == "BP")
#按照PValue从低到高排序[升序]
GOBP_dataset <- arrange(GOBP_dataset,desc(GOBP_dataset[, 4]))
GOBP_dataset <- head(GOBP_dataset, 10)
rownames(GOBP_dataset) = 1:nrow(GOBP_dataset)
GOBP_dataset$order=factor(rev(as.integer(rownames(GOBP_dataset))),labels = rev(GOBP_dataset$Description))
#气泡图#
KEGGplot <- ggplot(GOBP_dataset, aes(y = order, x = GeneRatio)) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(color = expression(pvalue), size = "Count",
       x = "Gene Ratio", y = "", title = "GOBP Pathway") +
  theme_bw() +
  theme(axis.text = element_text(size = 16, color = "black"),  # 调整坐标轴刻度字体大小和颜色
        axis.title = element_text(size = 16, color = "black"),  # 调整坐标轴标题字体大小和颜色
        plot.title = element_text(size = 16, color = "black"),  # 调整标题字体大小和颜色
        legend.text = element_text(size = 13),  # 调整图例文本的大小
        legend.title = element_text(size = 13),  # 调整图例标题的大小
        legend.key.size = unit(0.8, "cm")) +  # 调整图例符号的大小
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +  # 设置纵轴标签换行，宽度为20个字符
  scale_x_continuous(breaks = c(0.042, 0.046, 0.050))  # 这里根据需要设置具体的三个标尺值

# 打印图形
print(KEGGplot)
# 保存为 PDF 文件
ggsave("results/GOBP_pathway_enrichment.pdf", plot = KEGGplot, width = 7, height = 6)

# GO_CC
GOBP_dataset <- GO_dataset %>% filter(ONTOLOGY == "CC")
#按照PValue从低到高排序[升序]
GOBP_dataset <- arrange(GOBP_dataset,desc(GOBP_dataset[, 4]))
GOBP_dataset <- head(GOBP_dataset, 10)
rownames(GOBP_dataset) = 1:nrow(GOBP_dataset)
GOBP_dataset$order=factor(rev(as.integer(rownames(GOBP_dataset))),labels = rev(GOBP_dataset$Description))
#气泡图#
KEGGplot <- ggplot(GOBP_dataset, aes(y = order, x = GeneRatio)) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(color = expression(pvalue), size = "Count",
       x = "Gene Ratio", y = "", title = "GOCC Pathway") +
  theme_bw() +
  theme(axis.text = element_text(size = 16, color = "black"),  # 调整坐标轴刻度字体大小和颜色
        axis.title = element_text(size = 16, color = "black"),  # 调整坐标轴标题字体大小和颜色
        plot.title = element_text(size = 16, color = "black"),  # 调整标题字体大小和颜色
        legend.text = element_text(size = 13),  # 调整图例文本的大小
        legend.title = element_text(size = 13),  # 调整图例标题的大小
        legend.key.size = unit(0.8, "cm")) +  # 调整图例符号的大小
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) + # 设置纵轴标签换行，宽度为20个字符
 scale_x_continuous(breaks = c(0.03, 0.04, 0.05))  # 这里根据需要设置具体的三个标尺值

# 打印图形
print(KEGGplot)
# 保存为 PDF 文件
ggsave("results/GOCC_pathway_enrichment.pdf", plot = KEGGplot, width = 7, height = 6)

# GO_MF
GOBP_dataset <- GO_dataset %>% filter(ONTOLOGY == "MF")
#按照PValue从低到高排序[升序]
GOBP_dataset <- arrange(GOBP_dataset,desc(GOBP_dataset[, 4]))
GOBP_dataset <- head(GOBP_dataset, 10)
rownames(GOBP_dataset) = 1:nrow(GOBP_dataset)
GOBP_dataset$order=factor(rev(as.integer(rownames(GOBP_dataset))),labels = rev(GOBP_dataset$Description))
# 气泡图
KEGGplot <- ggplot(GOBP_dataset, aes(y = order, x = GeneRatio)) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(color = expression(pvalue), size = "Count",
       x = "Gene Ratio", y = "", title = "GOMF Pathway") +
  theme_bw() +
  theme(axis.text = element_text(size = 16, color = "black"),  # 调整坐标轴刻度字体大小和颜色
        axis.title = element_text(size = 16, color = "black"),  # 调整坐标轴标题字体大小和颜色
        plot.title = element_text(size = 16, color = "black"),  # 调整标题字体大小和颜色
        legend.text = element_text(size = 13),  # 调整图例文本的大小
        legend.title = element_text(size = 13),  # 调整图例标题的大小
        legend.key.size = unit(0.8, "cm")) +  # 调整图例符号的大小
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +  # 设置纵轴标签换行，宽度为20个字符
  scale_x_continuous(breaks = c(0.025, 0.035, 0.045))  # 这里根据需要设置具体的三个标尺值

# 显示图形
print(KEGGplot)
# 保存为 PDF 文件
ggsave("results/GOMF_pathway_enrichment.pdf", plot = KEGGplot, width = 7.5, height = 6)

#####################################################################################################################
#KEGG
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
write.csv(KEGG, 'results/KEGG.csv', row.names = F)
##气泡图
KEGG_dataset <- read.csv('results/KEGG.csv')
# 假设 GeneRatio 列的格式是 "a/b"
KEGG_dataset$GeneRatio <- sapply(KEGG_dataset$GeneRatio, function(x) {
  # 将分数转换为小数
  if (grepl("/", x)) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    return(parts[1] / parts[2])  # 计算小数
  } else {
    return(as.numeric(x))  # 如果不是分数，直接转换为数值
  }
})
#按照PValue从低到高排序[升序]
KEGG_dataset <- arrange(KEGG_dataset,desc(KEGG_dataset[, 5]))
KEGG_dataset <- head(KEGG_dataset, 10)
rownames(KEGG_dataset) = 1:nrow(KEGG_dataset)
KEGG_dataset$order=factor(rev(as.integer(rownames(KEGG_dataset))),labels = rev(KEGG_dataset$Description))
#气泡图#
KEGGplot <- ggplot(KEGG_dataset, aes(y = order, x = GeneRatio)) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(color = expression(pvalue), size = "Count",
       x = "Gene Ratio", y = "", title = "KEGG Pathway") +
  theme_bw() +
  theme(axis.text = element_text(size = 16, color = "black"),  # 调整坐标轴刻度字体大小和颜色
        axis.title = element_text(size = 16, color = "black"),  # 调整坐标轴标题字体大小和颜色
        plot.title = element_text(size = 16, color = "black"),  # 调整标题字体大小和颜色
        legend.text = element_text(size = 13),  # 调整图例文本的大小
        legend.title = element_text(size = 13),  # 调整图例标题的大小
        legend.key.size = unit(0.8, "cm")) +  # 调整图例符号的大小
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) + # 设置纵轴标签换行，宽度为20个字符
  scale_x_continuous(breaks = c(0.02, 0.04, 0.06))  # 这里根据需要设置具体的三个标尺值

# 打印图形
print(KEGGplot)
# 保存为 PDF 文件
ggsave("results/KEGG_pathway_enrichment.pdf", plot = KEGGplot, width = 7, height = 6)


# #GSEA
# names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
# info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
# GSEA_input <- info_merge$Log2FoldChange
# names(GSEA_input) = info_merge$ENTREZID
# GSEA_input = sort(GSEA_input, decreasing = TRUE)
# GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.05)#GSEA富集分析
# 
# write.csv(GSEA_KEGG, paste0(routine_fig2, 'GSEA.csv'), row.names = F)