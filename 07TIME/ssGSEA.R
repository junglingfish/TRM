setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results/')

library(ggpubr)
library(stringr)
library(tidyHeatmap)
library(tidyverse)
library(RColorBrewer)
library(linkET)

data <- read.table(paste0(routine_data, 'results/geo_tumor_exp_cleaned.txt'), check.names = F)
# data <- t(data)
riskscore <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)

gene <- read.csv('28_immunecells.csv')
# 使用 split 函数按细胞类型创建列表
geneset_cell_list <- split(gene$Metagene, gene$Cell.type)

# # load(paste0(routine_fig9, 'tumor_immunecell_data.rdata'))
# geneset_cell <- read.csv(paste0(routine_fig9, '16_cellgenes.CSV'))
# # geneset_cell <- load(paste0(routine_fig9, 'ssGSEA28.Rdata'))
# 
# 
# # 将矩阵转为列表格式，并去除空字符串
# geneset_cell_list <- lapply(colnames(geneset_cell), function(colname) {
#   list_content <- geneset_cell[[colname]]
#   # 使用 Filter 函数去除空字符串
#   list_content <- Filter(function(x) x != "", list_content)
#   return(list_content)
# })
# 
# # 将每个元素命名为对应的列名
# names(geneset_cell_list) <- colnames(geneset_cell)

library(GSVA)
data <- as.matrix(data)

gsvaP <- ssgseaParam(
  exprData = data,
  geneSets = geneset_cell_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)

#进行gsva分析
re <- gsva(gsvaP)        #注意表达谱exp载入后需转化为matrix，前面已转换
re <- t(re)
write.csv(re, 'results/ssgsea_result.csv')

im_ssgsea <- read.csv('results/ssgsea_result.csv', check.names = F)

ssgsea_long <- im_ssgsea %>% 
  pivot_longer(- ID,names_to = "cell_type",values_to = "Score")
head(ssgsea_long)

#筛选关键基因
genes <- colnames(riskscore)[3:10]

genes_expr <- as.data.frame(t(data[rownames(data) %in% genes,]))
genes_expr <- genes_expr[match(im_ssgsea$ID,rownames(genes_expr)),]
identical(im_ssgsea$ID,rownames(genes_expr))
## [1] TRUE

rownames(riskscore) <- riskscore$Samples
# 按行名顺序对 riskscore 进行重新排序，使其与 gene_expr 的行名一致
riskscore_aligned <- riskscore[rownames(genes_expr), "riskscore", drop = FALSE]
# 将 riskScore 列添加到 gene_expr 作为最后一列
genes_expr$riskscore <- riskscore_aligned$riskscore

cor_res <- correlate(genes_expr, im_ssgsea[,-1],method = "spearman")

# 热图1
qcorrplot(cor_res) +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

# 整理数据
df_r <- cor_res$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1, names_to = "cell_type", values_to = "correlation")

df_p <- cor_res$p %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1, names_to = "cell_type", values_to = "pvalue")

df_cor <- df_r %>% 
  left_join(df_p, by = c("gene", "cell_type")) %>% 
  mutate(stars = cut(pvalue, breaks = c(-Inf, 0.05, 0.01, 0.001, Inf), right = FALSE, labels = c("***", "**", "*", " ")))

# 设定基因顺序，使两张热图的基因顺序一致
gene_order <- unique(df_cor$gene)
# 反转基因顺序
df_cor$gene <- factor(df_cor$gene, levels = rev(c(setdiff(gene_order, "riskscore"), "riskscore")))

# 绘制第二张热图
library(ggplot2)

ggplot(df_cor, aes(cell_type, gene)) +
  geom_tile(aes(fill = correlation)) +
  geom_text(aes(label = stars), color = "black", size = 4) +
  scale_fill_gradient2(
    low = '#67B26F', high = '#F2AA9D', mid = 'white',
    limit = c(-1, 1), name = paste0("*    p < 0.05\n\n**  p < 0.01\n\n*** p < 0.001\n\nCorrelation")
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks.y = element_blank(),
    panel.background = element_blank()
  )



##7基因关系图 "OAS2", "TRIM22", "PSMB9", "BATF2", "HLA-F", "TAP1", "DDX60"
# riskscore 18 10
df_riskScore_scatter <- im_ssgsea %>% 
  mutate(riskscore = genes_expr[,"riskscore"],.before = 1) %>% 
  pivot_longer(-c(1,2),names_to = "cell_type",values_to = "score")

ggplot(df_riskScore_scatter, aes(riskscore,score))+
  geom_point()+
  geom_smooth(method = "lm",color="blue")+
  stat_cor(method = "spearman",color="red")+
  facet_wrap(~cell_type,scales = "free_y",ncol = 7)