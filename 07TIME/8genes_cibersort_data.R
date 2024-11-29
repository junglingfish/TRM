setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

risk <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)
keygene <- colnames(risk)[3:10]
exp <- read.table(paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt'), check.names = F)
cibersort <- read.csv("results/geo_cibersort_result.csv", row.names = 1)
cibersort <- cibersort[, -((ncol(cibersort)-2):ncol(cibersort))]

exp <- as.matrix(exp)
exp <- t(exp)
# keygene <- rownames(gene)


library(dplyr)
library(tibble)

# 假设 exp 是基因表达矩阵，keygene 是包含基因名称的向量，cibersort 是免疫浸润情况的数据框
# 示例：
# exp <- matrix(rnorm(100), nrow = 10, dimnames = list(paste0("Sample", 1:10), paste0("Gene", 1:10)))
# keygene <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
# cibersort <- data.frame(Sample = paste0("Sample", 1:10), Value = rnorm(10))

# 假设 exp 和 cibersort 已经定义，并且 keygene 是基因名的字符向量
for (gene in keygene) {
  # 获取该基因的表达数据，并确保转换为数值型
  gene_expression <- as.numeric(exp[, gene, drop = FALSE])
  
  # 检查是否有有效的数值
  if (length(gene_expression) == 0 || all(is.na(gene_expression))) {
    warning(paste("Gene", gene, "has no valid expression data."))
    next
  }
  
  # 计算中位数
  median_expression <- median(gene_expression, na.rm = TRUE)
  
  # 按照中位数划分高低表达组
  group <- ifelse(gene_expression > median_expression, "High", "Low")
  
  # 将分组信息转换为数据框
  group_df <- data.frame(group = group, row.names = rownames(exp))
  
  # 添加分组到 cibersort 中，注意行名匹配
  cibersort_with_group <- cbind(cibersort[rownames(cibersort) %in% rownames(group_df), ], 
                                group = group_df[rownames(group_df) %in% rownames(cibersort), ])
  
  # 输出 CSV 文件
  write.csv(cibersort_with_group, file = paste0("results/gene/", gene, ".csv"), row.names = TRUE)
}


