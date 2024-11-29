setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

library(ConsensusClusterPlus)
# library(doMPI) #rank设定4以上需要用这个包
library(dplyr)
library(survival)
library(survminer)
library(tidyverse)
library(glmnet)

dat <- read.table(paste0(routine_data, 'GEO/GSE53625/GSE53625_tumor_count_loged.txt'), header = TRUE, row.names = 1, check.names=FALSE)
# dat <- read.table(paste0(routine_result, 'tumor_expr_batch.txt'), header = TRUE, row.names = 1, check.names=FALSE)



gene <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)
genemarker <- colnames(gene)[3:10]

# 假设 dat 是你的数据矩阵，genemarker 是包含基因名的向量
filtered_dat <- dat[rownames(dat) %in% genemarker, ]

df <- filtered_dat

#sweep函数减去中位数进行标准化
exprSet = sweep(df,1, apply(df,1,median,na.rm=T))
par(mfrow = c(1,2))
boxplot(df[,1:20],main = "before")
boxplot(exprSet[,1:20],main = "after")

exprSet <- as.matrix(exprSet)

#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(exprSet, maxK=10, reps=500, pItem=0.8, pFeature=1, title='ConsensusCluster', clusterAlg="km", distance="euclidean",
                               corUse = 'everything',  innerLinkage = "ward.D2", seed=2222, plot="pdf")

#3 类
# 假设 results[[2]][["consensusClass"]] 是你的结果
consensus_class <- results[[2]][["consensusClass"]]



# 将结果转换为数据框
df <- data.frame(Sample = names(consensus_class), Cluster = as.numeric(consensus_class))

# 保存为 CSV 文件
write.csv(df, "results/consensus_class_results.csv", row.names = FALSE)


#########################################################################################################
##KM
# dat <- read.table(paste0(routine_result, 'tumor_expr_batch.txt'), header = TRUE, row.names = 1)
group <- read.csv('results/consensus_class_results.csv',row.names = 1)
# group <- read.csv(paste0(routine_sfig1, 'NMF_group.csv'),row.names = 1)
clinic <- read.csv(paste0(routine_data, 'clinic/GEO_clinic_data.csv'),row.names = 1)
# clinic <- read.csv(paste0(routine_data, 'clinic/all_clinic_data.csv'),row.names = 1)

##合并clinic和group数据
# 获取行名
rownames_group <- rownames(group)

# 筛选行名并进行修改
new_rownames <- sapply(rownames_group, function(x) {
  if (startsWith(x, "TCGA") && substr(x, 15, 15) == "1") {
    # 保留前12个字符并将 '.' 替换为 '-'
    return(gsub("\\.", "-", substr(x, 1, 12)))
  } else if (startsWith(x, "ec")) {
    # 保留 '_' 之前的所有内容
    return(sub("_.*", "", x))
  } else {
    return(NA)  # 删除其他行
  }
})

# 过滤掉 NA 的行名
valid_rows <- !is.na(new_rownames)
group <- group[valid_rows, , drop = FALSE]

# 设置新的行名
rownames(group) <- new_rownames[valid_rows]

# 将两个矩阵转换为数据框
group_df <- as.data.frame(group)
clinic_df <- as.data.frame(clinic)

# 添加行名作为一列
group_df$Samples <- rownames(group_df)
clinic_df$Samples <- rownames(clinic_df)

# 合并两个数据框，基于样本名
merged_df <- inner_join(group_df, clinic_df, by = "Samples")

merged_risk <- inner_join(gene, group_df, by = "Samples")

# 如果需要，将合并后的数据框转换回矩阵
merged_matrix <- as.matrix(merged_df)

write.csv(merged_risk, 'results/cluster_risk.csv', row.names = F)

# 设置合并后的行名
rownames(merged_matrix) <- merged_matrix[, "Samples"]  # 第一列为行名
merged_matrix <- merged_matrix[, -which(colnames(merged_matrix) == "Samples")]  # 删除原行名列
merged_matrix <- as.data.frame(merged_matrix)

write.csv(merged_matrix, 'results/cluster_clinic_data.csv')

#KM分析
fit <- survfit(Surv(as.numeric(time), as.numeric(status)) ~ Cluster, data = merged_matrix)
#KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 5,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = TRUE,          # 显示风险表
                       risk.table.height = 0.25,   # 风险表的高度
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = "Kaplan-Meier Curve for OS", # 大标题
                       legend.labs = c('C1', 'C2'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       risk.table.xaxis = TRUE     # 确保风险表的 x 轴对齐
)
lasso_KM

# 保存 KM 生存曲线为 PDF
ggsave("results/Consensus/cluster2_KMplot.pdf", plot = lasso_KM$plot, width = 6, height = 5)
