setwd('D:/ZJU-FISH/doctor/TRM/04model/')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_3 = 'D:/ZJU-FISH/doctor/TRM/03keygene/'

# tcga.dat<-read.delim(paste0(routine_data, 'GEO/GSE53625/GSE53625_count_loged.txt'),sep='\t',header = T,row.names = 1,check.names = F)
# tcga.group<-read.delim(paste0(routine_data, 'GEO/GSE53625/GSE53625_count_group.txt'), sep='\t',header = T)
# tcga.TRM.score<-read.delim(paste0(routine_3, 'results/tcga.TRM.score.txt'),sep='\t',header = T)
# tcga.cli<-read.csv(paste0(routine_data, 'clinic/GEO_clinic_data.csv'), check.names = F)

# tcga.module.dat<- tcga.dat[rownames(tcga.cox.fit),tcga.cli$Samples]
# rownames(tcga.module.dat)=gsub('-','__',rownames(tcga.module.dat))

#################################################################################################################
### LASSO回归分析 ----
# 假设 uniresult 已经是单因素 Cox 回归分析中 p 值 < 0.05 的结果
# 提取包含显著基因的列名
# selected_genes <- uniresult$gene
selected_genes <- t(tcga.module.dat)

# 提取各样本的生存状态、生存时间，以及单因素 Cox 回归分析中 p 值 < 0.05 的基因在各样本中的表达情况
# unigene <- tcga.cli[, c("time", "status", selected_genes)]
# unigene_filtered <- unigene
# # 假设 unigene 是一个数据框，并且有一个名为 time 的列
# if ("time" %in% colnames(unigene)) {
#   # 查找 time 列中值为 0 的行索引
#   zero_time_indices <- which(unigene$time == 0)
#
#   # 删除这些行
#   unigene_filtered <- unigene[-zero_time_indices, ]
#
#   # 打印结果以确认
#   print(unigene_filtered)
# } else {
#   stop("The 'time' column does not exist in the data frame 'unigene'.")
# }
y <- data.matrix(Surv(time = tcga.cli$time, event = tcga.cli$status))
# 获取数据框的列数
# n <- ncol(unigene_filtered)
# 选择倒数 24 列的内容
exp <- selected_genes
fit <- glmnet(exp, y, family = 'cox', alpha = 1)
plot(fit, xvar = 'lambda')

best_seed <- NA
best_genes <- NULL
results <- list()  # 用于存储符合条件的seed值及对应的显著基因数量

for (seed in 1:5000) {
  set.seed(seed)
  
  # 交叉验证
  lasso_fit <- cv.glmnet(as.matrix(exp), y, family = 'cox', type.measure = 'deviance')
  
  # 筛选变量
  coefficient <- coef(lasso_fit, s = lasso_fit$lambda.min)
  Active.Index <- which(as.numeric(coefficient) != 0)
  sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
  
  # 检查显著基因数量是否在1到10之间
  num_genes <- length(sig_gene_multi_cox)
  # if (num_genes > 0 && num_genes <= 10) 
  if (num_genes == 6) 
    {
    results[[length(results) + 1]] <- list(seed = seed, num_genes = num_genes)
  }
}

# 输出符合条件的seed值及对应的显著基因数量
cat("Seeds with 1-10 significant genes:\n")
for (result in results) {
  cat("Seed:", result$seed, "Number of significant genes:", result$num_genes, "\n")
}



set.seed(551)
#交叉验证
lasso_fit <- cv.glmnet(as.matrix(exp), y, family = 'cox', type.measure = 'deviance')
plot(lasso_fit, label = T)
#筛选变量
coefficient <- coef(lasso_fit, s=lasso_fit$lambda.min)
Active.Index <- which(as.numeric(coefficient) != 0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
sig_gene_multi_cox#返回基因

#计算风险评分
# 将稀疏矩阵转换为普通矩阵
coefficient_matrix <- as.matrix(coefficient)
# 将矩阵转换为数据框
coefficient_df <- as.data.frame(coefficient_matrix)
# 保留第一列中非零的行
filtered_coefficient_df <- coefficient_df[coefficient_df[, 1] != 0, ]

######riskScore 二分绘制KM##########
# 筛选出在 sig_gene_multi_cox 中的基因，并确保这些基因在表达矩阵中存在
common_genes <- intersect(sig_gene_multi_cox, colnames(combined_matrix_1))
# 提取这些基因的表达数据
selected_expression_matrix <- combined_matrix_1[, common_genes]
# 提取这些基因对应的系数
gene_coefficients <- filtered_coefficient_df

# 计算每个样本的风险评分
risk_scores <- rowSums(sweep(selected_expression_matrix, 2, gene_coefficients, FUN = "*"))
# 创建结果矩阵
risk_scores_matrix <- data.frame(Sample = rownames(combined_matrix_1), RiskScore = risk_scores)
riskScore <- risk_scores_matrix
# 检查行名是否一致
identical(rownames(riskScore), rownames(combined_matrix_extended))
unigene <- combined_matrix_1[, c("time", "status", common_genes)]

# 将 riskScore 附加到 combined_matrix_extended 后面
riskScore_cli <- cbind(unigene, riskScore = riskScore$RiskScore)

write.csv(riskScore_cli, paste0(routine_fig3, 'riskscore_GEO_WGCNA.csv'))