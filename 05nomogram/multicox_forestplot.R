setwd('D:/ZJU-FISH/doctor/TRM/05nomogram')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'

library(survival)
library(survminer)
library(readxl)
library(tidyverse)
library(glmnet)
library(survival)
library(survminer)
library(tableone)  
library(forestplot)

#################################################################################################################
### 多因素 Cox 回归分析 ----
unigene <- read.csv('results/GEO_clinic_data_cleaned_multi.csv', row.names = 1)

# 进行多因素 Cox 回归分析
multicox <- coxph(Surv(time = time, event = status) ~ ., data = unigene)
multisum <- summary(multicox)
pfilter_s <- 1
# 提取所有基因的多因素 Cox 回归分析结果至 multiresult 对象中
Factor <- colnames(unigene)[1:(ncol(unigene)-2)]
HR <- multisum$coefficients[, 2]
L95CI <- multisum$conf.int[, 3]
H95CI <- multisum$conf.int[, 4]
pvalue <- multisum$coefficients[, 5]

# 创建结果数据框
multiresult <- data.frame(
  Factor = Factor,
  HR = round(HR, 6),
  L95CI = round(L95CI, 6),
  H95CI = round(H95CI, 6),
  pvalue = round(pvalue, 6)
)
# 呈现所有的基因
multiresult <- multiresult[multiresult$pvalue < pfilter_s,]
# 保存多因素 Cox 回归分析结果
write.csv(multiresult, file = 'results/multi_cox_allfactor_GEO.csv', row.names = FALSE)


#####################################################################################################
##多因素cox回归森林图
data <- read.csv('results/multi_cox_allfactor_GEO.csv')

# 计算 log10 值
data$log10HR <- log10(data$HR)
data$log10L95CI <- log10(data$L95CI)
data$log10H95CI <- log10(data$H95CI)

# 使用临时变量保留原始的 log10L95CI 和 log10H95CI 值
data$original_log10L95CI <- data$log10L95CI
data$original_log10H95CI <- data$log10H95CI

# 比较正负并处理绝对值的大小
data$log10L95CI <- with(data, ifelse(
  (original_log10L95CI > 0 & original_log10H95CI > 0) | (original_log10L95CI < 0 & original_log10H95CI < 0), 
  ifelse(abs(original_log10L95CI) > abs(original_log10H95CI), 
         formatC(original_log10H95CI, format = "f", digits = 4), 
         formatC(original_log10L95CI, format = "f", digits = 4)),
  formatC(original_log10H95CI, format = "f", digits = 4)
))

data$log10H95CI <- with(data, ifelse(
  (original_log10L95CI > 0 & original_log10H95CI > 0) | (original_log10L95CI < 0 & original_log10H95CI < 0), 
  ifelse(abs(original_log10L95CI) > abs(original_log10H95CI), 
         formatC(original_log10L95CI, format = "f", digits = 4), 
         formatC(original_log10H95CI, format = "f", digits = 4)),
  formatC(original_log10L95CI, format = "f", digits = 4)
))

# 格式化 HR 和置信区间到四位小数
data$log10HR <- formatC(data$log10HR, format = "f", digits = 4)

# 确保 log10L95CI 和 log10H95CI 都已经计算好
data$`log10HR(95% CI)` <- with(data, 
                               ifelse(
                                 log10L95CI < 0 & log10H95CI < 0,  # 检查是否全为负数
                                 paste0(
                                   formatC(log10HR, format = "f", digits = 4),  # 格式化 log10HR
                                   " (", 
                                   formatC(log10H95CI, format = "f", digits = 4),  # 直接使用 log10H95CI
                                   ", ", 
                                   formatC(log10L95CI, format = "f", digits = 4),  # 直接使用 log10L95CI
                                   ")"
                                 ),
                                 paste0(
                                   formatC(log10HR, format = "f", digits = 4),  # 格式化 log10HR
                                   " (", 
                                   ifelse(log10L95CI < log10H95CI, 
                                          formatC(log10L95CI, format = "f", digits = 4),  # 小的值
                                          formatC(log10H95CI, format = "f", digits = 4)),  # 否则放大的值
                                   ", ", 
                                   ifelse(log10L95CI < log10H95CI, 
                                          formatC(log10H95CI, format = "f", digits = 4),  # 大的值
                                          formatC(log10L95CI, format = "f", digits = 4)),  # 否则放小的值
                                   ")"
                                 )
                               )
)

# 检查结果
print(data$`log10HR(95% CI)`)


# 修改 pvalue 列
# 使用 ifelse 修改 pvalue 列，并确保格式为小数
data$pvalue <- ifelse(data$pvalue < 0.0001, 
                      '<0.0001', 
                      sprintf("%.4f", data$pvalue))
# 提取列名
column_names <- colnames(data)

# 将列名转换为一个矩阵行，并且设置列名为空，防止列名冲突
column_names_row <- matrix(column_names, nrow = 1)
colnames(column_names_row) <- NULL  # 确保没有列名

# 将列名行和原始矩阵转换为相同的数据类型
# 将 data 转换为 matrix
data_matrix <- as.matrix(data)

# 合并列名行和原始矩阵
data_with_colnames <- rbind(column_names_row, data_matrix)

result1 <- data_with_colnames

# 绘制森林图
fig5 <- forestplot(result1[,c(1,11,5)], 
                   mean=result1[,6],   
                   lower=result1[,9],  
                   upper=result1[,10], 
                   zero=0,           
                   boxsize=0.4,      
                   graph.pos="right",
                   hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                   "2" = gpar(lty=2),
                                   "7"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                   graphwidth=unit(0.25, "npc"),
                   xlab="Hazard Ratio",
                   xticks=c(-0.1, 0, 0.65),
                   is.summary=c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                   txt_gp=fpTxtGp(
                     label = gpar(cex = 1),
                     ticks = gpar(cex = 1), 
                     xlab = gpar(cex = 1.2), 
                     title = gpar(cex = 2)
                   ),
                   lwd.zero=2,
                   lwd.ci=2.5,
                   lwd.xaxis=3, 
                   lty.ci=1.5,
                   ci.vertices=TRUE,
                   ci.vertices.height=0.2, 
                   clip=c(0, 0.4),
                   ineheight=unit(6, 'mm'), 
                   line.margin=unit(6, 'mm'),
                   colgap=unit(14, 'mm'),
                   fn.ci_norm="fpDrawDiamondCI", 
                   col=fpColors(box='#ca7465', 
                                lines='black', 
                                zero="#7ac6ce")) 


print(fig5)  # 打印图形

# 保存森林图为 PDF 文件
pdf('results/multi_forestplot_GEO.pdf', width = 7.5, height = 4)  # 设置文件名和尺寸
print(fig5)  # 打印图形
dev.off()  # 关闭 PDF 设备
