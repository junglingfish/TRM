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



################################################################################################################
########单因素cox回归模型
# combined_matrix <- read.csv(paste0(routine_fig5, 'all_clinic_data.csv'))
combined_matrix_extended <- read.csv('results/GEO_clinic_riskscore.csv')

#数据清洗
# 假设 combined_matrix_extended 是你的数据框
# 1. 删除包含 '--' 的行
combined_matrix_extended <- combined_matrix_extended %>%
  filter(!grepl("--", pathologicT))

# 2. 删除 Nx 的行
combined_matrix_extended <- combined_matrix_extended %>%
  filter(pathologicN != "Nx")

# 3. 对 pathologicT 进行赋值
combined_matrix_extended <- combined_matrix_extended %>%
  mutate(pathologicT = case_when(
    pathologicT == "T1" ~ 1,
    pathologicT == "T2" ~ 2,
    pathologicT == "T3" ~ 3,
    pathologicT == "T4" ~ 4,
    TRUE ~ NA_real_
  ))

# 4. 对 pathologicN 进行赋值
combined_matrix_extended <- combined_matrix_extended %>%
  mutate(pathologicN = case_when(
    pathologicN == "N0" ~ 0,
    pathologicN == "N1" ~ 1,
    pathologicN == "N2" ~ 2,
    pathologicN == "N3" ~ 3,
    TRUE ~ NA_real_
  ))

# 4. 对 pathologicN 进行赋值
combined_matrix_extended <- combined_matrix_extended %>%
  mutate(Tumor.grade = case_when(
    Tumor.grade == "well" ~ 1,
    Tumor.grade == "moderately" ~ 2,
    Tumor.grade == "poorly" ~ 3,
    TRUE ~ NA_real_
  ))

# # 5. 对 pathologicStage 进行规范化
# combined_matrix_extended <- combined_matrix_extended %>%
#   mutate(pathologicStage = case_when(
#     pathologicStage == "I" | grepl("^Stage I", pathologicStage) ~ 1,
#     pathologicStage == "II" | grepl("^Stage II", pathologicStage) ~ 2,
#     pathologicStage == "III" | grepl("^Stage III", pathologicStage) ~ 3,
#     grepl("^Stage IV", pathologicStage) ~ 4,
#     TRUE ~ NA_real_
#   ))
write.csv(combined_matrix_extended, 'results/GEO_clinic_data_cleaned.csv', row.names = F)


combined_matrix_extended <- read.csv('results/GEO_clinic_data_cleaned.csv', check.names = F, row.names = 1)
# 设置 p 值的阈值（可以根据需要调整）
pfilter <- 1
# 新建空白数据框
uniresult <- data.frame()  
# 提取目标变量
time <- combined_matrix_extended$time
status <- combined_matrix_extended$status
# 提取倒数6列作为因变量
variables <- combined_matrix_extended[, c(10, 14, 15, 16, 17, 18, 19)]


# 使用 for 循环对倒数6列中的每一列进行单因素 COX 分析
for (i in colnames(variables)) {   
  # 进行单因素 COX 回归分析
  unicox <- coxph(Surv(time = time, event = status) ~ variables[[i]], data = combined_matrix_extended)  
  unisum <- summary(unicox)   
  
  # 提取 p 值
  pvalue <- round(unisum$coefficients[, 5], 5) 
  
  # 检查 p 值是否为向量，进行相应处理
  if (length(pvalue) > 1) {
    pvalue <- pvalue[1]  # 只提取第一个 p 值
  }
  
  # 进行结果筛选
  if (pvalue < pfilter) { 
    uniresult <- rbind(uniresult,
                       data.frame(
                         Factor = i,
                         HR = round(unisum$coefficients[, 2], 5),
                         L95CI = round(unisum$conf.int[, 3], 5),
                         H95CI = round(unisum$conf.int[, 4], 5),
                         pvalue = pvalue
                       ))
  }
}

# 保存单因素 COX 回归分析结果
write.csv(uniresult, file = 'results/1_cox_GEO.csv', row.names = FALSE)

#####################################################################################################
##单因素cox回归森林图
data <- read.csv('results/1_cox_GEO.csv')

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


# # 确保 data_with_colnames 有正确的行数
# n_rows <- nrow(result1)
# 
# # 更新 hrzl_lines 的行号以确保在数据范围内
# hrzl_lines_list <- list()
# if (n_rows >= 1) {
#   hrzl_lines_list[["1"]] <- gpar(lty=1, lwd=2)
# }
# if (n_rows >= 2) {
#   hrzl_lines_list[["2"]] <- gpar(lty=2)
# }
# if (n_rows >= 25) {
#   hrzl_lines_list[["25"]] <- gpar(lwd=2, lty=1, columns=c(1:4))
# }

# ##设置字体
# library(showtext)
# font_add("Times New Roman", "C:/Users/junglingfish/AppData/Local/Microsoft/Windows/Fonts/Times New Roman.ttf")  # 确保路径是正确的
# showtext_auto()

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
                                   "9"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                   graphwidth=unit(0.25, "npc"),
                   xlab="Hazard Ratio",
                   xticks=c(-0.35, 0, 0.7),
                   is.summary=c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
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

# 显示森林图
print(fig5)

# 保存森林图为 PDF 文件
pdf('results/1_cox_forestplot_GEO.pdf', width = 7.5, height = 4.5)  # 设置文件名和尺寸
print(fig5)  # 打印图形
dev.off()  # 关闭 PDF 设备