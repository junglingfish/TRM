library(survival)
library(survminer)
library(readxl)
library(tidyverse)
library(glmnet)

setwd('D:/ZJU-FISH/doctor/TRM/04model/')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_2 = 'D:/ZJU-FISH/doctor/TRM/03keygene/'
dir.create('results/clinic_KM')

data <- read.csv('GEO_clinic_riskscore.csv', row.names = 1,  check.names = F)

# 1. 删除包含 '--' 的行
data_cleaned <- data[!apply(data, 1, function(row) any(grepl('--', row))), ]

# 2. 删除 Nx 的行
data_cleaned <- data_cleaned[data_cleaned$pathologicN != 'Nx', ]

# # 3. 对 pathologicStage 进行规范化
# data_cleaned$pathologicStage <- ifelse(
#   data_cleaned$pathologicStage %in% c('I', 'Stage IA', 'Stage IB'), 1,
#   ifelse(
#     data_cleaned$pathologicStage %in% c('II', 'Stage II', 'Stage IIA', 'Stage IIB'), 2,
#     ifelse(
#       data_cleaned$pathologicStage %in% c('III', 'Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC'), 3,
#       4
#     )
#   )
# )

data <- data_cleaned


##################################################################################################
### Age < 60
# 筛选 age 列中小于 60 的行
age_60 <- data[data$age < 60, ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=age_60)

# pdf("results/clinic_KM/Age_60_1.pdf"), width = 6, height = 4.8)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for Age < 60") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig1 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))
# lasso_KM

ggsave("results/clinic_KM/Age_60_1.pdf", plot = fig1, width = 6, height = 4.5)


##################################################################################################
### Age >= 60
# 筛选 age 列中大于等于 60 的行
age_60 <- data[data$age >= 60, ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=age_60)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for Age ≥ 60") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig2 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/Age_60_2.pdf", plot = fig2, width = 6, height = 4.5)


##################################################################################################
### Male
# 筛选 age 列中大于等于 60 的行
male <- data[data$gender == 'male', ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=male)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for Male") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig3 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/Male.pdf", plot = fig3, width = 6, height = 4.5)


##################################################################################################
### Female
# 筛选 age 列中大于等于 60 的行
female <- data[data$gender == 'female', ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=female)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for Female") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig4 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/Female.pdf", plot = fig4, width = 6, height = 4.5)


##################################################################################################
### Stage 1
# 筛选 age 列中大于等于 60 的行
stage1 <- data[data$pathologicStage == 1, ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=stage1)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for Stage I") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig5 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/Stage I.pdf", plot = fig5, width = 6, height = 4.5)


##################################################################################################
### Stage 2
# 筛选 age 列中大于等于 60 的行
stage2 <- data[data$pathologicStage == 2, ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=stage2)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for Stage II") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig6 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/Stage II.pdf", plot = fig6, width = 6, height = 4.5)


##################################################################################################
### Stage 3-4
# 筛选 age 列中大于等于 60 的行
stage34 <- data[data$pathologicStage >= 3, ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=stage34)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for Stage III-IV") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig7 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/Stage III-IV.pdf", plot = fig7, width = 6, height = 4.5)


##################################################################################################
### T1-2
# 筛选 pathologicT 列为 T1 和 T2 的行
T12 <- data[data$pathologicT %in% c("T1", "T2"), ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=T12)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for T1-2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig8 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/T1-2.pdf", plot = fig8, width = 6, height = 4.5)


##################################################################################################
### T3-4
# 筛选 pathologicT 列为 T3 和 T4 的行
T34 <- data[data$pathologicT %in% c("T3", "T4"), ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=T34)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for T3-4") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig9 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/T3-4.pdf", plot = fig9, width = 6, height = 4.5)


##################################################################################################
### N0
# 筛选 pathologicN 列为 N0
N0 <- data[data$pathologicN %in% c("N0"), ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=N0)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for N0") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig10 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/N0.pdf", plot = fig10, width = 6, height = 4.5)



##################################################################################################
### N1,2,3
# 筛选 pathologicN 列为 N1,2,3
N123 <- data[data$pathologicN %in% c("N1", "N2", "N3"), ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=N123)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for N1-3") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig11 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/N1-3.pdf", plot = fig11, width = 6, height = 4.5)


##################################################################################################
### Tumor.grade
well <- data[data$Tumor.grade %in% c("well"), ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=well)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for well") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig12 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/well.pdf", plot = fig12, width = 6, height = 4.5)


##################################################################################################
### Tumor.grade
moderately <- data[data$Tumor.grade %in% c("moderately"), ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=moderately)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for moderately") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig13 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/moderately.pdf", plot = fig13, width = 6, height = 4.5)


##################################################################################################
### Tumor.grade
poorly <- data[data$Tumor.grade %in% c("poorly"), ]
#KM分析
fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=poorly)
# KM生存曲线
lasso_KM <- ggsurvplot(fit,
                       pval = TRUE,                # 显示 p 值
                       pval.coord = c(0, 0.2),     # p 值位置坐标
                       pval.size = 7,              # p 值字体大小
                       conf.int = TRUE,            # 显示生存率的 95% CI
                       risk.table = FALSE,         # 不显示风险表
                       palette = c("#ed0000", "#00468b"), # 对调颜色
                       title = NULL,               # 不在这里设置标题
                       legend.labs = c('high' , 'low'), # 标签对调
                       legend.title = "",          # 改图例名称
                       surv.median.line = "hv",    # 中位生存期线
                       xlab = "Time (years)"        # 设置 x 轴标题
)

# 添加标题并设置格式
lasso_KM$plot <- lasso_KM$plot + 
  ggtitle("Survival Proportions for poorly") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 23),  # 调整标题大小
        axis.title.x = element_text(size = 20),             # x轴标题加粗并调整大小
        axis.title.y = element_text(size = 20),             # y轴标题加粗并调整大小
        legend.text = element_text(size = 15))              # 图例文本加粗并调整大小

# 显示 p 值加粗
fig14 <- lasso_KM$plot + 
  theme(plot.subtitle = element_text(face = "bold"))

# # 显示生存曲线
# lasso_KM
ggsave("results/clinic_KM/poorly.pdf", plot = fig14, width = 6, height = 4.5)



fig<-ggarrange(fig1,fig2,fig3,fig4,fig5,
               fig6,fig7,fig8,fig9,fig10,
               fig11,fig12,fig13,fig14,
               nrow = 5,ncol = 3,labels = '',
               heights = c(1, 1, 1, 1, 1), # Adjusts spacing between rows
               widths = c(1, 1, 1)         # Adjusts spacing between columns
               )
fig
ggsave('results/clinic_KMplot.pdf',fig,height = 28,width = 18)
