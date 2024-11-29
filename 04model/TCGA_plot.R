setwd('D:/ZJU-FISH/doctor/TRM/04model/')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_2 = 'D:/ZJU-FISH/doctor/TRM/03keygene/'

##### ROC curve all data
# 加载必要的包
library(survival)
library(timeROC)
library(pROC)
library(survminer)

# 加载必要的包
library(pROC)
library(survival)
library(timeROC)

# 读取数据
riskScore_cli <- read.delim('results/tcga_risk.txt', sep='\t',header = T,check.names = F)

# 定义时间点（以年为单位）
time_2_years <- 1
time_3_years <- 3
time_5_years <- 5
# 创建生存状态矩阵
survival_status_2_years <- ifelse(riskScore_cli$time > time_2_years, 0, riskScore_cli$status)
survival_status_3_years <- ifelse(riskScore_cli$time > time_3_years, 0, riskScore_cli$status)
survival_status_5_years <- ifelse(riskScore_cli$time > time_5_years, 0, riskScore_cli$status)
# 进一步调整生存状态
survival_status_2_years[survival_status_2_years == 1] <- 1  # 死亡标记为1
survival_status_2_years[survival_status_2_years == 0 & riskScore_cli$time <= time_2_years] <- 0  # 存活标记为0
survival_status_3_years[survival_status_3_years == 1] <- 1
survival_status_3_years[survival_status_3_years == 0 & riskScore_cli$time <= time_3_years] <- 0
survival_status_5_years[survival_status_5_years == 1] <- 1
survival_status_5_years[survival_status_5_years == 0 & riskScore_cli$time <= time_5_years] <- 0

# 将结果转化为数据框
survival_status_2_years_df <- data.frame(Survival_Status = survival_status_2_years)
survival_status_3_years_df <- data.frame(Survival_Status = survival_status_3_years)
survival_status_5_years_df <- data.frame(Survival_Status = survival_status_5_years)


# 计算 ROC 曲线
roc_2year <- roc(survival_status_2_years ~ riskScore_cli$riskscore, 
                 data = riskScore_cli)
roc_3year <- roc(survival_status_3_years ~ riskScore_cli$riskscore, 
                 data = riskScore_cli)
roc_5year <- roc(survival_status_5_years ~ riskScore_cli$riskscore, 
                 data = riskScore_cli)

# 使用 smooth 函数平滑 ROC 曲线
roc_2year_smoothed <- smooth(roc_2year)
roc_3year_smoothed <- smooth(roc_3year)
roc_5year_smoothed <- smooth(roc_5year)

# 计算 AUC 和 95% 置信区间
auc_2year <- auc(roc_2year_smoothed)
auc_3year <- auc(roc_3year_smoothed)
auc_5year <- auc(roc_5year_smoothed)

ci_2year <- ci.auc(roc_2year_smoothed)
ci_3year <- ci.auc(roc_3year_smoothed)
ci_5year <- ci.auc(roc_5year_smoothed)

# 设置绘图区域
plot(roc_2year_smoothed, col = "#00468b", lwd = 2, main = "Time-dependent ROC Curves",
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     cex.axis = 1.25, cex.lab = 1.25, font.lab = 2, font.axis = 2, xaxt = 'n')
# 添加背景网格
grid(col = "lightgray")  # 设置网格颜色
plot(roc_3year_smoothed, add = TRUE, col = "#ed0000", lwd = 2)
plot(roc_5year_smoothed, add = TRUE, col = "#42b540", lwd = 2)

# 添加图例
legend(x = 0.5, y = 0.15,
       legend = c(paste0("1-Year AUC: ", round(auc_2year, 3), " (95% CI: ", round(ci_2year[1], 3), "-", round(ci_2year[2], 3), ")"),
                  paste0("3-Year AUC: ", round(auc_3year, 3), " (95% CI: ", round(ci_3year[1], 3), "-", round(ci_3year[2], 3), ")"),
                  paste0("5-Year AUC: ", round(auc_5year, 3), " (95% CI: ", round(ci_5year[1], 3), "-", round(ci_5year[2], 3), ")")),
       col = c("#00468b", "#ed0000", "#42b540"),
       lwd = 2,
       cex = 1,        # 图例字体大小
       text.font = 2,
       bty = "n")    # 图例字体加粗

# 修改 x 轴的刻度和标签
axis(1, at = seq(0, 1, by = 0.2), labels = round(1 - seq(0, 1, by = 0.2), 2), cex.axis = 1.25, font.axis = 2, line = 1)
roc_plot <- recordPlot()
# 打开 PDF 设备
pdf("results/TCGA_ROCplot.pdf", width = 6, height = 5.25)
# 重放已记录的绘图
replayPlot(roc_plot)
# 关闭 PDF 设备
dev.off()


fit <- survfit(Surv(time, as.numeric(status)) ~ Risk, data=riskScore_cli)

# KM生存曲线
lasso_KM <- ggsurvplot(
  fit,
  pval = TRUE,                # 显示 p 值
  pval.coord = c(0, 0.2),     # p 值位置坐标
  pval.size = 5,              # p 值字体大小
  conf.int = TRUE,            # 显示生存率的 95% CI
  risk.table = TRUE,          # 显示风险表
  risk.table.height = 0.25,   # 风险表的高度
  palette = c("#ed0000", "#00468b"), # 对调颜色
  title = "Kaplan-Meier Curve for OS", # 大标题
  legend.labs = c('high', 'low'), # 标签对调
  legend.title = "",          # 改图例名称
  surv.median.line = "hv",    # 中位生存期线
  risk.table.xaxis = TRUE     # 确保风险表的 x 轴对齐
)

p1=lasso_KM$plot+theme_bw()+
  theme(axis.text.y=element_text(family="Times",face="plain"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.title = element_text(family="Times",face="plain"),
        legend.text = element_text(family="Times",face="plain"))
p2=lasso_KM$table+theme_bw()+
  theme(axis.text.y=element_text(family="Times",face="plain"),
        plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
        plot.title=element_blank(),
        legend.position=c(1,1), 
        legend.justification=c(1,1),
        legend.title = element_text(family="Times",face="plain"),
        legend.text = element_text(family="Times",face="plain"))

g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
ggsave(filename = 'results/TCGA_KMplot.pdf', plot = g2, width = 6, height = 5.25, units = "in")



####riskplot
###三图联合绘制
#################################################################################################################
# 加载必要的包
data_frame <- riskScore_cli
data_frame$status <- as.character(data_frame$status)


# 将数据框按照风险评分排序
data_frame <- data_frame[order(data_frame$riskscore), ]
common_gene <- colnames(data_frame)[3:10] # 假设基因列从第3列到第11列
# 对基因表达数据进行Z-score标准化
data_frame[, common_gene] <- apply(data_frame[, common_gene], 2, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
})

# # 将标准化后的数据限制在 -3 到 3 之间
# data_frame[, common_gene] <- apply(data_frame[, common_gene], 2, function(x) {
#   pmax(pmin(x, 2), -2)  # 超过范围的值设为边缘值
# })

# 定义x轴范围
x_range <- c(1, nrow(data_frame))

# 1. 风险评分分布图
risk_plot <- ggplot(data_frame, aes(x = seq_along(riskscore), y = riskscore, color = Risk)) +
  geom_point() +
  labs(y = "Risk Score") +
  theme_minimal() +
  scale_color_manual(values = c("Low" = "#34499d", "High" = "#e92428")) +
  theme(axis.text.x = element_blank(), # 隐藏x轴的文本
        axis.ticks.x = element_line(color = "black"), # 显示x轴刻度线
        axis.text.y = element_text(color = "black"), # 显示y轴刻度线的文本
        axis.ticks.y = element_line(color = "black"), # 显示y轴刻度线
        axis.title.x = element_blank(), # 隐藏x轴标题
        panel.grid = element_blank(), # 隐藏网格线
        axis.line = element_line(color = "black", size = 1)) + # 加粗黑色轴线
  scale_x_continuous(limits = x_range) # 设置x轴范围

# 2. 生存状态图
surv_plot <- ggplot(data_frame, aes(x = seq_along(time), y = time, color = factor(status))) +
  geom_point() +
  labs(y = "Survival Time") +
  theme_minimal() +
  scale_color_manual(values = c("1" = "#e92428", "0" = "#34499d"), labels = c("1" = "Dead", "0" = "Alive")) +
  theme(axis.text.x = element_blank(), # 隐藏x轴的文本
        axis.ticks.x = element_line(color = "black"), # 显示x轴刻度线
        axis.text.y = element_text(color = "black"), # 显示y轴刻度线的文本
        axis.ticks.y = element_line(color = "black"), # 显示y轴刻度线
        axis.title.x = element_blank(), # 隐藏x轴标题
        panel.grid = element_blank(), # 隐藏网格线
        axis.line = element_line(color = "black", size = 1)) + # 加粗黑色轴线
  scale_x_continuous(limits = x_range) # 设置x轴范围

# 3. 基因表达热图
# 首先，重塑数据框以适合绘制热图
gene_data <- data_frame[, c("Risk", common_gene, "Samples")]
gene_data_melt <- melt(gene_data, id.vars = c("Risk", "Samples"))

# 将样本顺序排序并应用于热图数据
gene_data_melt$Sample <- factor(gene_data_melt$Samples, levels = data_frame$Samples)
gene_data_melt$variable <- factor(gene_data_melt$variable, levels = common_gene)

# 绘制热图
heatmap_plot <- ggplot(gene_data_melt, aes(x = Sample, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#34499d", "white", "#e92428")) +
  labs(x = "Samples", y = "Genes") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), # 隐藏 x 轴上的文本
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

# 确保所有图的x轴范围一致
x_range <- c(1, nrow(data_frame))

# 更新每个图以使用相同的x轴范围
risk_plot <- risk_plot + coord_cartesian(xlim = x_range)
surv_plot <- surv_plot + coord_cartesian(xlim = x_range)
heatmap_plot <- heatmap_plot + coord_cartesian(xlim = x_range)

# 将三部分组合在一起
p3 <- grid.arrange(risk_plot, surv_plot, heatmap_plot, ncol = 1, heights = c(2, 2, 2))
ggsave(filename = 'results/TCGA_riskplot.pdf', plot = p3, width = 6, height = 8, units = "in")
