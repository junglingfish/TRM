setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

########################################################################
#riskscore
exp <- read.table(paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt'))
# exp <- read.table(paste0(routine_result, 'tumor_exp_cleaned.txt'), check.names = F)
# exp <- t(exp)
group <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)
rownames(group) <- group$Samples
exp <- t(exp)
# gene <- c('CD274', 'CTLA4', 'HAVCR2', 'LAG3', 'PDCD1', 'PDCD1LG2', 'TIGIT', 'SIGLEC15')
gene <- read.table('ICB.txt', check.names = F)
gene$V1 <- gsub("^HLA.*", "", gene$V1)

# 筛选 exp 中的列，保留在 gene 中出现的列
filtered_exp <- exp[, colnames(exp) %in% gene$V1]

# 创建一个新的数据框，将 group 的列添加到 filtered_exp 中
# 首先，确保 group 按照 exp 的行名顺序排列
group_ordered <- group[match(rownames(filtered_exp), rownames(group)), , drop = FALSE]

# 合并 filtered_exp 和 group
result <- data.frame(filtered_exp, group = group_ordered$Risk)

library(ggplot2)
library(ggpubr)
library(reshape2)

# 假设 result 是你的数据框
result_long <- melt(result, id.vars = "group")

# 确保 group 的值没有缺失
result_long <- result_long[!is.na(result_long$group), ]

# 绘制箱线图
p <- ggplot(result_long, aes(x = variable, y = value, color = group)) +
  geom_boxplot(fill = NA, linewidth = 0.8, outlier.colour = NA, 
               position = position_dodge(width = 0.7),  
               width = 0.4) +  
  scale_color_manual(values = c("High" = "#e92428", "Low" = "#34499d")) +
  labs(x = "", y = "Gene expression") +
  scale_x_discrete(limits = unique(result_long$variable)) +  # 添加 x 轴标尺
  scale_y_continuous(breaks = seq(0, max(result_long$value, na.rm = TRUE), by = 5)) +  # 添加 y 轴标尺
  theme_minimal() +
  theme(panel.grid = element_blank(),  # 取消背景网格线
        axis.title.y = element_text(size = 20),  # 调整 y 轴标题大小
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15, color = "black"),   # 调整 x 轴标尺大小
        axis.text.y = element_text(size = 15, color = "black"),  # 调整 y 轴标尺大小
        axis.line = element_line(colour = "black", size = 1), 
        axis.ticks = element_line(size = 1), # 调整坐标轴刻度线的长度
        axis.ticks.length = unit(0.2, "cm"))  # 设置坐标轴标尺颜色为黑色

# 添加显著性标记
p <- p + geom_signif(comparisons = list(c("High", "Low")),
                     map_signif_level = TRUE,
                     position = position_dodge(width = 0.7),
                     test = "t.test",  # 可以选择 "wilcox.test" 或其他测试
                     size = 0.8,
                     tip_length = 0.02,  # 设置标记线的长度
                     textsize = 8)       # 设置显著性标记的字体大小

# # 添加小点
# p <- p + geom_jitter(aes(color = group), 
#                      position = position_jitterdodge(dodge.width = 0.7), 
#                      size = 1.5, alpha = 0.6)  # 设置点的大小和透明度

# 添加误差线，确保与箱线图对齐
p <- p + stat_boxplot(geom = "errorbar",
                      aes(ymin = after_stat(ymax), color = group), 
                      position = position_dodge(width = 0.7),  # 确保误差线对齐
                      width = 0.3, size = 0.8) +  # 调整 width 与箱线图一致
  stat_boxplot(geom = "errorbar",
               aes(ymax = after_stat(ymin), color = group), 
               position = position_dodge(width = 0.7),  # 确保误差线对齐
               width = 0.3, size = 0.8)  # 调整 width 与箱线图一致

# 显示图形
print(p)

# 使用 ggsave() 函数保存图形
ggsave(filename = "results/ICB_cluster_boxplot.pdf",  plot = p, width = 24, height = 7, units = "in")

library(dplyr)
library(broom)

# 计算每个基因的 p 值
p_values <- result_long %>%
  group_by(variable) %>%
  summarize(p_value = wilcox.test(value ~ group)$p.value) %>%
  ungroup()

# 处理 p_values 数据框，添加显著性标记
filtered_p_values <- p_values %>%
  filter(p_value <= 0.05) %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "NS"
  ),
  p_value = ifelse(p_value < 0.001, "<0.001", p_value))

# 查看并保存结果
write.table(filtered_p_values, 'results/riskscore_p.txt', quote = F, sep = '\t', row.names = F)

