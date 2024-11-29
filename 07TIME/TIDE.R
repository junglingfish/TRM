setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
routine_7 = 'D:/ZJU-FISH/doctor/TRM/07TIME/results/gene/'
# dir.create('results/TIME/TIDE/')

tide <- read.csv('TIDE.csv', row.names = 1, check.names = F)
riskscore_group <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'), sep='\t',header = T,check.names = F)
rownames(riskscore_group) <- riskscore_group$Samples

# 获取原始行名
original_row_names <- rownames(tide)

# 初始化一个逻辑向量，标记要保留的行
to_keep <- rep(FALSE, length(original_row_names))

# # 处理以 'TCGA' 开头的行名
# for (i in seq_along(original_row_names)) {
#   if (grepl("^TCGA", original_row_names[i])) {
#     if (substr(original_row_names[i], 15, 15) == "1") {
#       to_keep[i] <- TRUE  # 保留第15位为 '1' 的行
#       original_row_names[i] <- gsub("\\.", "-", substr(original_row_names[i], 1, 12))  # 替换 . 为 -
#     }
#   }
# }

# 处理以 'ec' 开头的行名
for (i in seq_along(original_row_names)) {
  if (grepl("^ec", original_row_names[i])) {
    to_keep[i] <- TRUE
    original_row_names[i] <- sub("_.*", "", original_row_names[i])  # 保留 '_' 之前的内容
  }
}

# 从tide矩阵中删除不需要的行
tide <- tide[to_keep, ]

# 更新行名
rownames(tide) <- original_row_names[to_keep]

# 获取tide的行名
tide_row_names <- rownames(tide)

# 将riskscore_group按照tide的行名顺序排列
riskscore_group_ordered <- riskscore_group$Risk[match(tide_row_names, rownames(riskscore_group))]
# 将riskscore_group添加到tide矩阵的最后一列
tide <- cbind(tide, Risk = riskscore_group_ordered)

# 创建一个向量存储文件名
files <- colnames(riskscore_group)[3:10]
group_data_list <- list()  # 用于存储group列

# 获取tide的行名
tide_row_names <- rownames(tide)

# 循环读取每个文件并提取group列
for (file in files) {
  # 读取CSV文件
  df <- read.csv(paste0(routine_7, file, '.csv'), row.names = 1, check.names = FALSE)
  
  # 提取group列并根据tide的行名排序
  group_data <- df$group[match(tide_row_names, rownames(df))]
  
  # 将提取的group列存入列表
  group_data_list[[file]] <- group_data
}

# 将所有group列合并为一个数据框，按列合并
combined_group_data <- do.call(cbind, group_data_list)

# 修改列名为对应的变量名
colnames(combined_group_data) <- files

# 检查结果
head(combined_group_data)

# 将tide和combined_group_data按列合并
final_data <- cbind(tide, combined_group_data)
write.csv(final_data, 'results/TIDE_group.csv')


tide_matrix <- read.csv('results/TIDE_group.csv', row.names = 1, check.names = F)


# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)


tide_data <- tide_matrix %>%
  select(1, 2:10) %>%  # 根据需要选择列
  pivot_longer(cols = -1, names_to = "Group", values_to = "Status")  # 转化为长格式


# 确保按原顺序排列
tide_data$Group <- factor(tide_data$Group, levels = colnames(tide_matrix)[2:10])

# 绘制箱线图
p <- ggplot(tide_data, aes(x = Group, y = TIDE, color = Status)) +
  geom_boxplot(fill = NA, linewidth = 1,  
               position = position_dodge(width = 0.7),  
               width = 0.5) + # 设置边框的粗细
  scale_color_manual(values = c("High" = "#e92428", "Low" = "#34499d")) +
  labs(x = "Group", y = "TIDE") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold", color = 'black'),
    axis.text.y = element_text(size = 16, face = "bold", colour = 'black'),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    panel.grid.major = element_blank(),  # 删除主网格线
    panel.grid.minor = element_blank(),   # 删除次网格线
    axis.line = element_line(size = 1.25, color = "black"),  # 添加坐标轴线
    axis.ticks = element_line(size = 1),  # 添加刻度线
    axis.ticks.length = unit(0.25, "cm")  # 调节刻度线长度
  )

# 添加误差线，确保与箱线图对齐
p <- p + stat_boxplot(geom = "errorbar",
                      aes(ymin = after_stat(ymax), color = Status), 
                      position = position_dodge(width = 0.7),  # 确保误差线对齐
                      width = 0.3, size = 1) +  # 调整 width 与箱线图一致
  stat_boxplot(geom = "errorbar",
               aes(ymax = after_stat(ymin), color = Status), 
               position = position_dodge(width = 0.7),  # 确保误差线对齐
               width = 0.3, size = 1)  # 调整 width 与箱线图一致

# 显示图形
print(p)

ggsave("results/TIDE_boxplot.pdf", plot = p, width = 10, height = 6)

library(dplyr)
library(broom)

# 计算每个基因的 p 值
p_values <- tide_data %>%
  group_by(Group) %>%
  summarize(p_value = wilcox.test(TIDE ~ Status)$p.value) %>%
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
write.table(filtered_p_values, 'results/TIDE_p.txt', quote = F, sep = '\t', row.names = F)
