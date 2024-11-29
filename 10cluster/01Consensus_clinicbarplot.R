# 加载必要的库
library(ggplot2)
library(dplyr)

setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

data <- read.csv('results/cluster_clinic_data.csv', row.names = 1)

# 假设数据存储在 box_data 中，将其转换为数据框方便数据处理
box_data <- as.data.frame(data)

# 1. 删除包含 '--' 的行
box_data_cleaned <- box_data[!apply(box_data, 1, function(row) any(grepl('--', row))), ]

# 2. 删除 Nx 的行
box_data_cleaned <- box_data_cleaned[box_data_cleaned$pathologicN != 'Nx', ]

# 3. 对 pathologicStage 进行规范化
box_data_cleaned$pathologicStage <- ifelse(
  box_data_cleaned$pathologicStage %in% c('I', 'Stage IA', 'Stage IB'), 1,
  ifelse(
    box_data_cleaned$pathologicStage %in% c('II', 'Stage II', 'Stage IIA', 'Stage IIB'), 2,
    ifelse(
      box_data_cleaned$pathologicStage %in% c('III', 'Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC'), 3,
      4
    )
  )
)

library(ggplot2)
library(dplyr)

#####################################################################################################
# gender
# 计算每个 Cluster 中 male 和 female 的占比
gender_summary <- box_data_cleaned %>%
  group_by(Cluster, gender) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# 绘制柱状图
p <- ggplot(gender_summary, aes(x = factor(Cluster), y = percentage, fill = gender)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.5) +  # 添加边框和调节粗细
  scale_fill_manual(values = c("male" = "#3636b3", "female" = "#b23535")) +
  labs(x = "Cluster", y = "Percentage (%)", fill = "Gender") +  # 设置坐标轴标题
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # 设置坐标轴线的颜色和粗细
    axis.ticks = element_line(color = "black", size = 1),  # 设置刻度线的颜色和粗细
    axis.title.x = element_text(size = 18, color = "black"),  # 设置横坐标轴标题的字体颜色和大小
    axis.title.y = element_text(size = 18, color = "black"),  # 设置纵坐标轴标题的字体颜色和大小
    axis.text.x = element_text(size = 16, color = "black"),    # 设置横坐标轴刻度标签的字体颜色和大小
    axis.text.y = element_text(size = 16, color = "black"),    # 设置纵坐标轴刻度标签的字体颜色和大小
    panel.grid.major = element_blank(),           # 移除主网格线
    panel.grid.minor = element_blank()            # 移除次网格线
  ) +
  ggtitle("Gender Distribution by Cluster")

# 保存图形为 PDF
ggsave("results/Consensus/gender_distribution_by_cluster.pdf", plot = p, width = 5, height = 6)

#####################################################################################################
#pathologicT
# 计算每个 Cluster 中 pathologicT 类别的占比
t_cluster_summary <- box_data_cleaned %>%
  group_by(Cluster, pathologicT) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# 绘制柱状图
p_t_cluster <- ggplot(t_cluster_summary, aes(x = factor(Cluster), y = percentage, fill = pathologicT)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.5) +  # 添加边框和调节粗细
  scale_fill_manual(values = c("T1" = "#3636b3", "T2" = "#36b3b3", "T3" = "#73b235", "T4" = "#b23535")) +  # 设置颜色
  labs(x = "Cluster", y = "Percentage (%)", fill = "T Stage") +  # 设置坐标轴标题
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # 设置坐标轴线的颜色和粗细
    axis.ticks = element_line(color = "black", size = 1),  # 设置刻度线的颜色和粗细
    axis.title.x = element_text(size = 18, color = "black"),  # 设置横坐标轴标题的字体颜色和大小
    axis.title.y = element_text(size = 18, color = "black"),  # 设置纵坐标轴标题的字体颜色和大小
    axis.text.x = element_text(size = 16, color = "black"),    # 设置横坐标轴刻度标签的字体颜色和大小
    axis.text.y = element_text(size = 16, color = "black"),    # 设置纵坐标轴刻度标签的字体颜色和大小
    panel.grid.major = element_blank(),         # 移除主网格线
    panel.grid.minor = element_blank()          # 移除次网格线
  ) +
  ggtitle("Pathologic T Stage Distribution by Cluster")

# 保存图形为 PDF
ggsave("results/Consensus/pathologic_t_distribution_by_cluster.pdf", plot = p_t_cluster, width = 5, height = 6)

#####################################################################################################
#pathologicN
# 计算每个 Cluster 中 pathologicN 类别的占比
n_cluster_summary <- box_data_cleaned %>%
  group_by(Cluster, pathologicN) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# 绘制柱状图
p_n_cluster <- ggplot(n_cluster_summary, aes(x = factor(Cluster), y = percentage, fill = pathologicN)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.5) +  # 添加边框和调节粗细
  scale_fill_manual(values = c("N0" = "#3636b3", "N1" = "#36b3b3", "N2" = "#73b235", "N3" = "#b23535")) +  # 设置颜色
  labs(x = "Cluster", y = "Percentage (%)", fill = "N Stage") +  # 设置坐标轴标题
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # 设置坐标轴线的颜色和粗细
    axis.ticks = element_line(color = "black", size = 1),  # 设置刻度线的颜色和粗细
    axis.title.x = element_text(size = 18, color = "black"),  # 设置横坐标轴标题的字体颜色和大小
    axis.title.y = element_text(size = 18, color = "black"),  # 设置纵坐标轴标题的字体颜色和大小
    axis.text.x = element_text(size = 16, color = "black"),    # 设置横坐标轴刻度标签的字体颜色和大小
    axis.text.y = element_text(size = 16, color = "black"),    # 设置纵坐标轴刻度标签的字体颜色和大小
    panel.grid.major = element_blank(),         # 移除主网格线
    panel.grid.minor = element_blank()          # 移除次网格线
  ) +
  ggtitle("Pathologic N Stage Distribution by Cluster")

# 保存图形为 PDF
ggsave("results/Consensus/pathologic_n_distribution_by_cluster.pdf", plot = p_n_cluster, width = 5, height = 6)


#####################################################################################################
#pathologicStage
# 假设 box_data_cleaned 是你的数据框

# 计算每个 Cluster 中 pathologicStage 类别的占比
stage_cluster_summary <- box_data_cleaned %>%
  group_by(Cluster, pathologicStage) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# 绘制柱状图
p_stage_cluster <- ggplot(stage_cluster_summary, aes(x = factor(Cluster), y = percentage, fill = factor(pathologicStage))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.5) +  # 添加边框和调节粗细
  scale_fill_manual(values = c("1" = "#3636b3", "2" = "#36b3b3", "3" = "#73b235", "4" = "#b23535")) +  # 设置颜色
  labs(x = "Cluster", y = "Percentage (%)", fill = "TNM Stage") +  # 设置坐标轴标题
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # 设置坐标轴线的颜色和粗细
    axis.ticks = element_line(color = "black", size = 1),  # 设置刻度线的颜色和粗细
    axis.title.x = element_text(size = 18, color = "black"),  # 设置横坐标轴标题的字体颜色和大小
    axis.title.y = element_text(size = 18, color = "black"),  # 设置纵坐标轴标题的字体颜色和大小
    axis.text.x = element_text(size = 16, color = "black"),    # 设置横坐标轴刻度标签的字体颜色和大小
    axis.text.y = element_text(size = 16, color = "black"),    # 设置纵坐标轴刻度标签的字体颜色和大小
    panel.grid.major = element_blank(),         # 移除主网格线
    panel.grid.minor = element_blank()          # 移除次网格线
  ) +
  ggtitle("Pathologic Stage Distribution by Cluster")

# 保存图形为 PDF
ggsave("results/Consensus/pathologic_stage_distribution_by_cluster.pdf", plot = p_stage_cluster, width = 5, height = 6)

#####################################################################################################
#Tumor.grade
# 假设 box_data_cleaned 是你的数据框
# 计算每个 Cluster 中 Tumor.grade 类别的占比
stage_cluster_summary <- box_data_cleaned %>%
  group_by(Cluster, Tumor.grade) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# 绘制柱状图
p_stage_cluster <- ggplot(stage_cluster_summary, aes(x = factor(Cluster), y = percentage, fill = factor(Tumor.grade, levels = c("poorly", "moderately", "well")))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.5) +  # 添加边框和调节粗细
  scale_fill_manual(values = c("poorly" = "#b23535", "moderately" = "#73b235", "well" = "#3636b3")) +  # 设置颜色
  labs(x = "Cluster", y = "Percentage (%)", fill = "Tumor Grade") +  # 设置坐标轴标题
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # 设置坐标轴线的颜色和粗细
    axis.ticks = element_line(color = "black", size = 1),  # 设置刻度线的颜色和粗细
    axis.title.x = element_text(size = 18, color = "black"),  # 设置横坐标轴标题的字体颜色和大小
    axis.title.y = element_text(size = 18, color = "black"),  # 设置纵坐标轴标题的字体颜色和大小
    axis.text.x = element_text(size = 16, color = "black"),    # 设置横坐标轴刻度标签的字体颜色和大小
    axis.text.y = element_text(size = 16, color = "black"),    # 设置纵坐标轴刻度标签的字体颜色和大小
    panel.grid.major = element_blank(),         # 移除主网格线
    panel.grid.minor = element_blank()          # 移除次网格线
  ) +
  ggtitle("Tumor Grade Distribution by Cluster")


# 保存图形为 PDF
ggsave("results/Consensus/Tumorgrade_distribution_by_cluster.pdf", plot = p_stage_cluster, width = 5, height = 6)


#####################################################################################################
#age
# 假设 box_data_cleaned 是你的数据框
# 筛选出 Cluster 为 1 和 2 的数据
age_data <- box_data_cleaned %>%
  filter(Cluster %in% c(1, 2))

# 绘制箱线图
p_age <- ggplot(age_data, aes(x = factor(Cluster), y = age, fill = factor(Cluster))) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +  # 绘制箱线图，并设置异常值的颜色和形状
  labs(x = "Cluster", y = "Age", fill = "Cluster") +  # 设置坐标轴标题
  scale_fill_manual(values = c("1" = "#3636b3", "2" = "#b23535")) +  # 设置填充颜色
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # 设置坐标轴线的颜色和粗细
    axis.ticks = element_line(color = "black", size = 1),  # 设置刻度线的颜色和粗细
    axis.title.x = element_text(size = 18, color = "black"),  # 设置横坐标轴标题的字体颜色和大小
    axis.title.y = element_text(size = 18, color = "black"),  # 设置纵坐标轴标题的字体颜色和大小
    axis.text.x = element_text(size = 16, color = "black"),    # 设置横坐标轴刻度标签的字体颜色和大小
    axis.text.y = element_text(size = 16, color = "black"),    # 设置纵坐标轴刻度标签的字体颜色和大小
    panel.grid.major = element_blank(),         # 移除主网格线
    panel.grid.minor = element_blank()          # 移除次网格线
  ) +
  ggtitle("Age Distribution by Cluster")

# 保存图形为 PDF
ggsave("results/Consensus/age_distribution_by_cluster.pdf", plot = p_age, width = 4, height = 6)


###cibersort

