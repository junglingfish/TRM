setwd('D:/ZJU-FISH/doctor/TRM/05nomogram')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'

library(dplyr)
library(tidyr)


clinic <- read.csv(paste0(routine_data, 'clinic/GEO_clinic_data.csv'))
riskscore <- read.delim(paste0(routine_4, 'results/gse53625.risk.txt'),sep='\t',header = T)

clinic_updated<-merge(riskscore[3:15],clinic,by='Samples')

box_data <- clinic_updated

# 假设数据存储在 box_data 中，将其转换为数据框方便数据处理
box_data <- as.data.frame(box_data)

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

# # 4. 删除 pathologicStage 列中为 NA 的行（对应 Stage IV）
# box_data_cleaned <- box_data_cleaned[!is.na(box_data_cleaned$pathologicStage), ]

box_data <- box_data_cleaned

write.csv(box_data, 'results/GEO_clinic_riskscore.csv')

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
         
# 根据年龄划分组
box_data$age_group <- ifelse(box_data$age < 60, "<60", "≥60")

# 根据性别划分组
# 这里可以直接使用box_data中的gender列

# 根据pathologicT分期进行分组
box_data$pathologicT_group <- ifelse(box_data$pathologicT %in% c("T1", "T2"), "T1&T2", "T3&T4")

# 根据pathologicN分期进行分组
box_data$pathologicN_group <- ifelse(box_data$pathologicN == "N0", "N0", "N1,2,3")

# 修改pathologicStage分组
box_data$pathologicStage <- as.character(box_data$pathologicStage)
box_data$pathologicStage[box_data$pathologicStage == "1"] <- "Stage I"
box_data$pathologicStage[box_data$pathologicStage == "2"] <- "Stage II"
box_data$pathologicStage[box_data$pathologicStage %in% c("3", "4")] <- "Stage III and IV"

box_data <- as.data.frame(box_data)

# 汇总数据
box_data_long <- box_data %>%
  dplyr::select(riskscore, age_group, gender, pathologicStage, pathologicT_group, pathologicN_group, Tumor.grade) %>%
  tidyr::pivot_longer(cols = c(age_group, gender, pathologicStage, pathologicT_group, pathologicN_group, Tumor.grade), 
                      names_to = "Group", 
                      values_to = "GroupValue")


# 设置 GroupValue 列的因子水平顺序
box_data_long$GroupValue <- factor(box_data_long$GroupValue, 
                                   levels = c("<60", "≥60", "male", "female", 
                                              "T1&T2", "T3&T4", 
                                              "N0", "N1,2,3", 
                                              "Stage I", "Stage II", "Stage III and IV", 
                                              "well", "moderately", "poorly"))

# 绘制箱线图
p <- ggplot(box_data_long, aes(x = GroupValue, y = riskscore, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               color = "black", size = 1) +  # 添加黑色加粗边框
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 0, alpha = 0) +  # 删除散点
  stat_boxplot(geom = "errorbar",
               aes(ymin = after_stat(ymax)), 
               position = position_dodge(0.8),  # 确保误差线对齐
               width = 0.4, size = 0.6, color = "black") +  # 设置为黑色
  stat_boxplot(geom = "errorbar",
               aes(ymax = after_stat(ymin)), 
               position = position_dodge(0.8),  # 确保误差线对齐
               width = 0.4, size = 0.6, color = "black") +  # 设置为黑色
  scale_fill_manual(values = c("age_group" = "#27447C", 
                               "gender" = "#4871B3", 
                               "pathologicT_group" = "#E73C36",
                               "pathologicN_group" = "#991F22", 
                               "Tumor.grade" = "#B88640", 
                               "pathologicStage" = "#168676")) +
  labs(x = '', y = 'Risk Score', title = '') +  # 隐藏 x 轴和图标题
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # 隐藏 x 轴标题
    plot.title = element_blank(),  # 隐藏整张图标题
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, color = "black"),  # 加粗字体
    axis.text.y = element_text(size = 20, color = "black"),  # 加粗字体
    axis.title.y = element_text(size = 30),  # 加粗字体
    axis.ticks = element_line(color = "black", size = 0.7),  # 添加黑色刻度线
    axis.ticks.length = unit(0.2, "cm"),  # 设置刻度线长度
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),   # 移除次要网格线
    panel.border = element_blank()  # 移除所有边框
  ) +
  # 手动添加左侧和下侧边框
  geom_hline(yintercept = -Inf, color = "black", size = 1.5) +  # 下边框
  geom_vline(xintercept = -Inf, color = "black", size = 1.5)  # 左边框

# 添加显著性标注，设置 y_position 以确保标注可见
p <- p + geom_signif(comparisons = list(c("<60", "≥60"),
                                        c("male", "female"),
                                        c("T1&T2", "T3&T4"),
                                        c("N0", "N1,2,3"), 
                                        c("Stage I", "Stage II"),
                                        c("Stage I", "Stage III and IV"),
                                        c("Stage II", "Stage III and IV"),
                                        c("well", "moderately"),
                                        c("moderately", "poorly"),
                                        c("well", "poorly")),
                     map_signif_level = function(p) {
                       if (p <= 0.001) {
                         return("***")
                       } else if (p <= 0.01) {
                         return("**")
                       } else if (p <= 0.15) {
                         return("*")
                       } else {
                         return("NS.")
                       }
                     },
                     y_position = c(2.7, 2.7, 2.7, 2.7, 2.7, 3, 2.7, 2.7, 2.7, 3), textsize = 7, size = 0.8)  # 根据数据范围调整这些值

# 显示箱线图
print(p)

# 保存箱线图为 PDF 文件
ggsave('results/GEO_factor_boxplot.pdf', plot = p, width = 10, height = 8, dpi = 300)
