setwd('D:/ZJU-FISH/doctor/TRM/07TIME')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results')

library(tidyr)
library(ggpubr)
library(tibble)


# 设置文件路径
input_folder <- "results/gene/"  # CSV 文件所在文件夹
output_folder <- "results/gene/"  # 输出 PDF 文件夹

# 检查输出文件夹是否存在，不存在则创建
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# 遍历每个 CSV 文件
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
for (file_path in csv_files) {
  # 读取数据
  data <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  # 提取感兴趣的列
  result <- data[, c("T.cells.gamma.delta", "NK.cells.resting", 
                     "NK.cells.activated", "Monocytes", 
                     "Macrophages.M2", "Neutrophils", "group"), drop = FALSE]
  
  # 将行名转换为 sample 列
  result <- result %>% rownames_to_column("sample")
  
  # 重整数据格式
  b <- gather(result, key = "CIBERSORT", value = "Proportion", -c(group, sample))
  
  # 绘制箱线图
  p <- ggboxplot(b, x = "CIBERSORT", y = "Proportion",
                 fill = NULL,
                 color = "group",
                 palette = c("Low" = "#34499d", "High" = "#e92428"),
                 size = 0.6,
                 width = 0.5) +
    stat_compare_means(aes(group = group),
                       method = "wilcox.test",
                       label = "p.signif",
                       size = 8,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", ""))) +
    geom_boxplot(aes(color = group), 
                 position = position_dodge(0.8), 
                 size = 0.6,
                 width = 0.5, 
                 outlier.colour = NA) +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 20))
  
  # 获取文件名并保存 PDF
  output_filename <- paste0(output_folder, basename(file_path), ".pdf")
  ggsave(filename = output_filename, plot = p, device = "pdf", width = 7, height = 6)
  
  print(paste("Saved:", output_filename))
}