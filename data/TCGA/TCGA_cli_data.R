library(readr)
library(tidyverse)

setwd('D:/ZJU-FISH/doctor/model/data/TCGA/')
routine = 'D:/ZJU-FISH/doctor/model/data/TCGA/'

##########################################################XML
##XML
library(TCGAbiolinks)
# 下载XML临床数据
query1 <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "bcr xml"
)
GDCdownload(query1)

# 解析patient信息,指定clinical.info
clinical.patient.xml <- GDCprepare_clinic(query1, clinical.info = "patient")
dim(clinical.patient.xml) 
#提取样品ID
bcr_patient_barcode=clinical.patient.xml$bcr_patient_barcode
#提取年龄
days_to_birth=clinical.patient.xml$days_to_birth
#提取出生年份
days_to_last_followup=clinical.patient.xml$days_to_last_followup
#提取性别
days_to_death=clinical.patient.xml$days_to_death
#提取诊断时间
day_of_form_completion=clinical.patient.xml$day_of_form_completion
#提取生存时间
month_of_form_completion=clinical.patient.xml$month_of_form_completion
#提取生存状态
year_of_form_completion=clinical.patient.xml$year_of_form_completion
xml_clinical=cbind(bcr_patient_barcode,
                   days_to_birth,
                   days_to_last_followup,
                   days_to_death,
                   day_of_form_completion,
                   month_of_form_completion,
                   year_of_form_completion)
# 将矩阵导出为 .txt 文件
# write.table(xml_clinical, file = "clinical_patient.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.csv(xml_clinical,file = "clinical_xml_clinical.csv",quote = F)
# clinical.followup.xml <- GDCprepare_clinic(query1, clinical.info = "follow_up")

##################################################################################################
###整合clinical矩阵
clinical=read_tsv("clinical.cart.2024-09-13/clinical.tsv")
#提取样品ID
ID=clinical$case_submitter_id
#提取年龄
age=clinical$age_at_index
#提取出生年份
year=clinical$year_of_birth
#提取性别
gender=clinical$gender
#提取诊断时间
diagnosis=clinical$days_to_diagnosis
#提取生存时间
time=clinical$days_to_death
#提取生存状态
status=clinical$vital_status
#提取TMN分期
pathologicT=clinical$ajcc_pathologic_t
pathologicM=clinical$ajcc_pathologic_m
pathologicN=clinical$ajcc_pathologic_n
#提取stage分期
pathologicStage=clinical$ajcc_pathologic_stage

TCGA_clinical=cbind(ID,
                    age,
                    year,
                    gender,
                    diagnosis,
                    time,
                    status,
                    pathologicT,
                    pathologicM,
                    pathologicN,
                    pathologicStage)

TCGA_clinical=as.data.frame(TCGA_clinical)
#删除重复ID
duplicated(TCGA_clinical$ID)
TCGA_final<-TCGA_clinical[!duplicated(TCGA_clinical$ID),]
#导出文件
rownames(TCGA_final)=TCGA_final$ID
TCGA_final=TCGA_final[,2:ncol(TCGA_final)]
write.csv(TCGA_final,file = "TCGA_clinical.csv",quote = F)


#################################################################################################
##整合clinical数据
clinical=read.csv("TCGA_clinical.csv")
# 读取导出的 .txt 文件
clinical_matrix <- read.csv("clinical_xml_clinical.csv")
# 找到重合的内容
overlap <- intersect(clinical_matrix[, "bcr_patient_barcode"], clinical[, "case_submitter_id"])
# 在 clinical_matrix 中保留仅包含重合内容的行
clinical_matrix_filtered <- clinical_matrix[clinical_matrix[, "bcr_patient_barcode"] %in% overlap, ]
# merge
# 使用 match 获取顺序索引
sorted_clinical_matrix <- clinical_matrix[match(clinical$case_submitter_id, clinical_matrix$bcr_patient_barcode), ]
# 合并 days_to_last_followup 和 days_to_death，互补填充 NA
sorted_clinical_matrix$time <- ifelse(
  is.na(sorted_clinical_matrix$days_to_last_followup), 
  sorted_clinical_matrix$days_to_death, 
  sorted_clinical_matrix$days_to_last_followup
)
# 确保 'time' 列在两个矩阵中存在
if ("time" %in% colnames(sorted_clinical_matrix) & "time" %in% colnames(clinical)) {
  # 替换 clinical_matrix 中的 time 列
  clinical$time <- sorted_clinical_matrix$time
} else {
  stop("确保两个矩阵中都包含 'time' 列")
}

sorted_clinical_matrix_selected <- sorted_clinical_matrix[, (ncol(sorted_clinical_matrix)-3):(ncol(sorted_clinical_matrix)-1)]
# 合并选定的列
# 假设我们想要按列合并这两个子矩阵
merged_matrix <- cbind(clinical, sorted_clinical_matrix_selected)
# 将可能的非数值数据转换为数值型
merged_matrix$year <- as.numeric(as.character(merged_matrix$year))
merged_matrix$age <- as.numeric(as.character(merged_matrix$age))
merged_matrix$year_of_form_completion <- as.numeric(as.character(merged_matrix$year_of_form_completion))
merged_matrix$month_of_form_completion <- as.numeric(as.character(merged_matrix$month_of_form_completion))
merged_matrix$day_of_form_completion <- as.numeric(as.character(merged_matrix$day_of_form_completion))
origin_matrix <- merged_matrix
# 计算起始日期和终止日期，并填充空缺的 time 列
# 首先，处理 time 列中的 NA 值
merged_matrix$time <- with(merged_matrix, {
  # 生成起始日期：年 + 年龄的日期
  start_date <- as.Date(paste(year + age, "01-01", sep = "-"))
  # 生成终止日期：基于 year_of_form_completion、month_of_form_completion 和 day_of_form_completion
  end_date <- as.Date(paste(year_of_form_completion, month_of_form_completion, day_of_form_completion, sep = "-"))
  # 计算天数差
  time_diff <- as.numeric(difftime(end_date, start_date, units = "days"))
  # 用计算得到的天数填充 NA 值
  ifelse(is.na(time), time_diff, time)
})

write.csv(merged_matrix,file = "TCGA_clinical_95_final.csv",quote = F)
