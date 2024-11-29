setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'

dir.create('results')

#???ð?
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)

pFilter=0.001                  #pvalue?Ĺ???????
expFile= paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt')          #?????????ļ?
ClusterFile='results/cluster_clinic_data.csv'      #???͵Ľ????ļ?

allDrugs=c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")

#??ȡ?????????ļ?,???????ݽ??д???
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
# attributes(rt) <- attributes(rt)[-which(names(attributes(rt)) == "dim")]
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
dimnames=list(rownames(rt),colnames(rt))
data=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

# #ɾ????????Ʒ
# group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2","1",group)
# data=data[,group==0]
# data=t(data)
# rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
# data=avereps(data)
# data=t(data)

#??ȡ???ͽ????ļ?
ClusterRT=read.csv(ClusterFile, row.names=1)

colnames(ClusterRT)[colnames(ClusterRT) == "Risk"] <- "Cluster"

ClusterRT$Cluster <- as.factor(ClusterRT$Cluster)
# data <- as.matrix(data)

# 创建一个空的列表来存储所有的图形
figures <- list()

for(drug in allDrugs){
  # 预测药物敏感性
  senstivity = pRRopheticPredict(data, drug, selection = 1)
  senstivity = senstivity[senstivity != "NaN"]
  
  # 交集提取样本数据
  sameSample = intersect(row.names(ClusterRT), names(senstivity))
  Cluster = ClusterRT[sameSample, "Cluster", drop = F]
  senstivity = senstivity[sameSample]
  rt = cbind(Cluster, senstivity)
  
  # 设置比较类型
  type = levels(factor(rt[,"Cluster"]))
  comp = combn(type, 2)
  my_comparisons = list()
  for(i in 1:ncol(comp)){
    my_comparisons[[i]] <- comp[,i]
  }
  
  # 获取p值
  if(length(levels(factor(rt[,"Cluster"]))) > 2){
    test = kruskal.test(senstivity ~ Cluster, data = rt)
  } else {
    test = wilcox.test(senstivity ~ Cluster, data = rt)
  }
  pvalue = test$p.value
  if(pvalue < pFilter){
    # 绘制箱型图并调整显著性标注位置
    boxplot = ggplot(data = rt, aes(x = Cluster, y = senstivity, fill = Cluster)) + 
      scale_fill_manual(values = c("#ed0000", "#00468b")) +
      geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                  linewidth = 0.8, color = "black") +
      geom_boxplot(notch = TRUE, outlier.size = -1, 
                   color = "black", linewidth = 0.8, alpha = 0.7) +
      geom_point(shape = 21, size = 2, 
                 position = position_jitterdodge(), 
                 color = "black", alpha = 1) +
      theme_bw() + 
      ylab(paste0(drug, " senstivity (IC50)")) +
      xlab('Cluster') +
      theme(axis.text.x = element_text(size = 12, color = "black"),
            axis.ticks = element_line(linewidth = 0.2, color = "black"),
            axis.ticks.length = unit(0.2, "cm"),
            legend.position = "none",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 12)) + 
      stat_compare_means(comparisons = my_comparisons)  # 调整显著性标注位置
    
    # # 使用expand_limits扩展y轴范围
    # boxplot = boxplot + expand_limits(y = max(rt$senstivity) + 0.2)
    # 
    # 将每个药物的图形存入列表
    fig_name = paste0("fig", length(figures) + 1)
    figures[[fig_name]] = boxplot
  }
}


# 将所有图形放在一个大图上，4行7列
final_plot = ggarrange(plotlist = figures, ncol = 7, nrow = 4)

# 保存最终的大图
pdf(file = "results/drug/AllDrugSensitivity.pdf", width = 24, height = 14)
print(final_plot)
dev.off()

