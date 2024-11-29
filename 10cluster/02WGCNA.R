library(WGCNA)
library(stringr)

options(stringsAsFactors = FALSE) 
enableWGCNAThreads() ## 打开多线程

setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results/WGCNA')


exprMat <- "GSE53625_tumor_count_loged.txt"

# 官方推荐 "signed" 或 "signed hybrid"
type = "unsigned"

# 相关性计算
# 官方推荐 biweight mid-correlation & bicor
# corType: pearson or bicor
# 为与原文档一致，故未修改
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

##导入数据
dataExpr <- read.table(paste0(routine_data, 'results/GEO_tumor_exp_cleaned.txt'), header = TRUE, row.names = 1, check.names=FALSE)
# dataExpr <- read.table(paste0(routine_data_GEO, 'GSE53625_tumor_count_loged.txt'), header = TRUE, row.names = 1, check.names=FALSE)


## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                0),]
# max(quantile(m.mad, probs=seq(0, 1, 0.75))[2],0.01)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

# ## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

##2离群，delete
dataExpr <- dataExpr[!rownames(dataExpr) %in% c("ec112"), ]
dataExpr <- as.data.frame(dataExpr)

## 软阈值的筛选原则是使构建的网络更符合无标度网络特征
powers = c(c(1:20), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType = type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")

power = sft$powerEstimate
power

##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.4,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)

# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。
table(net$colors)

### 层级聚类树展示各个模块
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

### 模块之间相关性图
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)

# ##### 可视化基因网络 (TOM plot)
# # 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# # 否则需要再计算一遍，比较耗费时间
# # TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
# load(net$TOMFiles[1], verbose=T)
# 
# TOM <- as.matrix(TOM)
# 
# dissTOM = 1-TOM
# # Transform dissTOM with a power to make moderately strong
# # connections more visible in the heatmap
# plotTOM = dissTOM^7
# # Set diagonal to NA for a nicer plot
# diag(plotTOM) = NA
# # Call the plot function
# 
# # 这一部分特别耗时，行列同时做层级聚类
# TOMplot(plotTOM, net$dendrograms, moduleColors,
#     main = "Network heatmap plot, all genes")

#### 关联表型数据
trait <- "results/WGCNA/group_for_WGCNA.csv"
exp_tumor <- read.csv(trait, row.names = 1)
# 2. 修改列名
# 2.1 对以 'TCGA' 开头的列，保留第15位字符为 '1' 的列，修改列名
tcga_columns <- rownames(exp_tumor)[grep("^TCGA", rownames(exp_tumor))]
tcga_columns_filtered <- tcga_columns[substr(tcga_columns, 15, 15) == "1"]
new_tcga_columns <- substr(tcga_columns_filtered, 1, 12)
new_tcga_columns <- gsub("\\.", "-", new_tcga_columns)
# 2.2 对以 'ec' 开头的列，保留 '_' 之前的内容
ec_columns <- rownames(exp_tumor)[grep("^ec", rownames(exp_tumor))]
new_ec_columns <- sub("_.*$", "", ec_columns)
# 合并筛选后的列名，并修改列名
filtered_columns <- c(tcga_columns_filtered, ec_columns)
filtered_names <- c(new_tcga_columns, new_ec_columns)
names(filtered_names) <- filtered_columns
# 生成新的 exp 矩阵
exp_filtered <- exp_tumor[filtered_columns, ]
exp_filtered <- as.data.frame(exp_filtered)
rownames(exp_filtered) <- filtered_names


# 读入表型数据，不是必须的
# if(trait != "") {
# traitData <- read.csv(file=trait, row.names=1)
traitData <- exp_filtered
sampleName = rownames(dataExpr)
traitData = traitData[match(sampleName, rownames(traitData)), ]
# }

# traitData <- as.factor(traitData)

### 模块与表型数据关联
# if (corType=="pearsoon") {
modTraitCor = cor(MEs_col, traitData, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
# } else {
#   modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
#   modTraitCor = modTraitCorP$bicor
#   modTraitP   = modTraitCorP$p
# }
## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 3), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData),
               yLabels = colnames(MEs_col),
               cex.lab = 1,
               ySymbols = colnames(MEs_col), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 1, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# names(dataExpr)[moduleColors=="pink"]
names(dataExpr)[moduleColors=="green"]
write.csv(names(dataExpr)[moduleColors=="green"], 'results/WGCNA/green_module_gene.csv')
# #差异基因
# dif_gene <- read.csv(paste0(routine_sfig1, "cluster.DESeq2.select.csv"), row.names = 1)
# # final_genes <- intersect(rownames(dif_gene), names(dataExpr)[moduleColors=="pink"])
# final_genes <- intersect(rownames(dif_gene), names(dataExpr)[moduleColors=="pink"])

###genes
weight =as.data.frame(traitData$Cluster1)
names(weight) ="weight"

modNames = substring(names(MEs_col), 3)

# 计算MM的P值
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# 计算性状和基因表达量之间的相关性（GS）
geneTraitSignificance = as.data.frame(cor(dataExpr, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

# module = "pink"
module = 'green'
column = match(module, modNames)
moduleGenes = moduleColors==module
pdf("results/WGCNA/Module membership vs gene significance_Green.pdf",width = 7, height=7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, 
                                            column]),abs(geneTraitSignificance[moduleGenes, 1]), xlab = 
                     paste("Module Membership in", module, "module"), ylab = "Gene 
  significance for body weight", main = paste("Module membership 
  vs gene significance"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 
                     1.2, col = 'black')
dev.off()

##
GS_MM <- cbind(Gene = names(dataExpr)[moduleColors=="black"], MM = geneModuleMembership[moduleGenes, column], GS = abs(geneTraitSignificance[moduleGenes, 1]))
GS_MM <- as.data.frame(GS_MM)
# 确保 MM 和 GS 列为数值类型
GS_MM$MM <- as.numeric(GS_MM$MM)
GS_MM$GS <- as.numeric(GS_MM$GS)
filtered_GSMM <- GS_MM[abs(GS_MM[, "MM"]) > 0.8 & abs(GS_MM[, "GS"]) > 0.2, ]
write.csv(filtered_GSMM, "results/WGCNA/GS_MM_keygenes_Green.csv")

# ####gene signifaction
# names(dataExpr)
# # 返回brown模块中所有ID
# # names(dataExpr)[moduleColors=="pink"]
# names(dataExpr)[moduleColors=="pink"]
# # # 导入注释文件
# # annot = read.csv(file = "GeneAnnotation.csv");
# # dim(annot)
# # names(annot)
# probes = names(dataExpr)
# # probes2annot = match(probes, annot$substanceBXH)
# # # 统计没有注释到的基因
# # sum(is.na(probes2annot))
# # 创建数据集，包含探测ID ，
# geneInfo0 = data.frame(substanceBXH = probes, 
#                        # geneSymbol = annot$gene_symbol[probes2annot], 
#                        # LocusLinkID = annot$LocusLinkID[probes2annot], 
#                        moduleColor = moduleColors, 
#                        geneTraitSignificance, GSPvalue)
# # 通过显著性对模块进行排序
# modOrder = order(-abs(cor(MEs, weight, use = "p")))
# # 添加模块成员
# for (mod in 1:ncol(geneModuleMembership))
# {
#   oldNames = names(geneInfo0)
#   geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, 
#                                                          modOrder[mod]],
#                          MMPvalue[, modOrder[mod]]);
#   names(geneInfo0) = c(oldNames, paste("MM.", 
#                                        modNames[modOrder[mod]], sep=""),
#                        paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# # 对基因进行排序
# geneOrder = order(geneInfo0$moduleColor, abs(geneInfo0$GS.weight))
# geneInfo = geneInfo0[geneOrder, ]
# # 导出
# write.csv(geneInfo, file = "geneInfo.csv")
# 
# ###GSMM基因 merge DEGs
# gsmm <- read.csv('GS_MM_keygenes.csv', row.names = 1)
# deg <- read.csv('cluster.DESeq2.select.csv', row.names = 1)
# 
# gsmm <- rownames(gsmm)
# deg <- rownames(deg)
# finalgene <- intersect(gsmm, deg)
# 
# write.csv(finalgene, 'finalkeygene.csv')
