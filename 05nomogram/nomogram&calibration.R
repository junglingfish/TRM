setwd('D:/ZJU-FISH/doctor/TRM/05nomogram')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'

library(survival)
library(rms)
library(regplot)
library(ggplot2)

data <- read.csv('results/GEO_clinic_data_cleaned.csv', row.names = 1)

# 创建 datadist 对象并指定数据分布
dd <- datadist(data)
options(datadist = "dd")

# data$gender = as.factor(data$gender)
# data$pathologicT = as.factor(data$pathologicT)
# data$pathologicN = as.factor(data$pathologicN)
data$pathologicStage = as.factor(data$pathologicStage)
# data$Tumor.grade = as.factor((data$Tumor.grade))

coxfit <- cph(Surv(time, status) ~ riskscore + pathologicStage + age,
              data = data, x=T,y=T,surv = T
)

p1 <- regplot(coxfit#对观测2的六个指标在列线图上进行计分展示
        ,observation = data[10,] #也可以不展示
        ,obscol = "#e92428"
        #预测3年和5年的死亡风险，此处单位是day
        ,title="Nomogram for ESCC"
        ,failtime = c(2, 3, 5)
        ,prfail = TRUE #cox回归中需要TRUE
        ,showP = T #是否展示统计学差异
        ,droplines = F#观测2示例计分是否画线
        #,colors = mg_colors[1:3] #用前面自己定义的颜色
        #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
        #,interval="confidence"
        #,rank="decreasing"
        #,clickable=T
        ,points=TRUE)

ggsave('results/ESCC_nomogram.pdf',p1$plot,height = 8,width = 14)


#校准曲线
library(rms)

# 确保所有生存时间大于0
data <- data[data$time > 0, ]

# 调整模型，确保x=T, y=T
f2 <- psm(Surv(time, status) ~ riskscore + pathologicStage + pathologicN + pathologicT + gender + age, data = data, x = TRUE, y = TRUE, dist = 'lognormal')

# 构建校准曲线
cal2 <- calibrate(f2, cmethod='KM', method="boot", u=1, m=25, B=1000)   # 2年
cal3 <- calibrate(f2, cmethod='KM', method="boot", u=3, m=25, B=1000)   # 3年
cal5 <- calibrate(f2, cmethod='KM', method="boot", u=5, m=25, B=1000)   # 5年

# 打开PDF设备
pdf("results/calibration_compare.pdf", width = 8, height = 8)

# 设置边距并绘制黑色边框
par(mar = c(5, 5, 4, 2) + 0.1)  # 调整边距大小，增加空间

# 绘制2年校准曲线
plot(cal2, lwd = 2, lty = 0, errbar.col = c("#2166AC"),
     bty = "l", # 只画左边和下边框
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = expression(bold("Nomogram-predicted OS (%)")), 
     ylab = expression(bold("Observed OS (%)")),
     col = c("#2166AC"),
     cex.lab = 1.6, cex.axis = 1.6, cex.main = 1.2, cex.sub = 0.6)
lines(cal2[, c('mean.predicted', "KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)

# 绘制3年校准曲线
plot(cal3, lwd = 2, lty = 0, errbar.col = c("#B2182B"),
     xlim = c(0, 1), ylim = c(0, 1), col = c("#B2182B"), add = TRUE)
lines(cal3[, c('mean.predicted', "KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

# 绘制5年校准曲线
plot(cal5, lwd = 2, lty = 0, errbar.col = c("#e18d2d"),
     xlim = c(0, 1), ylim = c(0, 1), col = c("#e18d2d"), add = TRUE)
lines(cal5[, c('mean.predicted', "KM")],
      type = 'b', lwd = 1, col = c("#e18d2d"), pch = 16)

# 添加 y=x 参考线
abline(0, 1, lwd = 2, lty = 3, col = c("#224444"))

# 添加图例
legend("bottomright", # 图例的位置
       legend = c("1-year", "3-year", "5-year"), # 图例文字
       col = c("#2166AC", "#B2182B", "#e18d2d"), # 图例线的颜色
       lwd = 2, # 图例中线的粗细
       cex = 1.6, # 图例字体大小
       bty = "n") # 不显示图例边框

box(col = "black", lwd = 2)      # 添加黑色边框，线宽设为2
# 可选：添加图标题
# title(main = "Calibration Curves for OS at Different Time Points", cex.main = 1.5)

# 关闭PDF设备
dev.off()

#决策曲线
library(rmda)
data$status <- ifelse(data$status == "0", 0, 1) #保证你的输入data全为数值型，这里将Disease的状态定义为0，1。CT为0，disease为1
####简单模型---单个自变量
simple_penk<-decision_curve(status ~ riskscore,
                            data = data,family = binomial(link ='logit'),
                            thresholds = seq(0,1, by = 0.01),
                            confidence.intervals= 0.95,
                            study.design = 'case-control',
                            population.prevalence = 0.3)
simple_gcgr<-decision_curve(status ~ pathologicStage,
                            data = data,family = binomial(link ='logit'),
                            thresholds = seq(0,1, by = 0.01),
                            confidence.intervals= 0.95,
                            study.design = 'case-control',
                            population.prevalence = 0.3)
simple_ag<-decision_curve(status ~ age,
                            data = data,family = binomial(link ='logit'),
                            thresholds = seq(0,1, by = 0.01),
                            confidence.intervals= 0.95,
                            study.design = 'case-control',
                            population.prevalence = 0.3)
####复杂模型---多个自变量
complex<-decision_curve(status ~ age + riskscore + pathologicStage,
                        data = data,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence = 0.3)

#####创建list将所有模型画在一张图中
List<- list(simple_penk,simple_gcgr,simple_ag, complex)
plot_decision_curve(List,
                    curve.names=c('Riskscore','TNM Stage', 'Age', 'Nomogram'),
                    cost.benefit.axis =T,col= c("#0072B2", "#E69F00", "#009E73", "#D55E33"),
                    confidence.intervals=F, 
                    standardize = T)


##AUC

library("timeROC")
cox_auc=data
# 将gender列转换为数值
cox_auc$gender <- ifelse(cox_auc$gender == "male", 1, 0)

cox_auc$Riskscore=as.numeric(cox_auc$riskscore)
cox_auc$pathologicT=as.numeric(as.factor(cox_auc$pathologicT))
cox_auc$pathologicN=as.numeric(as.factor(cox_auc$pathologicN))
cox_auc$Tumor.grade=as.numeric(as.factor(cox_auc$Tumor.grade))
cox_auc$Stage=as.numeric(as.factor(cox_auc$pathologicStage))
cox_auc$gender=as.numeric(as.factor(cox_auc$gender))
cox_auc$age=as.numeric(as.factor(cox_auc$age))

ROC.DSST.age=timeROC(T=cox_auc$time,
                     delta=cox_auc$status,
                     marker=cox_auc$age,
                     cause=1,weighting="marginal",
                     times=1:5,
                     iid=TRUE)
ROC.DSST.gender=timeROC(T=cox_auc$time,
                        delta=cox_auc$status,
                        marker=cox_auc$gender,
                        cause=1,weighting="marginal",
                        times=1:5,
                        iid=TRUE)
ROC.DSST.pathologicT=timeROC(T=cox_auc$time,
                         delta=cox_auc$status,
                         marker=cox_auc$pathologicT,
                         cause=1,weighting="marginal",
                         times=1:5,
                         iid=TRUE)

ROC.DSST.pathologicN=timeROC(T=cox_auc$time,
                         delta=cox_auc$status,
                         marker=cox_auc$pathologicN,
                         cause=1,weighting="marginal",
                         times=1:5,
                         iid=TRUE)

ROC.DSST.Tumor.grade=timeROC(T=cox_auc$time,
                         delta=cox_auc$status,
                         marker=cox_auc$Tumor.grade,
                         cause=1,weighting="marginal",
                         times=1:5,
                         iid=TRUE)

ROC.DSST.pathologicStage=timeROC(T=cox_auc$time,
                       delta=cox_auc$status,
                       marker=cox_auc$Stage,
                       cause=1,weighting="marginal",
                       times=1:5,
                       iid=TRUE)

ROC.DSST.Risk=timeROC(T=cox_auc$time,
                      delta=cox_auc$status,
                      marker=cox_auc$riskscore,
                      cause=1,weighting="marginal",
                      times=1:5,
                      iid=TRUE)
ROC.DSST.Nomo<-timeROC(T=cox_auc$time,
                       delta=cox_auc$status,
                       marker=cox_auc$Riskscore,
                       other_markers=as.matrix(cox_auc[,c("pathologicT","pathologicN")]),
                       cause=1,
                       weighting="cox",
                       times=1:5,
                       iid=F)
ROC.DSST.age$AUC
ROC.DSST.gender$AUC
ROC.DSST.pathologicT$AUC
ROC.DSST.pathologicN$AUC
ROC.DSST.Tumor.grade$AUC
ROC.DSST.pathologicStage$AUC
ROC.DSST.Risk$AUC
ROC.DSST.Nomo$AUC

mg_colors=ggsci::pal_lancet()(9)
pdf('results/AUC.pdf',height = 7,width = 7)
plotAUCcurve(ROC.DSST.Nomo,conf.int=F,col=mg_colors[1])
plotAUCcurve(ROC.DSST.Risk,conf.int=F,col=mg_colors[2],add=TRUE)
plotAUCcurve(ROC.DSST.pathologicStage,conf.int=F,col=mg_colors[3],add=TRUE)
plotAUCcurve(ROC.DSST.age,conf.int=F,col=mg_colors[4],add=TRUE)
plotAUCcurve(ROC.DSST.pathologicT,conf.int=F,col=mg_colors[5],add=TRUE)
plotAUCcurve(ROC.DSST.pathologicN,conf.int=F,col=mg_colors[6],add=TRUE)
plotAUCcurve(ROC.DSST.Tumor.grade,conf.int=F,col=mg_colors[7],add=TRUE)
plotAUCcurve(ROC.DSST.gender,conf.int=F,col=mg_colors[8],add=TRUE)

legend("topright",c("Nomogram","RiskScore", 'Tumor grade', 'TNM Stage',"Age",'T Stage','N Stage','Gender')
       ,col=mg_colors[c(1:8)],lty=1,lwd=2)

dev.off()