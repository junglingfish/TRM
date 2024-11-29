# install.packages("Rtsne")
# install.packages("ggplot2")

setwd('D:/ZJU-FISH/doctor/TRM/10cluster')
routine = 'D:/ZJU-FISH/doctor/TRM/'
routine_data = 'D:/ZJU-FISH/doctor/TRM/data/'
routine_4 = 'D:/ZJU-FISH/doctor/TRM/04model/'
dir.create('results/PCA/')

#???ð?
library(Rtsne)
library(ggplot2)
library(ggpubr)

clusterFile="results/cluster_risk.csv"     #???͵Ľ????ļ?


#??ȡ?????ļ?,??ȡ????
rt=read.csv(clusterFile)
rt$Cluster <- as.factor((rt$Cluster))
data=rt[c(3:(ncol(rt)-6))]

#PCA????
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
PCA = data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], rt[,c("Risk","Cluster")])
#???Ʒ??յ?PCAͼ
pdf(file="results/PCA/PCA.risk.pdf", width=5.5, height=4.5)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Risk)) +
  scale_colour_manual(name="Risk",  values =c("red", "blue"))+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()
#???Ʒ??͵?PCAͼ
bioCol=c("#0066FF","#FF9900","#00DB00","#FF0000","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt$Cluster))]
pdf(file="results/PCA/PCA.cluster.pdf", width=5.5, height=4.5)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
  scale_colour_manual(name="Cluster",  values =bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


#t-SNE????
tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
tsne=data.frame(tSNE1=tsneOut$Y[,1], tSNE2=tsneOut$Y[,2], rt[,c("Risk","Cluster")])	
#???Ʒ??յ?tSNEͼ
pdf(file="results/PCA/tSNE.risk.pdf", width=5.5, height=4.5)       #???????????ļ?
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = Risk)) +
  scale_colour_manual(name="Risk",  values =c("red", "blue"))+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()
#???Ʒ??͵?tSNEͼ
pdf(file="results/PCA/tSNE.Cluster.pdf", width=5.5, height=4.5)       #???????????ļ?
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = Cluster)) +
  scale_colour_manual(name="Cluster",  values =bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


####PCA Fig
pca <- PCA
Fig1a.taxa.pca <- ggplot(pca,aes(PC1,PC2))+
  stat_ellipse(aes(fill=pca$Risk), geom = 'polygon', color = NA, alpha = 0.2)+
  geom_point(size=2,aes(col=Risk,shape=Cluster))+ 
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   legend.position = "left")

pca$Group = paste(pca$Risk,"|", pca$Cluster,sep = "")
Fig1a.taxa.pca

Fig1a.taxa.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=Group, fill=Risk, linetype=Cluster),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  theme_bw() +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),
        legend.position = "none") + 
  labs(fill="")
Fig1a.taxa.pc1.density

Fig1a.taxa.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=Group, fill=Risk, linetype=Cluster),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  theme_bw() +
  theme(axis.title.y = element_blank(),  # 隐藏纵轴标题
        axis.text.y = element_blank()) + # 隐藏纵轴刻度数值
  labs(fill="") + 
  coord_flip()
Fig1a.taxa.pc2.density

figs1a <- ggarrange(Fig1a.taxa.pca,Fig1a.taxa.pc2.density,nrow = 1,ncol = 2,widths = c(2,0.75))
fig <- ggarrange(Fig1a.taxa.pc1.density, figs1a,nrow = 2,ncol = 1,heights = c(0.5, 2))
fig
