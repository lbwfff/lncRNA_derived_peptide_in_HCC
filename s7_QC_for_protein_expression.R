#####################################
#Normalization and QC needed before biostatistics
#PCA look at TMT batches, LBF status, pep PCA alone, PLSDA
#lncPep intensity in LBF data is too sparse, which may have made it difficult to be biometrically significant

load('LBF_protein_exp.RData')

protein<-as.data.frame(protein);rownames(protein)<-protein$Protein;protein<-protein[,-1]

f<-function(x) sum(is.na(x))
rownum<-apply(protein,1,f)
protein<-protein[rownum < ncol(protein)*0.95 ,] # Remove too many missing lines first, then fill them in

protein <- t(apply(protein, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

lncexp<-protein[grep('sPep',rownames(protein)),] 

#PCA
library(ggrepel)

pca1 <- prcomp(t(protein),center = TRUE,scale. = TRUE)

df1 <- pca1$x 
df1 <- as.data.frame(df1) 

summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

library(ggplot2)
library(MetBrewer)
group<-substr(rownames(df1),6,6)

pdf("./Figure4/LBF_protein_PCA.pdf",width = 6,height = 6)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 2)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 2)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off() 

#No need to manually Normalize, misstate has already done it

#lncpep
pca1 <- prcomp(t(lncexp),center = TRUE,scale. = TRUE)
df1 <- pca1$x 
df1 <- as.data.frame(df1) 
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

library(ggplot2)
library(MetBrewer)
group<-substr(rownames(df1),6,6)

pdf("./Figure4/LBF_lncPep_PCA.pdf",width = 6,height = 6)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 2)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 2)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off()

#PLSDA
library(mixOmics)
plsda<-plsda(t(lncexp),group,ncomp = 2)
plsda<-plotIndiv(plsda, ind.names = group, ellipse = TRUE, legend =TRUE)
plsda<-plsda[["df"]]

pdf("./Figure4/LBF_lncPep_plsda.pdf",width = 6,height = 6)
ggplot(data = plsda,aes(x = x,y = y,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 
  geom_point(size = 3.5)+
  labs(x = 'PLSDA comp1',y = 'PLSDA comp2',color = "Condition",title = "PLSDA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 2)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 2)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off() #没问题

##TMT
load('TMT_exp.RData')

protein<-as.data.frame(exparray);rownames(protein)<-protein$Protein;protein<-protein[,-1]

f<-function(x) sum(is.na(x))
rownum<-apply(protein,1,f)
protein<-protein[rownum < ncol(protein)*0.95 ,] #先去除过于稀疏的，然后再补全

protein <- t(apply(protein, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

#PCA looks at T_vs_P, and in addition also looks at batch effects
anno<-read.table('./clinic/captac_anno.txt',header = T)

library(tidyr)
long_anno <- pivot_longer(anno, cols = c(X127N, X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131N,X131C), names_to = "Variable", values_to = "Value")

pca1 <- prcomp(t(protein),center = TRUE,scale. = TRUE)

df1 <- pca1$x 
df1 <- as.data.frame(df1) 

summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")


anno<-long_anno[match(rownames(df1),long_anno$Value),]
group<-substr(anno$Value,1,1)


pdf("./Figure4/TMT_protein_PCA.pdf",width = 6,height = 6)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 2)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 2)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off()

group<-anno$AnalyticalSample

pdf("./Figure4/TMT_protein_PCA_batch.pdf",width = 6,height = 6)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 33)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 33)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off() #TMT data has a batch effect and needs correction

#limma

protein <- limma::removeBatchEffect(protein, batch=anno$AnalyticalSample, 
                                    group=substr(anno$Value,1,1))

pca1 <- prcomp(t(protein),center = TRUE,scale. = TRUE)
df1 <- pca1$x 
df1 <- as.data.frame(df1) 
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

group<-substr(anno$Value,1,1)

pdf("./Figure4/TMT_protein_PCA_removebatch.pdf",width = 6,height = 6)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 2)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 2)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off()

group<-anno$AnalyticalSample

pdf("./Figure4/TMT_protein_PCA_removebatch_batch.pdf",width = 6,height = 6)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 33)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 33)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off() 

#lncPep
lncexp<-protein[grep('sPep',rownames(protein)),] 

pca1 <- prcomp(t(lncexp),center = TRUE,scale. = TRUE)
df1 <- pca1$x 
df1 <- as.data.frame(df1) 
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

library(ggplot2)
library(MetBrewer)
group<-substr(anno$Value,1,1)

pdf("./Figure4/TMT_lncPep_PCA.pdf",width = 6,height = 6)
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 2)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 2)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off()

#PLSDA
library(mixOmics)
plsda<-plsda(t(lncexp),group,ncomp = 2)
plsda<-plotIndiv(plsda, ind.names = group, ellipse = TRUE, legend =TRUE)
plsda<-plsda[["df"]]

pdf("./Figure4/TMT_lncPep_plsda.pdf",width = 6,height = 6)
ggplot(data = plsda,aes(x = x,y = y,color = group))+
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  geom_point(size = 3.5)+
  labs(x = 'PLSDA comp1',y = 'PLSDA comp2',color = "Condition",title = "PLSDA Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = rev(met.brewer("Troy", 2)))+
  scale_colour_manual(values = rev(met.brewer("Troy", 2)))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))+
  theme(aspect.ratio=1)
dev.off()


