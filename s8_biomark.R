###########################################################


load('LBF_protein_exp.RData')

protein<-as.data.frame(protein);rownames(protein)<-protein$Protein;protein<-protein[,-1]

f<-function(x) sum(is.na(x))
rownum<-apply(protein,1,f)
protein<-protein[rownum < ncol(protein)*0.95 ,]

protein <- t(apply(protein, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

lncexp<-protein[grep('sPep',rownames(protein)),] 
protein<-protein[!(rownames(protein) %in% rownames(lncexp)),]

group<-substr(colnames(lncexp),6,6)
#
library(glmnet)
set.seed(23725)

cv<-cv.glmnet(as.matrix(t(protein)),group,nfold=10,family='binomial')
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(protein)),group,family = "binomial")
plot(fit,xvar = "lambda",label = TRUE, las=1)

pred_class <-predict(cv, as.matrix(t(protein)), s = "lambda.min", type = "response") #class

cv<-cv.glmnet(as.matrix(t(lncexp)),group,nfold=10,family='binomial')
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(lncexp)),group,family = "binomial")
plot(fit,xvar = "lambda",label = TRUE, las=1)

pred_class2 <-predict(cv, as.matrix(t(lncexp)), s = "lambda.min", type = "response") #class


modelA <- glm(ifelse(group=='T',1,0) ~ pred_class, family = "binomial" )
modelB <- glm(ifelse(group=='T',1,0) ~ pred_class + pred_class2, family = "binomial" )
modelC <- glm(ifelse(group=='T',1,0) ~ pred_class2, family = "binomial" )

anova(modelA, modelB, test = "LRT") 
anova(modelC, modelB, test = "LRT") 


rocA <- roc(ifelse(group=='T',1,0),  predict(modelA, type="response"))
rocB <- roc(ifelse(group=='T',1,0), predict(modelB, type="response"))
rocC <- roc(ifelse(group=='T',1,0),  predict(modelC, type="response"))

g <- ggroc(list(Protein=rocA,merge=rocB,lncPep=rocC),size=0.8)+theme_bw(base_size=18)+
  theme(panel.grid.major =element_blank(), 
        panel.background=element_rect(size =1.1,fill='transparent', color='black'),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        legend.title = element_blank())+
  theme(legend.key = element_blank(),
        legend.title = element_blank())+
  geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
  scale_x_reverse(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1)+
  coord_fixed(ratio = 0.8)+
  labs(x="Sensitivity", y = "Specificity")+
  scale_colour_manual(values=met.brewer("Derain", 5,type ='discrete')[c(1,2,4)])

pdf("./biomark/LBF_classfier_ROC.pdf",width = 6,height = 6)
g
dev.off()

##

load('TMT_exp.RData')

protein<-as.data.frame(exparray);rownames(protein)<-protein$Protein;protein<-protein[,-1]

f<-function(x) sum(is.na(x))
rownum<-apply(protein,1,f)
protein<-protein[rownum < ncol(protein)*0.95 ,] 

protein <- t(apply(protein, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

anno<-read.table('./clinic/captac_anno.txt',header = T)

library(tidyr)
long_anno <- pivot_longer(anno, cols = c(X127N, X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131N,X131C), names_to = "Variable", values_to = "Value")

pca1 <- prcomp(t(protein),center = TRUE,scale. = TRUE)
df1 <- pca1$x 
df1 <- as.data.frame(df1) 

anno<-long_anno[match(rownames(df1),long_anno$Value),]
group<-substr(anno$Value,1,1)

protein <- limma::removeBatchEffect(protein, batch=anno$AnalyticalSample, 
                                    group=substr(anno$Value,1,1))
lncexp<-protein[grep('sPep',rownames(protein)),] 

protein<-protein[!(rownames(protein) %in% rownames(lncexp)),]

set.seed(23725)

cv<-cv.glmnet(as.matrix(t(protein)),group,nfold=10,family='binomial')
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(protein)),group,family = "binomial")
plot(fit,xvar = "lambda",label = TRUE, las=1)

pred_class <-predict(cv, as.matrix(t(protein)), s = "lambda.min", type = "response") #class

cv<-cv.glmnet(as.matrix(t(lncexp)),group,nfold=10,family='binomial')
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(lncexp)),group,family = "binomial")
plot(fit,xvar = "lambda",label = TRUE, las=1)

pred_class2 <-predict(cv, as.matrix(t(lncexp)), s = "lambda.min", type = "response") #class


modelA <- glm(ifelse(group=='T',1,0) ~ pred_class, family = "binomial" )
modelB <- glm(ifelse(group=='T',1,0) ~ pred_class + pred_class2, family = "binomial" )
modelC <- glm(ifelse(group=='T',1,0) ~ pred_class2, family = "binomial" )

anova(modelA, modelB, test = "LRT") 
anova(modelC, modelB, test = "LRT") 

rocA <- roc(ifelse(group=='T',1,0),  predict(modelA, type="response"))
rocB <- roc(ifelse(group=='T',1,0), predict(modelB, type="response"))
rocC <- roc(ifelse(group=='T',1,0),  predict(modelC, type="response"))

g <- ggroc(list(Protein=rocA,merge=rocB,lncPep=rocC),size=0.8)+theme_bw(base_size=18)+
  theme(panel.grid.major =element_blank(), 
        panel.background=element_rect(size =1.1,fill='transparent', color='black'),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        legend.title = element_blank())+
  theme(legend.key = element_blank(),
    legend.title = element_blank())+
  geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
  scale_x_reverse(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1)+
  coord_fixed(ratio = 0.8)+
  labs(x="Sensitivity", y = "Specificity")+
  scale_colour_manual(values=met.brewer("Derain", 5,type ='discrete')[c(1,2,4)])

pdf("./biomark/TMT_classfier_ROC.pdf",width = 6,height = 6)
g
dev.off()

#####################################
#####################################
#Prognostic models


load('LBF_protein_exp.RData')

#lncpep
lncexp<-protein[grep('sPep',protein$Protein),] 

f<-function(x) sum(is.na(x))
rownum<-apply(lncexp,1,f)
lncexp<-lncexp[rownum < ncol(lncexp)*0.95 ,] 

lncexp <- t(apply(lncexp, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

lncexp<-as.data.frame(lncexp)
colnames(lncexp)<-substr(colnames(lncexp),1,6)
rownames(lncexp)<-lncexp$Protei
lncexp<-lncexp[,-1]

#protein
load('LBF_protein_exp.RData')
protein<-protein[-grep('sPep',protein$Protein),] 

f<-function(x) sum(is.na(x))
rownum<-apply(protein,1,f)
protein<-protein[rownum < ncol(protein)*0.8 ,] 

protein <- t(apply(protein, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

protein<-as.data.frame(protein)
colnames(protein)<-substr(colnames(protein),1,6)
rownames(protein)<-protein$Protei
protein<-protein[,-1]

#meta
meta<-data.table::fread('clinic/iprox_meta.csv')
meta<-meta[meta$`Profiling data QC status`=='Pass',]
meta<-meta[meta$`Phospho- proteome datab`==1,]

meta<-data.frame(sample=c(meta$`Case No.`),
                 OS=c(meta$`Died of recurrence`),
                 OStime=c(meta$`Total follow up period (m)`),
                 rec=c(meta$`Cancer recurrence`),
                 rectime=c(meta$`Disease free survival`))
meta$sample<-paste0(meta$sample,'_T')

lncexp<-lncexp[,match(meta$sample,colnames(lncexp))]
protein<-protein[,match(meta$sample,colnames(protein))]

lncexp <- as.data.frame(
  lapply(lncexp, function(x) {
    if (is.factor(x) || is.character(x)) {
      as.numeric(as.character(x))} else {x}}))
protein <- as.data.frame(
  lapply(protein, function(x) {
    if (is.factor(x) || is.character(x)) {
      as.numeric(as.character(x))} else {x}}))

# mad<-apply(lncexp,1,mad)
# lncexp<-lncexp[mad>1,]

mad<-apply(protein,1,mad)
protein<-protein[mad>2,]

#
library(glmnet)

net<-data.frame(time=c(meta$OStime)*30,status=c(meta$OS))
cv<-cv.glmnet(as.matrix(t(protein)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(protein)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class <-predict(cv, as.matrix(t(protein)), s = "lambda.min", type = "response") #class

cv<-cv.glmnet(as.matrix(t(lncexp)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(lncexp)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class2 <-predict(cv, as.matrix(t(lncexp)), s = "lambda.min", type = "response") #class

library(survival)
library(timeROC)
library(survivalROC)

modelA <- coxph(Surv(time, status) ~ pred_class, data = net)
modelB <- coxph(Surv(time, status) ~ pred_class + pred_class2, data = net)
modelC <- coxph(Surv(time, status) ~ pred_class2, data = net)

anova(modelA, modelB, test = "LRT") 
anova(modelC, modelB, test = "LRT") 

net$pred_class3<- predict(modelB, type="lp")

p<-plottimeroc(net)

library(patchwork)
pdf("biomark/LBF_timeROC.pdf",width = 15,height = 4)
wrap_plots(p,nrow=1, guides="keep") 
dev.off()

#

net<-data.frame(time=c(meta$rectime)*30,status=c(meta$rec))
cv<-cv.glmnet(as.matrix(t(protein)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(protein)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class <-predict(cv, as.matrix(t(protein)), s = "lambda.min", type = "response") #class

cv<-cv.glmnet(as.matrix(t(lncexp)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(lncexp)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class2 <-predict(cv, as.matrix(t(lncexp)), s = "lambda.min", type = "response") #class


modelA <- coxph(Surv(time, status) ~ pred_class, data = net) 
modelB <- coxph(Surv(time, status) ~ pred_class + pred_class2, data = net)
modelC <- coxph(Surv(time, status) ~ pred_class2, data = net)

anova(modelA, modelB, test = "LRT") 
anova(modelC, modelB, test = "LRT") 

net$pred_class3<- predict(modelB, type="lp")

p<-plottimeroc(net)

library(patchwork)
pdf("biomark/LBF_timeROC_rec.pdf",width = 15,height = 4)
wrap_plots(p,nrow=1, guides="keep") 
dev.off()


#########

load('TMT_exp.RData')

protein<-as.data.frame(exparray);rownames(protein)<-protein$Protein;protein<-protein[,-1]

f<-function(x) sum(is.na(x))
rownum<-apply(protein,1,f)
protein<-protein[rownum < ncol(protein)*0.95 ,] 

protein <- t(apply(protein, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

anno<-read.table('./clinic/captac_anno.txt',header = T)
long_anno <- pivot_longer(anno, cols = c(X127N, X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131N,X131C), names_to = "Variable", values_to = "Value")
anno<-long_anno[match(colnames(protein),long_anno$Value),]
protein <- limma::removeBatchEffect(protein, batch=anno$AnalyticalSample, 
                                    group=substr(anno$Value,1,1))

lncexp<-protein[grep('sPep',rownames(protein)),]
protein<-protein[-grep('sPep',rownames(protein)),]

# mad<-apply(lncexp,1,mad)
# lncexp<-lncexp[mad>0,]

mad<-apply(protein,1,mad)
protein<-protein[mad>0.2,]

meta<-data.table::fread('clinic/captac_meta.csv')

meta<-data.frame(sample=c(substr(meta$`Tumor (T) sample ID`,2,99)),
                 OS=c(meta$`Survial  (1, dead; 0, alive)`),
                 OStime=c(meta$`Overall survial (month)`),
                 rec=c(meta$`Recurrence  (1, yes; 0, no)`),
                 rectime=c(meta$`Recurrence-free survival (month)`)) 
meta$sample<-paste0('T',meta$sample)

lncexp<-lncexp[,match(meta$sample,colnames(lncexp))]
protein<-protein[,match(meta$sample,colnames(protein))]

#

net<-data.frame(time=c(meta$OStime)*30,status=c(meta$OS))
cv<-cv.glmnet(as.matrix(t(protein)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(protein)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class <-predict(cv, as.matrix(t(protein)), s = "lambda.min", type = "response") #class

cv<-cv.glmnet(as.matrix(t(lncexp)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(lncexp)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class2 <-predict(cv, as.matrix(t(lncexp)), s = "lambda.min", type = "response") #class

modelA <- coxph(Surv(time, status) ~ pred_class, data = net)
modelB <- coxph(Surv(time, status) ~ pred_class + pred_class2, data = net)
modelC <- coxph(Surv(time, status) ~ pred_class2, data = net)

anova(modelA, modelB, test = "LRT") 
anova(modelC, modelB, test = "LRT") 

net$pred_class3<- predict(modelB, type="lp")

p<-plottimeroc(net)

library(patchwork)
pdf("biomark/TMT_timeROC.pdf",width = 15,height = 4)
wrap_plots(p,nrow=1, guides="keep") 
dev.off()

#

net<-data.frame(time=c(meta$rectime)*30,status=c(meta$rec))
cv<-cv.glmnet(as.matrix(t(protein)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(protein)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class <-predict(cv, as.matrix(t(protein)), s = "lambda.min", type = "response") #class

cv<-cv.glmnet(as.matrix(t(lncexp)),as.matrix(net),nfold=10,family='cox',alpha  = 0.5,lambda.min.ratio = 1e-5)
fit.train <- cv$glmnet.fit
plot(cv,las=1)

fit <- glmnet(as.matrix(t(lncexp)),as.matrix(net),family = "cox",alpha  = 0.5)
plot(fit,xvar = "lambda",label = TRUE, las=1)

net$pred_class2 <-predict(cv, as.matrix(t(lncexp)), s = "lambda.min", type = "response") #class


modelA <- coxph(Surv(time, status) ~ pred_class, data = net) #为什么复发数据这么差，
modelB <- coxph(Surv(time, status) ~ pred_class + pred_class2, data = net)
modelC <- coxph(Surv(time, status) ~ pred_class2, data = net)

anova(modelA, modelB, test = "LRT") 
anova(modelC, modelB, test = "LRT") 

net$pred_class3<- predict(modelB, type="lp")

p<-plottimeroc(net)

library(patchwork)
pdf("biomark/TMT_timeROC_rec.pdf",width = 15,height = 4)
wrap_plots(p,nrow=1, guides="keep") 
dev.off()


