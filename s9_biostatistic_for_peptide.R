######################################
#Here are the main biostatistics codes

load('LBF_protein_exp.RData')

lncexp<-protein[grep('sPep',protein$Protein),] 


# myindex<-TNSMD::generate_index('./writeindex_2.csv','other',0,'sPep')
# 
# f <-function(x) unlist(strsplit(x['name'],'[|]'))[2]
# myindex$ID<-apply(myindex,1,f)
# f <-function(x) unlist(strsplit(x['name'],'[|]'))[3]
# myindex$ORFname<-apply(myindex,1,f)
# 
# f <-function(x) paste0(unlist(strsplit(x['ORFname'],' '))[2],' ',unlist(strsplit(x['ORFname'],' '))[3],' ',unlist(strsplit(x['ORFname'],' '))[4])
# myindex$ORFname<-apply(myindex,1,f)

f <-function(x) unlist(strsplit(x['Protein'],'[|]'))[2]
lncexp$Protein<-apply(lncexp,1,f)
match<-myindex[match(lncexp$Protein,myindex$ID),]

lncexp<-as.data.frame(lncexp)
rownames(lncexp)<-match$ORFname
lncexp<-lncexp[,-1]

meta<-data.table::fread('clinic/iprox_meta.csv')
meta<-meta[meta$`Profiling data QC status`=='Pass',]
meta<-meta[meta$`Phospho- proteome datab`==1,]

meta<-data.frame(sample=c(meta$`Case No.`),
                 OS=c(meta$`Died of recurrence`),
                 OStime=c(meta$`Total follow up period (m)`),
                 rec=c(meta$`Cancer recurrence`),
                 rectime=c(meta$`Disease free survival`))

colnames(lncexp)<-substr(colnames(lncexp),1,6)

biostat<-function(exparray,meta){
  or<-data.frame()
  os<-data.frame()
  rec<-data.frame()
  
  exparray<-exparray[,grep(paste0(meta$sample,collapse = '|'),colnames(exparray))]
  
  for (i in 1:nrow(exparray)){
    
    exp<-as.data.frame(t(exparray[i,]))
    exp$group<-c(NA)
    colnames(exp)[1]<-c('intensity')
    exp<-exp[!is.na(exp$intensity),]
    exp$intensity<-c(scale(exp$intensity))
    
    exp$group<-ifelse(grepl('T',rownames(exp)),1,0)
    
    if (nrow(exp[exp$group==1,])>=nrow(meta)*0.025 & nrow(exp[exp$group==0,])>=nrow(meta)*0.025) {
      
      #计算肿瘤和非肿瘤组逻辑回归的OR值
      
      fit.full<-glm(group~.,
                    data=exp,family = binomial())
      
      fit.result<-summary(fit.full)
      
      result<-as.data.frame(exp(cbind(OR = coef(fit.full), 
                                      confint(fit.full))))
      result$p<-as.numeric(fit.result$coefficients[,4])
      result<-data.frame(result[-1,c(1,4,2,3)]) 
      result$Var<-rownames(exparray)[i]
      colnames(result)<-c("OR","Pvalue","OR_1","OR_2",'lncPep')
      
      or<-rbind(or,result)
      
    } else{
      print(paste0('Skip logistic regression for ',rownames(exparray)[i]))
    }
    if (nrow(exp[exp$group==1,])>=nrow(meta)*0.05) {
      
      coxarray<-exp[exp$group==1,]
      match<-meta[match(rownames(coxarray),paste0(meta$sample,'_T')),]
      
      coxarray$OStime<-match$OStime
      coxarray$OS<-match$OS
      coxarray$rectime<-match$rectime
      coxarray$rec<-match$rec
      
      coxarray<-coxarray[!is.na(coxarray$OStime),]
      
      library("survival")
      library("survminer")
      #生存
      suv.cox<-coxph(Surv(OStime,OS)~intensity, data=coxarray) 
      
      result2<-data.frame(OR=c(summary(suv.cox)$coef[,'exp(coef)']),
                          Pvalue=c(summary(suv.cox)$coef[, "Pr(>|z|)"]),
                          OR_1=c(summary(suv.cox)$conf.int[, "lower .95"]),
                          OR_2=c(summary(suv.cox)$conf.int[, "upper .95"]),
                          lncPep=c(rownames(exparray)[i]))
      
      
      #复发
      res.cox<-coxph(Surv(rectime,rec)~intensity, data=coxarray) 
      
      result3<-data.frame(OR=c(summary(res.cox)$coef[,'exp(coef)']),
                          Pvalue=c(summary(res.cox)$coef[, "Pr(>|z|)"]),
                          OR_1=c(summary(res.cox)$conf.int[, "lower .95"]),
                          OR_2=c(summary(res.cox)$conf.int[, "upper .95"]),
                          lncPep=c(rownames(exparray)[i]))
      
      os<-rbind(os,result2)
      rec<-rbind(rec,result3)
      
    } else{
      print(paste0('Skip Cox regression for ',rownames(exparray)[i]))
    }
    
  }
  
  return(list(or,os,rec))
  
}

LBF_unimpu<-biostat(lncexp,meta)

########################################
#minimum value imputation

load('LBF_protein_exp.RData')

lncexp<-protein[grep('sPep',protein$Protein),] 

f <-function(x) unlist(strsplit(x['Protein'],'[|]'))[2]
lncexp$Protein<-apply(lncexp,1,f)
match<-myindex[match(lncexp$Protein,myindex$ID),]

lncexp<-as.data.frame(lncexp)
rownames(lncexp)<-match$ORFname
lncexp<-lncexp[,-1]

f<-function(x) sum(is.na(x))
rownum<-apply(lncexp,1,f)
lncexp<-lncexp[rownum < ncol(lncexp)*0.95 ,] 

lncexp <- t(apply(lncexp, 1, function(row) {
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

lncexp<-as.data.frame(lncexp)
colnames(lncexp)<-substr(colnames(lncexp),1,6)

LBF_minimpu<-biostat(lncexp,meta)

save(LBF_unimpu,LBF_minimpu,file = 'Figure4/LBF_suv.RData')
##############################################

load('TMT_exp.RData') 

protein<-as.data.frame(exparray);rownames(exparray)<-exparray$Protein;exparray<-exparray[,-1]

anno<-read.table('./clinic/captac_anno.txt',header = T)
long_anno <- pivot_longer(anno, cols = c(X127N, X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131N,X131C), names_to = "Variable", values_to = "Value")
anno<-long_anno[match(colnames(exparray),long_anno$Value),]
exparray <- limma::removeBatchEffect(exparray, batch=anno$AnalyticalSample, 
                                    group=substr(anno$Value,1,1)) 

lncexp<-exparray[grep('sPep',rownames(exparray)),]
match<-myindex[match(rownames(lncexp),myindex$ID),]
lncexp<-as.data.frame(lncexp)
rownames(lncexp)<-match$ORFname

meta<-data.table::fread('clinic/captac_meta.csv')

meta<-data.frame(sample=c(substr(meta$`Tumor (T) sample ID`,2,99)),
                 OS=c(meta$`Survial  (1, dead; 0, alive)`),
                 OStime=c(meta$`Overall survial (month)`),
                 rec=c(meta$`Recurrence  (1, yes; 0, no)`),
                 rectime=c(meta$`Recurrence-free survival (month)`)) 

biostat<-function(exparray,meta){
  or<-data.frame()
  os<-data.frame()
  rec<-data.frame()
  
  for (i in 1:nrow(exparray)){
    
    exp<-as.data.frame(t(exparray[i,]))
    exp$group<-c(NA)
    colnames(exp)[1]<-c('intensity')
    exp<-exp[!is.na(exp$intensity),]
    exp$intensity<-c(scale(exp$intensity))
    
    exp$group<-ifelse(grepl('T',rownames(exp)),1,0)
    
    if (nrow(exp[exp$group==1,])>=nrow(meta)*0.025 & nrow(exp[exp$group==0,])>=nrow(meta)*0.025) {
      
      fit.full<-glm(group~.,
                    data=exp,family = binomial())
      fit.result<-summary(fit.full)

      if (fit.result$coefficients[2,'Pr(>|z|)']<0.99) {
      
      coefs <- coef(fit.full)
      confint_values <- confint(fit.full) 
      result <-as.data.frame(exp(cbind(OR = coef(fit.full), confint(fit.full))))
        
      result$p<-as.numeric(fit.result$coefficients[,4])
      result<-data.frame(result[-1,c(1,4,2,3)]) 
      result$Var<-rownames(exparray)[i]
      colnames(result)<-c("OR","Pvalue","OR_1","OR_2",'lncPep')
      
      or<-rbind(or,result) } else{print(paste0('Skip logistic regression for ',rownames(exparray)[i]))}
      
    } else{
      print(paste0('Skip logistic regression for ',rownames(exparray)[i]))
    }
    if (nrow(exp[exp$group==1,])>=nrow(meta)*0.05) {
      
      coxarray<-exp[exp$group==1,]
      match<-meta[match(rownames(coxarray),paste0('T',meta$sample)),]
      
      coxarray$OStime<-match$OStime
      coxarray$OS<-match$OS
      coxarray$rectime<-match$rectime
      coxarray$rec<-match$rec
      
      coxarray<-coxarray[!is.na(coxarray$OStime),]
      
      library("survival")
      library("survminer")
      #生存
      suv.cox<-coxph(Surv(OStime,OS)~intensity, data=coxarray) 
      
      result2<-data.frame(OR=c(summary(suv.cox)$coef[,'exp(coef)']),
                          Pvalue=c(summary(suv.cox)$coef[, "Pr(>|z|)"]),
                          OR_1=c(summary(suv.cox)$conf.int[, "lower .95"]),
                          OR_2=c(summary(suv.cox)$conf.int[, "upper .95"]),
                          lncPep=c(rownames(exparray)[i]))
      
      
      #复发
      res.cox<-coxph(Surv(rectime,rec)~intensity, data=coxarray) 
      
      result3<-data.frame(OR=c(summary(res.cox)$coef[,'exp(coef)']),
                          Pvalue=c(summary(res.cox)$coef[, "Pr(>|z|)"]),
                          OR_1=c(summary(res.cox)$conf.int[, "lower .95"]),
                          OR_2=c(summary(res.cox)$conf.int[, "upper .95"]),
                          lncPep=c(rownames(exparray)[i]))
      
      os<-rbind(os,result2)
      rec<-rbind(rec,result3)
      
    } else{
      print(paste0('Skip Cox regression for ',rownames(exparray)[i]))
    }
    
  }
  
  output<-list(or=or,
               os=os,
               rec=rec)
  
  return(output)
  
}

TMT_nonimp<-biostat(lncexp,meta)

#minimum value imputation for TMT data

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
match<-myindex[match(rownames(lncexp),myindex$ID),]
lncexp<-as.data.frame(lncexp)
rownames(lncexp)<-match$ORFname

meta<-data.table::fread('clinic/captac_meta.csv')

meta<-data.frame(sample=c(substr(meta$`Tumor (T) sample ID`,2,99)),
                 OS=c(meta$`Survial  (1, dead; 0, alive)`),
                 OStime=c(meta$`Overall survial (month)`),
                 rec=c(meta$`Recurrence  (1, yes; 0, no)`),
                 rectime=c(meta$`Recurrence-free survival (month)`)) #sampleID的命名规则不一样，有点懒得改了

TMT_minimp<-biostat(lncexp,meta)


###############################################################
#heatmap
heatarray<-data.frame(lnPep=unique(c(LBF_TvsP$lncPep,LBF_REC$lncPep,TMT_TvsP$lncPep,TMT_REC$lncPep)),
                      LBF_OR=c(NA), LBF_OR2=c(NA),
                      TMT_OR=c(NA), TMT_OR2=c(NA),
                      LBF_HR3=c(NA), LBF_HR4=c(NA),
                      TMT_HR3=c(NA), TMT_HR4=c(NA),
                      LBF_HR5=c(NA),LBF_HR6=c(NA),
                      TMT_HR5=c(NA), TMT_HR6=c(NA)) 

heatarray_p<-data.frame(lnPep=unique(c(LBF_TvsP$lncPep,LBF_REC$lncPep,TMT_TvsP$lncPep,TMT_REC$lncPep)),
                      LBF_OR=c(NA), LBF_OR2=c(NA),
                      TMT_OR=c(NA), TMT_OR2=c(NA),
                      LBF_HR3=c(NA), LBF_HR4=c(NA),
                      TMT_HR3=c(NA), TMT_HR4=c(NA),
                      LBF_HR5=c(NA),LBF_HR6=c(NA),
                      TMT_HR5=c(NA), TMT_HR6=c(NA))

{
match<-LBF_unimpu[[1]][match(heatarray$lnPep,LBF_unimpu[[1]]$lncPep),]
heatarray$LBF_OR<-match$OR
heatarray_p$LBF_OR<-match$Pvalue

match<-LBF_minimpu[[1]][match(heatarray$lnPep,LBF_minimpu[[1]]$lncPep),]
heatarray$LBF_OR2<-match$OR
heatarray_p$LBF_OR2<-match$Pvalue

match<-LBF_unimpu[[2]][match(heatarray$lnPep,LBF_unimpu[[2]]$lncPep),]
heatarray$LBF_HR3<-match$OR
heatarray_p$LBF_HR3<-match$Pvalue

match<-LBF_minimpu[[2]][match(heatarray$lnPep,LBF_minimpu[[2]]$lncPep),]
heatarray$LBF_HR4<-match$OR
heatarray_p$LBF_HR4<-match$Pvalue

match<-LBF_unimpu[[3]][match(heatarray$lnPep,LBF_unimpu[[3]]$lncPep),]
heatarray$LBF_HR5<-match$OR
heatarray_p$LBF_HR5<-match$Pvalue

match<-LBF_minimpu[[3]][match(heatarray$lnPep,LBF_minimpu[[3]]$lncPep),]
heatarray$LBF_HR6<-match$OR
heatarray_p$LBF_HR6<-match$Pvalue

match<-TMT_nonimp[[1]][match(heatarray$lnPep,TMT_nonimp[[1]]$lncPep),]
heatarray$TMT_OR<-match$OR
heatarray_p$TMT_OR<-match$Pvalue

match<-TMT_minimp[[1]][match(heatarray$lnPep,TMT_minimp[[1]]$lncPep),]
heatarray$TMT_OR2<-match$OR
heatarray_p$TMT_OR2<-match$Pvalue

match<-TMT_nonimp[[2]][match(heatarray$lnPep,TMT_nonimp[[2]]$lncPep),]
heatarray$TMT_HR3<-match$OR
heatarray_p$TMT_HR3<-match$Pvalue

match<-TMT_minimp[[2]][match(heatarray$lnPep,TMT_minimp[[2]]$lncPep),]
heatarray$TMT_HR4<-match$OR
heatarray_p$TMT_HR4<-match$Pvalue

match<-TMT_nonimp[[3]][match(heatarray$lnPep,TMT_nonimp[[3]]$lncPep),]
heatarray$TMT_HR5<-match$OR
heatarray_p$TMT_HR5<-match$Pvalue

match<-TMT_minimp[[3]][match(heatarray$lnPep,TMT_minimp[[3]]$lncPep),]
heatarray$TMT_HR6<-match$OR
heatarray_p$TMT_HR6<-match$Pvalue
}

library(ComplexHeatmap)
library(circlize)
library(MetBrewer)

rownames(heatarray)<-heatarray$lnPep
heatarray<-heatarray[,-1]
rownames(heatarray_p)<-heatarray_p$lnPep
heatarray_p<-heatarray_p[,-1]

col_fun1 = colorRamp2(c(0,1,5), c('#22B5AF','white','#F57F17'))

pdf('./Figure4/test_heatmap.pdf',width = 10,height = 15)
draw(
  Heatmap(heatarray,na_col = "grey",col = col_fun1,
          name = "OR/HR value",
          column_title = NULL,
          cluster_columns = F,
          cluster_rows =F,
          show_row_dend = T,
          show_column_names =T,
          show_heatmap_legend = T,
          row_names_side = "left",
          width = ncol(heatarray)*unit(4, "mm"), 
          height = nrow(heatarray)*unit(4, "mm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(heatarray_p[i, j] < 0.05 & !is.na(heatarray_p[i, j]))
              grid.text('*', x, y, gp = gpar(fontsize = 10))
          }
  )
)
dev.off() 


or_sig<-unique(c(rownames(heatarray)[heatarray_p$LBF_OR<0.05],
                 rownames(heatarray)[heatarray_p$LBF_OR2<0.05],
                 rownames(heatarray)[heatarray_p$TMT_OR<0.05],
                 rownames(heatarray)[heatarray_p$TMT_OR2<0.05]))
suv_sig<-unique(c(c(rownames(heatarray)[heatarray_p$LBF_HR3<0.05],
                    rownames(heatarray)[heatarray_p$LBF_HR4<0.05],
                    rownames(heatarray)[heatarray_p$TMT_HR3<0.05],
                    rownames(heatarray)[heatarray_p$TMT_HR4<0.05])))
rec_sig<-unique(c(c(rownames(heatarray)[heatarray_p$LBF_HR5<0.05],
                    rownames(heatarray)[heatarray_p$LBF_HR6<0.05],
                    rownames(heatarray)[heatarray_p$TMT_HR5<0.05],
                    rownames(heatarray)[heatarray_p$TMT_HR6<0.05])))

library(ggvenn)
library(MetBrewer)

venn<-list('Cancer-related'=or_sig[!is.na(or_sig)],
           'Survival-related'=suv_sig[!is.na(suv_sig)],
           'Recurrence-related'=rec_sig[!is.na(rec_sig)])

pdf("./Figure4/VENN.pdf",width = 4,height = 4)
ggvenn(venn,fill_color = met.brewer('Hiroshige',n=3),
       stroke_size = 0.5, set_name_size = 4)
dev.off()

##############################################################
#Box lines and forest plots

orplot<-rbind(LBF_unimpu[[1]],LBF_minimpu[[1]],TMT_nonimp[[1]],TMT_minimp[[1]])
orplot$source<-c(rep('Jiang et al 2019',nrow(LBF_unimpu[[1]])),
                  rep('Jiang et al 2019',nrow(LBF_minimpu[[1]])),
                  rep('Gao et al 2019',nrow(TMT_nonimp[[1]])),
                  rep('Gao et al 2019',nrow(TMT_minimp[[1]])))

orplot$type<-c(rep('Remove missing values',nrow(LBF_unimpu[[1]])),
               rep('minimum value imputation',nrow(LBF_minimpu[[1]])),
               rep('Remove missing values',nrow(TMT_nonimp[[1]])),
               rep('minimum value imputation',nrow(TMT_minimp[[1]])))

orplot<-orplot[orplot$lncPep %in% c('PPIAP79 lncORF 1','SEPTIN7P8 lncORF 1'),]

p<-list()

p[[1]]<-
ggplot(orplot[orplot$source=='Jiang et al 2019',], aes(y=reorder(lncPep,OR), x=OR,  colour=type)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,5))+
  theme(legend.position = "bottom")+
  theme(aspect.ratio=0.5)

p[[2]]<-
ggplot(orplot[orplot$source!='Jiang et al 2019',], aes(y=reorder(lncPep,OR), x=OR,  colour=type)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  # coord_cartesian(xlim=c(0.1,5))+
  theme(legend.position = "bottom")+
  theme(aspect.ratio=0.5)

pdf("./Figure4/twolncpep_forest.pdf",width = 8,height = 6)
wrap_plots(p,nrow=2, guides="collect") 
dev.off()


{
  ##四个表达矩阵
  
  #
  load('LBF_protein_exp.RData')
  lncexp<-protein[grep('sPep',protein$Protein),]
  f <-function(x) unlist(strsplit(x['Protein'],'[|]'))[2]
  lncexp$Protein<-apply(lncexp,1,f)
  match<-myindex[match(lncexp$Protein,myindex$ID),]
  lncexp<-as.data.frame(lncexp)
  rownames(lncexp)<-match$ORFname
  lbf_unimpexp<-lncexp[,-1]
  
  #
  load('LBF_protein_exp.RData')
  lncexp<-protein[grep('sPep',protein$Protein),] 
  f <-function(x) unlist(strsplit(x['Protein'],'[|]'))[2]
  lncexp$Protein<-apply(lncexp,1,f)
  match<-myindex[match(lncexp$Protein,myindex$ID),]
  lncexp<-as.data.frame(lncexp)
  rownames(lncexp)<-match$ORFname
  lncexp<-lncexp[,-1]
  f<-function(x) sum(is.na(x))
  rownum<-apply(lncexp,1,f)
  lncexp<-lncexp[rownum < ncol(lncexp)*0.95 ,] 
  lncexp <- t(apply(lncexp, 1, function(row) {
    row[is.na(row)] <- min(row, na.rm = TRUE)
    return(row)}))
  lncexp<-as.data.frame(lncexp)
  colnames(lncexp)<-substr(colnames(lncexp),1,6)
  lbf_impexp<-lncexp
  
  #
  
  load('TMT_exp.RData') 
  protein<-as.data.frame(exparray);rownames(exparray)<-exparray$Protein;exparray<-exparray[,-1]
  anno<-read.table('./clinic/captac_anno.txt',header = T)
  long_anno <- pivot_longer(anno, cols = c(X127N, X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131N,X131C), names_to = "Variable", values_to = "Value")
  anno<-long_anno[match(colnames(exparray),long_anno$Value),]
  exparray <- limma::removeBatchEffect(exparray, batch=anno$AnalyticalSample, 
                                       group=substr(anno$Value,1,1))
  lncexp<-exparray[grep('sPep',rownames(exparray)),]
  match<-myindex[match(rownames(lncexp),myindex$ID),]
  lncexp<-as.data.frame(lncexp)
  rownames(lncexp)<-match$ORFname
  tmt_unimpexp<-lncexp
  
  #
  load('TMT_exp.RData')
  protein<-as.data.frame(exparray);rownames(protein)<-protein$Protein;protein<-protein[,-1]
  f<-function(x) sum(is.na(x))
  rownum<-apply(protein,1,f)
  protein<-protein[rownum < ncol(protein)*0.95 ,] #先去除过于稀疏的，然后再补全
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
  match<-myindex[match(rownames(lncexp),myindex$ID),]
  lncexp<-as.data.frame(lncexp)
  rownames(lncexp)<-match$ORFname
  tmt_impexp<-lncexp
}

p<-list()

lbf<-data.frame(intensity=as.numeric(c(lbf_impexp['PPIAP79 lncORF 1',],lbf_unimpexp['PPIAP79 lncORF 1',])),
                group=c(substr(colnames(lbf_impexp),6,6),substr(colnames(lbf_unimpexp),6,6)),
                type=c(rep('minimum value imputation',ncol(lbf_impexp)),rep('Remove missing values',ncol(lbf_unimpexp))))
lbf<-lbf[!is.na(lbf$intensity),]

p[[1]]<-
ggplot(lbf, 
       aes(x = group, y = intensity)) +
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
  labs(y="log2 (intensity)")+ 
  theme_classic(base_size = 10)+ 
  scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
  theme_classic(base_size = 10)+ 
  theme(legend.position = "none")+
  xlab(NULL)+
  theme(aspect.ratio=2)+
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_line(color = "gray", size = 0.25))+
  facet_wrap(~ type,scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 0.95, label.y = 28.1)

tmt<-data.frame(intensity=as.numeric(c(tmt_impexp['PPIAP79 lncORF 1',],tmt_unimpexp['PPIAP79 lncORF 1',])),
                group=c(substr(colnames(tmt_impexp),1,1),substr(colnames(tmt_unimpexp),1,1)),
                type=c(rep('minimum value imputation',ncol(tmt_impexp)),rep('Remove missing values',ncol(tmt_unimpexp))))
tmt<-tmt[!is.na(tmt$intensity),]

p[[2]]<-
ggplot(tmt, 
       aes(x = group, y = intensity)) +
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
  labs(y="log2 (intensity)")+ 
  theme_classic(base_size = 10)+ 
  scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
  theme_classic(base_size = 10)+ 
  theme(legend.position = "none")+
  xlab(NULL)+
  theme(aspect.ratio=2)+
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_line(color = "gray", size = 0.25))+
  facet_wrap(~ type,scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 0.95, label.y = 19)

library(patchwork)
pdf("./Figure4/RPKM_PPIAP79_o1_box.pdf",width = 8,height = 8)
wrap_plots(p,nrow=2, guides="collect") 
dev.off()

p<-list()

lbf<-data.frame(intensity=as.numeric(c(lbf_impexp['SEPTIN7P8 lncORF 1',],lbf_unimpexp['SEPTIN7P8 lncORF 1',])),
                group=c(substr(colnames(lbf_impexp),6,6),substr(colnames(lbf_unimpexp),6,6)),
                type=c(rep('minimum value imputation',ncol(lbf_impexp)),rep('Remove missing values',ncol(lbf_unimpexp))))
lbf<-lbf[!is.na(lbf$intensity),]

p[[1]]<-
  ggplot(lbf, 
         aes(x = group, y = intensity)) +
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
  labs(y="log2 (intensity)")+ 
  theme_classic(base_size = 10)+ 
  scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
  theme_classic(base_size = 10)+ 
  theme(legend.position = "none")+
  xlab(NULL)+
  theme(aspect.ratio=2)+
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_line(color = "gray", size = 0.25))+
  facet_wrap(~ type,scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 0.95, label.y = 26.2)

tmt<-data.frame(intensity=as.numeric(c(tmt_impexp['SEPTIN7P8 lncORF 1',],tmt_unimpexp['SEPTIN7P8 lncORF 1',])),
                group=c(substr(colnames(tmt_impexp),1,1),substr(colnames(tmt_unimpexp),1,1)),
                type=c(rep('minimum value imputation',ncol(tmt_impexp)),rep('Remove missing values',ncol(tmt_unimpexp))))
tmt<-tmt[!is.na(tmt$intensity),]

p[[2]]<-
  ggplot(tmt, aes(x = group, y = intensity)) +
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
  labs(y="log2 (intensity)")+ 
  theme_classic(base_size = 10)+ 
  scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
  theme_classic(base_size = 10)+ 
  theme(legend.position = "none")+
  xlab(NULL)+
  theme(aspect.ratio=2)+
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_line(color = "gray", size = 0.25))+
  facet_wrap(~ type,scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 0.95, label.y = 18.6) 


pdf("./Figure4/RPKM_SEPTIN7P8_o1_box.pdf",width = 8,height = 8)
wrap_plots(p,nrow=2, guides="collect") 
dev.off()


###############################################################
#Other genes including ：
#recurrence: POTEKP lncORF 2, HSPD1P4 lncORF 1 
#Survival: HNRNPA1P36 lncORF 1

tmt<-data.frame(intensity=as.numeric(c(tmt_impexp['POTEKP lncORF 2',],tmt_unimpexp['POTEKP lncORF 2',])),
                group=c(substr(colnames(tmt_impexp),1,1),substr(colnames(tmt_unimpexp),1,1)),
                sample=c(colnames(tmt_impexp),colnames(tmt_unimpexp)),
                type=c(rep('minimum value imputation',ncol(tmt_impexp)),rep('Remove missing values',ncol(tmt_unimpexp))))
tmt<-tmt[!is.na(tmt$intensity),]

p<-list()

p[[1]]<-
ggplot(tmt, aes(x = group, y = intensity)) +
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
  labs(y="log2 (intensity)")+ 
  theme_classic(base_size = 10)+ 
  scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
  theme_classic(base_size = 10)+ 
  theme(legend.position = "none")+
  xlab(NULL)+
  theme(aspect.ratio=2)+
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_line(color = "gray", size = 0.25))+
  facet_wrap(~ type,scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 0.9, label.y = 23.2) 

meta<-data.table::fread('clinic/captac_meta.csv')

meta<-data.frame(sample=c(substr(meta$`Tumor (T) sample ID`,2,99)),
                 OS=c(meta$`Survial  (1, dead; 0, alive)`),
                 OStime=c(meta$`Overall survial (month)`),
                 rec=c(meta$`Recurrence  (1, yes; 0, no)`),
                 rectime=c(meta$`Recurrence-free survival (month)`))
tmt<-tmt[tmt$group=='T',]
meta<-meta[match(substr(tmt$sample,2,99),meta$sample),]

tmt$rec<-meta$rec
tmt$rectime<-meta$rectime

library(survminer)
library(survival)
library(ggpubr)
library(gridExtra)

matt<-tmt[tmt$type!='minimum value imputation',]
more.med.exp.index<-which(matt$intensity>=quantile(matt$intensity)[4]) 
less.med.exp.index<-which(matt$intensity< quantile(matt$intensity)[2])
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(rectime,rec) ~ status, data = matt)
s.diff<-survdiff(Surv(rectime,rec) ~ status, data = matt)

p[[2]]<-ggsurvplot(s.fit,
                        data=matt,
                        palette="Pastel1",
                        pval = TRUE,
                        pval.method = TRUE,
                        conf.int = TRUE,
                        xlab = 'Time (Month)',ylab='Disease-free probability',
                        ggtheme = theme_survminer(),
                        surv.median.line = 'hv')

p[[2]]<-p[[2]]$plot+ theme(aspect.ratio=1)

pdf("./Figure4/POTEKPORF2_REC_nonimp.pdf",width = 4,height = 8)
wrap_plots(p,nrow=2, guides="keep")  
dev.off()

matt<-tmt[tmt$type=='minimum value imputation',]
more.med.exp.index<-which(matt$intensity>=quantile(matt$intensity)[4]) 
less.med.exp.index<-which(matt$intensity< quantile(matt$intensity)[2])
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(rectime,rec) ~ status, data = matt)
s.diff<-survdiff(Surv(rectime,rec) ~ status, data = matt)

p[[2]]<-ggsurvplot(s.fit,
                   data=matt,
                   palette="Pastel1",
                   pval = TRUE,
                   pval.method = TRUE,
                   conf.int = TRUE,
                   xlab = 'Time (Month)',ylab='Disease-free probability',
                   ggtheme = theme_survminer(),
                   surv.median.line = 'hv')

p[[2]]<-p[[2]]$plot+ theme(aspect.ratio=1)

pdf("./Figure4/POTEKPORF2_REC_minimp.pdf",width = 4,height = 8)
wrap_plots(p,nrow=2, guides="keep") 
dev.off()

tmt<-data.frame(intensity=as.numeric(c(tmt_impexp['HSPD1P4 lncORF 1',],tmt_unimpexp['HSPD1P4 lncORF 1',])),
                group=c(substr(colnames(tmt_impexp),1,1),substr(colnames(tmt_unimpexp),1,1)),
                sample=c(colnames(tmt_impexp),colnames(tmt_unimpexp)),
                type=c(rep('minimum value imputation',ncol(tmt_impexp)),rep('Remove missing values',ncol(tmt_unimpexp))))
tmt<-tmt[!is.na(tmt$intensity),]

p<-list()

p[[1]]<-
  ggplot(tmt, aes(x = group, y = intensity)) +
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
  labs(y="log2 (intensity)")+ 
  theme_classic(base_size = 10)+ 
  scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
  theme_classic(base_size = 10)+ 
  theme(legend.position = "none")+
  xlab(NULL)+
  theme(aspect.ratio=2)+
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_line(color = "gray", size = 0.25))+
  facet_wrap(~ type,scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 0.9, label.y = 15.2) 

meta<-data.table::fread('clinic/captac_meta.csv')

meta<-data.frame(sample=c(substr(meta$`Tumor (T) sample ID`,2,99)),
                 OS=c(meta$`Survial  (1, dead; 0, alive)`),
                 OStime=c(meta$`Overall survial (month)`),
                 rec=c(meta$`Recurrence  (1, yes; 0, no)`),
                 rectime=c(meta$`Recurrence-free survival (month)`))
tmt<-tmt[tmt$group=='T',]
meta<-meta[match(substr(tmt$sample,2,99),meta$sample),]

tmt$rec<-meta$rec
tmt$rectime<-meta$rectime

matt<-tmt[tmt$type!='minimum value imputation',]
more.med.exp.index<-which(matt$intensity>=quantile(matt$intensity)[4]) 
less.med.exp.index<-which(matt$intensity< quantile(matt$intensity)[2])
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(rectime,rec) ~ status, data = matt)
s.diff<-survdiff(Surv(rectime,rec) ~ status, data = matt)

p[[2]]<-ggsurvplot(s.fit,
                   data=matt,
                   palette="Pastel1",
                   pval = TRUE,
                   pval.method = TRUE,
                   conf.int = TRUE,
                   xlab = 'Time (Month)',ylab='Disease-free probability',
                   ggtheme = theme_survminer(),
                   surv.median.line = 'hv')

p[[2]]<-p[[2]]$plot+ theme(aspect.ratio=1)

pdf("./Figure4/HSPD1P4O1_nonimp.pdf",width = 4,height = 8)
wrap_plots(p,nrow=2, guides="keep")  
dev.off()

matt<-tmt[tmt$type=='minimum value imputation',]
more.med.exp.index<-which(matt$intensity>quantile(matt$intensity)[4]) 
less.med.exp.index<-which(matt$intensity< quantile(matt$intensity)[2])
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(rectime,rec) ~ status, data = matt)
s.diff<-survdiff(Surv(rectime,rec) ~ status, data = matt)

p[[2]]<-ggsurvplot(s.fit,
                   data=matt,
                   palette="Pastel1",
                   pval = TRUE,
                   pval.method = TRUE,
                   conf.int = TRUE,
                   xlab = 'Time (Month)',ylab='Disease-free probability',
                   ggtheme = theme_survminer(),
                   surv.median.line = 'hv')

p[[2]]<-p[[2]]$plot+ theme(aspect.ratio=1)

pdf("./Figure4/HSPD1P4O1_minimp.pdf",width = 4,height = 8)
wrap_plots(p,nrow=2, guides="keep")  
dev.off()

#

tmt<-data.frame(intensity=as.numeric(c(tmt_impexp['HNRNPA1P36 lncORF 1',],tmt_unimpexp['HNRNPA1P36 lncORF 1',])),
                group=c(substr(colnames(tmt_impexp),1,1),substr(colnames(tmt_unimpexp),1,1)),
                sample=c(colnames(tmt_impexp),colnames(tmt_unimpexp)),
                type=c(rep('minimum value imputation',ncol(tmt_impexp)),rep('Remove missing values',ncol(tmt_unimpexp))))
tmt<-tmt[!is.na(tmt$intensity),]

p<-list()

p[[1]]<-
  ggplot(tmt, aes(x = group, y = intensity)) +
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,aes(fill=group),show.legend = F)+
  labs(y="log2 (intensity)")+ 
  theme_classic(base_size = 10)+ 
  scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
  theme_classic(base_size = 10)+ 
  theme(legend.position = "none")+
  xlab(NULL)+
  theme(aspect.ratio=2)+
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_line(color = "gray", size = 0.25))+
  facet_wrap(~ type,scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 0.9, label.y = 16.8) 

meta<-data.table::fread('clinic/captac_meta.csv')

meta<-data.frame(sample=c(substr(meta$`Tumor (T) sample ID`,2,99)),
                 OS=c(meta$`Survial  (1, dead; 0, alive)`),
                 OStime=c(meta$`Overall survial (month)`),
                 rec=c(meta$`Recurrence  (1, yes; 0, no)`),
                 rectime=c(meta$`Recurrence-free survival (month)`))
tmt<-tmt[tmt$group=='T',]
meta<-meta[match(substr(tmt$sample,2,99),meta$sample),]

tmt$OS<-meta$OS
tmt$OStime<-meta$OStime

matt<-tmt[tmt$type=='minimum value imputation',]
more.med.exp.index<-which(matt$intensity>=quantile(matt$intensity)[4])
less.med.exp.index<-which(matt$intensity< quantile(matt$intensity)[2])
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OStime,OS) ~ status, data = matt)
s.diff<-survdiff(Surv(OStime,OS) ~ status, data = matt)

p[[2]]<-ggsurvplot(s.fit,
                   data=matt,
                   palette="Pastel1",
                   pval = TRUE,
                   pval.method = TRUE,
                   conf.int = TRUE,
                   xlab = 'Time (Month)',ylab='Survival probability',
                   ggtheme = theme_survminer(),
                   surv.median.line = 'hv')

p[[2]]<-p[[2]]$plot+ theme(aspect.ratio=1)

pdf("./Figure4/HNRNPA1P36O1_SUV_minimp.pdf",width = 4,height = 8)
wrap_plots(p,nrow=2, guides="keep")  
dev.off()

matt<-tmt[tmt$type!='minimum value imputation',]
more.med.exp.index<-which(matt$intensity>=quantile(matt$intensity)[4]) 
less.med.exp.index<-which(matt$intensity< quantile(matt$intensity)[2])
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OStime,OS) ~ status, data = matt)
s.diff<-survdiff(Surv(OStime,OS) ~ status, data = matt)

p[[2]]<-ggsurvplot(s.fit,
                   data=matt,
                   palette="Pastel1",
                   pval = TRUE,
                   pval.method = TRUE,
                   conf.int = TRUE,
                   xlab = 'Time (Month)',ylab='Survival probability',
                   ggtheme = theme_survminer(),
                   surv.median.line = 'hv')

p[[2]]<-p[[2]]$plot+ theme(aspect.ratio=1)

pdf("./Figure4/HNRNPA1P36O1_SUV_nonimp.pdf",width = 4,height = 8)
wrap_plots(p,nrow=2, guides="keep")  
dev.off()


###############################################################
#forest plot

p<-list()

LBF_TvsP<-rbind(LBF_unimpu[[1]],LBF_minimpu[[1]])
LBF_TvsP$group<-c(rep('Remove missing values',nrow(LBF_unimpu[[1]])),rep('minimum value imputation',nrow(LBF_minimpu[[1]])))

p[[1]]<-
  ggplot(LBF_TvsP, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,6))+
  theme(aspect.ratio=length(unique(LBF_TvsP$lncPep))/6)+
  theme(legend.position = "bottom")

LBF_SUV<-rbind(LBF_unimpu[[2]],LBF_minimpu[[2]])
LBF_SUV$group<-c(rep('Remove missing values',nrow(LBF_unimpu[[2]])),rep('minimum value imputation',nrow(LBF_minimpu[[2]])))
LBF_SUV<-LBF_SUV[ !(is.na(LBF_SUV$OR_2)| LBF_SUV$OR_2=='Inf'),]

p[[2]]<-
  ggplot(LBF_SUV, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Hazard Ratio per SD",
       title= paste0("Cox regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,6))+
  theme(aspect.ratio=length(unique(LBF_SUV$lncPep))/6)+
  theme(legend.position = "bottom")

LBF_REC<-rbind(LBF_unimpu[[3]],LBF_minimpu[[3]])
LBF_REC$group<-c(rep('Remove missing values',nrow(LBF_unimpu[[3]])),rep('minimum value imputation',nrow(LBF_minimpu[[3]])))
LBF_REC<-LBF_REC[ !(is.na(LBF_REC$OR_2)| LBF_REC$OR_2=='Inf'),]

p[[3]]<-
  ggplot(LBF_REC, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Hazard Ratio per SD",
       title= paste0("Cox regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,6))+
  theme(aspect.ratio=length(unique(LBF_REC$lncPep))/6)+
  theme(legend.position = "bottom")

library(patchwork)

pdf("./Figure4/LBF_test.pdf",width = 4,height = 12)
p[[1]]/p[[2]]/p[[3]]+plot_layout(heights = c(6,5,6))
dev.off()

#

p<-list()

TMT_TvsP<-rbind(TMT_nonimp[[1]],TMT_minimp[[1]])
TMT_TvsP$group<-c(rep('Remove missing values',nrow(TMT_nonimp[[1]])),rep('minimum value imputation',nrow(TMT_minimp[[1]])))

p[[1]]<-
  ggplot(TMT_TvsP, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,6))+
  theme(aspect.ratio=length(unique(TMT_TvsP$lncPep))/6)+
  theme(legend.position = "bottom")

TMT_SUV<-rbind(TMT_nonimp[[2]],TMT_minimp[[2]])
TMT_SUV$group<-c(rep('Remove missing values',nrow(TMT_nonimp[[2]])),rep('minimum value imputation',nrow(TMT_minimp[[2]])))
TMT_SUV<-TMT_SUV[ !(is.na(TMT_SUV$OR_2)| TMT_SUV$OR_2=='Inf'),]

p[[2]]<-
  ggplot(TMT_SUV, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Hazard Ratio per SD",
       title= paste0("Cox regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,6))+
  theme(aspect.ratio=length(unique(TMT_SUV$lncPep))/6)+
  theme(legend.position = "bottom")

TMT_REC<-rbind(TMT_nonimp[[3]],TMT_minimp[[3]])
TMT_REC$group<-c(rep('Remove missing values',nrow(TMT_nonimp[[3]])),rep('minimum value imputation',nrow(TMT_minimp[[3]])))
TMT_REC<-TMT_REC[ !(is.na(TMT_REC$OR_2)| TMT_REC$OR_2=='Inf'),]

p[[3]]<-
  ggplot(TMT_REC, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Hazard Ratio per SD",
       title= paste0("Cox regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,6))+
  theme(aspect.ratio=length(unique(TMT_REC$lncPep))/6)+
  theme(legend.position = "bottom")

pdf("./Figure4/TMT_test.pdf",width = 4,height = 40)
p[[1]]/p[[2]]/p[[3]]+plot_layout(heights = c(length(unique(TMT_TvsP$lncPep)),
                                             length(unique(TMT_SUV$lncPep)),
                                             length(unique(TMT_REC$lncPep))))
dev.off()



TvsP<-rbind(LBF[[1]],TMT[[1]])
TvsP$group<-c(rep('Jiang et al 2019',6),rep('Gao et al 2019',81))
TvsP<-TvsP[!is.na(TvsP$OR_2),]

pdf("./Figure4/T_vs_P.pdf",width = 8,height = 16)
ggplot(TvsP, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  theme(legend.position = "none")+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,5))+
  theme(aspect.ratio=2.5)
dev.off()

retain<-table(TvsP$lncPep)
retain<-retain[retain==2]
retain2<-TvsP$lncPep[TvsP$Pvalue<0.05]
retain<-c(names(retain),retain2)

retain<-TvsP[TvsP$lncPep %in% retain,]

p<-list()

p[[1]]<-
# pdf("./Figure4/T_vs_P_filter.pdf",width = 8,height = 8)
ggplot(retain, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  theme(legend.position = "none")+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,5))
# dev.off()


SUV<-rbind(LBF[[2]],TMT[[2]])
SUV$group<-c(rep('Jiang et al 2019',6),rep('Gao et al 2019',65))
SUV<-SUV[!is.na(SUV$OR_2),]

# pdf("./Figure4/SUV.pdf",width = 8,height = 16)
ggplot(SUV, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  theme(legend.position = "none")+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Hazard Ratio per SD",
       title= paste0("Cox regression (95% CI)") )+
  coord_cartesian(xlim=c(0.4,3.5))+
  theme(aspect.ratio=2.5)
# dev.off()

retain<-table(SUV$lncPep)
retain<-retain[retain==2]
retain2<-SUV$lncPep[SUV$Pvalue<0.05]
retain<-c(names(retain),retain2)
retain<-SUV[SUV$lncPep %in% retain,]

p[[2]]<-
# pdf("./Figure4/SUV_filter.pdf",width = 6,height = 6)
ggplot(retain, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  theme(legend.position = "none")+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,5))
# dev.off()

REC<-rbind(LBF[[3]],TMT[[3]])
REC$group<-c(rep('Jiang et al 2019',6),rep('Gao et al 2019',65))
REC<-REC[!is.na(REC$OR_2),]

pdf("./Figure4/REC.pdf",width = 8,height = 16)
ggplot(REC, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  theme(legend.position = "none")+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Hazard Ratio per SD",
       title= paste0("Cox regression (95% CI)") )+
  coord_cartesian(xlim=c(0.4,3.5))+
  theme(aspect.ratio=2.5)
dev.off()

retain<-table(REC$lncPep)
retain<-retain[retain==2]
retain2<-REC$lncPep[REC$Pvalue<0.05]
retain<-c(names(retain),retain2)
retain<-REC[REC$lncPep %in% retain,]

p[[3]]<-
# pdf("./Figure4/REC_filter.pdf",width = 6,height = 6)
ggplot(retain, aes(y=reorder(lncPep,OR), x=OR,  colour=group)) + 
  geom_errorbarh(aes(xmin=OR_1, xmax=OR_2), height=.3) +
  geom_point(size=2)+
  theme(legend.position = "none")+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds ratio per SD",
       title= paste0("Logistic regression (95% CI)") )+
  coord_cartesian(xlim=c(0.1,5))
# dev.off()

#########################################
library(patchwork)


combined_plot <- p[[1]] + p[[2]] + p[[3]] + plot_layout(ncol = 1, heights = c(14, 4, 3))

pdf("./Figure4/forest_test.pdf",width = 10,height = 20)
combined_plot
dev.off() 

#########################################

#venn
or_sig<-TvsP$lncPep[TvsP$Pvalue<0.05]
suv_sig<-SUV$lncPep[SUV$Pvalue<0.05]
rec_sig<-REC$lncPep[REC$Pvalue<0.05]

library(ggvenn)
library(MetBrewer)

venn<-list('Cancer-related'=unique(or_sig),
           'Survival-related'=unique(suv_sig),
           'Recurrence-related'=unique(rec_sig))

pdf("./Figure4/VENN.pdf",width = 4,height = 4)
ggvenn(venn,fill_color = met.brewer('Hiroshige',n=3),
               stroke_size = 0.5, set_name_size = 4)
dev.off()

suv_sig[suv_sig %in% or_sig]

match$sequence[match$ORFname=='ACTBP4 lncORF 1']
match$sequence[match$ORFname=='CES1P1 lncORF 6']
match$sequence[match$ORFname=='PPIAP79 lncORF 1']



