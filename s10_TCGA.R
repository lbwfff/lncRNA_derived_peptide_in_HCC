##############################

library(data.table)

meta<-fread('TCGA/TCGA.LIHC.sampleMap_LIHC_clinicalMatrix')
tcgasuv<-fread('TCGA/survival_LIHC_survival.txt')

# tcgaexp<-fread('TCGA/TCGA.LIHC.sampleMap_HiSeqV2_percentile.gz')
load('TCGA/lihcexp.RData')
tcgaexp<-lihcexp
# tcgaexp<-tcgaexp[!is.na(tcgaexp$`TCGA-FV-A495-01`),]
tcgaexp<-as.data.frame(tcgaexp)
rownames(tcgaexp)<-tcgaexp$sample
tcgaexp<-tcgaexp[,-1]

# gene_anno<-data.frame(rtracklayer::import.gff('../biodata/gencode.v46.annotation.gtf',format = 'gtf'))
# gene_anno<-gene_anno[gene_anno$type=='transcript',]
# geneid<-gene_anno[gene_anno$gene_name %in% sankytable$genename,]

table(substr(rownames(tcgaexp),1,15) %in% substr(geneid$gene_id,1,15))

lncexp<-tcgaexp[(substr(rownames(tcgaexp),1,15) %in% substr(geneid$gene_id,1,15)),]
geneid<-geneid[match(substr(rownames(lncexp),1,15),substr(geneid$gene_id,1,15)),]
rownames(lncexp)<-geneid$gene_name


#######################################################

library("survival")
library("survminer")

output<-data.frame()

for (i in 1:nrow(lncexp)){
  
  #逻辑回归
  exp<-as.data.frame(t(lncexp[i,]))
  exp$group<-substr(rownames(exp),14,15)
  exp$tumor<-ifelse(exp$group==11,'P','T')
  exp$forlr<-ifelse(exp$tumor=='P',0,1)
  
  colnames(exp)[1]<-c('TPM')
  if ( nrow(exp[exp$TPM>(-9),]) > 20  ) {
  
  forlr<-exp[,c(1,4)]
  fit.full<-glm(forlr~.,
                data=forlr,family = binomial())
  fit.result<-summary(fit.full)
  
  coefs <- coef(fit.full)
  confint_values <- confint(fit.full) 
  
  result <-as.data.frame(exp(cbind(OR = coef(fit.full), confint(fit.full))))
  result$p<-as.numeric(fit.result$coefficients[,4])
  result<-data.frame(result[-1,c(1,4,2,3)]) 
  colnames(result)<-c("OR","Pvalue","OR_1","OR_2")
  
  #Cox
  exp<-exp[exp$tumor=='T',]
  suvinf<-tcgasuv[match(rownames(exp),tcgasuv$sample),]
  
  exp<-cbind(exp,suvinf[,c(3:10)])
  
  expos<-exp[!is.na(exp$OS),]
  expdss<-exp[!is.na(exp$DSS),]
  expdfi<-exp[!is.na(exp$DFI),]
  exppfi<-exp[!is.na(exp$PFI),]
  
 
  os<-coxph(Surv(OS.time,OS)~TPM, data=expos) 
  dss<-coxph(Surv(DSS.time,DSS)~TPM, data=expos) 
  dfi<-coxph(Surv(DFI.time,DFI)~TPM, data=expos) 
  pfi<-coxph(Surv(PFI.time,PFI)~TPM, data=expos) 
  
  result2<-data.frame(OSHR=c(summary(os)$coef[,'exp(coef)']),
                      OSPvalue=c(summary(os)$coef[, "Pr(>|z|)"]),
                      OSHR_1=c(summary(os)$conf.int[, "lower .95"]),
                      OSHR_2=c(summary(os)$conf.int[, "upper .95"]),
                      
                      DSSHR=c(summary(dss)$coef[,'exp(coef)']),
                      DSSPvalue=c(summary(dss)$coef[, "Pr(>|z|)"]),
                      DSSHR_1=c(summary(dss)$conf.int[, "lower .95"]),
                      DSSHR_2=c(summary(dss)$conf.int[, "upper .95"]),
                      
                      DFIHR=c(summary(dfi)$coef[,'exp(coef)']),
                      DFIPvalue=c(summary(dfi)$coef[, "Pr(>|z|)"]),
                      DFIHR_1=c(summary(dfi)$conf.int[, "lower .95"]),
                      DFIHR_2=c(summary(dfi)$conf.int[, "upper .95"]),
                      
                      PFIHR=c(summary(pfi)$coef[,'exp(coef)']),
                      PFIPvalue=c(summary(pfi)$coef[, "Pr(>|z|)"]),
                      PFIHR_1=c(summary(pfi)$conf.int[, "lower .95"]),
                      PFIHR_2=c(summary(pfi)$conf.int[, "upper .95"]),
                      
                      lncRNA=c(rownames(lncexp)[i]))
  
  all<-cbind(result,result2)
  
  output<-rbind(output,all) } else{ print('pass') }
}

#######################################################
#heatmap

mspepgenelist<-rownames(heatarray)
mspepgenelist<-sub(" .*", "", mspepgenelist)

tcgamatch<-output[match(mspepgenelist,output$lncRNA),]

heatarray_sup<-cbind(heatarray,tcgamatch[,c(1,5,9,13,17)])
heatarrayp_sup<-cbind(heatarray_p,tcgamatch[,c(2,6,10,14,18)])

pdf('./Figure4/heatmap_sup.pdf',width = 16,height = 15)
draw(
  Heatmap(heatarray_sup,na_col = "grey",col = col_fun1,
          name = "OR/HR value",
          column_title = NULL,
          cluster_columns = F,
          cluster_rows =F,
          show_row_dend = T,
          show_column_names =T,
          show_heatmap_legend = T,
          row_names_side = "left",
          width = ncol(heatarray_sup)*unit(4, "mm"), 
          height = nrow(heatarray_sup)*unit(4, "mm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(heatarrayp_sup[i, j] < 0.05 & !is.na(heatarrayp_sup[i, j]))
              grid.text('*', x, y, gp = gpar(fontsize = 10))
          }
  )
)
dev.off() 

###############
#Visualize the two prognostic peptides

{
  exp<-as.data.frame(t(lncexp[which(rownames(lncexp)=='HNRNPA1P36'),]))
  exp$group<-substr(rownames(exp),14,15)
  exp$tumor<-ifelse(exp$group==11,'P','T')
  exp$forlr<-ifelse(exp$tumor=='P',0,1)
  colnames(exp)[1]<-'TPM'
  
  p<-list()
  
  p[[1]]<-
    ggplot(exp, aes(x = tumor, y = TPM)) +
    stat_boxplot(geom = "errorbar",width=0.2)+
    geom_boxplot(outlier.shape = 1,aes(fill=tumor),show.legend = F)+
    labs(y="log2 (TPM)")+ 
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
    stat_compare_means(method = "wilcox.test", label.x = 0.9, label.y = -1) 
  
  exp<-exp[exp$tumor=='T',]
  suvinf<-tcgasuv[match(rownames(exp),tcgasuv$sample),]
  
  exp<-cbind(exp,suvinf[,c(3:10)])
  
  expos<-exp[!is.na(exp$OS),]
  expdss<-exp[!is.na(exp$DSS),]
  expdfi<-exp[!is.na(exp$DFI),]
  exppfi<-exp[!is.na(exp$PFI),]
  
  library(survminer)
  library(survival)
  library(ggpubr)
  library(gridExtra)
  
  matt<-exppfi
  more.med.exp.index<-which(matt$TPM>=quantile(matt$TPM)[4]) 
  less.med.exp.index<-which(matt$TPM<= quantile(matt$TPM)[2])
  matt$status<-NA
  matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
  matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')
  
  s.fit<-survfit(Surv(PFI.time,PFI) ~ status, data = matt)
  s.diff<-survdiff(Surv(PFI.time,PFI) ~ status, data = matt)
  
  p[[2]]<-ggsurvplot(s.fit,
                     data=matt,
                     palette="Pastel1",
                     pval = TRUE,
                     pval.method = TRUE,
                     conf.int = TRUE,
                     xlab = 'Time (Month)',ylab='Progression-free probability',
                     ggtheme = theme_survminer(),
                     surv.median.line = 'hv')
  
  p[[2]]<-p[[2]]$plot+ theme(aspect.ratio=1)
  
  # pdf("./Figure4/HNRNPA1P36_TCGA.pdf",width = 6,height = 4)
  # wrap_plots(p,nrow=1, guides="keep")  
  # dev.off()
  # 
  
  exp<-as.data.frame(t(lncexp[which(rownames(lncexp)=='POTEKP'),]))
  exp$group<-substr(rownames(exp),14,15)
  exp$tumor<-ifelse(exp$group==11,'P','T')
  exp$forlr<-ifelse(exp$tumor=='P',0,1)
  colnames(exp)[1]<-'TPM'
  
  # p<-list()
  
  p[[3]]<-
    ggplot(exp, aes(x = tumor, y = TPM)) +
    stat_boxplot(geom = "errorbar",width=0.2)+
    geom_boxplot(outlier.shape = 1,aes(fill=tumor),show.legend = F)+
    labs(y="log2 (TPM)")+ 
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
    stat_compare_means(method = "wilcox.test", label.x = 0.9, label.y = -1) 
  
  exp<-exp[exp$tumor=='T',]
  suvinf<-tcgasuv[match(rownames(exp),tcgasuv$sample),]
  
  exp<-cbind(exp,suvinf[,c(3:10)])
  
  expos<-exp[!is.na(exp$OS),]
  expdss<-exp[!is.na(exp$DSS),]
  expdfi<-exp[!is.na(exp$DFI),]
  exppfi<-exp[!is.na(exp$PFI),]
  
  library(survminer)
  library(survival)
  library(ggpubr)
  library(gridExtra)
  
  matt<-expos
  more.med.exp.index<-which(matt$TPM>=quantile(matt$TPM)[4]) 
  less.med.exp.index<-which(matt$TPM<= quantile(matt$TPM)[2])
  matt$status<-NA
  matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
  matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')
  
  s.fit<-survfit(Surv(OS.time,OS) ~ status, data = matt)
  s.diff<-survdiff(Surv(OS.time,OS) ~ status, data = matt)
  
  p[[4]]<-ggsurvplot(s.fit,
                     data=matt,
                     palette="Pastel1",
                     pval = TRUE,
                     pval.method = TRUE,
                     conf.int = TRUE,
                     xlab = 'Time (Month)',ylab='Survival probability',
                     ggtheme = theme_survminer(),
                     surv.median.line = 'hv')
  
  p[[4]]<-p[[4]]$plot+ theme(aspect.ratio=1)
  
  matt<-expdss
  more.med.exp.index<-which(matt$TPM>=quantile(matt$TPM)[4]) 
  less.med.exp.index<-which(matt$TPM<= quantile(matt$TPM)[2])
  matt$status<-NA
  matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
  matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')
  
  s.fit<-survfit(Surv(DSS.time,DSS) ~ status, data = matt)
  s.diff<-survdiff(Surv(DSS.time,DSS) ~ status, data = matt)
  
  p[[5]]<-ggsurvplot(s.fit,
                     data=matt,
                     palette="Pastel1",
                     pval = TRUE,
                     pval.method = TRUE,
                     conf.int = TRUE,
                     xlab = 'Time (Month)',ylab='Disease-Specific survival probability',
                     ggtheme = theme_survminer(),
                     surv.median.line = 'hv')
  
  p[[5]]<-p[[5]]$plot+ theme(aspect.ratio=1)
  
  
  matt<-exppfi
  more.med.exp.index<-which(matt$TPM>=quantile(matt$TPM)[4]) 
  less.med.exp.index<-which(matt$TPM<= quantile(matt$TPM)[2])
  matt$status<-NA
  matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
  matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')
  
  s.fit<-survfit(Surv(PFI.time,PFI) ~ status, data = matt)
  s.diff<-survdiff(Surv(PFI.time,PFI) ~ status, data = matt)
  
  p[[6]]<-ggsurvplot(s.fit,
                     data=matt,
                     palette="Pastel1",
                     pval = TRUE,
                     pval.method = TRUE,
                     conf.int = TRUE,
                     xlab = 'Time (Month)',ylab='Progression-free probability',
                     ggtheme = theme_survminer(),
                     surv.median.line = 'hv')
  
  p[[6]]<-p[[6]]$plot+ theme(aspect.ratio=1)
  
  library(patchwork)
  pdf("./Figure4/TCGA_suv.pdf",width = 30,height = 5)
  wrap_plots(p,nrow=1, guides="keep")  
  dev.off()
  
  
}

# OS（Overall Survival）
# DSS（Disease-Specific Survival）
# DFI（Disease-Free Interval）
# PFI（Progression-Free Interval）






