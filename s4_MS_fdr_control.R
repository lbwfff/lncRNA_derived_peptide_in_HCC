source('source.R')

###############################################
#Compare the number of mass spectrometry spectra obtained by MS under different methods.

stat<-function(path){

list<-list.files(path)
list<-list[grep('psms.tsv',list)]

f <-function(x) unlist(strsplit(x,'_percolator_'))[1]
list<-sapply(list,f)
list<-as.character(list)
list<-unique(list)

cache<-data.frame()

for (i in list){
  
  target<-paste0(path,i,'_percolator_target_psms.tsv')
  decoy<-paste0(path,i,'_percolator_decoy_psms.tsv')
  target<-data.table::fread(target)
  decoy<-data.table::fread(decoy)
  
  target<-target[grep('tr',target$proteinIds),]
  target<-target[grep('HUMAN',target$proteinIds,invert=T),]
  decoy<-decoy[grep('tr',decoy$proteinIds),]
  decoy<-decoy[grep('HUMAN',decoy$proteinIds,invert=T),]

  globe<-data.frame(score=c(target$score,decoy$score),
                    group=c(rep('target',length(target$score)),rep('decoy',length(decoy$score))))
  
  hold<-get_threshold(globe)
  
  poscontrol<-target[target$score>hold,]
  
  stat<-data.frame(num=c(nrow(poscontrol)),sample=c(i))
  cache<-rbind(cache,stat)
}

return(cache)
}

group<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_group_index1/01/')
group$group<-c('Group_FDR_index1')

group2<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_group_index2/01/')
group2$group<-c('Group_FDR_index2')

group3<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_group_index3/01/')
group3$group<-c('Group_FDR_index3')

twostep<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_twopass2_index1/01/')
twostep$group<-c('Two_step_index1')

twostep2<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_twopass2_index2/01/')
twostep2$group<-c('Two_step_index2') 

twostep3<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_twopass2_index3/01/')
twostep3$group<-c('Two_step_index3')


plot<-rbind(
  group,group2,group3,twostep,twostep2,twostep3
)

library('ggplot2')
library('ggsignif')
library('MetBrewer')

compaired <- list(c("Two_step_index2","Group_FDR_index2"),
                  c("Two_step_index1","Two_step_index2"),
                  c("Group_FDR_index2","Group_FDR_index1"))

pdf("TMT_SC.pdf",width = 6,height = 6)
ggplot(data=plot,aes_string(x='group',y='num'))+   
  geom_boxplot(aes(color=group))+
  geom_signif(comparisons = compaired,step_increase = 0.12,vjust =0.2,
              map_signif_level = T,test = wilcox.test)+
  scale_color_manual(values=met.brewer("Egypt",n=6,type="continuous"))+
  scale_fill_manual(values=met.brewer("Egypt",n=6,type="continuous"))+
  theme_bw() + 
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab('Spectra Count')+
  xlab(NULL)+
  theme(aspect.ratio=1)
dev.off()

######

########################################


stat<-function(path){
  
  list<-list.files(path,full.names=T,all.files=T,recursive=T)
  list<-list[grep('psms.tsv',list)]
  
  f <-function(x) unlist(strsplit(x,'_percolator_'))[1]
  list<-sapply(list,f)
  list<-as.character(list)
  list<-unique(list)
  
  cache<-data.frame()
  
  for (i in list){
    
    target<-paste0(i,'_percolator_target_psms.tsv')
    decoy<-paste0(i,'_percolator_decoy_psms.tsv')
    target<-data.table::fread(target)
    decoy<-data.table::fread(decoy)
    
    target<-target[grep('tr',target$proteinIds),]
    target<-target[grep('HUMAN',target$proteinIds,invert=T),]
    decoy<-decoy[grep('tr',decoy$proteinIds),]
    decoy<-decoy[grep('HUMAN',decoy$proteinIds,invert=T),]
    
    globe<-data.frame(score=c(target$score,decoy$score),
                      group=c(rep('target',length(target$score)),rep('decoy',length(decoy$score))))
    
    hold<-get_threshold(globe)
    
    poscontrol<-target[target$score>hold,]
    
    stat<-data.frame(num=c(nrow(poscontrol)),sample=c(i))
    cache<-rbind(cache,stat)
  }
  
  return(cache)
}

group<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_group_index1/')
group$group<-c('Group_FDR_index1')

group2<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_group_index2/')
group2$group<-c('Group_FDR_index2')

group3<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_group_index3/')
group3$group<-c('Group_FDR_index3')

group4<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_twopass2_index1/')
group4$group<-c('Two_step_index1')

group5<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_twopass2_index2/')
group5$group<-c('Two_step_index2')

group6<-stat('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_twopass2_index3/')
group6$group<-c('Two_step_index3')


plot<-rbind(
  group,group2,group3,group4,group5,group6
)

library('ggplot2')
library('ggsignif')
library('MetBrewer')

compaired <- list(c("Two_step_index2","Group_FDR_index2"),
                  c("Two_step_index3","Two_step_index2"),
                  c("Group_FDR_index3","Group_FDR_index2"))

pdf("LBF_SC.pdf",width = 6,height = 6)
ggplot(data=plot,aes_string(x='group',y='num'))+   
  geom_boxplot(aes(color=group))+
  geom_signif(comparisons = compaired,step_increase = 0.12,vjust =0.2,
              map_signif_level = T,test = wilcox.test)+
  scale_color_manual(values=met.brewer("Egypt",n=6,type="continuous"))+
  scale_fill_manual(values=met.brewer("Egypt",n=6,type="continuous"))+
  theme_bw() + 
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab('Spectra Count')+
  xlab(NULL)+
  theme(aspect.ratio=1)
dev.off()

##########################################
#Cumulative Histogram

# setwd('/scratch/lb4489/project/hcc/R_code/')

group1<-data.table::fread('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_group_index1/01/peptide.tsv')
group1$group<-c('Group_FDR_index1')
group2<-data.table::fread('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_group_index2/01/peptide.tsv')
group2$group<-c('Group_FDR_index2')
group3<-data.table::fread('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_group_index3/01/peptide.tsv')
group3$group<-c('Group_FDR_index3')

group4<-data.table::fread('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_twopass2_index1/01/peptide.tsv')
group4$group<-c('Two_step_index1')
group5<-data.table::fread('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_twopass2_index2/01/peptide.tsv')
group5$group<-c('Two_step_index2')
group6<-data.table::fread('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/tmt_twopass2_index3/01/peptide.tsv')
group6$group<-c('Two_step_index3')

plot<-rbind(
  group1,
  group2,
  group3,
  group4,
  group5,group6
)

plot<-plot[grep('tr',plot$Protein),]
plot<-plot[grep('HUMAN',plot$`Mapped Proteins`,invert=T),]

mapping<-data.frame(seq=unique(plot$Peptide),
                    type=c(NA))

pattern2 <- Biostrings::readAAStringSet('/scratch/lb4489/bioindex/uniprot_human_nonolnc.fasta')
mapping$match<-NA
mapping$snpmatch<-NA

library(Biostrings)

for (i in 1:nrow(mapping)){
  
  pattern1 <- Biostrings::AAString(mapping$seq[i])
  
  match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=0)
  match2<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)
  
  
  mapping$match[i]<-sum(match)
  mapping$snpmatch[i]<-sum(match2) 
  
}

mapping$type<-ifelse(mapping$match==0 & mapping$snpmatch==0,'Unique',
                     ifelse(mapping$match==0 & mapping$snpmatch>0,'SAP','Miss-cleavage'))

mapping<-mapping[match(plot$Peptide,mapping$seq),]
plot$type<-mapping$type

library('ggplot2')
library('ggsignif')
library('MetBrewer')

plot$type<-factor(plot$type,levels = c('Unique','SAP','Miss-cleavage'))

pdf("TMT_peptide_bar.pdf",width = 8,height = 6)

ggplot(data=plot,aes(x=group,fill = type))+   
  geom_bar( position = "stack") +
  scale_color_manual(values=met.brewer("Hokusai3",n=3,type="continuous"))+
  scale_fill_manual(values=met.brewer("Hokusai3",n=3,type="continuous"))+
  theme_bw() + 
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab('Peptide Count')+
  xlab(NULL)+
  theme(aspect.ratio=1)

dev.off()


####

mergelbfpep<-function(path){
  dir<-list.dirs(path)[-1]
  dir<-dir[-grep('MSBooster',dir)]
  
  cache<-data.frame()
  for (i in dir){
    peptide<-data.table::fread(paste0(i,'/peptide.tsv'))
    cache<-rbind(cache,peptide)
  }
  return(cache)
}

group1<-mergelbfpep('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_group_index1/')
group1$group<-c('Group_FDR_index1')
group2<-mergelbfpep('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_group_index2/')
group2$group<-c('Group_FDR_index2')
group3<-mergelbfpep('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_group_index3/')
group3$group<-c('Group_FDR_index3')

group4<-mergelbfpep('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_twopass2_index1/')
group4$group<-c('Two_step_index1')
group5<-mergelbfpep('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_twopass2_index2/')
group5$group<-c('Two_step_index2')
group6<-mergelbfpep('/scratch/lb4489/project/hcc/iprox/fragpipe/newfig2/lbf_twopass2_index3/')
group6$group<-c('Two_step_index3')

plot<-rbind(
  group1,
  group2,
  group3,
  group4,group5,group6
)

plot<-plot[grep('tr',plot$Protein),]
plot<-plot[grep('HUMAN',plot$`Mapped Proteins`,invert=T),]

mapping<-data.frame(seq=unique(plot$Peptide),
                    match=c(NA),
                    snpmatch=c(NA))

for (i in 1:nrow(mapping)){
  
  pattern1 <- Biostrings::AAString(mapping$seq[i])
  
  match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=0)
  match2<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)
  
  
  mapping$match[i]<-sum(match)
  mapping$snpmatch[i]<-sum(match2) 
  
}

mapping$type<-ifelse(mapping$match==0 & mapping$snpmatch==0,'Unique',
                     ifelse(mapping$match==0 & mapping$snpmatch>0,'SAP','Miss-cleavage'))

mapping<-mapping[match(plot$Peptide,mapping$seq),]
plot$type<-mapping$type

plot$type<-factor(plot$type,levels = c('Unique','SAP','Miss-cleavage'))

pdf("LBF_peptide_bar.pdf",width = 8,height = 6)

ggplot(data=plot,aes(x=group,fill = type))+   
  geom_bar( position = "stack") +
  scale_color_manual(values=met.brewer("Hokusai3",n=3,type="continuous"))+
  scale_fill_manual(values=met.brewer("Hokusai3",n=3,type="continuous"))+
  theme_bw() + 
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab('Peptide Count')+
  xlab(NULL)+
  theme(aspect.ratio=1)

dev.off()


#########################################
#Representative FDR control chart

setwd('/scratch/lb4489/project/hcc/R_code/')

p1<-plot_fdr('../captac/fragpipe/groupfdr_2/01/20171001_HF_ZHW_total_liver01_F01_percolator_target_psms.tsv',
         '../captac/fragpipe/groupfdr_2/01/20171001_HF_ZHW_total_liver01_F01_percolator_decoy_psms.tsv','sp')

p2<-plot_fdr('../captac/fragpipe/twostep/sec_index2/01/20171001_HF_ZHW_total_liver01_F01_sub_percolator_target_psms.tsv',
             '../captac/fragpipe/twostep/sec_index2/01/20171001_HF_ZHW_total_liver01_F01_sub_percolator_decoy_psms.tsv','sp')


p3<-plot_fdr('../iprox/fragpipe/onestep_2/L001_P/CNHPP_HCC_LC_profiling_L001_P_F1_percolator_target_psms.tsv',
             '../iprox/fragpipe/onestep_2/L001_P/CNHPP_HCC_LC_profiling_L001_P_F1_percolator_decoy_psms.tsv','sp')

p4<-plot_fdr('../iprox/fragpipe/twostep/sec_index2/L001_P/CNHPP_HCC_LC_profiling_L001_P_F1_sub_percolator_target_psms.tsv',
             '../iprox/fragpipe/twostep/sec_index2/L001_P/CNHPP_HCC_LC_profiling_L001_P_F1_sub_percolator_decoy_psms.tsv','sp')

plot<-list(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],
           p2[[1]],p2[[2]],p2[[3]],p2[[4]],p2[[5]],p2[[6]],
           p3[[1]],p3[[2]],p3[[3]],p3[[4]],p3[[5]],p3[[6]],
           p4[[1]],p4[[2]],p4[[3]],p4[[4]],p4[[5]],p4[[6]])


library(patchwork)
pdf("FDR_plot.pdf",width = 12,height = 24)
wrap_plots(plot,nrow=8, guides="collect") 
dev.off()


wrap_plots(p,nrow=2, guides="collect") 

