
######################
#Comparison of CDS and lncORF
plot<-rbind(codeall[codeall$ORF_type=='annotated',],
            result_lnc)
plot$group<-ifelse(plot$ORF_type=='annotated','CDS','lncRNA-ORF')

library('ggplot2')
p<-list()

p[[1]]<-
ggplot(plot, aes(-log10(pval_combined), after_stat(density), colour = group)) +
  geom_density(adjust = 2, size=1)+     
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
  coord_cartesian(xlim=c(0,100))+
  xlab('Combinde Pvalue')

plot$tricer_phase_score<-as.numeric(plot$tricer_phase_score)
p[[2]]<-
ggplot(plot[!is.na(plot$tricer_phase_score)], aes(tricer_phase_score, after_stat(density), colour = group)) +
  geom_density(adjust = 2, size=1)+     
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))  +
  xlab('Phase score')

plot$read_density<-as.numeric(plot$read_density)
p[[3]]<-
ggplot(plot[!is.na(plot$read_density)], aes(log2(read_density), after_stat(density), colour = group)) +
  geom_density(adjust = 2, size=1)+     
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
  xlab('Read density')



library(patchwork)
pdf("./demo/Figure1.pdf",width = 15,height = 5)
wrap_plots(p,nrow=1, guides="collect") 
dev.off()

#################################################
#有点想把这部分重写了，不然数量关系会特别乱
cpc<-data.table::fread('./demo/ribocode_all/output.txt.txt')
gensouce<-cpc[grep('ENST',cpc$`#ID`),]
nonsouce<-cpc[-grep('ENST',cpc$`#ID`),]

{
  tishlnc<-rbind(tish[grep('pseudo|lncRNA',tish$GeneType),],
                tish[tish$GeneType=='',])
  tricerlnc<-tricer[grep('pseudo|lncRNA|assumed_protein_coding',tricer$gene_type),]
  codelnc<-codeall[grepl('pseudo|lncRNA|None',codeall$gene_type) & codeall$code==1,]
  
}


# alllnc<-codeall[!(codeall$gene_type %in% c('artifact','protein_coding','TEC','IG_C_gene','TR_C_gene',
#                                            'IG_V_gene','ribozyme','TR_C_gene','snoRNA','scaRNA','TEC','miRNA',
#                                            'misc_RNA','snRNA','vault_RNA')),]

# gennon<-alllnc$transcript_id[alllnc$gene_type!='None' & alllnc$save<(1)]
# genco<-alllnc$transcript_id[alllnc$gene_type!='None' & alllnc$save>0]
# gencoover<-alllnc$transcript_id[alllnc$gene_type!='None' & alllnc$save>1]
# gennon<-gennon[!(gennon %in% c(genco,gencoover))]
# 
# nonnon<-alllnc$transcript_id[alllnc$gene_type=='None' & alllnc$save<(1)]
# nonco<-alllnc$transcript_id[alllnc$gene_type=='None' & alllnc$save>0]
# noncoover<-alllnc$transcript_id[alllnc$gene_type=='None' & alllnc$save>1]
# nonnon<-nonnon[!(nonnon %in% c(nonco,noncoover))]

gencodelnc<-c(unique(tishlnc$Tid)[grep('ENST',unique(tishlnc$Tid))],
              unique(codelnc$transcript_id)[grep('ENST',unique(codelnc$transcript_id))],
              unique(tricerlnc$transcript_id)[grep('ENST',unique(tricerlnc$transcript_id))])
gencodelncove<-names(table(gencodelnc)[table(gencodelnc)>1])

noncodelnc<-c(unique(tishlnc$Tid)[-grep('ENST',unique(tishlnc$Tid))],
              unique(codelnc$transcript_id)[-grep('ENST',unique(codelnc$transcript_id))],
              unique(tricerlnc$transcript_id)[-grep('ENST',unique(tricerlnc$transcript_id))])
noncodelncove<-names(table(noncodelnc)[table(noncodelnc)>1])

gtffile<-as.data.frame(rtracklayer::import.gff('./data_prep/merge.gtf',format = 'gtf')) 
gtffile<-gtffile[!duplicated(gtffile$transcript_id),]
gtffile<-gtffile[grepl('pseudo|lncRNA',gtffile$gene_type) | grepl('NONHSAT',gtffile$transcript_id),]
alllnc<-gtffile$transcript_id
alllnc<-alllnc[!(alllnc %in% c(gencodelnc,gencodelncove,noncodelnc,noncodelncove))]
genco<-alllnc[-grep('NONHSAT',alllnc)]
nonco<-alllnc[grep('NONHSAT',alllnc)]

plot1<-cpc[cpc$`#ID` %in% c(genco,gencodelnc,gencodelncove)]
plot1$group<-ifelse(plot1$`#ID` %in% gencodelncove,'over',
                    ifelse(plot1$`#ID` %in% gencodelnc,'with','without'))

p1<-wilcox.test(plot1$coding_probability[plot1$group %in% c('with')],
                plot1$coding_probability[plot1$group %in% c('without')])[["p.value"]]
p2<-wilcox.test(plot1$coding_probability[plot1$group %in% c('over')],
                plot1$coding_probability[plot1$group %in% c('without')])[["p.value"]]

p<-list()

p[[1]]<-
ggplot(plot1, aes(coding_probability, colour = group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2)+   
  theme_bw(base_size = 18) +
  theme(legend.position='right',      
        panel.grid = element_blank()) +
  labs(x="Coding probability", y="Cumulative distribution") + 
  scale_color_manual(values=met.brewer("Hokusai1", 7)[c(6,5,4)])+
  scale_x_continuous(expand = c(0,0))+  
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1) +
  annotate('text',x=0.6,y=0.4,label=paste0('Wilcox P = ',signif(p1,3)),colour = met.brewer("Hokusai1", 7)[c(5)],size=4)+
  annotate('text',x=0.6,y=0.3,label=paste0('Wilcox P = ',signif(p2,3)),colour = met.brewer("Hokusai1", 7)[c(6)],size=4)

plot2<-cpc[cpc$`#ID` %in% c(nonco,noncodelnc,noncodelncove)]
plot2$group<-ifelse(plot2$`#ID` %in% noncodelncove,'over',
                    ifelse(plot2$`#ID` %in% noncodelnc,'with','without'))

p1<-wilcox.test(plot2$coding_probability[plot2$group %in% c('with')],
                plot2$coding_probability[plot2$group %in% c('without')])[["p.value"]]
p2<-wilcox.test(plot2$coding_probability[plot2$group %in% c('over')],
                plot2$coding_probability[plot2$group %in% c('without')])[["p.value"]]

p[[2]]<-
ggplot(plot2, aes(coding_probability, colour = group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2)+   
  theme_bw(base_size = 18) +
  theme(legend.position='right',         
        panel.grid = element_blank()) +
  labs(x="Coding probability", y="Cumulative distribution") + 
  scale_color_manual(values=met.brewer("Hokusai1", 7)[c(6,5,4)])+
  scale_x_continuous(expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1) +
  annotate('text',x=0.6,y=0.4,label=paste0('wilcox P = ',signif(p1,3)),colour = met.brewer("Hokusai1", 7)[c(5)],size=4)+
  annotate('text',x=0.6,y=0.3,label=paste0('wilcox P = ',signif(p2,3)),colour = met.brewer("Hokusai1", 7)[c(6)],size=4)

library(patchwork)
pdf("./qc/Figure1_CPC2.pdf",width = 10,height = 5)
wrap_plots(p,nrow=1, guides="collect") 
dev.off()

################################################
#
gtf<-rbind(data.table::fread('demo/ribocode_all/ribocode_allseq_longest.gtf'),
           data.table::fread('demo/ribocode_all/ribocode_allseq.gtf'))

# gtf$V9<-gsub('gene_id','gene_esmble',gtf$V9)
# gtf$V9<-gsub('orf_id','gene_id',gtf$V9)

gtf<-gtf[!duplicated(paste0(gtf$V3,'_',gtf$V9)),]
gtf$V3<-gsub('ORF','CDS',gtf$V3)

gtf<-gtf[gtf$V3=='exon',]

bed<-data.frame(V1=c(gtf$V1),V2=c(gtf$V4),V3=c(gtf$V5),
                V4=c(paste0('exon_',1:nrow(gtf))),V5=c(0),V6=c(gtf$V7))

# write.table(bed[bed$V6=='+',],file = 'allexon_pos.bed',sep = '\t',quote = F,col.names = F,row.names = F)
# write.table(bed[bed$V6=='-',],file = 'allexon_neg.bed',sep = '\t',quote = F,col.names = F,row.names = F)

pos<-rbind(
  cbind(data.table::fread('./conservation/pos_1.tab'),
           data.table::fread('./conservation/pos_2.tab'),
           data.table::fread('./conservation/pos_3.tab')),
  cbind(data.table::fread('./conservation/neg_1.tab'),
        data.table::fread('./conservation/neg_1.tab'),
        data.table::fread('./conservation/neg_1.tab'))
)


pos<-pos[match(bed$V4,pos$V1),]

#
colnames(pos)[c(4, 10, 16)]<-c('frame1','frame2','frame3')
pos$frame1<-as.numeric(pos$frame1);pos$frame2<-as.numeric(pos$frame2);pos$frame3<-as.numeric(pos$frame3)
pos[is.na(pos)]<-0

library(progress)

pb <- progress_bar$new(total = nrow(pos))

sorted_values <- apply(pos[, c(4, 10, 16)], 1, function(values) {
  pb$tick()
  sort(values)
})

result <- data.frame(t(sorted_values))
colnames(result) <- c("Min", "Median", "Max")

gtf <- cbind(gtf, result)

#Merge and then calculate the average score
f <-function(x) unlist(strsplit(x['V9'],';'))[1]
gtf$orf_ID<-apply(gtf,1,f)
f <-function(x) unlist(strsplit(x['orf_ID'],'[ ]'))[2]
gtf$orf_ID<-apply(gtf,1,f)
gtf$orf_ID<-gsub('["]','',gtf$orf_ID)
gtf$length<-pos$V2

gtf<-gtf[!duplicated(paste0(gtf$orf_ID,'_',gtf$V4,'_',gtf$V5)),]
plot<-gtf[,c(13,14,10,11,12)]

plot$Min<-as.numeric(plot$Min);plot$Median<-as.numeric(plot$Median);plot$Max<-as.numeric(plot$Max)
plot$length<-plot$length+1

library(dplyr)

plot2 <- plot %>%
  group_by(orf_ID) %>%
  summarise(across(length:Max, sum))

plot2$diff<-(plot2$Max-((plot2$Min+plot2$Median)/2))/plot2$length
plot2$maxave<-plot2$Max/plot2$length

#plot

alllnc<-codeall[!(codeall$gene_type %in% c('artifact','protein_coding','TEC','IG_C_gene','TR_C_gene')),]
cache<-plot2[match(alllnc$ORF_ID,plot2$orf_ID),]
alllnc$diff<-cache$diff
alllnc$maxqve<-cache$maxave

plot1<-alllnc[alllnc$gene_type!='None' ,]
plot1$group<-ifelse(plot1$save > 1,'with','without')

p1<-wilcox.test(plot1$diff[plot1$group %in% c('with')],
                plot1$diff[plot1$group %in% c('without')])[["p.value"]]

p<-list()

p[[1]]<-
  ggplot(plot1, aes(diff, colour = group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2)+   
  theme_bw(base_size = 18) +
  theme(legend.position='right',         
        panel.grid = element_blank()) +
  labs(x="Phylocsf Score Differences", y="Cumulative distribution") + 
  scale_color_manual(values=met.brewer("Hokusai1", 7)[c(6,5)])+
  scale_x_continuous(expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1) +
  coord_cartesian(ylim=c(0.5,1))+
  annotate('text',x=8,y=0.6,label=paste0('wilcox P = ',signif(p1,3)),colour = met.brewer("Hokusai1", 7)[c(5)],size=4)

p1<-wilcox.test(plot1$maxqve[plot1$group %in% c('with')],
                plot1$maxqve[plot1$group %in% c('without')])[["p.value"]]

p[[2]]<-
  ggplot(plot1, aes(maxqve, colour = group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2)+  
  theme_bw(base_size = 18) +
  theme(legend.position='right',         
        panel.grid = element_blank()) +
  labs(x="Max Phylocsf Score", y="Cumulative distribution") + 
  scale_color_manual(values=met.brewer("Hokusai1", 7)[c(6,5)])+
  scale_x_continuous(expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1) +
  annotate('text',x=0.6,y=0.4,label=paste0('wilcox P = ',signif(p1,3)),colour = met.brewer("Hokusai1", 7)[c(5)],size=4)


plot2<-alllnc[alllnc$gene_type=='None',]
plot2$group<-ifelse(plot2$save>1,'with','without')

p1<-wilcox.test(plot2$diff[plot2$group %in% c('with')],
                plot2$diff[plot2$group %in% c('without')])[["p.value"]]

p[[3]]<-
  ggplot(plot2, aes(diff, colour = group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2)+  
  theme_bw(base_size = 18) +
  theme(legend.position='right',        
        panel.grid = element_blank()) +
  labs(x="Phylocsf Score Differences", y="Cumulative distribution") + 
  scale_color_manual(values=met.brewer("Hokusai1", 7)[c(6,5)])+
  scale_x_continuous(expand = c(0,0))+  
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1) +
  coord_cartesian(ylim=c(0.5,1))+
  annotate('text',x=12,y=0.6,label=paste0('wilcox P = ',signif(p1,3)),colour = met.brewer("Hokusai1", 7)[c(5)],size=4)

p1<-wilcox.test(plot2$maxqve[plot2$group %in% c('with')],
                plot2$maxqve[plot2$group %in% c('without')])[["p.value"]]

p[[4]]<-
  ggplot(plot2, aes(maxqve, colour = group)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2)+  
  theme_bw(base_size = 18) +
  theme(legend.position='right',          #c(.25,.75),
        panel.grid = element_blank()) +
  labs(x="Max Phylocsf Score", y="Cumulative distribution") + 
  scale_color_manual(values=met.brewer("Hokusai1", 7)[c(6,5)])+
  scale_x_continuous(expand = c(0,0))+  
  scale_y_continuous(expand = c(0,0))+
  theme(aspect.ratio=1) +
  annotate('text',x=0.6,y=0.4,label=paste0('wilcox P = ',signif(p1,3)),colour = met.brewer("Hokusai1", 7)[c(5)],size=4)

#

library(patchwork)
pdf("./demo/FigureSX_Phylocsf.pdf",width = 10,height = 10)
wrap_plots(p,nrow=2, guides="collect") 
dev.off()

######################################################

#######################
#remove duplicate
library(dplyr)
dup <- merge %>%
  group_by(aaseq) %>%
  filter(n() > 1) %>%
  ungroup()

merge2<-merge[!(merge$aaseq %in% dup$aaseq),]

new_dup <- dup %>%
  group_by(aaseq) %>%
  summarize(
    name = paste(na.omit(name), collapse = ";"),  
    genename = paste(na.omit(genename), collapse = ";"),
    geneid = paste(na.omit(geneid), collapse = ";"),
    genetype = paste(na.omit(genetype), collapse = ";"),
    transcript_id = paste(na.omit(transcript_id), collapse = ";"),
    tishp = ifelse(all(is.na(tishp)), NA, paste(na.omit(tishp), collapse = ";")) ,
    codep = ifelse(all(is.na(codep)), NA, paste(na.omit(codep), collapse = ";"))  ,
    code = sum(code, na.rm = TRUE)  ,
    tish = sum(tish, na.rm = TRUE)  ,
    tricer_phase_score = max(all(is.na(tricer_phase_score)), NA, paste(na.omit(tricer_phase_score), collapse = ";"))  ,
    read_density = max(all(is.na(read_density)), NA, paste(na.omit(read_density), collapse = ";"))  ,
    tricer = sum(tricer, na.rm = TRUE)  ,
    save = sum(save, na.rm = TRUE)  
  ) %>%
  ungroup()

new_dup$tricer[new_dup$tricer>1]<-1
new_dup$save<-(new_dup$code+new_dup$tish+new_dup$tricer)
new_dup<-new_dup[,colnames(merge2)]

merge2<-rbind(merge2,new_dup)

#######################

#venn

library(ggvenn)
library(MetBrewer)

# venn<-list(Ribocode=alllnc$ORF_ID[alllnc$code>0],
#            RiboTish=alllnc$ORF_ID[alllnc$tish>0],
#            RiboTricer=alllnc$ORF_ID[alllnc$tricer>0])

alllnc<-codeall[!(codeall$gene_type %in% c('artifact','protein_coding','TEC','IG_C_gene','TR_C_gene',
                                           'IG_V_gene','ribozyme','TR_C_gene','snoRNA','scaRNA','TEC','miRNA',
                                           'misc_RNA','snRNA','vault_RNA')),]

venn<-list(Ribocode=unique(merge2$aaseq[merge2$code>0]),
           RiboTish=unique(merge2$aaseq[merge2$tish>0]),
           RiboTricer=unique(merge2$aaseq[merge2$tricer>0]))

pdf("./demo/venn_3tools.pdf",width = 5,height = 4)
ggvenn(venn,fill_color = met.brewer('Hiroshige',n=5),
          stroke_size = 0.5, set_name_size = 4)
dev.off()

#Length distribution
plot<-data.frame(length=c(alllnc$ORF_length[alllnc$code>0],alllnc$ORF_length[alllnc$tish>0],alllnc$ORF_length[alllnc$tricer>0]),
                 group=c(rep('Ribocode',length(alllnc$ORF_length[alllnc$code>0])),
                         rep('RiboTish',length(alllnc$ORF_length[alllnc$tish>0])),
                         rep('RiboTricer',length(alllnc$ORF_length[alllnc$tricer>0]))))

pdf("./demo/length_3tools.pdf",width = 6,height = 5)
ggplot(plot, aes(length, after_stat(density), colour = group)) +
  geom_density(adjust = 1, size=1)+ 
  xlim(0,300)+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF",'#96345A74'))   
dev.off()


################################
#And finally, two examples
# ENSG00000270170.2_196942717_196943016_99
# ENSG00000175701.11_110212525_110212355_56 

tricer<-data.table::fread('./demo/merge.txt')

ENSG00000270170<-tricer[tricer$gene_id=='ENSG00000270170.2',]
ENSG00000175701<-tricer[tricer$gene_id=='ENSG00000175701.11' & tricer$ORF_type=='annotated',]

write.csv(ENSG00000270170,file = 'ENSG00000270170_plot.csv')
write.csv(ENSG00000175701,file = 'ENSG00000175701_plot.csv')

library(ggpubr)
appleplot('KRASIM','ENSG00000270170_plot.csv')
appleplot('MPM','ENSG00000175701_plot.csv')

