
#The match was re-run before the Sankey diagram to give a more appropriate comparison of the peptide sources

tmtpep<-data.table::fread('./MS_result/TMT_combined_peptide.tsv')
tmtpep<-tmtpep[!duplicated(tmtpep$Peptide),]
tmtpep<-tmtpep[-grep('HUMAN',tmtpep$Protein),]
tmtpep<-tmtpep[-grep('HUMAN',tmtpep$`Mapped Proteins`),] 

load('MS_result/TMT_quant.RData')
data<-data[!duplicated(data$Peptide),]
data<-data[grep('sPep',data$ProteinID),]

tmtpep<-tmtpep[tmtpep$Peptide %in% data$Peptide,]
tmtpep<-tmtpep[,-c(20:30)]

# 
lbfpep<-data.table::fread('./MS_result/LBF_combined_peptide.tsv')
lbfpep<-lbfpep[!duplicated(lbfpep$`Peptide Sequence`),]
lbfpep<-lbfpep[-grep('HUMAN',lbfpep$Protein),]
lbfpep<-lbfpep[-grep('HUMAN',lbfpep$`Mapped Proteins`),]

load('MS_result/LBF_quant.RData')
quant<-quant[!duplicated(quant$PeptideSequence),]
quant<-quant[grep('sPep',quant$ProteinName),]

lbfpep<-lbfpep[lbfpep$`Peptide Sequence` %in% gsub("\\[.*?\\]|[a-z]", "", quant$PeptideSequence),]
length(unique(lbfpep$Protein))
###################################

sankytable<-data.frame(pepgroup=unique(c(lbfpep$`Protein ID`,tmtpep$`Protein ID`)),
                       seq=c(NA))

for (i in 1:nrow(sankytable)){
  pep1<-lbfpep$`Peptide Sequence`[lbfpep$`Protein ID`==sankytable$pepgroup[i]]
    pep2<-tmtpep$Peptide[tmtpep$`Protein ID`==sankytable$pepgroup[i]]
    peplist<-unique(c(pep1,pep2))
    sankytable$seq[i]<-paste0(peplist,collapse = '&')
}

map<-merge[merge$save>0,]
sankytable$orflist<-c(NA)
sankytable$orindexlist<-c(NA)
sankytable$ribosource<-c(NA)

library(TNSMD)
myindex<-generate_index('./writeindex_2.csv','other',0,'sPep')

f <-function(x) unlist(strsplit(x['name'],'[|]'))[2]
myindex$ID<-apply(myindex,1,f)
f <-function(x) unlist(strsplit(x['name'],'[ ]'))[2]
myindex$genename<-apply(myindex,1,f)

for (i in 1:nrow(sankytable)){
  patterns <- unlist(strsplit(sankytable$seq[i],'&')) 
  matches_list <- lapply(patterns, function(pat) grep(pat, map$aaseq))
  common_matches <- Reduce(intersect, matches_list)
  
  source<-merge[merge$name %in% c(map$name[common_matches]),]
  code<-ifelse(sum(source$code)>0,T,F)
  tish<-ifelse(sum(source$tish)>0,T,F)
  tricer<-ifelse(sum(source$tricer)>0,T,F)
  source<-c('RiboCode'[code],'Ribo-TISH'[tish],'ribotricer'[tricer])
  
  sankytable$orflist[i]<-paste0(map$name[common_matches],collapse = ',')
  sankytable$ribosource[i]<-paste0(source,collapse = ' & ')
  
  patterns <- unlist(strsplit(sankytable$seq[i],'&')) 
  matches_list <- lapply(patterns, function(pat) grep(pat, myindex$sequence))
  common_matches <- Reduce(intersect, matches_list)
  sankytable$orindexlist[i]<-paste0(myindex$ID[common_matches],collapse = ',')

}

sankytable$tpart<-c(NA)
sankytable$tpartsource<-c(NA)

for (i in 1:nrow(sankytable)){
  patterns <- unlist(strsplit(sankytable$seq[i],'&')) 
  matches_list <- lapply(patterns, function(pat) grep(pat, part3$pepseq))
  common_matches <- Reduce(intersect, matches_list)
  
  source<-unique(part3$source[common_matches])
  
  sankytable$tpart[i]<-paste0(part3$genename[common_matches],collapse = ',')
  
  sankytable$tpartsource[i]<-paste0(source,collapse = ',')
}

# LBF<-paste0(lbfpep$`Protein ID`,lbfpep$`Mapped Proteins`,collapse = ',')
# LBF<-paste0(unlist(strsplit(LBF,'tr')),collapse = ',')
# LBF<-paste0(unlist(strsplit(LBF,'[|]')),collapse = ',')
# LBF<-unlist(strsplit(LBF,','))
# LBF<-unique(LBF)
# LBF<-LBF[!(LBF %in% c('',' '))]
# 
# TMT<-paste0(tmtpep$Protein,tmtpep$`Mapped Proteins`,collapse = ',')
# TMT<-paste0(unlist(strsplit(TMT,'tr')),collapse = ',')
# TMT<-paste0(unlist(strsplit(TMT,'[|]')),collapse = ',')
# TMT<-unlist(strsplit(TMT,','))
# TMT<-unique(TMT)
# TMT<-TMT[!(TMT %in% c('',' '))]
# 
# plot3<-data.frame(lncORF=c(LBF,TMT),
#                   from=c(rep('LBF',length(LBF)),rep('TMT',length(TMT))))

# f <-function(x) unlist(strsplit(x['name'],'[|]'))[2]
# myindex$ID<-apply(myindex,1,f)
# 
# match<-myindex[match(plot3$lncORF,myindex$ID),]
# plot3$name<-match$ORFname
# plot3$aaseq<-match$sequence

#
# cache<-alllnc[alllnc$save!=0,]
# 
# plot3$ribo<-ifelse(plot3$aaseq %in% cache$AAseq,1,0)
# plot3$genecode<-ifelse(plot3$aaseq %in% genelnc$V22,1,0) #None of the peptides from GENCODE
# plot3$uniprot<-ifelse(plot3$aaseq %in% filter$seq,1,0)

#

# plot3$ribocode<-ifelse(plot3$aaseq %in% cache$AAseq[cache$code>0],1,0)
# plot3$ribotish<-ifelse(plot3$aaseq %in% cache$AAseq[cache$tish>0],1,0)
# plot3$ribotricer<-ifelse(plot3$aaseq %in% cache$AAseq[cache$tricer>0],1,0)
# 
# gene_anno<-data.frame(aaseq=c(cache$AAseq,filter$seq),
#                       genetype=c(cache$gene_type,filter$type))
# 
# gene_anno$adjtype<-ifelse(gene_anno$genetype=='lncRNA','lncRNA',
#                           ifelse(gene_anno$genetype=='None','NONCODE sourced lncRNA',
#                                  ifelse(grepl('pseudo',gene_anno$genetype),'Pseudogene','Others')))
# 
# gene_anno<-gene_anno[match(plot3$aaseq,gene_anno$aaseq),]
# 
# plot3$genetype<-c(gene_anno$adjtype)

##############################################
#visualization

sankytable$mssource<-ifelse(sankytable$pepgroup %in% lbfpep$`Protein ID` & sankytable$pepgroup %in% tmtpep$`Protein ID`,
                            'Jiang et al 2019 & Gao et al 2019', ifelse(sankytable$pepgroup %in% lbfpep$`Protein ID`,'Jiang et al 2019','Gao et al 2019'))

match<-myindex[match(sankytable$pepgroup,myindex$ID),]

sankytable$genename<-match$genename

# gene_anno<-data.frame(rtracklayer::import.gff('../biodata/gencode.v46.annotation.gtf',format = 'gtf'))
# gene_anno<-gene_anno[gene_anno$type=='gene',]
gene_match<-gene_anno[match(sankytable$genename,gene_anno$gene_name),]
sankytable$genetype<-ifelse(is.na(gene_match$gene_type),'NONCODE',ifelse(
  gene_anno$gene_type=='lncRNA','lncRNA','pseudogene'
))

library(ggsankey)
library(ggplot2)
library(dplyr) 

plot3<-sankytable
plot3$Index<-ifelse(plot3$orflist !='' & plot3$tpartsource!='' ,'Ribo-seq & UniProt',
                    ifelse(plot3$orflist !='','Ribo-seq','UniProt'))

# plot3$tools<-ifelse(plot3$ribocode ==1 & plot3$ribotish==1 &plot3$ribotricer==1 ,'RiboCode & RiboTish & RiboTricer',
#                     ifelse(plot3$ribocode ==1 & plot3$ribotish==1 ,'RiboCode & RiboTish',
#                            ifelse(plot3$ribotish==1 & plot3$ribotricer==1,'RiboTish & RiboTricer',
#                                   ifelse(plot3$ribocode ==1,'RiboCode',
#                                          ifelse(plot3$ribotish==1,'RiboTish',
#                                                 ifelse(plot3$ribotricer==1,'RiboTricer','UniProt'))))))

# plot3$From<-ifelse(plot3$from == 'LBF','Jiang et al 2019','Gao et al 2019')

plot3$ribosource[plot3$ribosource=='']<-c('UniProt')

df <- plot3 %>%
  make_long(mssource, Index, ribosource,genetype) 

pdf("./Figure3/sankey.pdf",width = 6,height = 4)
ggplot(df, aes(x = x, next_x = next_x, 
               node = node, next_node = next_node,
               fill = factor(node),label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 2.5, color = 1, fill = "white") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 16) +xlab(NULL)+
  theme(legend.position = "none")+
  theme(aspect.ratio=0.5)
dev.off()


# cache<-alllnc[alllnc$save!=0,]
# cache$MS<-ifelse(cache$AAseq %in% plot3$aaseq,'MS','Ribo-seq')
# compare<-cache[,c('gene_type','chrom','start_codon','strand','ORF_length','ORF_tstart','ORF_tstop',
#                   'Psites_coverage_frame0','Psites_frame0_RPKM','pval_combined','tishp','tricer_phase_score',
#                   'read_density','code','tish','tricer','diff','maxqve','MS')]
# 
# compare<-compare[!compare$gene_type %in% c('miRNA','misc_RNA','scaRNA','snoRNA','snRNA','vault_RNA'),]
# compare$gene_type<-ifelse(grepl('pseudoge',compare$gene_type),'pseudogene',compare$gene_type)
# compare$gene_type[compare$gene_type=='None']<-c('NONCODE lncRNA')
# 
# compare$tishp[is.na(compare$tishp)]<-1
# compare$tricer_phase_score[is.na(compare$tricer_phase_score)]<-0
# compare$read_density[is.na(compare$read_density)]<-0
# 
# compare$Psites_coverage_frame0<- as.numeric(gsub("%", "", compare$Psites_coverage_frame0)) / 100
# compare$tricer_phase_score<-as.numeric(compare$tricer_phase_score)
# compare$read_density<-as.numeric(compare$read_density)

msorflist<-paste0(sankytable$orindexlist,collapse = ',')
msorflist<-unlist(strsplit(msorflist,','))
msorflist<-unique(msorflist)
msorflist<-msorflist[!(msorflist %in% c('',' '))]

mltable<-data.frame(aaseq=c(myindex$sequence[myindex$ID %in% msorflist],myindex$sequence[!(myindex$ID %in% msorflist)]),
                    group=c(rep('MS peptide',length(myindex$sequence[myindex$ID %in% msorflist])),
                            rep('Ribo-seq peptide',length(myindex$sequence[!(myindex$ID %in% msorflist)])) ))

library(Biostrings)
library(seqinr)

data(aaindex)

stats  <- do.call(rbind, sapply(mltable$aaseq, FUN = function(x) {
  
  result <- AAstat(s2c(x), plot = FALSE)
  length<-nchar(x)
  compo <- as.vector(result$Compo)  
  prop <- unlist(result$Prop) 
  pi <- result$Pi  
  
  data.frame(
    A = compo[1]/length, C = compo[2]/length, D = compo[3]/length, E = compo[4]/length, F = compo[5]/length, 
    G = compo[6]/length, H = compo[7]/length, I = compo[8]/length, K = compo[9]/length, L = compo[10]/length, 
    M = compo[11]/length, N = compo[12]/length, P = compo[13]/length, Q = compo[14]/length, R = compo[15]/length, 
    S = compo[16]/length, T = compo[17]/length, V = compo[18]/length, W = compo[19]/length, Y = compo[20]/length, 
    Tiny = prop["Tiny"], Small = prop["Small"], Aliphatic = prop["Aliphatic"], 
    Aromatic = prop["Aromatic"], Non_polar = prop["Non.polar"], Polar = prop["Polar"], 
    Charged = prop["Charged"], Basic = prop["Basic"], Acidic = prop["Acidic"], 
    Pi = pi
  )
}, simplify = FALSE))

prop_mean <- t(sapply(mltable$aaseq, FUN=function(x) c(lapply(aaindex, FUN=function(y) mean(y$I[aaa(s2c(x))], na.rm=TRUE)), Length=getLength(x), PMW=pmw(s2c(x)), PI=computePI(s2c(x)), unlist(stats["Prop",x]))))

mlarray<-cbind(stats,prop_mean)

# save(mlarray,file = './Figure3/MLarray.RData')
load('./Figure3/MLarray.RData')

library(caret)
# dummies <- dummyVars(MS ~ ., data = compare)
# mlarray<-predict(dummies, newdata = compare)
# mlarray<-as.data.frame(mlarray)
mlarray$label<-ifelse(mltable$group=='MS peptide',1,0)
mlarray$label<-factor(mlarray$label)

mlarray <- data.frame(lapply(mlarray, function(column) {if (is.list(column)) {return(sapply(column, function(x) if(is.list(x)) unlist(x) else x))} else { return(column)}}))

library(glmnet)
cv<-cv.glmnet(as.matrix(mlarray[,-578]),as.numeric(mlarray$label),nfold=10,family='binomial')

fit.train <- cv$glmnet.fit
pdf(file = "./Figure3/lasso1.pdf",width =8,height = 6)
plot(cv,las=1)
dev.off()

fit <- glmnet(as.matrix(mlarray[,-578]),as.numeric(mlarray$label),family = "gaussian")
pdf(file = "./Figure3/lasso2.pdf",width =8,height = 6)
plot(fit,xvar = "lambda",label = TRUE, las=1)
dev.off()

coef.min = coef(cv, s = "lambda.min")  
active.min = which(coef.min@i != 0) #
lasso_feature <- coef.min@Dimnames[[1]][coef.min@i+1] 

mlfilter<-mlarray[,colnames(mlarray) %in% c(lasso_feature,'label')]

#Needs downsampling, do 10 downsamples, then 50 iterations of Boruta each, for a total of 500 iterations

library(Boruta)

borutaresult<-data.frame()

for(i in 1:10) {

set.seed(20240926+i)
  
mlarray_down <- downSample(x = mlfilter[, -ncol(mlfilter)],
                         y = mlfilter$label)

Var.Selec<-Boruta(Class~., data= mlarray_down)

impmerge<-attStats(Var.Selec)
impmerge$runs<-i
impmerge$variable<-rownames(impmerge)

# imparray<-Var.Selec[["ImpHistory"]]
# imparray<-reshape2::melt(imparray)
# imparray$runs<-i

borutaresult<-rbind(borutaresult,impmerge)

}

# borutaresult$type<-ifelse(grepl('chrom',borutaresult$variable),'chrom',
#                           ifelse(grepl('gene_type',borutaresult$variable),'gene type',
#                                  ifelse(grepl('start_codon',borutaresult$variable),'start codon',
#                                         ifelse(borutaresult$variable %in% c('diff','maxqve'),'conservation',
#                                                ifelse(borutaresult$variable %in% c('ORF_tstart','ORF_tstop','strand+','strand-'),'location',
#                                                       ifelse(borutaresult$variable %in% c('tish','code','tricer'),'tools','other'
#                                                       ))))))


p<-ggplot(borutaresult, 
          aes(x = reorder(variable,medianImp), y = medianImp)) +
    stat_boxplot(geom = "errorbar",width=0.2)+
    geom_boxplot(outlier.shape = 1,aes(fill=variable),show.legend = F)+
    # geom_boxplot(width=0.8,aes(fill=group),colour='black',alpha = 1,outlier.shape = NA)+
    labs(y="Median Importance (Z-Score)")+ 
    theme_classic(base_size = 10)+ 
    scale_fill_manual(values=MetBrewer::met.brewer("Hokusai1", 24))+
    theme_classic(base_size = 10)+ 
    theme(legend.position = "none")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank())+
    coord_flip()+
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
          panel.grid.minor = element_line(color = "gray", size = 0.25))


pdf("./Figure3/boruta_aaindex.pdf",width = 4,height = 8)
p
dev.off()

#Length has an effect, consistent with previous reports

##################################################
#box plot
library(ggpubr)

p<-list()

# compare[1:5,]

# p[[1]]<-
# ggplot(compare[compare$gene_type!='NONCODE lncRNA'], 
#        aes(x = MS, y = log2(Psites_frame0_RPKM))) +
#   stat_boxplot(geom = "errorbar",width=0.2)+
#   geom_boxplot(outlier.shape = NA,aes(fill=MS),show.legend = F)+
#   labs(y="log2 (Frame0 RPKM)")+ 
#   theme_classic(base_size = 10)+ 
#   scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
#   theme_classic(base_size = 10)+ 
#   theme(legend.position = "none")+
#   xlab(NULL)+
#   theme(aspect.ratio=2)+
#   theme(axis.title.x = element_text(size = 16, color = "black"),
#         axis.title.y = element_text(size = 16, color = "black"),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         legend.text = element_text(size = 14, color = "black"),
#         legend.background = element_rect(fill = "white"),
#         legend.key = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = NA, color = "black"),
#         panel.grid.minor = element_line(color = "gray", size = 0.25))+
#   coord_cartesian(ylim=c(-5,10))+
#   stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 10)
# 
# p[[2]]<-
# ggplot(compare[compare$gene_type!='NONCODE lncRNA'], 
#        aes(x = MS, y = Psites_coverage_frame0)) +
#   stat_boxplot(geom = "errorbar",width=0.2)+
#   geom_boxplot(outlier.shape = NA,aes(fill=MS),show.legend = F)+
#   labs(y="Frame0 coverage")+ 
#   theme_classic(base_size = 10)+ 
#   scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
#   theme_classic(base_size = 10)+ 
#   theme(legend.position = "none")+
#   xlab(NULL)+
#   theme(aspect.ratio=2)+
#   theme(axis.title.x = element_text(size = 16, color = "black"),
#         axis.title.y = element_text(size = 16, color = "black"),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         legend.text = element_text(size = 14, color = "black"),
#         legend.background = element_rect(fill = "white"),
#         legend.key = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = NA, color = "black"),
#         panel.grid.minor = element_line(color = "gray", size = 0.25))+
#   coord_cartesian(ylim=c(0,1))+
#   stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1)
# 
# p[[3]]<-
# ggplot(compare[compare$gene_type!='NONCODE lncRNA'], 
#        aes(x = MS, y = -log10(pval_combined))) +
#   stat_boxplot(geom = "errorbar",width=0.2)+
#   geom_boxplot(outlier.shape = NA,aes(fill=MS),show.legend = F)+
#   labs(y="-log10 (combined pvalue)")+ 
#   theme_classic(base_size = 10)+ 
#   scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
#   theme_classic(base_size = 10)+ 
#   theme(legend.position = "none")+
#   xlab(NULL)+
#   theme(aspect.ratio=2)+
#   theme(axis.title.x = element_text(size = 16, color = "black"),
#         axis.title.y = element_text(size = 16, color = "black"),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         legend.text = element_text(size = 14, color = "black"),
#         legend.background = element_rect(fill = "white"),
#         legend.key = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = NA, color = "black"),
#         panel.grid.minor = element_line(color = "gray", size = 0.25))+
#   coord_cartesian(ylim=c(0,10))+
#   stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 10)
# 
# library(patchwork)
# pdf("./Figure3/RPKM_GENCODE.pdf",width = 15,height = 5)
# wrap_plots(p,nrow=1, guides="collect") 
# dev.off()

p<-list()

plot<-mlfilter
plot$group<-ifelse(plot$label==1,'Ms peptide','Ribo peptide')

 p[[1]]<-
   ggplot(plot, 
          aes(x = group, y = Length)) +
   stat_boxplot(geom = "errorbar",width=0.2)+
   geom_boxplot(outlier.shape = NA,aes(fill=group),show.legend = F)+
   labs(y="Peptide length")+ 
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
   coord_cartesian(ylim=c(0,400))+
   stat_compare_means(method = "wilcox.test",label.y = 370)

 p[[2]]<-
   ggplot(plot, 
          aes(x = group, y = GEOR030105)) +
   stat_boxplot(geom = "errorbar",width=0.2)+
   geom_boxplot(outlier.shape = NA,aes(fill=group),show.legend = F)+
   labs(y="AA index: GEOR030105")+ 
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
   coord_cartesian(ylim=c(0.9,1.15))+
   stat_compare_means(method = "wilcox.test",label.y = 1.13)
 
 p[[3]]<-
   ggplot(plot, 
          aes(x = group, y = GEOR030101)) +
   stat_boxplot(geom = "errorbar",width=0.2)+
   geom_boxplot(outlier.shape = NA,aes(fill=group),show.legend = F)+
   labs(y="AA index: GEOR030101")+ 
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
   coord_cartesian(ylim=c(0.94,1.08))+
   stat_compare_means(method = "wilcox.test",label.y = 1.06)
 
 p[[4]]<-
   ggplot(plot, 
          aes(x = group, y = RICJ880107)) +
   stat_boxplot(geom = "errorbar",width=0.2)+
   geom_boxplot(outlier.shape = NA,aes(fill=group),show.legend = F)+
   labs(y="AA index: RICJ880107")+ 
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
   coord_cartesian(ylim=c(0.7,1.5))+
   stat_compare_means(method = "wilcox.test",label.y = 1.5)
 

 p[[5]]<-  
  ggplot(plot, 
         aes(x = group, y = K)) +
    stat_boxplot(geom = "errorbar",width=0.2)+
    geom_boxplot(outlier.shape = NA,aes(fill=group),show.legend = F)+
    labs(y="AA ratio: K")+ 
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
    coord_cartesian(ylim=c(0,0.15))+
    stat_compare_means(method = "wilcox.test",label.y =0.15)

  library(patchwork)
  pdf("./Figure3/AAindex.pdf",width = 15,height = 10)
  wrap_plots(p,nrow=2, guides="collect") 
  dev.off()

# Modeled separately based on these elements: length,aaindex,all together
# Physicochemical factors may be another factor affecting the detection of lncRNA-encoded peptides
  
  ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       sampling = "down")

  
  set.seed(5627)
  
  mlarray$label<-factor(ifelse(mlarray$label==1,'M','R'),levels = c('M','R'))
  
  set.seed(618)
  inTraining <- createDataPartition(mlarray$label, p = .7, list = FALSE)
  training <- mlarray[ inTraining,]
  testing  <- mlarray[-inTraining,]
  
  trainset1<-training[,colnames(training) %in% c(lasso_feature[2:5],'label')]
  trainset2<-training[,colnames(training) %in% c(lasso_feature[6:24],'label')]
  trainset3<-training[,colnames(training) %in% c(lasso_feature[25],'label')]
  trainset4<-training[,colnames(training) %in% c(lasso_feature,'label')]
  
  down_inside <- train(label ~ ., data = trainset1, 
                        method = "knn",
                        metric = "ROC",tuneLength = 4,
                        trControl = ctrl)

  library('pROC')
  probs = predict(down_inside,testing,type = "prob")
  ROC1 = roc(response = testing$label,
               predictor = probs$M,
               levels = levels(testing$label))
  
  down_inside <- train(label ~ ., data = trainset3, 
                       method = "rf",
                       metric = "ROC",tuneLength = 4,
                       trControl = ctrl)
  
  probs = predict(down_inside,testing,type = "prob")
  ROC2 = roc(response = testing$label,
            predictor = probs$M,
            levels = levels(testing$label))
  
  down_inside <- train(label ~ ., data = trainset2, 
                       method = "rf",
                       metric = "ROC",tuneLength = 4,
                       trControl = ctrl)
  
  probs = predict(down_inside,testing,type = "prob")
  ROC3 = roc(response = testing$label,
            predictor = probs$M,
            levels = levels(testing$label)) 
  
  down_inside <- train(label ~ ., data = trainset4, 
                       method = "rf",
                       metric = "ROC",tuneLength = 4,
                       trControl = ctrl)
  
  probs = predict(down_inside,testing,type = "prob")
  ROC4 = roc(response = testing$label,
            predictor = probs$M,
            levels = levels(testing$label))
  
  g <- ggroc(list(AAratio=ROC1,Length=ROC2,AAindex=ROC3,ALL=ROC4),size=0.8)
  
  pdf("Figure3/ROC.pdf",width = 8,height = 6)
  g+theme_bw(base_size=18)+
    theme(panel.grid.major =element_blank(), 
          panel.background=element_rect(size =1.1,fill='transparent', color='black'),
          panel.grid.minor = element_blank(),panel.border = element_blank(),
          legend.title = element_blank())+
    geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
    scale_x_reverse(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme(aspect.ratio=1)+
    coord_fixed(ratio = 0.8)+
    scale_colour_manual(values=MetBrewer::met.brewer("Derain", 5,type ='discrete')[c(1,2,4,5)])
  dev.off()
  
# p<-list()
# 
# p[[1]]<-
#   ggplot(compare[compare$gene_type=='NONCODE lncRNA'], 
#          aes(x = MS, y = log2(Psites_frame0_RPKM))) +
#   stat_boxplot(geom = "errorbar",width=0.2)+
#   geom_boxplot(outlier.shape = NA,aes(fill=MS),show.legend = F)+
#   labs(y="log2 (Frame0 RPKM)")+ 
#   theme_classic(base_size = 10)+ 
#   scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
#   theme_classic(base_size = 10)+ 
#   theme(legend.position = "none")+
#   xlab(NULL)+
#   theme(aspect.ratio=2)+
#   theme(axis.title.x = element_text(size = 16, color = "black"),
#         axis.title.y = element_text(size = 16, color = "black"),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         legend.text = element_text(size = 14, color = "black"),
#         legend.background = element_rect(fill = "white"),
#         legend.key = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = NA, color = "black"),
#         panel.grid.minor = element_line(color = "gray", size = 0.25))+
#   coord_cartesian(ylim=c(-5,10))+
#   stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 10)
# 
# p[[2]]<-
#   ggplot(compare[compare$gene_type=='NONCODE lncRNA'], 
#          aes(x = MS, y = Psites_coverage_frame0)) +
#   stat_boxplot(geom = "errorbar",width=0.2)+
#   geom_boxplot(outlier.shape = NA,aes(fill=MS),show.legend = F)+
#   labs(y="Frame0 coverage")+ 
#   theme_classic(base_size = 10)+ 
#   scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
#   theme_classic(base_size = 10)+ 
#   theme(legend.position = "none")+
#   xlab(NULL)+
#   theme(aspect.ratio=2)+
#   theme(axis.title.x = element_text(size = 16, color = "black"),
#         axis.title.y = element_text(size = 16, color = "black"),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         legend.text = element_text(size = 14, color = "black"),
#         legend.background = element_rect(fill = "white"),
#         legend.key = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = NA, color = "black"),
#         panel.grid.minor = element_line(color = "gray", size = 0.25))+
#   coord_cartesian(ylim=c(0,1))+
#   stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1)
# 
# p[[3]]<-
#   ggplot(compare[compare$gene_type=='NONCODE lncRNA'], 
#          aes(x = MS, y = -log10(pval_combined))) +
#   stat_boxplot(geom = "errorbar",width=0.2)+
#   geom_boxplot(outlier.shape = NA,aes(fill=MS),show.legend = F)+
#   labs(y="-log10 (combined pvalue)")+ 
#   theme_classic(base_size = 10)+ 
#   scale_fill_manual(values=MetBrewer::met.brewer("Cassatt1", 2))+
#   theme_classic(base_size = 10)+ 
#   theme(legend.position = "none")+
#   xlab(NULL)+
#   theme(aspect.ratio=2)+
#   theme(axis.title.x = element_text(size = 16, color = "black"),
#         axis.title.y = element_text(size = 16, color = "black"),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         legend.text = element_text(size = 14, color = "black"),
#         legend.background = element_rect(fill = "white"),
#         legend.key = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = NA, color = "black"),
#         panel.grid.minor = element_line(color = "gray", size = 0.25))+
#   coord_cartesian(ylim=c(0,10))+
#   stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 10)
# 
# library(patchwork)
# pdf("./Figure3/RPKM_NONCODE.pdf",width = 15,height = 5)
# wrap_plots(p,nrow=1, guides="collect") 
# dev.off()

#
# p<-list()
# 
# p[[1]]<-
# ggplot(compare[compare$gene_type!='NONCODE lncRNA'], aes(ORF_length, after_stat(density), colour = MS)) +
#   geom_density(adjust = 2, size=1)+     
#   theme_classic(base_size = 18)+
#   theme(axis.text.x = element_text(size=12))+
#   theme(aspect.ratio=0.7)+
#   scale_color_manual(values = met.brewer("Cassatt1", 2))+
#   xlab('ORF length')+ylab('Density')
# 
# p[[2]]<-
# ggplot(compare[compare$gene_type!='NONCODE lncRNA'], aes(diff, after_stat(density), colour = MS)) +
#   geom_density(adjust = 2, size=1)+     
#   theme_classic(base_size = 18)+
#   theme(axis.text.x = element_text(size=12))+
#   theme(aspect.ratio=0.7)+
#   scale_color_manual(values = met.brewer("Cassatt1", 2))+
#   xlab('PhyloCSF score difference (per base)')+ylab('Density')
# 
# #这个差值怎么可以的是负值呢？怎么理解？翻一下代码
# 
# p[[3]]<-
# ggplot(compare[compare$gene_type!='NONCODE lncRNA'], aes(maxqve, after_stat(density), colour = MS)) +
#   geom_density(adjust = 2, size=1)+     
#   theme_classic(base_size = 18)+
#   theme(axis.text.x = element_text(size=12))+
#   theme(aspect.ratio=0.7)+
#   scale_color_manual(values = met.brewer("Cassatt1", 2))+
#   xlab('Max PhyloCSF score (per base)')+ylab('Density')
# 
# pdf("./Figure3/ORF_length_GENCODE.pdf",width = 15,height = 5)
# wrap_plots(p,nrow=1, guides="collect") 
# dev.off()
# 
# p<-list()
# 
# p[[1]]<-
#   ggplot(compare[compare$gene_type=='NONCODE lncRNA'], aes(ORF_length, after_stat(density), colour = MS)) +
#   geom_density(adjust = 2, size=1)+     
#   theme_classic(base_size = 18)+
#   theme(axis.text.x = element_text(size=12))+
#   theme(aspect.ratio=0.7)+
#   scale_color_manual(values = met.brewer("Cassatt1", 2))+
#   xlab('ORF length')+ylab('Density')
# 
# p[[2]]<-
#   ggplot(compare[compare$gene_type=='NONCODE lncRNA'], aes(diff, after_stat(density), colour = MS)) +
#   geom_density(adjust = 2, size=1)+     
#   theme_classic(base_size = 18)+
#   theme(axis.text.x = element_text(size=12))+
#   theme(aspect.ratio=0.7)+
#   scale_color_manual(values = met.brewer("Cassatt1", 2))+
#   xlab('PhyloCSF score difference (per base)')+ylab('Density')
# 
# 
# p[[3]]<-
#   ggplot(compare[compare$gene_type=='NONCODE lncRNA'], aes(maxqve, after_stat(density), colour = MS)) +
#   geom_density(adjust = 2, size=1)+     
#   theme_classic(base_size = 18)+
#   theme(axis.text.x = element_text(size=12))+
#   theme(aspect.ratio=0.7)+
#   scale_color_manual(values = met.brewer("Cassatt1", 2))+
#   xlab('Max PhyloCSF score (per base)')+ylab('Density')
# 
# pdf("./Figure3/ORF_length_NONCODE.pdf",width = 15,height = 5)
# wrap_plots(p,nrow=1, guides="collect") 
# dev.off()

# 
# orfbed<-data.frame(V1=c(cache$chrom),
#                    V2=c(cache$ORF_gstart),
#                    V3=c(cache$ORF_gstop),
#                    V4=c(cache$ORF_ID),
#                    V5=c(0),
#                    V6=c(cache$strand))
# 
# riboorf<-orfbed[compare$MS=='Ribo-seq',]
# msorf<-orfbed[compare$MS!='Ribo-seq',]
# 
# write.table(riboorf,file = 'RiboseqORF.bed',quote = F,sep = '\t',row.names = F,col.names = F)
# write.table(msorf,file = 'MSORF.bed',quote = F,sep = '\t',row.names = F,col.names = F)
alllength<-read.csv('./data_prep/alltranslength.csv')
alllength<-alllength[alllength$transid %in% cache$transcript_id,]

orfdis<-data.frame()

for (i in 1:nrow(cache)){
  dis<-data.frame(
    position=c(round( (cache$ORF_tstart[i]/alllength$length[which(alllength$transid==cache$transcript_id[i])])*100,0):round( (cache$ORF_tstop[i]/alllength$length[which(alllength$transid==cache$transcript_id[i])])*100,0)),
    group=c(NA)
  )
  dis$group<-cache$MS[i]
  orfdis<-rbind(orfdis,dis)
}

pdf("./Figure3/ORF_position.pdf",width = 6,height = 5)
ggplot(orfdis, aes(position, after_stat(density), colour = group)) +
  geom_density(adjust = 2, size=1)+     
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = met.brewer("Cassatt1", 2))+
  xlab('Relative position')+ylab('Density')
dev.off()

mean(alllength$length[alllength$transid %in% cache$transcript_id[cache$MS=='MS']])
mean(alllength$length[alllength$transid %in% cache$transcript_id[cache$MS!='MS']])

match<-alllength[match(cache$transcript_id,alllength$transid),]
cache$translength<-match$length

ggplot(cache, aes(log2(translength), after_stat(density), colour = MS)) +
  geom_density(adjust = 2, size=1)+     
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = met.brewer("Cassatt1", 2))+
  xlab('Transcript length')+ylab('Density')


################################################
#Missing values

load('LBF_protein_exp.RData')

lncexp<-protein[grep('sPep',protein$Protein),]
f<-function(x) sum(!is.na(x))
rownum<-apply(lncexp,1,f)
f <-function(x) unlist(strsplit(x['Protein'],'[|]'))[2]
lbfpep<-apply(lncexp,1,f)

load('TMT_exp.RData')
lncexp<-exparray[grep('sPep',exparray$Protein),]
f<-function(x) sum(!is.na(x))
rownum2<-apply(lncexp,1,f)
tmtpep<-as.character(lncexp$Protein)

missrate<-data.frame(rate=c(rownum/224,rownum2/330),
                     group=c(rep('Jiang et al 2019',16),rep('Gao et al 2019',90)))

pdf("./Figure3/Detected.pdf",width = 6,height = 4)
ggplot(missrate, aes(rate, after_stat(density), colour = group)) +
  geom_density(adjust = 2, size=1)+     
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(size=12))+
  theme(aspect.ratio=0.7)+
  scale_color_manual(values = met.brewer("Derain", 2))+
  xlab('Detected ratio')+ylab('Density')
dev.off()

#venn其实包含在sankey了

# library(ggvenn)
# library(MetBrewer)
# 
# venn<-list('Jiang et al 2019'=unique(lbfpep),
#            'Gao et al 2019'=unique(tmtpep))
# 
# pdf("./Figure3/VENN.pdf",width = 4,height = 4)
# ggvenn(venn,fill_color = met.brewer('Derain',n=2),
#        stroke_size = 0.5, set_name_size = 4)
# dev.off()

#################
write.csv(cache,file = 'alllnc_inf.csv')
write.csv(codeall[codeall$ORF_type=='annotated',],file = 'allcds_inf.csv')

load('getmotif.RData')
kozak<-getmotif[!is.na(getmotif$kozak),]

library(ggseqlogo)
library(ggplot2)

p<-list()

p[[1]]<-
ggseqlogo( kozak$kozak[kozak$group=='Ribo-seq'] )+
  coord_cartesian(ylim=c(0,0.2))
p[[2]]<-
ggseqlogo( kozak$kozak[kozak$group=='MS'] )+
  coord_cartesian(ylim=c(0,0.2))
p[[3]]<-
ggseqlogo( kozak$kozak[kozak$group=='CDS'] )+
  coord_cartesian(ylim=c(0,0.3))

library(patchwork)
pdf("./Figure3/kozak.pdf",width = 8,height = 8)
wrap_plots(p,nrow=3, guides="collect") 
dev.off()





