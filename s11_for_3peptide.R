#############################
#Ribo-seq visualization, MS visualization, single-gene GSEA
#PPIAP79 lncORF 1,sPep_25500
#POTEKP lncORF 2,sPep_25465
#HNRNPA1P36 lncORF 1,sPep_4148

#ENSG00000214853.5_85000127_84999849_92 ， SIYQEKFDDENFILK HTGPGILSMANAGPGTNVSQFFICTAK
#ENSG00000204434.6_131626007_131627266_419， 
#CPEALFQPCFLGMESCGIHKTTFNSIVKSDVDIR,
#LCYVALDSEQEMAMAASSSSVEKSYELPDGQVITIGNER
#ENSG00000231942.3_81807772_81808599_275, GFAFVTFHDHDSMDKTVIQK,HYTMNGHNCEVR

match$sequence[match$ORFname=='PPIAP79 lncORF 1']
match$sequence[match$ORFname=='POTEKP lncORF 2']
match$sequence[match$ORFname=='HNRNPA1P36 lncORF 1']

load('MS_result/TMT_quant.RData')
data[1:5,]

data<-data[!duplicated(data$Peptide),]

pote<-read.csv('Figure5/sPep_25465psm.csv',header = F)
hnrn<-read.csv('Figure5/sPep_4148psm.csv',header = F)

#coverage

library(ggplot2)

######

ppiaseq<-match$sequence[match$ORFname=='PPIAP79 lncORF 1']
ppiapsm<-read.csv('Figure5/sPep_25500psm.csv',header = F)

long_string <- ppiaseq

short_strings <- unique(ppiapsm$V3)

special_strings <- c("HTGPGILSMANAGPGTNVSQFFICTAK",'SIYQEKFDDENFILK')

coverage<-data.frame()

df <- data.frame(
  position = 1:nchar(long_string),
  char = strsplit(long_string, "")[[1]],
  covered = rep(0, nchar(long_string)),
  special = rep(0, nchar(long_string)),
  lncpep = c('PPIAP79 lncORF 1')
)

for (short_string in short_strings) {
  start_pos <- gregexpr(short_string, long_string)[[1]]
  for (pos in start_pos) {
    if (pos > 0) { 
      end_pos <- pos + nchar(short_string) - 1
      df$covered[pos:end_pos] <- 1 
    }
  }
}
for (special_string in special_strings) {
  special_pos <- gregexpr(special_string, long_string)[[1]]
  for (pos in special_pos) {
    if (pos > 0) {
      end_pos <- pos + nchar(special_string) - 1
      df$special[pos:end_pos] <- 1
    }
  }
}

coverage<-rbind(coverage,df)

#

ppiaseq<-match$sequence[match$ORFname=='POTEKP lncORF 2']
ppiapsm<-read.csv('Figure5/sPep_25465psm.csv',header = F)

long_string <- ppiaseq

short_strings <- unique(ppiapsm$V3)

special_strings <- c('CPEALFQPCFLGMESCGIHKTTFNSIVK','LCYVALDSEQEMAMAASSSSVEK',
                    'CPEALFQPCFLGMESCGIHKTTFNSIVKSDVDIR','LCYVALDSEQEMAMAASSSSVEKSYELPDGQVITIGNER')

df <- data.frame(
  position = 1:nchar(long_string),
  char = strsplit(long_string, "")[[1]],
  covered = rep(0, nchar(long_string)),
  special = rep(0, nchar(long_string)),
  lncpep = c('POTEKP lncORF 2')
)


for (short_string in short_strings) {
  start_pos <- gregexpr(short_string, long_string)[[1]]
  for (pos in start_pos) {
    if (pos > 0) { 
      end_pos <- pos + nchar(short_string) - 1
      df$covered[pos:end_pos] <- 1 
    }
  }
}
for (special_string in special_strings) {
  special_pos <- gregexpr(special_string, long_string)[[1]]
  for (pos in special_pos) {
    if (pos > 0) {
      end_pos <- pos + nchar(special_string) - 1
      df$special[pos:end_pos] <- 1
    }
  }
}
coverage<-rbind(coverage,df)

#

ppiaseq<-match$sequence[match$ORFname=='HNRNPA1P36 lncORF 1']
ppiapsm<-read.csv('Figure5/sPep_4148psm.csv',header = F)

long_string <- ppiaseq

short_strings <- unique(ppiapsm$V3)

special_strings <- c('GFAFVTFHDHDSMDKTVIQK','GFAFVTFHDHDSMDK','HYTMNGHNCEVR')

df <- data.frame(
  position = 1:nchar(long_string),
  char = strsplit(long_string, "")[[1]],
  covered = rep(0, nchar(long_string)),
  special = rep(0, nchar(long_string)),
  lncpep = c('HNRNPA1P36 lncORF 1')
)


for (short_string in short_strings) {
  start_pos <- gregexpr(short_string, long_string)[[1]]
  for (pos in start_pos) {
    if (pos > 0) { 
      end_pos <- pos + nchar(short_string) - 1
      df$covered[pos:end_pos] <- 1 
    }
  }
}
for (special_string in special_strings) {
  special_pos <- gregexpr(special_string, long_string)[[1]]
  for (pos in special_pos) {
    if (pos > 0) {
      end_pos <- pos + nchar(special_string) - 1
      df$special[pos:end_pos] <- 1
    }
  }
}
coverage<-rbind(coverage,df)


pdf("./Figure5/3pep_covered.pdf",width = 12,height = 4)
ggplot(coverage, aes(x = position, y = 1, fill = interaction(covered, special))) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("0.0" = "lightblue",     
                               "1.0" = "purple",          
                               "1.1" = "red"),     
                    labels = c("Not Covered", "Peptide covered", "Unique peptide covered"),
                    name = "Coverage Status") +
  theme_minimal() +
  labs(x = "Position", y = "", title = "Coverage of Substrings in Long String with Special Characters") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~ lncpep,scales = "free_x") +
  theme(aspect.ratio=0.1)
dev.off()

###########################################

library(tidyr)

load('TMT_exp.RData')

protein<-as.data.frame(exparray);rownames(protein)<-protein$Protein;protein<-protein[,-1]

anno<-read.table('./clinic/captac_anno.txt',header = T)
long_anno <- pivot_longer(anno, cols = c(X127N, X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131N,X131C), names_to = "Variable", values_to = "Value")
anno<-long_anno[match(colnames(protein),long_anno$Value),]
protein <- limma::removeBatchEffect(protein, batch=anno$AnalyticalSample, 
                                     group=substr(anno$Value,1,1))


#PPIAP79 lncORF 1,sPep_25500
#POTEKP lncORF 2,sPep_25465
#HNRNPA1P36 lncORF 1,sPep_4148

cor<-t(protein)
cor<-cor[grep('P',rownames(cor)),]
PPIAP79<-cor[!is.na(cor[,'sPep_25500']),]
f<-function(x) sum(is.na(x))
colnum<-apply(PPIAP79,2,f)
PPIAP79<-PPIAP79[,colnum < 1]  #nrow(PPIAP79)*0.1
PPIAP79[is.na(PPIAP79)]<-0

cortest<-t(cor(PPIAP79[,'sPep_25500'],PPIAP79,method=c('spearman')))
cortest<-as.data.frame(cortest)
cortest$uniprotname<-rownames(cortest)

library('org.Hs.eg.db')
library(clusterProfiler)

genelist<-bitr(cortest$uniprotname, fromType = "UNIPROT",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)
cortest<-cortest[cortest$uniprotname %in% genelist$UNIPROT,]
genelist<-genelist[match(cortest$uniprotname,genelist$UNIPROT),]
cortest$entre<-genelist$ENTREZID
cortest<-cortest[order(cortest$V1,decreasing = T),]
cortest<-cortest[!duplicated(cortest$entre),]

gsea0<-as.numeric(cortest$V1)
names(gsea0)<-cortest$entre
edo2 <- clusterProfiler::gseGO(gsea0,OrgDb= 'org.Hs.eg.db',ont='ALL',pvalueCutoff = 0.1,eps=0)

p<-list()

p[[1]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='BP'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:BP")

p[[2]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='MF'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:MF")

p[[3]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='CC'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:CC")

# library(patchwork)
# pdf("./Figure5/PPIAP79_lncORF_GSEA_onlyP.pdf",width = 22,height = 6)
# wrap_plots(p,nrow=1) 
# dev.off()

#
cor<-t(protein)
cor<-cor[grep('P',rownames(cor)),]
PPIAP79<-cor[!is.na(cor[,'sPep_25465']),]
f<-function(x) sum(is.na(x))
colnum<-apply(PPIAP79,2,f)
PPIAP79<-PPIAP79[,colnum < 1]  #nrow(PPIAP79)*0.1
PPIAP79[is.na(PPIAP79)]<-0

cortest<-t(cor(PPIAP79[,'sPep_25465'],PPIAP79,method='spearman'))
cortest<-as.data.frame(cortest)
cortest$uniprotname<-rownames(cortest)
genelist<-bitr(cortest$uniprotname, fromType = "UNIPROT",
               toType = c( "ENTREZID",'SYMBOL'),
               OrgDb = org.Hs.eg.db)
cortest<-cortest[cortest$uniprotname %in% genelist$UNIPROT,]
genelist<-genelist[match(cortest$uniprotname,genelist$UNIPROT),]
cortest$entre<-genelist$ENTREZID
cortest$SYMBOL<-genelist$SYMBOL
cortest<-cortest[order(cortest$V1,decreasing = T),]
cortest<-cortest[!duplicated(cortest$entre),]
gsea0<-as.numeric(cortest$V1)
names(gsea0)<-cortest$entre
edo2 <- clusterProfiler::gseGO(gsea0,OrgDb= 'org.Hs.eg.db',ont='ALL',pvalueCutoff = 0.05,eps=0)

# p<-list()

p[[4]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='BP'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:BP")

p[[5]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='MF'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:MF")

p[[6]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='CC'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:CC")

# library(patchwork)
# pdf("./Figure5/POTEKP_lncORF2_GSEA_onlyP.pdf",width = 22,height = 6)
# wrap_plots(p,nrow=1) 
# dev.off()

##HNRNPA1P36 lncORF 1,sPep_4148
cor<-t(protein)
cor<-cor[grep('P',rownames(cor)),]
PPIAP79<-cor[!is.na(cor[,'sPep_4148']),]
f<-function(x) sum(is.na(x))
colnum<-apply(PPIAP79,2,f)
PPIAP79<-PPIAP79[,colnum < nrow(PPIAP79)*0.1] 
PPIAP79[is.na(PPIAP79)]<-0

cortest<-t(cor(PPIAP79[,'sPep_4148'],PPIAP79,method=c('spearman')))
cortest<-as.data.frame(cortest)
cortest$uniprotname<-rownames(cortest)
genelist<-bitr(cortest$uniprotname, fromType = "UNIPROT",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)
cortest<-cortest[cortest$uniprotname %in% genelist$UNIPROT,]
genelist<-genelist[match(cortest$uniprotname,genelist$UNIPROT),]
cortest$entre<-genelist$ENTREZID
cortest<-cortest[order(cortest$V1,decreasing = T),]
cortest<-cortest[!duplicated(cortest$entre),]
gsea0<-as.numeric(cortest$V1)
names(gsea0)<-cortest$entre
edo2 <- clusterProfiler::gseGO(gsea0,OrgDb= 'org.Hs.eg.db',ont='ALL',pvalueCutoff = 0.1,eps=0)

# p<-list()

p[[7]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='BP'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:BP")

p[[8]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='MF'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:MF")

p[[9]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='CC'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:CC")

library(patchwork)
pdf("./Figure5/lncORF1_GSEA_onlyP.pdf",width = 22,height = 18)
wrap_plots(p,nrow=3)
dev.off()


load('LBF_protein_exp.RData')
protein<-as.data.frame(protein);rownames(protein)<-protein$Protein;protein<-protein[,-1]
cor<-t(protein)
cor<-cor[grep('P',rownames(cor)),]
PPIAP79<-cor[!is.na(cor[,'tr|sPep_25500|sPep_25500']),]
f<-function(x) sum(is.na(x))
colnum<-apply(PPIAP79,2,f)
PPIAP79<-PPIAP79[,colnum < 1]  #nrow(PPIAP79)*0.1
PPIAP79[is.na(PPIAP79)]<-0
cortest<-t(cor(PPIAP79[,'tr|sPep_25500|sPep_25500'],PPIAP79))
cortest<-as.data.frame(cortest)
cortest$ID<-rownames(cortest)

f <-function(x) unlist(strsplit(x['ID'],'[|]'))[2]
cortest$uniprotname<-apply(cortest,1,f)

genelist<-bitr(cortest$uniprotname, fromType = "UNIPROT",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)
cortest<-cortest[cortest$uniprotname %in% genelist$UNIPROT,]
genelist<-genelist[match(cortest$uniprotname,genelist$UNIPROT),]
cortest$entre<-genelist$ENTREZID
cortest<-cortest[order(cortest$V1,decreasing = T),]
cortest<-cortest[!duplicated(cortest$entre),]

gsea0<-as.numeric(cortest$V1)
names(gsea0)<-cortest$entre
edo2 <- clusterProfiler::gseGO(gsea0,OrgDb= 'org.Hs.eg.db',ont='ALL',pvalueCutoff = 0.1,eps=0)

p<-list()

p[[1]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='BP'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:BP")

p[[2]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='MF'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:MF")

p[[3]]<-
  GseaVis::gseaNb(object = edo2,
                  geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='CC'][1:5],
                  subPlot = 2,
                  termWidth = 35,
                  # legend.position = c(0.8,0.8),
                  addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
                  pvalX = 0.05,pvalY = 0.05,pDigit=3)+
  labs(title="GO:CC")

library(patchwork)
pdf("./Figure5/PPIAP79_lncORF_GSEA_LBF.pdf",width = 22,height = 6)
wrap_plots(p,nrow=1) 
dev.off()


#######################
load('TMT_exp.RData')

protein<-as.data.frame(exparray);rownames(protein)<-protein$Protein;protein<-protein[,-1]

anno<-read.table('./clinic/captac_anno.txt',header = T)
long_anno <- pivot_longer(anno, cols = c(X127N, X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131N,X131C), names_to = "Variable", values_to = "Value")
anno<-long_anno[match(colnames(protein),long_anno$Value),]
protein <- limma::removeBatchEffect(protein, batch=anno$AnalyticalSample, 
                                    group=substr(anno$Value,1,1))
#PPIAP79 lncORF 1,sPep_25500
#POTEKP lncORF 2,sPep_25465
#HNRNPA1P36 lncORF 1,sPep_4148
protein<-as.data.frame(protein)
protein<-protein[,grep('P',colnames(protein))]
protein<-protein[rownames(protein) %in% c('sPep_25500','sPep_25465','sPep_4148'),]
cor<-as.data.frame( t(protein))
colnames(cor) <-c('POTEKP lncORF 2','PPIAP79 lncORF 1','HNRNPA1P36 lncORF 1')
cor<-cor[!is.na(cor$`PPIAP79 lncORF 1`) & !is.na(cor$`HNRNPA1P36 lncORF 1`),]

rvalue<-cor(cor)
pvalue<-ggcorrplot::cor_pmat(cor)

library('corrplot')

pdf("./Figure5/3pep_corplot.pdf",width = 6,height = 6)
corrplot(corr =rvalue, p.mat = pvalue,method = "circle",
         type="upper",col=rev(COL2('RdBu', 200)),
         tl.pos="lt", tl.cex=1, tl.col="black",tl.srt = 45,tl.offset=0.5,
         insig="label_sig",sig.level = c(.001, .01, .05),
         pch.cex = 0.8,pch.col = "black") 
corrplot(corr = rvalue, method = "number",
         type="lower",add=TRUE,col=rev(COL2('RdBu', 200)),
         tl.pos = "n",cl.pos = "n",diag=FALSE,
         number.digits = 3,number.cex = 1,number.font = NULL)
dev.off() 

