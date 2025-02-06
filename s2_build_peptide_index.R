
##########################
#Comparison of ORF results for predictions from different tools

tish<-rbind(data.table::fread('./demo/ribotish_liver_sep.txt'),#ribotish result
            data.table::fread('./demo/ribotish_liver_longest.txt'))
tish<-tish[!duplicated(paste0(tish$Gid,'_',tish$GenomePos,'_',tish$AALen,'_',tish$Seq)),]

code<-rbind(data.table::fread('demo/ribocode_all/ribocode_allseq.txt'), #Ribocode result
              data.table::fread('demo/ribocode_all/ribocode_allseq_longest.txt'))

code<-code[!duplicated(paste0(code$ORF_ID,'_',code$AAseq)),]

library(Biostrings)
seq <- DNAStringSet(tish$Seq)
aaseq <- translate(seq)
tish$aaseq<-gsub("\\.","",paste(aaseq))

tricer<-data.table::fread('./demo/merge.txt')  #ribotricer result
tricer<-tricer[,-18]
tricer<-tricer[order(tricer$phase_score),]
tricer<-tricer[!duplicated(tricer$ORF_ID),]
tricer<-tricer[tricer$phase_score!='phase_score',]

tricer<-tricer[!duplicated(paste0(tricer$gene_id,'_',tricer$length,'_',tricer$phase_score)),]

########################
#
table(code$gene_type)
table(tish$GeneType)
table(tricer$gene_type)

gencode<-data.table::fread('./Ribo-seq_ORFs.bed') #GENCODE Ribo-ORF
choth<-read.csv('Chothani_2022.csv') #Chothani 2022 paper
choth<-choth[-grep('CATG',choth$Gene.ID),]

library(ggvenn)
library(MetBrewer)

venn<-list(Ribocode=unique(substr(code$gene_id[code$gene_type=='lncRNA'],1,15)),
           RiboTish=unique(substr(tish$Gid[tish$GeneType=='lncRNA'],1,15)),
           RiboTricer=unique(substr(tricer$gene_id[tricer$gene_type=='lncRNA'],1,15)), 
           GENCODE=unique(gencode$V18[gencode$V20=='lncRNA']),
           Chothani=unique(choth$Gene.ID[choth$Gene.type=='lincRNA']),
           UniProt=unique(substr(filter$id[filter$type=='lncRNA'],1,15)))

p<-list()

p[[1]]<-ggvenn(venn,fill_color = met.brewer('Hiroshige',n=5),
          stroke_size = 0.5, set_name_size = 4)

p[[2]]<-ggvenn(venn[-4],fill_color = met.brewer('Hiroshige',n=5),
       stroke_size = 0.5, set_name_size = 4) 

p[[3]]<-ggvenn(venn[-c(4:5)],fill_color = met.brewer('Hiroshige',n=5),
               stroke_size = 0.5, set_name_size = 4)

library(patchwork)
pdf("./qc/venn.pdf",width = 15,height = 4)
wrap_plots(p,nrow=1, guides="collect") 
dev.off()

library(vegan)
library(pheatmap)

venn_matrix <- sapply(venn, function(x) {
  unique_genes <- unique(unlist(venn))
  as.numeric(unique_genes %in% x)
})

jaccard_dist <- vegdist(t(venn_matrix), method = "jaccard")

jaccard_corr <- 1 - as.matrix(jaccard_dist)

library('corrplot')
# diag(jaccard_corr) <- NA

pdf("./qc/corrplot.pdf",width = 7,height = 7)
corrplot(jaccard_corr, method = "color", 
         addCoef.col = "black", 
         tl.col = "black",  diag = FALSE,    
         tl.srt = 45) 
dev.off()

#################################################
#########################
#In this project we focus on lncRNA and pseudogene.

# #tricer
# 
# f <-function(x) unlist(strsplit(x['ORF_ID'],'_'))[2]
# tricer$start<-apply(tricer,1,f)
# f <-function(x) unlist(strsplit(x['ORF_ID'],'_'))[3]
# tricer$end<-apply(tricer,1,f)
# 
# pos<-tricer[tricer$strand=='+',]
# neg<-tricer[tricer$strand=='-',]
# 
# f <-function(x) unlist(strsplit(x['ORF_ID'],'_'))[3]
# neg$start<-apply(neg,1,f)
# f <-function(x) unlist(strsplit(x['ORF_ID'],'_'))[2]
# neg$end<-apply(neg,1,f)
# 
# collapse<-rbind(pos,neg)
# collapse$adjid<-paste0(collapse$transcript_id,'_',collapse$end)
# 
# aaseq<-data.table::fread('./demo/ntsequence.txt')
# aaseq<-aaseq[match(collapse$ORF_ID,aaseq$ORF_ID),]
# collapse$aaseq<-aaseq$sequence

#########################
# 
# codeall<-rbind(data.table::fread('demo/ribocode_all/ribocode_allseq.txt'),
#                data.table::fread('demo/ribocode_all/ribocode_allseq_longest.txt'))
# 
# codeall<-codeall[!duplicated(paste0(codeall$ORF_ID,'_',codeall$ORF_length)),]

#合并tish

codeall<-code

f <-function(x) unlist(strsplit(x['GenomePos'],'[:]'))[2]
tish$position<-apply(tish,1,f)
f <-function(x) unlist(strsplit(x['position'],'-'))[1]
tish$ORF_gstart<-apply(tish,1,f)
f <-function(x) unlist(strsplit(x['position'],'-'))[2]
tish$ORF_gstop<-apply(tish,1,f)

library(dplyr)

codeall <- codeall %>%
  mutate(
    temp = ORF_gstart, 
    ORF_gstart = if_else(strand == "-", ORF_gstop-1, ORF_gstart-1),  
    ORF_gstop = if_else(strand == "-", temp, ORF_gstop)  
  ) %>%
  select(-temp)

codeall$new_name<-paste0(codeall$transcript_id,'_',codeall$ORF_gstart,'_',codeall$ORF_gstop)
codealnc<-codeall[grepl('pseudo|lncRNA|None',codeall$gene_type),]
codealnc<-codealnc[!duplicated(codealnc$AAseq),]

tish$new_name<-paste0(tish$Tid,'_',tish$ORF_gstart,'_',tish$ORF_gstop)
tishlnc<-rbind(tish[grep('pseudo|lncRNA',tish$GeneType),],
               tish[tish$GeneType=='',])
tishlnc<-tishlnc[!duplicated(tishlnc$aaseq),]

merge<-data.frame(name=c(codealnc$new_name,tishlnc$new_name),
                  genename=c(codealnc$gene_name,tishlnc$Symbol),
                  geneid=c(codealnc$gene_id,tishlnc$Gid),
                  genetype=c(codealnc$gene_type,tishlnc$GeneType),
                  transcript_id=c(codealnc$transcript_id,tishlnc$Tid),
                  aaseq=c(codealnc$AAseq,sapply(tishlnc$aaseq, function(x) substr(x, 1, nchar(x) - 1)))
                  ) 
merge<-merge[!duplicated(merge$name),]

match<-tishlnc[match(merge$name,tishlnc$new_name),]
match2<-codealnc[match(merge$name,codealnc$new_name),]

merge$tishp<-match$RiboPvalue
merge$codep<-match2$pval_combined

merge$code<-ifelse(merge$codep<0.05,1,0)
merge$tish<-ifelse(merge$tishp<0.05,1,0)

#######
#Merge tricer
f <-function(x) unlist(strsplit(x['ORF_ID'],'_'))[2]
tricer$start<-apply(tricer,1,f)
f <-function(x) unlist(strsplit(x['ORF_ID'],'_'))[3]
tricer$end<-apply(tricer,1,f)

tricer$start<-as.numeric(tricer$start);tricer$end<-as.numeric(tricer$end)

tricer <- tricer %>%
  mutate(
    start = if_else(strand == "-", start-4, start-1),  
    end = if_else(strand == "-", end, end+3)  ) 

tricer$new_name<-paste0(tricer$transcript_id,'_',tricer$start,'_',tricer$end)

match<-tricer[match(merge$name,tricer$new_name),]
merge$tricer_phase_score<-match$phase_score
merge$read_density<-match$read_density
merge$tricer<-ifelse(merge$tricer_phase_score>0.428,1,0)

merge$code[is.na(merge$code)]<-c(0)
merge$tish[is.na(merge$tish)]<-c(0)
merge$tricer[is.na(merge$tricer)]<-c(0)

merge$save<-rowSums(merge[,c(9,10,13)])

result<-merge[merge$save>1,]

result_lnc<-result[!(result$gene_type %in% c('artifact','protein_coding','TEC','IG_C_gene','TR_C_gene')),]

#This is the result before collapse

collapse<-rbind(data.table::fread('./demo/ribocode_all/ribocode_allseq_collapsed.txt'),
                data.table::fread('./demo/ribocode_all/ribocode_allseq_longest_collapsed.txt'))
test<-collapse[collapse$ORF_ID %in% result_lnc$ORF_ID,] 

collapse<-result_lnc[result_lnc$ORF_ID %in% collapse$ORF_ID,] 

test<-collapse[!duplicated(collapse$AAseq),] 


#####################################
#writer aaindex

library('seqinr')
seq<-merge[,c('ORF_ID','AAseq')]

library(TNSMD)

other<-data.frame(genename=c(merge$gene_name),
                  orftype=c('lncORF'),
                  pepseq=c(merge$AAseq),
                  source=c('Ribo'))

genelnc<-gencode[gencode$V20 %in% c('lncRNA'),]
  
part3<-data.frame(genename=c(genelnc$V19,filter$genename),
         orftype=c('lncORF'),
         pepseq=c(genelnc$V22,filter$seq),
         source=c(rep('GENCODE',length(genelnc$V22)),rep('Uniprot',length(filter$seq))))  #GENCODE and uniprot-derived peptides


myindex<-rbind(other,part3)
myindex<-myindex[!duplicated(paste0(myindex$genename,'_',myindex$pepseq)),]
myindex<-myindex[order(myindex$genename),]

write.csv(myindex,'./writeindex.csv')

myindex<-generate_index('./writeindex.csv','other',0,'sPep')

# collapse$index_name<-other$name

characters<-myindex$sequence
list<-as.list(characters)
seqinr::write.fasta(sequences=list, names=myindex$name,
                    file.out='relncpep.fasta', open='w', nbchar=60)


#################################

##############################
#In addition to the above we construct two types of indexes, 
#one taking the results of the concatenation set 
#and the other the ORFfinder results of all transcripts with FPKM greater than 1

#index2
result<- merge[merge$save>0,]
nonover<-result
nonover<-nonover[!duplicated(paste0(nonover$gene_name,'_',nonover$AAseq)),]

library(TNSMD)

other<-data.frame(genename=c(nonover$genename),
                  orftype=c('lncORF'),
                  pepseq=c(nonover$aaseq),
                  source=c('Ribo'))

myindex<-rbind(other,part3)
myindex<-myindex[!duplicated(paste0(myindex$genename,'_',myindex$pepseq)),]
myindex<-myindex[order(myindex$genename),]

write.csv(myindex,'./writeindex_2_new.csv')

myindex<-generate_index('./writeindex_2_new.csv','other',0,'sPep')

characters<-myindex$sequence
list<-as.list(characters)
seqinr::write.fasta(sequences=list, names=myindex$name,
                    file.out='relncpep_2_Oct15.fasta', open='w', nbchar=60)

######
#index3
count<-data.table::fread('demo/feature.count') #ORFfinder results

kb <-count$Length/1000
count$whatever<-0

rpk<-count[,c(7:8)]/kb

fpkm<-t(t(rpk)/colSums(count[,c(7:8)])*10^6)
fpkm=as.data.frame(fpkm)
count$fpkm<-fpkm$`./merge.sorted.bam`
count<-count[count$fpkm>0.5,]

library(rtracklayer)

gtf <- data.frame(import.gff('../biodata/gencode.v46.annotation.gtf',format = 'gtf'))
gtf<-gtf[gtf$type=='transcript',]
gtf<-gtf[match(count$Geneid,gtf$transcript_id),]
count$genetype<-gtf$gene_type
count$genetype[is.na(count$genetype)]<-c('noncode')
count<-count[!(count$genetype %in% c('artifact','protein_coding','TEC','IG_C_gene','TR_C_gene',
                                     'vault_RNA','scaRNA','snoRNA','snRNA','misc_RNA',
                                     'IG_J_gene','Mt_rRNA','Mt_tRNA','ribozyme','rRNA_pseudogene')),]


allrna <- Biostrings::readDNAStringSet('transcripts_sequence.fa')
name<-allrna@ranges@NAMES

f <-function(x) unlist(strsplit(x,' '))[1]
name<-sapply(name,f)
name<-as.character(name)

allrna <-allrna [name %in% c(merge$transcript_id[merge$save>0])]

writeXStringSet(allrna, 'transcripts_lnc.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

allrna <- Biostrings::readAAStringSet('lncorffinder_7Oct.fasta') #会多挺多的

allran<-data.frame(name=c(allrna@ranges@NAMES),
                   seq=c(sequence = gsub("\\.","",paste(allrna))))
allran<-allran[!duplicated(allran$seq),]

f <-function(x) unlist(strsplit(x['name'],':'))[1]
allran$transid<-apply(allran,1,f)
f <-function(x) unlist(strsplit(x['transid'],'_'))[2]
allran$transid<-apply(allran,1,f)


other<-data.frame(genename=c(allran$transid),
                  orftype=c('lncORF'),
                  pepseq=c(allran$seq),
                  source=c('Ribo'))

myindex<-rbind(other,part3)
myindex<-myindex[!duplicated(paste0(myindex$genename,'_',myindex$pepseq)),]
myindex<-myindex[order(myindex$genename),]

write.csv(myindex,'./writeindex_3.csv')

myindex<-generate_index('./writeindex_3.csv','other',0,'sPep')

characters<-myindex$sequence
list<-as.list(characters)
seqinr::write.fasta(sequences=list, names=myindex$name,
                    file.out='relncpep_3.fasta', open='w', nbchar=60)

