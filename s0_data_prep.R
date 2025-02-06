####################################################
# To filter redundant transcripts we compared NONCODE annotations with GENCODE annotations using GffCompare
# Filtering NONCODE transcripts

noncode<-data.table::fread('./data_prep/NONCODEv6_human_hg38_lncRNA.gtf')
noncode$V9 <- gsub("\r", "", noncode$V9)
write.table(noncode,file = './data_prep/NONCODE_adj.gtf',sep = '\t',row.names = F,col.names = F,quote = F)

anno<-data.table::fread('data_prep/NONCODE_compare.NONCODE_adj.gtf.refmap')
anno<-data.table::fread('data_prep/NONCODE_compare.NONCODE_adj.gtf.tmap')

library(rtracklayer)

GENCODE <- as.data.frame(import.gff('data_prep/gencode.v46.annotation.gtf',format = 'gtf')) 
GENCODE<-GENCODE[GENCODE$type=='gene',]

match<-GENCODE[match(anno$ref_gene_id,GENCODE$gene_id),]
anno$ref_type<-match$gene_type
anno$ref_name<-match$gene_name

filter<-anno[anno$ref_type=='protein_coding',]

pool1<-filter[filter$class_code=='=',]
pool2<-anno[anno$class_code=='s',]
pool3<-filter[filter$num_exons=='1',]
filterlist<-rbind(pool1,pool2,pool3)

NONCODE <- as.data.frame(import.gff('data_prep/NONCODE_adj.gtf',format = 'gtf')) 
NONCODE<-NONCODE[!(NONCODE$transcript_id %in% filterlist$qry_id),]


#Adjusting the NONCODE annotation format

library(progress)
pb <- progress_bar$new(total = length(unique(NONCODE$gene_id)))

adjnoncode<-data.frame()

for (i in unique(NONCODE$gene_id)){
  
  pb$tick()
  
  inf<-NONCODE[NONCODE$gene_id==i,]
  gene<-data.frame(seqnames=c(unique(inf$seqnames)),start=c(min(inf$start)),end=c(max(inf$end)),
                   width=c(max(inf$end)-min(inf$start)+1),strand=c(unique(inf$strand)),
                   source=c('Cufflinks'),type=c('gene'),score=c(0),phase=c(NA),
                   gene_id=c(i),transcript_id=c(NA),
                   FPKM=c(0),exon_number=c(NA))
  gene<-rbind(gene,inf)
  adjnoncode<-rbind(adjnoncode,gene)
}

save(adjnoncode,file = 'data_prep/adjusted_NONCODE.RData')

#merge GENCODE and NONCODE annotation

GENCODE <- as.data.frame(import.gff('data_prep/gencode.v46.annotation.gtf',format = 'gtf')) 
GENCODE[1:5,]

merge<-data.frame(seqnames=c(GENCODE$seqnames,adjnoncode$seqnames),
                  start=c(GENCODE$start,adjnoncode$start),end=c(GENCODE$end,adjnoncode$end),
                  width=c(GENCODE$width,adjnoncode$width),strand=c(GENCODE$strand,adjnoncode$strand),
                  source=c(GENCODE$source,adjnoncode$source),type=c(GENCODE$type,adjnoncode$type),
                  score=c(NA),phase=c(NA),gene_id=c(GENCODE$gene_id,adjnoncode$gene_id),
                  gene_type=c(GENCODE$gene_type,rep(NA,nrow(adjnoncode))),gene_name=c(GENCODE$gene_name,adjnoncode$gene_id), 
                  level=c(GENCODE$level,rep(NA,nrow(adjnoncode))),tag=c(GENCODE$tag,rep(NA,nrow(adjnoncode))),
                  transcript_id=c(GENCODE$transcript_id,adjnoncode$transcript_id),transcript_type=c(GENCODE$transcript_type,rep(NA,nrow(adjnoncode))),
                  transcript_name=c(GENCODE$transcript_name,adjnoncode$transcript_id), 
                  exon_number=c(GENCODE$exon_number,adjnoncode$exon_number))

# export.gff(merge,'./data_prep/merge.gtf',format='gtf')


gtffile<-data.table::fread('./data_prep/merge.gtf',skip = 3)
gtffile$V2<-c(as.character(GENCODE$source),adjnoncode$source)
gtffile$V3<-c(as.character(GENCODE$type),adjnoncode$type)

gtffile<-gtffile[gtffile$V7!='.',]
table(gtffile$V7)

gtffile<-gtffile[-grep('alt',gtffile$V1),]
gtffile<-gtffile[-grep('random',gtffile$V1),]
gtffile<-gtffile[-grep('chrUn',gtffile$V1),]

write.table(gtffile,'data_prep/GENCODE_NONCODE_merged.gtf',sep = '\t',col.names = F,row.names = F,quote = F)

#####################################
#Obtaining uniprot-annotated lncRNA encoded peptide

library(Biostrings)

unport<-readAAStringSet('./data_prep/uniprot_human.fasta')
uniprot<-unport@ranges@NAMES
uniprot

f <-function(x) unlist(strsplit(x,'[|]'))[3]
uniprot<- sapply(uniprot, f)
uniprot<-as.character(uniprot)

f <-function(x) unlist(strsplit(x,'[=]'))[4]
uniprot<- sapply(uniprot, f)
uniprot<-as.character(uniprot)

f <-function(x) unlist(strsplit(x,'[ ]'))[1]
uniprot<- sapply(uniprot, f)
uniprot<-as.character(uniprot)

library(rtracklayer)

gtf <- as.data.frame(import.gff('../biodata/gencode.v46.annotation.gtf',format = 'gtf') ) #gencode annotation
gtf<-gtf[gtf$type=='gene',]

cache<-gtf[match(uniprot,gtf$gene_name),]

uniprot<-data.frame(title=c(unport@ranges@NAMES),
                    seq=gsub("\\.","",paste(unport)),
                    genename=c(uniprot),
                    type=c(cache$gene_type),
                    id=c(cache$gene_id))


filter<-uniprot[uniprot$type %in% c('lncRNA','processed_pseudogene','transcribed_processed_pseudogene','transcribed_unitary_pseudogene',
                                    'transcribed_unprocessed_pseudogene','translated_processed_pseudogene','unitary_pseudogene','unprocessed_pseudogene'),]

#list of uniprot annotation lncRNA-peptide

save(filter,file = 'uniprot_lncPep.RData')

unport<-readAAStringSet('./data_prep/uniprot_human.fasta')
unport<-unport[!(unport@ranges@NAMES %in% filter$title)]

writeXStringSet(unport, 'uniprot_human_nonolnc.fasta', append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")




