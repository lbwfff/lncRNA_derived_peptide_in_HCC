#################################
#Record the length of reads after manual filtering

samplelist<-list.files('./qc/')
samplelist<-samplelist[grep('out_qual',samplelist)]
samplelist<-data.frame(sample=c(samplelist),
                       length=c(NA))

samplelist$length<-c('27,28,29,30','27,28,29,30','27,28,29,30','27,28,29,30','27,28,29,30',
                     '29','29','29,30','28,29','29','29,30','28','29,30','29,29','28','28,29','28,29')

#Quality control results using RiboTish

library(jsonlite)

cache1<-data.frame()
cache2<-data.frame()
cache3<-data.frame()
cache4<-data.frame()

for (i in 1:nrow(samplelist)){
  inf<-read.table(paste0('./qc/',samplelist$sample[i]),fill = T,sep = '\t')
  sel<-unlist(strsplit(samplelist$length[i],'[,]'))
  
  length <- gsub("(\\d+):", '"\\1":', inf[1,1])
  length <- fromJSON(length)
  length<-length[names(length) %in% sel]
  length <- as.data.frame(t(as.data.frame(length)))
  length$length<-substr(rownames(length),2,99)
  length$sample<-samplelist$sample[i]
  cache1<-rbind(cache1,length)
  
  start <- gsub("(\\d+):", '"\\1":', inf[2,1])
  start <- fromJSON(start)
  start<-start[names(start) %in% sel]
  start <- as.data.frame(start)
  start$Sum<-rowSums(start)
  start$sample<-samplelist$sample[i]
  start<-start[,((ncol(start)-1):ncol(start))]
  start$order<-c(1:60)-40
  cache2<-rbind(cache2,start)
  
  stop <- gsub("(\\d+):", '"\\1":', inf[3,1])
  stop <- fromJSON(stop)
  stop<-stop[names(stop) %in% sel]
  stop <- as.data.frame(stop)
  stop$Sum<-rowSums(stop)
  stop$sample<-samplelist$sample[i]
  stop<-stop[,((ncol(stop)-1):ncol(stop))]
  stop$order<-c(1:60)-40
  cache3<-rbind(cache3,stop)
  
  lengthframe <- gsub("(\\d+):", '"\\1":', inf[4,1])
  lengthframe <- fromJSON(lengthframe)
  lengthframe<-lengthframe[names(lengthframe) %in% sel]
  lengthframe <- as.data.frame(t(as.data.frame(lengthframe)))
  lengthframe$length<-substr(rownames(lengthframe),2,99)
  lengthframe <- reshape2::melt(lengthframe)
  lengthframe$sample<-samplelist$sample[i]
  cache4<-rbind(cache4,lengthframe)

}

#The bedtools intersect function was used as a statistic for the distribution of different feature reads of mRNAs.

# data.frame(group=c(rep('CDS',166634039),rep('5UTR',6466945),rep('3UTR',1555301),rep('ncRNA',2694873)),
#            whatever=c(NA))

plot<-data.frame(group=c(rep('CDS',1666340),rep('5UTR',64669),rep('3UTR',15553),rep('ncRNA',26948)),
                 whatever=c(NA))


library(ggpie) 
library(MetBrewer)

pdf("./FigureS1/reads_distr.pdf",width =5,height = 5)
ggpie(data = plot, group_key = "group", count_type = "full",
      label_info = "ratio", label_type = "horizon",
      label_size = 4, label_pos = "out",label_gap=0.001 )+
  scale_fill_manual(values=met.brewer('Hiroshige',n=4))+
  labs(fill='Feature type')
dev.off()

rm(plot)

########################################
#plot
p<-list()

library('ggplot2')
library('MetBrewer')
library('reshape2')

plot<-cache1

p[[1]]<-
ggplot(plot, aes(x=length, weight = V1, fill = sample))+
  geom_bar( position = "stack")+ 
  theme_minimal()+theme_classic()+
  labs(y=NULL,x='Length')+
  theme_classic(base_size = 18)+  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))+
  theme(legend.text = element_text(size=8))+
  scale_fill_manual(values=met.brewer('Hiroshige',n=17))+
  theme(panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line()) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),  
    labels = scales::label_number(scale = 1e-6, suffix = " M") )+
  theme(legend.position = "none")+
  theme(aspect.ratio=2)

#
plot <- aggregate(value ~ variable, data = cache4, sum)
plot$group<-c('Frame 1','Frame 2','Frame 3')
plot$value<-plot$value/sum(plot$value)

p[[2]]<-
  ggplot(plot, aes(x=group, weight = value, fill = group))+
  geom_bar()+ 
  theme_minimal()+theme_classic()+
  labs(x=NULL,y=NULL)+
  theme_classic(base_size = 18)+  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))+
  theme(legend.text = element_text(size=8))+
  theme(panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line()) +
  scale_fill_manual(values=met.brewer('Hiroshige',n=4))+
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),  
    labels = scales::label_number() )+
  theme(aspect.ratio=2)

#

plot <- dcast(cache2, order ~ sample, value.var = "Sum")

plot<-data.frame(order=plot$order,
                 freq=rowSums(plot[,-1]))
plot$group<-rep(c('Frame 3','Frame 1','Frame 2'),20)

p[[3]]<-
ggplot(plot, aes(x=order, weight = freq, fill = group))+
  geom_bar()+ 
  theme_minimal()+theme_classic()+
  labs(x='Dist. from start codon (nt)',y=NULL)+
  theme_classic(base_size = 18)+  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))+
  theme(legend.text = element_text(size=8))+
  theme(panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line()) +
  scale_fill_manual(values=met.brewer('Hiroshige',n=4))+
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),  
    labels = scales::label_number(scale = 1e-3, suffix = " K") ) 

plot <- dcast(cache3, order ~ sample, value.var = "Sum")

plot<-data.frame(order=plot$order,
                 freq=rowSums(plot[,-1]))
plot$group<-rep(c('Frame 3','Frame 1','Frame 2'),20)

p[[4]]<-
ggplot(plot, aes(x=order, weight = freq, fill = group))+
  geom_bar()+ 
  theme_minimal()+theme_classic()+
  labs(x='Dist. from stop codon (nt)',y=NULL)+
  theme_classic(base_size = 18)+  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=0.9))+
  theme(panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line()) +
  theme(legend.text = element_text(size=8))+
  scale_fill_manual(values=met.brewer('Hiroshige',n=4))+
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),  
    labels = scales::label_number(scale = 1e-3, suffix = " K") ) 


library(patchwork)
pdf("./demo/qc/qc1.pdf",width = 16,height = 6)
wrap_plots(p,nrow=1, guides="collect") 
dev.off()


##########################################################

allgene<-data.table::fread('reads_dis/ALL_gene.bed')

library(rtracklayer)
gtf <-  as.data.frame(import.gff('../biodata/gencode.v46.annotation.gtf',format = 'gtf'))
gtf<-gtf[gtf$type=='transcript',]

table(substr(allgene$V4,1,15) %in% substr(gtf$transcript_id,1,15))

allnc<-allgene[substr(allgene$V4,1,15) %in% substr(gtf$transcript_id[gtf$gene_type!='protein_coding'],1,15) ]
almrna<-allgene[substr(allgene$V4,1,15) %in% substr(gtf$transcript_id[gtf$gene_type=='protein_coding'],1,15) ]

write.table(allnc,file = 'reads_dis/alllnc.bed',sep = '\t',quote = F,row.names = F,col.names = F)
write.table(almrna,file = 'reads_dis/allmrna.bed',sep = '\t',quote = F,row.names = F,col.names = F)

