#

lbf<-data.table::fread('/scratch/lb4489/project/hcc/iprox/fragpipe/allrun/msstats_filter.csv')

prot_psm<-lbf[grep('HUMAN',lbf$ProteinName),]
lnc_psm<-lbf[-grep('HUMAN',lbf$ProteinName),]

lnc_psm<-lnc_psm[-grep('HEVBR|BOVIN|CHICK|PIG|RABIT|SHEEP',lnc_psm$ProteinName),]

pepfilter<-data.frame(pep=c(unique(lnc_psm$PeptideSequence)),
                      match=c(NA),
                      sapmatch=c(NA))
pepfilter$pep<-gsub("\\[.*?\\]|[a-z]", "", pepfilter$pep)

pattern2 <- Biostrings::readAAStringSet('/scratch/lb4489/project/hcc/R_code/all_protein.fasta')

for (i in 1:nrow(pepfilter)){
  
  pattern1 <- Biostrings::AAString(pepfilter$pep[i])
  
  match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=0)
  match2<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)
  
  pepfilter$match[i]<-sum(match)
  pepfilter$sapmatch[i]<-sum(match2) 
  
}

lnc_psm_filter<-lnc_psm[gsub("\\[.*?\\]|[a-z]", "", lnc_psm$PeptideSequence) %in% pepfilter$pep[pepfilter$match==0 & pepfilter$sapmatch==0],]
length(unique(lnc_psm_filter$ProteinName))

quant<-rbind(prot_psm,lnc_psm_filter)

# save(quant,file = '/scratch/lb4489/project/hcc/R_code/LBF_quant.RData')

quant

match<-data.frame(run=c(unique(quant$Run)),
                  Fraction=c(NA))

f <- function(x) {
  parts <- unlist(strsplit(x['run'], '_'))  
  tail(parts, 1)  }
match$Fraction<-apply(match,1,f)
write.csv(match,file = '/scratch/lb4489/project/hcc/R_code/anno_LBF.csv')
match<-read.csv('/scratch/lb4489/project/hcc/R_code/anno_LBF.csv',header = F)
match<-match[match(quant$Run,match$V2),]
quant$Fraction<-match$V3

# save(quant,file = '/scratch/lb4489/project/hcc/R_code/LBF_quant.RData')
load('/scratch/lb4489/project/hcc/R_code/LBF_quant.RData')

processedData <- MSstats::dataProcess(quant)

save(processedData,file = '/scratch/lb4489/project/hcc/R_code/LBF_msstats.RData')

protein<-MSstats::quantification(processedData)

save(protein,file = '/scratch/lb4489/project/hcc/R_code/LBF_protein_exp.RData')

################################################################


pep<-data.table::fread('/scratch/lb4489/project/hcc/captac/fragpipe/all_run/tmt-report/abundance_peptide_MD.tsv')

nonadd<-data.frame(Peptide.Sequence=c(pep$Peptide),
                   Protein=c(pep$ProteinID))

nonadd<-nonadd[grep('sPep',nonadd$Protein),]

peplist<-nonadd$Peptide.Sequence
peplist<-unique(peplist)
  
peplist<-data.frame(seq=c(peplist),
                      match=c(NA))
  
pattern2 <- Biostrings::readAAStringSet('/scratch/lb4489/project/hcc/R_code/all_protein.fasta')


for (i in 1:nrow(peplist)){
    
    pattern1 <- Biostrings::AAString(peplist$seq[i])
    
    match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)
    
    peplist$match[i]<-sum(match)
    }
  
list2<-peplist$seq[peplist$match>0]
  
quant<-pep

quant<-quant[!(quant$Peptide %in% list2),]
  
data<-as.data.frame(quant)
data<-data[,c(3,4,10:ncol(data))]
  
data<-tidyr::gather(data, key = "sample",
                      value = "Intensity", -c('ProteinID', 'Peptide'),
                      na.rm = FALSE, convert = FALSE, factor_key = FALSE)
  
data<-as.data.frame(data)
data<-data[!is.na(data$Intensity),]

save(data,file = '/scratch/lb4489/project/hcc/R_code/TMT_quant.RData')
  
annation<-data.table::fread('/scratch/lb4489/project/hcc/captac/fragpipe/all_run/experiment_annotation.tsv',header = T)
annation$run<-rep(paste0('run',annation$plex))
annation$condition[annation$condition=='POOL']<-c('Norm')
annation$replicate<-c(1)
annation$replicate[annation$condition=='POOL']<-c('Norm')
  
input<-data.frame(ProteinName=c(data$ProteinID),PeptideSequence=c(data$Peptide),
                    Charge=c(NA),PSM=c(NA),Mixture=c('Mixture1'),TechRepMixture=c(1),
                    Run=c(NA),Channel=c(NA),BioReplicate=c(NA),
                    Condition=c(NA),Intensity=c(data$Intensity),sample=c(data$sample))
  
match<-annation[match(input$sample,annation$sample),]
  
input$Run<-match$run
input$Channel<-match$channel
input$BioReplicate<-match$replicate
input$Condition<-match$condition
input$Intensity<-2^input$Intensity
  
input$Charge<-c(1)
input$PSM<-input$PeptideSequence

quant.msstats <- MSstatsTMT::proteinSummarization(input,
                                                    method="msstats",
                                                    global_norm=TRUE,
                                                    reference_norm=TRUE,
                                                    remove_norm_channel = TRUE,
                                                    remove_empty_channel = TRUE)
  
protein<-quant.msstats[["ProteinLevelData"]]
  
exparray <- tidyr::spread(protein[,c(5,6,8)], key = Condition, value = Abundance)



