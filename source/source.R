plot_fdr <- function(target,decoy,label) {
  
  target<-read.table(target,header = T)
  decoy<-read.table(decoy,header = T)
  
  plot<-data.frame(
    score=c(target$score,decoy$score),
    group=c(rep('target',nrow(target)),rep('decoy',nrow(decoy)))
  )
  plot$group<-factor(plot$group,levels=c('target','decoy'))
  
  shres<-get_threshold(plot,0.01)
  
  p<-list()
  
  p[[1]]<-
    ggplot2::ggplot(plot, aes(score, fill = group, col = I("black"))) +
    geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
    labs(x = 'score', y = "",
         title = "Global") +
    scale_fill_manual(values = c("decoy" = "#FF9900","target" = "#009900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)+
    geom_vline(xintercept=c(max(shres)),lty=3,col="black",lwd=0.5)
  
  plot$pep<-runif(nrow(plot))
  p[[4]]<-
  ggplot(data = plot,aes(x=score,y=pep, colour=group))+
    geom_point(shape=21,size=4)+
    scale_colour_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)+
    labs(x = 'score', y = "")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    geom_vline(xintercept=c(max(shres)),lty=3,col="black",lwd=0.5)
  
  target_clas<-target[grep(label,target$proteinIds),]
  decoy_clas<-decoy[grep(label,decoy$proteinIds),]
  
  plot<-data.frame(
    score=c(target_clas$score,decoy_clas$score),
    group=c(rep('target',nrow(target_clas)),rep('decoy',nrow(decoy_clas)))
  )
  plot$group<-factor(plot$group,levels=c('target','decoy'))
  
  shres<-get_threshold(plot,0.01)
  
  p[[2]]<-
    ggplot2::ggplot(plot, aes(score, fill = group, col = I("black"))) +
    geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
    labs(x = 'score', y = "",
         title = "Canonical protein") +
    scale_fill_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)+
    geom_vline(xintercept=c(max(shres)),lty=3,col="black",lwd=0.5)
  
  plot$pep<-runif(nrow(plot))
  p[[5]]<-
    ggplot(data = plot,aes(x=score,y=pep, colour=group))+
    geom_point(shape=21,size=4)+
    scale_colour_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)+
    labs(x = 'score', y = "")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    geom_vline(xintercept=c(max(shres)),lty=3,col="black",lwd=0.5)
  
  target_clas<-target[-grep(label,target$proteinIds),]
  decoy_clas<-decoy[-grep(label,decoy$proteinIds),]
  
  plot<-data.frame(
    score=c(target_clas$score,decoy_clas$score),
    group=c(rep('target',nrow(target_clas)),rep('decoy',nrow(decoy_clas)))
  )
  plot$group<-factor(plot$group,levels=c('target','decoy'))
  
  shres<-get_threshold(plot,0.01)
  
  p[[3]]<-
    ggplot2::ggplot(plot, aes(score, fill = group, col = I("black"))) +
    geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
    labs(x = 'score', y = "",
         title = "sORF encoded peptide") +
    scale_fill_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)+
    geom_vline(xintercept=c(max(shres)),lty=3,col="black",lwd=0.5)
  
  plot$pep<-runif(nrow(plot))
  p[[6]]<-
    ggplot(data = plot,aes(x=score,y=pep, colour=group))+
    geom_point(shape=21,size=4)+
    scale_colour_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)+
    labs(x = 'score', y = "")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    geom_vline(xintercept=c(max(shres)),lty=3,col="black",lwd=0.5)
  
  return(p)
}

get_threshold <- function(data, threshold) {
  
  data<-data[order(data$score,decreasing = T),]
  
  filter<-data.frame(hold=quantile(data$score, probs = c(0,0.1,0.2,0.3,0.4,0.5,
                                                         0.6,0.7,0.8,0.9,1)),
                     fdr=c(NA))
  
  for (i in 1:nrow(filter)){
    more<-data[data$score>filter$hold[i],]
    a <- nrow(more) #[less$group=='target',]
    b <- nrow(more[more$group!='target',])
    filter$fdr[i] <- b/a
  }
  filter$fdr[filter$fdr=='NaN']<-0
  filter<-filter[order(filter$fdr),]
  
  insert_index <- findInterval(0.01, filter$fdr)
  
  data$fdr<-NA
  
  for(i in min(which(data$score<filter$hold[insert_index])):max(which(data$score>filter$hold[insert_index+1]))){
    
    less <- data[1:i,]
    a <- nrow(less) #[less$group=='target',]
    b <- nrow(less[less$group!='target',])
    data$fdr[i] <- b/a
    
    if (b/a > 0.01) {
      stop_index <- i
      break
    }
  }
  
  return(data$score[stop_index])
}
