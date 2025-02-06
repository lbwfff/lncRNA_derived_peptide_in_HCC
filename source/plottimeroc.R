
plottimeroc<-function(net) {

library(dplyr)
library(survivalROC)

roc_1yr <- survivalROC(
  Stime       = net$time, 
  status      = net$status,
  marker      = net$pred_class,
  span        = nrow(net)^-0.5,
  predict.time = 1*365
)

df_1yr <- data.frame(
  FP   = roc_1yr$FP,
  TP   = roc_1yr$TP,
  year = "1 year",                             
  AUC  = sprintf("%.3f", roc_1yr$AUC)           
)

roc_years <- c(2,3)

df_list <- list(df_1yr) 

for (i in roc_years) {
  roc_i <- survivalROC(
    Stime       = net$time, 
    status      = net$status,
    marker      = net$pred_class,
    span        = nrow(net)^-0.5,
    predict.time = i*365
  )

  df_i <- data.frame(
    FP   = roc_i$FP,
    TP   = roc_i$TP,
    year = paste0(i, " years"),
    AUC  = sprintf("%.3f", roc_i$AUC)
  )
  
  df_list[[length(df_list)+1]] <- df_i
}

df_all <- bind_rows(df_list)

df_all <- df_all %>%
  group_by(year) %>%
  mutate(
    label = paste0(
      year, 
      " (AUC=", unique(AUC), ")"
    )
  ) %>%
  ungroup()


myColors <- c("1 year (AUC=0.XX)"="#FFC1C1", 
              "3 years (AUC=0.XX)"="#63B8FF", 
              "5 years (AUC=0.XX)"="#FA8072")

p<-list()

p[[1]]<-
ggplot(df_all, aes(x=FP, y=TP, color=label)) +
  geom_line(size=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  coord_equal() +  
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title    = "ROC curve",
    x        = "False Positive Rate",
    y        = "True Positive Rate",
    color    = "Time & AUC"  
  ) +
  theme_bw(base_size=14) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))



roc_1yr <- survivalROC(
  Stime       = net$time, 
  status      = net$status,
  marker      = net$pred_class2,
  span        = nrow(net)^-0.5,
  predict.time = 1*365
)

df_1yr <- data.frame(
  FP   = roc_1yr$FP,
  TP   = roc_1yr$TP,
  year = "1 year",                             
  AUC  = sprintf("%.3f", roc_1yr$AUC)           
)

roc_years <- c(2,3)

df_list <- list(df_1yr) 

for (i in roc_years) {
  roc_i <- survivalROC(
    Stime       = net$time, 
    status      = net$status,
    marker      = net$pred_class2,
    span        = nrow(net)^-0.5,
    predict.time = i*365
  )
  
  df_i <- data.frame(
    FP   = roc_i$FP,
    TP   = roc_i$TP,
    year = paste0(i, " years"),
    AUC  = sprintf("%.3f", roc_i$AUC)
  )
  
  df_list[[length(df_list)+1]] <- df_i
}

df_all <- bind_rows(df_list)

df_all <- df_all %>%
  group_by(year) %>%
  mutate(
    label = paste0(
      year, 
      " (AUC=", unique(AUC), ")"
    )
  ) %>%
  ungroup()


myColors <- c("1 year (AUC=0.XX)"="#FFC1C1", 
              "3 years (AUC=0.XX)"="#63B8FF", 
              "5 years (AUC=0.XX)"="#FA8072")

p[[2]]<-
  ggplot(df_all, aes(x=FP, y=TP, color=label)) +
  geom_line(size=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  coord_equal() +  
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title    = "ROC curve",
    x        = "False Positive Rate",
    y        = "True Positive Rate",
    color    = "Time & AUC"  
  ) +
  theme_bw(base_size=14) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))




roc_1yr <- survivalROC(
  Stime       = net$time, 
  status      = net$status,
  marker      = net$pred_class3,
  span        = nrow(net)^-0.5,
  predict.time = 1*365
)

df_1yr <- data.frame(
  FP   = roc_1yr$FP,
  TP   = roc_1yr$TP,
  year = "1 year",                             
  AUC  = sprintf("%.3f", roc_1yr$AUC)           
)

roc_years <- c(2,3)

df_list <- list(df_1yr) 

for (i in roc_years) {
  roc_i <- survivalROC(
    Stime       = net$time, 
    status      = net$status,
    marker      = net$pred_class3,
    span        = nrow(net)^-0.5,
    predict.time = i*365
  )
  
  df_i <- data.frame(
    FP   = roc_i$FP,
    TP   = roc_i$TP,
    year = paste0(i, " years"),
    AUC  = sprintf("%.3f", roc_i$AUC)
  )
  
  df_list[[length(df_list)+1]] <- df_i
}

df_all <- bind_rows(df_list)

df_all <- df_all %>%
  group_by(year) %>%
  mutate(
    label = paste0(
      year, 
      " (AUC=", unique(AUC), ")"
    )
  ) %>%
  ungroup()


myColors <- c("1 year (AUC=0.XX)"="#FFC1C1", 
              "3 years (AUC=0.XX)"="#63B8FF", 
              "5 years (AUC=0.XX)"="#FA8072")

p[[3]]<-
  ggplot(df_all, aes(x=FP, y=TP, color=label)) +
  geom_line(size=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  coord_equal() + 
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title    = "ROC curve",
    x        = "False Positive Rate",
    y        = "True Positive Rate",
    color    = "Time & AUC"   
  ) +
  theme_bw(base_size=14) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))

return(p)

}



