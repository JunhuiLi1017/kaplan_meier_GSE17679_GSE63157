library(survminer)
library(survival)
library(ggplot2)
library(patchwork)

setwd("~/Dropbox (UMass Medical School)/Project/Junhui_Li_UMMS/Michael_Green/Shahid_Banday/20230309_kaplan_meier/")

##
##correlation with dataset with log2 transform
##see slide for more details
##-------------------------
##gse63157
##-------------------------
gse63157_raw <- read.table("ps_avgpres_gse63157geop85_huex10p_box1680745400-datagrabber-.txt",comment.char = "")
gse63157_raw[1:20,1:5]
gse63157_C1GALT1 <- gse63157_raw[c(1:17,which(gse63157_raw[,1] %in% "C1GALT1")),]
gse63157 <- as.data.frame(t(gse63157_C1GALT1))
dim(gse63157)

colnames(gse63157) <- c("probeset","age","age_at_enrollment","diagnosis","efs_time","efs_status","efs_event","experiment_date","id","os_time","os_status","primary_tumor_site","sex","study","survival_status","survival_time","tumor_content","C1GALT1_2989441","C1GALT1_2989442","C1GALT1_2989443","C1GALT1_2989444","C1GALT1_2989446","C1GALT1_2989447","C1GALT1_2989448")
dat <- gse63157[-c(1:2),]
dat_remove_nd <- dat

dat_remove_nd$age <- as.numeric(dat_remove_nd$age)
dat_remove_nd$efs_time <- as.numeric(dat_remove_nd$efs_time)
dat_remove_nd$os_time <- as.numeric(dat_remove_nd$os_time)
dat_remove_nd$C1GALT1_2989441_value <- as.numeric(dat_remove_nd$C1GALT1_2989441)
dat_remove_nd$C1GALT1_2989442_value <- as.numeric(dat_remove_nd$C1GALT1_2989442)
dat_remove_nd$C1GALT1_2989443_value <- as.numeric(dat_remove_nd$C1GALT1_2989443)
dat_remove_nd$C1GALT1_2989444_value <- as.numeric(dat_remove_nd$C1GALT1_2989444)
dat_remove_nd$C1GALT1_2989446_value <- as.numeric(dat_remove_nd$C1GALT1_2989446)
dat_remove_nd$C1GALT1_2989447_value <- as.numeric(dat_remove_nd$C1GALT1_2989447)
dat_remove_nd$C1GALT1_2989448_value <- as.numeric(dat_remove_nd$C1GALT1_2989448)


dat_remove_nd_sort <- dat_remove_nd[order(dat_remove_nd$C1GALT1_2989442_value, decreasing = FALSE),]
dat_remove_nd_sort$id <- 1:nrow(dat_remove_nd_sort)

dat_remove_nd_sort$os_status <- as.factor(dat_remove_nd_sort$os_status)
str(dat_remove_nd_sort)

## scan threshold
cutoff <- list()
cutoff["scan"]  <- list(96.9)
cutoff["median"]  <- list(81.8)
cutoff["average"]  <- list(87.1)
cutoff["first_quartile"]  <- list(54.4)
cutoff["last_quartile"]  <- list(109.9)
cutoff["first_vs_last_quartile"]  <- list(23)
cutoff["curtain"]  <- list(14)

#i=names(cutoff)[1]
for(i in names(cutoff)){
  cut_off <- cutoff[[i]]
  if(i=="first_vs_last_quartile" | i=="curtain"){
    dat_remove_nd_sort[(cut_off+1):nrow(dat_remove_nd_sort),"C1GALT1"] <- "high"
    dat_remove_nd_sort[1:cut_off,"C1GALT1"] <- "low"
  }else{
    dat_remove_nd_sort[dat_remove_nd_sort$C1GALT1_value > cut_off,"C1GALT1"] <- "high"
    dat_remove_nd_sort[dat_remove_nd_sort$C1GALT1_value <= cut_off,"C1GALT1"] <- "low"
    dat_remove_nd_sort$C1GALT1 <- as.factor(dat_remove_nd_sort$C1GALT1)
  }
  
  #Y = Surv(dat_remove_nd_sort$overall, dat_remove_nd_sort$status == "dead")
  #kmfit = survfit(Y ~ dat_remove_nd_sort$C1GALT1)
  #summary(kmfit, times = c(seq(0, 1000, by = 100)))
  
  fit <- survfit(Surv(overall, status== "dead") ~ C1GALT1,
                 data = dat_remove_nd_sort)
  
  p1 <- ggsurvplot(fit, data = dat_remove_nd_sort, risk.table = TRUE,pval = TRUE,surv.plot.height=0.85)+
    xlab("Follow up in months")+
    ylab("Overall survival probability")
  
  p2 <- ggplot(data=dat_remove_nd_sort,aes(x=id,y=C1GALT1_value,color=Event))+
    geom_point()+
    theme_classic()+
    geom_vline(xintercept = sum(dat_remove_nd_sort$C1GALT1=="low"), color="gray")+
    geom_hline(yintercept = cut_off, color="gray")+
    ylab("Expression")+
    xlab("sorted by expression")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "top")
  
  
  #(p1$plot/p1$table)|p2
  pdf(paste0("gse17679_2plots_cutoff_method_",i,".pdf"))
  plot_layout(widths = c(2, 1), heights = c(1, 1))
  print(p1$plot + p2 +  plot_layout(widths = c(3, 1)))
  dev.off()
  
  layout <- "
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  BBBBBBCC
  "
  
  pdf(paste0("gse17679_3plots_cutoff_method_",i,".pdf"))
  print(p1$plot + p1$table + p2 + 
          plot_layout(design = layout))
  dev.off()
}



















##-------------------------
##gse17679
##-------------------------
gse17679_raw <- read.table("ps_avgpres_gse17679geo117_u133p2_box1680745218-datagrabber-.txt",comment.char = "")
gse17679_raw[1:15,1:5]
gse17679_C1GALT1 <- gse17679_raw[c(1:10,which(gse17679_raw[,1] %in% "C1GALT1")),]
gse17679 <- as.data.frame(t(gse17679_C1GALT1))
gse17679[1:5,]
colnames(gse17679) <- c("hugo","age","efs","id","overall","sample_series_id","set","sex","state","status","C1GALT1")
dat <- gse17679[-c(1:2),]
dat[1:5,]
dat_remove_nd <- dat[!dat$efs %in% "nd",]

dat_remove_nd$age <- as.numeric(dat_remove_nd$age)
dat_remove_nd$efs <- as.numeric(dat_remove_nd$efs)
dat_remove_nd$overall <- as.numeric(dat_remove_nd$overall)
dat_remove_nd$C1GALT1_value <- as.numeric(dat_remove_nd$C1GALT1)
str(dat_remove_nd)
dat_remove_nd_sort <- dat_remove_nd[order(dat_remove_nd$C1GALT1_value, decreasing = FALSE),]
dat_remove_nd_sort$id <- 1:nrow(dat_remove_nd_sort)
dat_remove_nd_sort$Event[dat_remove_nd_sort$status=="dead"] <- "YES"
dat_remove_nd_sort$Event[dat_remove_nd_sort$status!="dead"] <- "NO"
dat_remove_nd_sort$Event <- as.factor(dat_remove_nd_sort$Event)
dat_remove_nd_sort$status
head(dat_remove_nd_sort)
str(dat_remove_nd_sort)
## scan threshold
cutoff <- list()
cutoff["scan"]  <- list(63.1)
cutoff["median"]  <- list(81.8)
cutoff["average"]  <- list(87.1)
cutoff["first_quartile"]  <- list(54.4)
cutoff["last_quartile"]  <- list(109.9)
cutoff["first_vs_last_quartile"]  <- list(23)
cutoff["curtain"]  <- list(14)

#i=names(cutoff)[1]
for(i in names(cutoff)){
  cut_off <- cutoff[[i]]
  if(i=="first_vs_last_quartile" | i=="curtain"){
    dat_remove_nd_sort[(cut_off+1):nrow(dat_remove_nd_sort),"C1GALT1"] <- "high"
    dat_remove_nd_sort[1:cut_off,"C1GALT1"] <- "low"
  }else{
    dat_remove_nd_sort[dat_remove_nd_sort$C1GALT1_value > cut_off,"C1GALT1"] <- "high"
    dat_remove_nd_sort[dat_remove_nd_sort$C1GALT1_value <= cut_off,"C1GALT1"] <- "low"
    dat_remove_nd_sort$C1GALT1 <- as.factor(dat_remove_nd_sort$C1GALT1)
  }
    
  #Y = Surv(dat_remove_nd_sort$overall, dat_remove_nd_sort$status == "dead")
  #kmfit = survfit(Y ~ dat_remove_nd_sort$C1GALT1)
  #summary(kmfit, times = c(seq(0, 1000, by = 100)))
  
  fit <- survfit(Surv(overall, status== "dead") ~ C1GALT1,
                 data = dat_remove_nd_sort)

  p1_overall <- ggsurvplot(fit, data = dat_remove_nd_sort, risk.table = TRUE,pval = TRUE,surv.plot.height=0.85)+
    xlab("Follow up in months")+
    ylab("Overall survival probability")
  
  p2 <- ggplot(data=dat_remove_nd_sort,aes(x=id,y=C1GALT1_value,color=Event))+
    geom_point()+
    theme_classic()+
    geom_vline(xintercept = sum(dat_remove_nd_sort$C1GALT1=="low"), color="gray")+
    geom_hline(yintercept = cut_off, color="gray")+
    ylab("Expression")+
    xlab("sorted by expression")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "top")
  
  
  #(p1$plot/p1$table)|p2
  pdf(paste0("gse17679_2plots_cutoff_method_",i,".overall.pdf"))
    plot_layout(widths = c(2, 1), heights = c(1, 1))
    print(p1_overall$plot + p2 +  plot_layout(widths = c(3, 1)))
  dev.off()
  
  layout <- "
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  AAAAAACC
  BBBBBBCC
  "
  
  pdf(paste0("gse17679_3plots_cutoff_method_",i,".overall.pdf"))
    print(p1_overall$plot + p1_overall$table + p2 + 
    plot_layout(design = layout))
  dev.off()
  
  
  fit_efs <- survfit(Surv(efs, status== "dead") ~ C1GALT1,
                 data = dat_remove_nd_sort)
  
  p1_efs <- ggsurvplot(fit.event.free, data = dat_remove_nd_sort, risk.table = TRUE,pval = TRUE,surv.plot.height=0.85)+
    xlab("Follow up in months")+
    ylab("Overall survival probability")
  
  
}

