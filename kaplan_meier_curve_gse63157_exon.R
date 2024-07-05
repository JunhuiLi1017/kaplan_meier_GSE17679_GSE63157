setwd("~/dropbox/Project/Junhui_Li_UMMS/Michael_Green/Shahid_Banday/20230309_kaplan_meier/result/gse63157_exon//")
library(survminer)
library(survival)
##------------------
## prepare the dataset
##------------------
method <- "var"
#method <- "iqr"
gse63157_exon_raw <- read.table(paste0("gse63157_exon_",method,".txt"),comment.char = "",header=T)
dat <- gse63157_exon_raw
# convert variable type
dat$efs_time <- as.numeric(dat$efs_time)
dat$os_time <- as.numeric(dat$os_time)
dat$C1GALT1_value <- as.numeric(dat$C1GALT1)
str(dat)

# sort dataframe based on C1GALT1_value columns
dat_sort <- dat[order(dat$C1GALT1_value, decreasing = FALSE),]
dat_sort$id <- 1:nrow(dat_sort)
head(dat_sort)
str(dat_sort)
cbind(dat_sort[,c(1,2,8)])
dat_sort$os_time <- dat_sort$os_time/30
dat_sort$efs_time <- dat_sort$efs_time/30

##-----------------------
# correlation plot
##-----------------------

dat_sort[1:5,]
dat_corr <- dat_sort
dat_corr$C1GALT1 <- as.numeric(dat_corr$C1GALT1)
dat_corr$C1GALT1_log2 <- log2(dat_corr$C1GALT1)
summary(dat_corr$C1GALT1)
dat_corr$GLI1 <- as.numeric(dat_corr$GLI1)
dat_corr$GLI1_log2 <- log2(dat_corr$GLI1)
summary(dat_corr$GLI1)
dat_corr$SMO <- as.numeric(dat_corr$SMO)
dat_corr$SMO_log2 <- log2(dat_corr$SMO)
summary(dat_corr$SMO)
str(dat_corr)


pc1 <- ggplot(data=dat_corr,aes(x=C1GALT1_log2,y=GLI1_log2))+
  geom_point()+
  theme_classic2()+
  stat_cor(method = "pearson", label.y = max(dat_corr$GLI1_log2)+1)

pc2 <- ggplot(data=dat_corr,aes(x=C1GALT1_log2,y=SMO_log2))+
  geom_point()+
  theme_classic2()+
  stat_cor(method = "pearson", label.y = max(dat_corr$SMO_log2)+1)

pdf(paste0("correlation_",method,"/gse63157_exon_C1GALT1_GLI1_correlation.pdf"))
print(pc1)
dev.off()
pdf(paste0("correlation_",method,"/gse63157_exon_C1GALT1_SMO_correlation.pdf"))
print(pc2)
dev.off()
# install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
dat_corr_gene <- dat_corr[,c("C1GALT1_log2","GLI1_log2","SMO_log2")]
pdf(paste0("correlation_",method,"/gse63157_exon_C1GALT1_GLI1_SMO_correlation.pdf"))
chart.Correlation(dat_corr_gene, histogram = TRUE, method = "pearson")
dev.off()

##-----------------------
# go through all cutoff methods
##-----------------------
## threshold for iqr method
# cutoff <- list()
# cutoff["scan"]  <- list(55)
# cutoff["median"]  <- list(42)
# cutoff["average"]  <- list(50)
# cutoff["first_quartile"]  <- list(21)
# cutoff["last_quartile"]  <- list(63)
# 
# cutoff_ef <- list()
# cutoff_ef["scan"]  <- list(49)
# cutoff_ef["median"]  <- list(42)
# cutoff_ef["average"]  <- list(50)
# cutoff_ef["first_quartile"]  <- list(21)
# cutoff_ef["last_quartile"]  <- list(63)

## threshold for var method
cutoff <- list()
cutoff["scan"]  <- list(48)
cutoff["median"]  <- list(42)
cutoff["average"]  <- list(51)
cutoff["first_quartile"]  <- list(21)
cutoff["last_quartile"]  <- list(63)

cutoff_ef <- list()
cutoff_ef["scan"]  <- list(47)
cutoff_ef["median"]  <- list(42)
cutoff_ef["average"]  <- list(51)
cutoff_ef["first_quartile"]  <- list(21)
cutoff_ef["last_quartile"]  <- list(63)

## overall
#i=names(cutoff)[5]
for(i in names(cutoff)){
  cut_off <- cutoff[[i]]
  dat_sort1 <- dat_sort
  if(i=="first_vs_last_quartile" | i=="curtain"){
    dat_sort1[(nrow(dat_sort1)-cut_off+1):nrow(dat_sort1),"C1GALT1"] <- "high"
    dat_sort1[1:cut_off,"C1GALT1"] <- "low"
    dat_sort1 <- dat_sort1[c(1:cut_off,(nrow(dat_sort1)-cut_off+1):nrow(dat_sort1)),]
  }else{
    dat_sort1[1:cut_off,"C1GALT1"] <- "low"
    dat_sort1[(cut_off+1):nrow(dat_sort1),"C1GALT1"] <- "high"
  }
  cut_off <- dat_sort1[cut_off,"C1GALT1_value"]
  dat_sort1$C1GALT1 <- as.factor(dat_sort1$C1GALT1)
  head(dat_sort1)
  fit_overall <- survfit(Surv(os_time, os_status== 1) ~ C1GALT1,
                         data = dat_sort1)
  
  p1_overall <- ggsurvplot(fit_overall,
                           legend.title = "",
                           palette = c("darkred", "darkblue"),
                           data = dat_sort1,
                           risk.table = TRUE,
                           pval = TRUE,
                           surv.plot.height=0.85)+
    xlab("Follow up in months")+
    ylab("Overall survival probability")
  dat_sort1$os_status <- as.factor(dat_sort1$os_status)
  p2 <- ggplot(data=dat_sort1,aes(x=id,y=C1GALT1_value,color=os_status))+
    geom_point()+
    theme_classic()+
    geom_vline(xintercept = sum(dat_sort1$C1GALT1=="low"), color="gray")+
    geom_hline(yintercept = cut_off, color="gray")+
    ylab("Expression")+
    xlab("sorted by expression")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "top")
  
  
  #(p1$plot/p1$table)|p2
  pdf(paste0("kaplan_meier_curve_",method,"/gse63157_trans_2plots_cutoff_method_",i,".overall.pdf"))
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
  pdf(paste0("kaplan_meier_curve_",method,"/gse63157_trans_3plots_cutoff_method_",i,".overall.pdf"))
  print(p1_overall$plot + p1_overall$table + p2 + 
          plot_layout(design = layout))
  dev.off()
}


## eventfree
#i=names(cutoff)[6]
for(i in names(cutoff_ef)){
  cut_off <- cutoff_ef[[i]]
  dat_sort1 <- dat_sort
  if(i=="first_vs_last_quartile" | i=="curtain"){
    dat_sort1[(nrow(dat_sort1)-cut_off+1):nrow(dat_sort1),"C1GALT1"] <- "high"
    dat_sort1[1:cut_off,"C1GALT1"] <- "low"
    dat_sort1 <- dat_sort1[c(1:cut_off,(nrow(dat_sort1)-cut_off+1):nrow(dat_sort1)),]
  }else{
    dat_sort1[1:cut_off,"C1GALT1"] <- "low"
    dat_sort1[(cut_off+1):nrow(dat_sort1),"C1GALT1"] <- "high"
  }
  cut_off <- dat_sort1[cut_off,"C1GALT1_value"]
  dat_sort1$C1GALT1 <- as.factor(dat_sort1$C1GALT1)
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
  fit_efs <- survfit(Surv(efs_time, efs_status== 1) ~ C1GALT1,
                     data = dat_sort1)
  
  p1_efs <- ggsurvplot(fit_efs, 
                       legend.title = "",
                       palette = c("darkred", "darkblue"),
                       data = dat_sort1, 
                       risk.table = TRUE,
                       pval = TRUE,
                       surv.plot.height=0.85)+
    xlab("Follow up in months")+
    ylab("eventfree survival probability")
  
  dat_sort1$os_status <- as.factor(dat_sort1$os_status)
  p2 <- ggplot(data=dat_sort1,aes(x=id,y=C1GALT1_value,color=os_status))+
    geom_point()+
    theme_classic()+
    geom_vline(xintercept = sum(dat_sort1$C1GALT1=="low"), color="gray")+
    geom_hline(yintercept = cut_off, color="gray")+
    ylab("Expression")+
    xlab("sorted by expression")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "top")
  
  #(p1$plot/p1$table)|p2
  pdf(paste0("kaplan_meier_curve_",method,"/gse63157_trans_2plots_cutoff_method_",i,".eventfree.pdf"))
  plot_layout(widths = c(2, 1), heights = c(1, 1))
  print(p1_efs$plot + p2 +  plot_layout(widths = c(3, 1)))
  dev.off()
  pdf(paste0("kaplan_meier_curve_",method,"/gse63157_trans_3plots_cutoff_method_",i,".eventfree.pdf"))
  print(p1_efs$plot + p1_efs$table + p2 + 
          plot_layout(design = layout))
  dev.off()
}


