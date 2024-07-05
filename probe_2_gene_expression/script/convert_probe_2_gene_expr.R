# library(devtools)
# install_github("emreg00/pepper", dependencies = T)

## convert probe to gene expression
# ?PEPPER::convert.probe.to.gene.expression

library(PEPPER)

setwd("~/Dropbox (UMass Medical School)/Project/Junhui_Li_UMMS/Michael_Green/Shahid_Banday/20230309_kaplan_meier/probe_2_gene_expression/")

## 1. read probe gene file(some probe id downloaded from R2 is not found in expression dataset)
## so we drop the probe id and gene id we prepared from R2, use probe id from expression dataset instead.
#probeid2geneid <- read.table("probeid_geneid.txt",sep="\t",header=T)

## 2. read probe expression file
## ---------exon---------------
gse63157_exon <- read.table("gse63157_exon.txt",header=F,comment.char = "")
dim(gse63157_exon)
gse63157_exon[1:25,1:5]

probeid2geneid <- gse63157_exon[gse63157_exon[,1] %in% c("GLI1","SMO","C1GALT1"),c(2,1)]
colnames(probeid2geneid) <- c("Probe","Gene")

gse63157_exon_3gene <- gse63157_exon[which(gse63157_exon[,1] %in% c("GLI1","SMO","C1GALT1")),]
colnames(gse63157_exon_3gene) <- gse63157_exon[1,]
rownames(gse63157_exon_3gene) <- gse63157_exon_3gene[,2]
gse63157_exon_3gene <- gse63157_exon_3gene[,-c(1:2)]
gse63157_exon_3gene[1:5,1:5]

## ---------transcript---------------
## no need to transform, since there is only on probe and gene.

## 3. convert probe expression to gene expression with all method

gse63157_exon_geneexpr_iqr <- PEPPER::convert.probe.to.gene.expression(gse63157_exon_3gene,probeid2geneid,selection.method = "iqr")
gse63157_exon_geneexpr_var <- PEPPER::convert.probe.to.gene.expression(gse63157_exon_3gene,probeid2geneid,selection.method = "var")
#gse63157_exon_geneexpr_med <- PEPPER::convert.probe.to.gene.expression(gse63157_exon_3gene,probeid2geneid,selection.method = "med")

genename <- rownames(gse63157_exon_geneexpr_iqr)
gse63157_exon_geneexpr_iqr <- cbind(matrix(rep(genename,2),3,2),gse63157_exon_geneexpr_iqr)
colnames(gse63157_exon_geneexpr_iqr)[1:2] <- rep("genename",2)
gse63157_exon_geneexpr_var <- cbind(matrix(rep(genename,2),3,2),gse63157_exon_geneexpr_var)
colnames(gse63157_exon_geneexpr_var)[1:2] <- rep("genename",2)



meta_data <- gse63157_exon[1:17,]
colnames(meta_data) <- colnames(gse63157_exon_geneexpr_iqr)

gse63157_exon_iqr <- t(rbind(meta_data,gse63157_exon_geneexpr_iqr))
gse63157_exon_var <- t(rbind(meta_data,gse63157_exon_geneexpr_var))
#gse63157_exon_med <- rbind(meta_data,gse63157_exon_geneexpr_med)

gse63157_exon_iqr1 <- gse63157_exon_iqr[,c(5:7,10:11,15:16,18:20)]
head(gse63157_exon_iqr1)
colnames(gse63157_exon_iqr1) <- c("efs_time","efs_status","efs_event",
                                  "os_time","os_status","survival_status",
                                  "survival_time","C1GALT1","GLI1","SMO")
gse63157_exon_iqr2 <- gse63157_exon_iqr1[-c(1:2),]

gse63157_exon_var1 <- gse63157_exon_var[,c(5:7,10:11,15:16,18:20)]
head(gse63157_exon_var1)
colnames(gse63157_exon_var1) <- colnames(gse63157_exon_iqr1)
gse63157_exon_var2 <- gse63157_exon_var1[-c(1:2),]

write.table(gse63157_exon_iqr2,"gse63157_exon_iqr.txt",sep="\t",
            quote=F)
write.table(gse63157_exon_var2,"gse63157_exon_var.txt",sep="\t",
            quote=F)



