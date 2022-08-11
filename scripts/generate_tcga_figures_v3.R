knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

`%notin%` <- Negate(`%in%`)
sdError <- function(x) sd(x)/sqrt(length(x))

# if (download_data){
#   download_files(pipeline_path)
# }
# 
# raw_tcga_filename=paste(pipeline_path,"inputs","LUAD_exprs.rd",sep="/")
# merged_tcga_filename=paste(pipeline_path,"inputs","merged_data_luad.csv",sep="/")
# 
# load(raw_tcga_filename)
# load(merged_tcga_filename)

setwd("G:/My Drive/Paper website/Scripts")

library(gplots)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(ggplot2)
library(Matrix.utils)
library(RColorBrewer)
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(corrplot)
library(reshape2)
library(pheatmap)
library("ggpubr")
library(reshape)
library(Rmisc)
library(gridExtra)

write_path <- "outputs/"
scaling_sigs <- c("ADH1B_CAF","FAP_CAF","LCAM")

#########################
# loading raw LUAD data and normalizing
load("inputs/LUAD_exprs.rd")
geneData_luad <- expr_mat
rm(expr_mat)

# loading raw LUSC data and normalizing
load("inputs/LUSC_exprs.rd")
geneData_lusc <- expr_mat
rm(expr_mat)

geneData <- cbind(geneData_luad, geneData_lusc)
rm(geneData_luad)
rm(geneData_lusc)

# reading gene lists
t1=read.table("inputs/gene_lists.txt",stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

# genes to be discarded
junk <- c(gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)
junk <- intersect(junk,rownames(geneData))

# dealing with duplicates

# dups <-  names(table(colnames(geneData))[table(colnames(geneData))==2])
geneData <- geneData[,!duplicated(colnames(geneData))]
# for (dup in dups){
  # dup_index <- grep(dup, colnames(geneData))
  # geneData[,dup] <- rowSums(geneData[,dup_index])
# }

# remove excluded genes from rna seq data
geneData <- geneData[!rownames(geneData)%in%junk,]

# normalize to total gene expression
norm_dat <- t(t(geneData)/rowSums(t(geneData)))
# regularize
reg <- 1E-8
# converting to log to normalize
norm_dat <- log2(norm_dat + reg)
length(colnames(norm_dat))
norm_dat <- norm_dat[,names(table(colnames(norm_dat)))]
length(colnames(norm_dat))

# clinical_data <- read.csv("G:/My Drive/Paper website/Scripts/inputs/clinical_luad.csv", row.names = 1)
clinical_data <- read.csv("G:/My Drive/Paper website/Scripts/inputs/clinical_luad_lusc.csv", row.names = 1)

tumor_patients_luad <- rownames(clinical_data)[clinical_data[,"sample_type_id"]%in%c("1")]

tumor_patients_luad <- tumor_patients_luad[tumor_patients_luad%in%colnames(norm_dat)]

norm_dat <- norm_dat[,tumor_patients_luad]

norm_dat_gene_scaled <- t(scale(t(norm_dat)))
rm(geneData)

#########################
# creating signature scores

LCAMhi <- c("LCAMhiPlasma","LCAMhiMoMac","LCAMhiTact")
LCAMlo <- c("LCAMloB","LCAMloAM","LCAMloDC2","LCAMloAZU1","LCAMloDC1")
sigs <- c(LCAMhi,LCAMlo,"ADH1B_CAF","FAP_CAF")

sig_matrix <- matrix(data = NA, nrow=length(sigs), ncol=ncol(norm_dat_gene_scaled))
rownames(sig_matrix) <- sigs
colnames(sig_matrix) <- colnames(norm_dat_gene_scaled)

for (num in 1:length(sigs)){
  sig_name <- sigs[num]
  sig <- gl[[sig_name]]
  sig <- intersect(sig,rownames(norm_dat_gene_scaled))
  if (length(sig) == 1){
    sig_norm_dat_score <- norm_dat_gene_scaled[sig,]
  }else{
    sig_norm_dat_score <- colMeans(norm_dat_gene_scaled[sig,])
  }
  sig_matrix[sig_name,] <- sig_norm_dat_score
}

# adding LCAM score
sig_names <- rownames(sig_matrix)
temp_row <- c(sig_names, "LCAM")
sig_matrix <- rbind(sig_matrix,NA)
rownames(sig_matrix) <- temp_row

LCAMhi <- c("LCAMhiPlasma","LCAMhiMoMac","LCAMhiTact")
LCAMlo <- c("LCAMloB","LCAMloAM","LCAMloDC2","LCAMloAZU1","LCAMloDC1")
sig_matrix["LCAM",] <- colMedians(sig_matrix[LCAMhi,])-colMedians(sig_matrix[LCAMlo,])

#########################
# merging with clincal data
sig_matrix <- t(sig_matrix)

clinical_data <- clinical_data[tumor_patients_luad,]

merged_data <- join.Matrix(sig_matrix,clinical_data, by.x=rownames(sig_matrix),
                           by.y=rownames(clinical_data), all.x = FALSE, all.y = FALSE)

# write.csv(merged_data, "G:/My Drive/Paper website/Scripts/inputs/merged_data_luad_lusc_v2.csv")

#########################
# figure 2G
# graphing patients expression of ADH1B CAF and FAP CAF gene signatures

merged_data <- read.csv("inputs/merged_data_luad_v2.csv", row.names = 1)

caf_genes <- c(gl$ADH1B_CAF, gl$FAP_CAF)
new_dat <- norm_dat_gene_scaled[rownames(norm_dat_gene_scaled)%in%caf_genes,]

sig_matrix <- t(sig_matrix)

caf_sig <- c("ADH1B_CAF","FAP_CAF")

# selecting luad only patients
new_dat <- new_dat[,colnames(new_dat)%in%rownames(merged_data)]

new_dat <- rbind(new_dat, sig_matrix[caf_sig,colnames(new_dat)])
new_dat <- new_dat[c(caf_genes, caf_sig),]

# # below line not needed if using merged_data_luad_v2.csv
# new_dat <- t(scale(t(new_dat)))

breaks = seq(-2.5,2.5,length.out=10)

new_dat <- t(new_dat)
new_dat <- new_dat[rownames(new_dat)%in%rownames(merged_data),]
merged_data <- merged_data[rownames(merged_data)%in%rownames(new_dat),]
new_dat <- cbind(new_dat, merged_data[,c("tumor_nuclei","sample_type_id")])

# selecting only primary tumor

# # below line not needed if using merged_data_luad_v2.csv
# new_dat <- new_dat[new_dat[,"sample_type_id"]%in%c("1"),]

new_dat[,"tumor_nuclei"] <- scale(new_dat[,"tumor_nuclei"])
new_dat <- t(new_dat)

for (sig in caf_sig){
  temp_dat <- new_dat[,order(new_dat[sig,])]
  
  temp_rownames <- rownames(temp_dat)
  temp_dat <- temp_dat[c(temp_rownames[1:(length(temp_rownames)-4)],"tumor_nuclei"),]
  
  pheatmap(temp_dat, cluster_rows=F, cluster_cols=F, scale = "none", show_colnames=F,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
           breaks = breaks)
}

rm(new_dat)
rm(temp_dat)

#########################
# figure S4D
# correlation of gene signatures with each other

all_genes <- c()
for (geneSet in gl){
  all_genes <- c(all_genes, geneSet)
}
all_genes <- intersect(all_genes,rownames(norm_dat_gene_scaled))

merged_data <- read.csv("inputs/merged_data_luad_v2.csv", row.names = 1)

# # # below 2 lines not needed if using merged_data_luad_v2.csv
# tumor_samps <- rownames(merged_data)[merged_data[,"sample_type_id"]%in%c("1")]
# norm_dat_gene_scaled <- norm_dat_gene_scaled[,tumor_samps]

gene_cors <- cor(t(norm_dat_gene_scaled[all_genes,rownames(merged_data)]))
na_genes <- !is.na(gene_cors[1,])
gene_cors <- as.matrix(gene_cors[na_genes,na_genes])
genes_of_int <- c(gl$Fib,gl$ADH1B_CAF,gl$FAP_CAF)
gene_cors <- gene_cors[genes_of_int,genes_of_int]

corrplot(gene_cors, method = 'color', order = "original", type = "upper", number.font = 1,
         tl.col = "black", tl.srt = 45, col=rev(brewer.pal(n=8, name="RdBu")))
rm(gene_cors)

#########################

# FIXXXXXXX #

# figure 2H
# histo pathology

# for (num in c("1","2","3")){

merged_data <- read.csv("inputs/merged_data_luad_lusc_v2.csv", row.names = 1)

score_of_interest <- c("ADH1B_CAF","FAP_CAF")
# selecting unique subtypes with > 10 patients
names_data <- c("Acinar predominant Adc","Papillary predominant Adc",
                "Solid predominant Adc","Micropapillary predom Adc","LUSC")

names_data <- c("Acinar predominant Adc","Papillary predominant Adc",
                "Solid predominant Adc","Micropapillary predom Adc")

merged_data[,"subtype"] <- as.character(merged_data[,"subtype"])
# renaming all unused subtypes as 'Other'
merged_data[merged_data[,"subtype"]%notin%names_data,"subtype"] <- "Other"
# selecting only primary tumor
merged_data <- merged_data[merged_data[,"sample_type_id"]%in%c("1"),]

#########
# merged_data <- merged_data[merged_data[,"tumor_stage_num"]%in%c(num),]
# nrow(merged_data)
#########

merged_data <- merged_data[,c(score_of_interest,"subtype")]

# removing other pathologies
merged_data <- merged_data[merged_data[,"subtype"]%in%names_data,]

# rescaling
# for (score in score_of_interest){
#   merged_data[,score] <- scale(merged_data[,score])
# }

melted_data <- melt(merged_data)
summarized_melted <- summarySE(melted_data, measurevar="value", groupvars=c("subtype","variable"))

cols <- c("khaki3","red")

g <- ggplot(summarized_melted, aes(x=factor(subtype, level = names_data), y=value, fill=variable)) +
  geom_errorbar(data=summarized_melted, mapping=aes(x=factor(subtype, level = names_data), ymin=value-se, ymax=value+se),
                size = 1.5, width = .6, position = "dodge", color = rep(cols,4)) +
  # ggtitle(paste("Stage",num)) +
  geom_point(size = 2, position = position_dodge(width = .6)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.ontop = TRUE
  ) +
  ylim(-1.5,1.5)
print(g)

for (path in names_data){
  n_subtypes <- length(names_data)
  temp_data <- merged_data[merged_data[,"subtype"]==path,]
  temp_p <- t.test(temp_data[,"ADH1B_CAF"],temp_data[,"FAP_CAF"])
  print(nrow(temp_data[temp_data[,"subtype"]==path,]))
  print(temp_p$p.value)
  # print(n_subtypes)
  # print(paste(path, "p value", p.adjust(temp_p$p.value,method="BH",n = n_subtypes)))
}
# }

#########################
# figure 2H
# pathology

# for (num in c("1","2","3")){

merged_data <- read.csv("inputs/merged_data_luad_lusc_v2.csv", row.names = 1)

score_of_interest <- c("ADH1B_CAF","FAP_CAF")

names_data <- c("LUAD","LUSC")

# selecting only primary tumor
merged_data <- merged_data[merged_data[,"sample_type_id"]%in%c("1"),]

# #########
# merged_data <- merged_data[merged_data[,"tumor_stage_num"]%in%c(num),]
# nrow(merged_data)
# #########

merged_data <- merged_data[,c(score_of_interest,"disease")]

# removing other pathologies
merged_data <- merged_data[!is.na(merged_data[,"disease"]),]

# rescaling
# for (score in score_of_interest){
#   merged_data[,score] <- scale(merged_data[,score])
# }

melted_data <- melt(merged_data)
summarized_melted <- summarySE(melted_data, measurevar="value", groupvars=c("disease","variable"))

cols <- c("khaki3","red")

g<-ggplot(summarized_melted, aes(x=factor(disease, level = names_data), y=value, fill=variable)) +
  geom_errorbar(data=summarized_melted, mapping=aes(x=factor(disease, level = names_data), ymin=value-se, ymax=value+se),
                size = 1.5, width = .6, position = "dodge", color = rep(cols,2)) +
  # ggtitle(paste("Stage",num)) +
  geom_point(size = 2, position = position_dodge(width = .6)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.ontop = TRUE
  ) +
  ylim(-0.5,0.5)
print(g)

for (path in names_data){
  n_subtypes <- length(unique(merged_data[,"disease"]))
  temp_data <- merged_data[merged_data[,"disease"]==path,]
  temp_p <- t.test(temp_data[,"ADH1B_CAF"],temp_data[,"FAP_CAF"])
  print(nrow(temp_data[temp_data[,"disease"]==path,]))
  print(temp_p$p.value)
  print(n_subtypes)
  print(paste(path, "p value", p.adjust(temp_p$p.value,method="BH",n = n_subtypes)))
}


temp_data <- merged_data[,]

temp_p <- t.test(temp_data[merged_data[,"disease"]=="LUAD","ADH1B_CAF"],temp_data[merged_data[,"disease"]=="LUSC","ADH1B_CAF"])
temp_p$p.value
temp_p <- t.test(temp_data[merged_data[,"disease"]=="LUAD","FAP_CAF"],temp_data[merged_data[,"disease"]=="LUSC","FAP_CAF"])
temp_p$p.value
# }

#########################
# figure 2H
# stages 2 and above are merged and considered 'high' stage

# troublshooting
path_subset <- "LUSC"

for (path_subset in c("LUAD","LUSC")){
  
merged_data <- read.csv("inputs/merged_data_luad_lusc_v2.csv", row.names = 1)
merged_data[,"tumor_stage_num"] <- as.integer(merged_data[,"tumor_stage_num"])
merged_data <- merged_data[merged_data[,"disease"]==path_subset,]

# selecting only primary tumor
merged_data <- merged_data[merged_data[,"sample_type_id"]%in%c("1"),]

score_of_interest <- c("ADH1B_CAF","FAP_CAF")
# for (score in score_of_interest){
#   merged_data[,score] <- scale(merged_data[,c(score)])
# }
merged_data <- merged_data[,c(score_of_interest,"tumor_stage_num")]
merged_data <- merged_data[!is.na(merged_data[,"tumor_stage_num"]),]

merged_data[merged_data[,"tumor_stage_num"] == 3,"tumor_stage_num"]  <- 2
merged_data[merged_data[,"tumor_stage_num"] == 4,"tumor_stage_num"]  <- 2

stages <- c(1,2)

merged_data[,"tumor_stage_num"] <- as.character(merged_data[,"tumor_stage_num"])

score_of_interest <- c("ADH1B_CAF","FAP_CAF")
# rescaling
# for (score in score_of_interest){
#   merged_data[,score] <- scale(merged_data[,score])
# }

melted_data <- melt(merged_data, id="tumor_stage_num")
summarized_melted <- summarySE(melted_data, measurevar="value", groupvars=c("tumor_stage_num","variable"))

cols <- c("khaki3","red")

# if adjusting number of stanges, "color = rep(cols,2))" must be changed
g<-ggplot(summarized_melted, aes(x=tumor_stage_num, y=value, fill=variable)) +
  geom_errorbar(data=summarized_melted, mapping=aes(x=tumor_stage_num, ymin=value-se, ymax=value+se),
                size = 1.5, width = .6, position = "dodge", color = rep(cols,2)) +
  geom_point(size = 2, position = position_dodge(width = .6)) +
  ylim(-0.5,0.5) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.ontop = TRUE
  )
print(g)

for (stage in c(1,2)){
  temp_p <- t.test(merged_data[merged_data[,"tumor_stage_num"] == stage,"ADH1B_CAF"],merged_data[merged_data[,"tumor_stage_num"] == stage,"FAP_CAF"])
  print(paste(path_subset, "p value", "stage", stage))
  print(temp_p$p.value)

  # print(paste(path_subset, "p value", "stage", stage, p.adjust(temp_p$p.value,method="BH",n = 2)))
}

for (caf in c("ADH1B_CAF","FAP_CAF")){
  temp_p <- t.test(merged_data[merged_data[,"tumor_stage_num"] == 2,caf],merged_data[merged_data[,"tumor_stage_num"] == 1,caf])
  print(paste(path_subset, "p value", "stage", caf))
  print(temp_p$p.value)
}
}

#########################
# figure 2I and 3G

merged_data <- read.csv("inputs/merged_data_luad_v2.csv", row.names = 1, stringsAsFactors = F)
sigs <- c("ADH1B_CAF","FAP_CAF","LCAM","delta")

# selecting only primary tumor

# # # below line not needed if using merged_data_luad_v2.csv
# merged_data <- merged_data[merged_data[,"sample_type_id"]%in%c("1"),]

merged_data <- merged_data[,c("ADH1B_CAF","FAP_CAF","LCAM","tumor_stage_num","TMB")]

# various residuals, etc
delta_caf <- merged_data[,"FAP_CAF"] - merged_data[,"ADH1B_CAF"]

delta_to_lcam_residuals <- lm(merged_data[,"LCAM"] ~ delta_caf)

merged_data[,"delta"] <- delta_caf
merged_data[,"delta"] <- scale(merged_data[,"delta"])

merged_data[,"TMB"] <- log2(as.numeric(merged_data[,"TMB"]+1))
delta_predict <- lm(merged_data[,"delta"] ~ merged_data[,"TMB"] + merged_data[,"tumor_stage_num"])

stats <- summary(lm(merged_data[,"delta"] ~ merged_data[,"TMB"]))

ggplot(merged_data, aes(x=delta, y=TMB)) +
  geom_point(size=3) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ggtitle(paste("R = ", round(sqrt(stats$r.squared),2)))

stats <- summary(delta_to_lcam_residuals)

ggplot(merged_data, aes(x=delta, y=LCAM)) +
  geom_point(size=3) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylim(-2,2) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3),limits = c(-3,3)) +
  ggtitle(paste("R = ", round(sqrt(stats$r.squared),2)))

ggplot(merged_data, aes(x=ADH1B_CAF, y=FAP_CAF)) +
  geom_point(size=3) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ggtitle(paste("R = ", round(sqrt(stats$r.squared),2)))

