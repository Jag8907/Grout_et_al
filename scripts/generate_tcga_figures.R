knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

`%notin%` <- Negate(`%in%`)
sdError <- function(x) sd(x)/sqrt(length(x))

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

#########################
#########################
#########################
#########################
# create download location for this LUAD
load("/inputs/LUAD_exprs.rd")
# add tcga description when downloading
# change to only 1 loading. merge clinical and raw tcga
#########################
#########################
#########################
#########################



#########################
# loading raw LUAD data and normaizing
load("/inputs/LUAD_exprs.rd")
geneData <- expr_mat
rm(expr_mat)

# reading 'junk' genes
t1=read.table("/inputs/excluded_genes.txt",stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

junk <- c(gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)
junk <- intersect(junk,rownames(geneData))

# reading signatures and creating a normlized dataset of genes
t1=read.table("/inputs/TCGA_sigs.txt",stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

dups <-  names(table(colnames(geneData))[table(colnames(geneData))==2])
for (dup in dups){
  dup_index <- grep(dup, colnames(geneData))
  geneData[,dup] <- rowSums(geneData[,dup_index])
}

reg <- 1E-8
geneData <- geneData[!rownames(geneData)%in%junk,]
norm_dat <- t(t(geneData)/rowSums(t(geneData)))
norm_dat <- log2(norm_dat + reg)
norm_dat <- norm_dat[,names(table(colnames(norm_dat)))]

norm_dat_gene_scaled <- t(scale(t(norm_dat)))
rm(geneData)

#########################
# removing 'junk' genes from LUAD data and creating signature scores

skips <- c("Junk", "RPs", "Ribosome", "IgVar", "Mito", "Metallo")
skipLength <- length(skips)

sig_matrix <- matrix(data = NA, nrow=(length(gl)-skipLength), ncol=ncol(norm_dat_gene_scaled))
rownames(sig_matrix) <- c(names(gl)[1:(length(gl)-skipLength)])
colnames(sig_matrix) <- colnames(norm_dat_gene_scaled)

for (num in 1:length(gl)){
  sig_name <- names(gl)[num]
  if  (sig_name%in%skips){
    next
  }
  sig <- gl[[num]]
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
# figure 3C

merged_data <- read.csv("/inputs/merged_data_luad.csv", row.names = 1)

caf_genes <- c(gl$ADH1B_CAF, gl$FAP_CAF)
new_dat <- norm_dat_gene_scaled[rownames(norm_dat_gene_scaled)%in%caf_genes,]

caf_sig <- c("ADH1B_CAF","FAP_CAF")
new_dat <- rbind(new_dat, sig_matrix[caf_sig,])
new_dat <- new_dat[c(caf_genes, caf_sig),]
new_dat <- t(scale(t(new_dat)))

breaks = seq(-2.5,2.5,length.out=10)

new_dat <- t(new_dat)
new_dat <- new_dat[rownames(new_dat)%in%rownames(merged_data),]
merged_data <- merged_data[rownames(merged_data)%in%rownames(new_dat),]
new_dat <- cbind(new_dat, merged_data[,c("tumor_nuclei","sample_type_id")])

# selecting only primary tumor
new_dat <- new_dat[new_dat[,"sample_type_id"]%in%c("1"),]
new_dat[,"tumor_nuclei"] <- scale(new_dat[,"tumor_nuclei"])
new_dat <- t(new_dat)

for (sig in caf_sig){
  new_dat <- new_dat[,order(new_dat[sig,])]
  
  temp_rownames <- rownames(new_dat)
  temp_dat <- new_dat[c(temp_rownames[1:(length(temp_rownames)-4)],"tumor_nuclei"),]
  
  pheatmap(temp_dat, cluster_rows=F, cluster_cols=F, scale = "none", show_colnames=F,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
           breaks = breaks)
}

rm(new_dat)

#########################
# figure S5D

all_genes <- c()
for (geneSet in gl){
  all_genes <- c(all_genes, geneSet)
}
all_genes <- intersect(all_genes,rownames(norm_dat_gene_scaled))

merged_data <- read.csv("/inputs/merged_data_luad.csv", row.names = 1)
tumor_samps <- rownames(merged_data)[merged_data[,"sample_type_id"]%in%c("1")]
norm_dat_gene_scaled <- norm_dat_gene_scaled[,tumor_samps]

gene_cors <- cor(t(norm_dat_gene_scaled[all_genes,]))
na_genes <- !is.na(gene_cors[1,])
gene_cors <- as.matrix(gene_cors[na_genes,na_genes])
genes_of_int <- c(gl$Fib,gl$ADH1B_CAF,gl$FAP_CAF)
gene_cors <- gene_cors[genes_of_int,genes_of_int]

corrplot(gene_cors, method = 'color', order = "original", type = "upper", number.font = 1,
         tl.col = "black", tl.srt = 45, col=rev(brewer.pal(n=8, name="RdBu")))
rm(gene_cors)

#########################
# figure 3D

merged_data <- read.csv("/inputs/merged_data_luad.csv", row.names = 1)

# selecting unique subtypes with > 10 patients
names_data <- c("Acinar predominant Adc","Papillary predominant Adc",
                "Solid predominant Adc","Micropapillary predom Adc")
merged_data[,"Pathology.Result.Updated.10.1"] <- as.character(merged_data[,"Pathology.Result.Updated.10.1"])
# renaming all unused subtypes as 'Other'
merged_data[merged_data[,"Pathology.Result.Updated.10.1"]%notin%names_data,"Pathology.Result.Updated.10.1"] <- "Other"
# selecting only primary tumor
merged_data <- merged_data[merged_data[,"sample_type_id"]%in%c("1"),]

score_of_interest <- c("ADH1B_CAF","FAP_CAF")
for (score in score_of_interest){
  merged_data[,score] <- scale(merged_data[,score])
}
merged_data <- merged_data[,c(score_of_interest,"Pathology.Result.Updated.10.1")]

# removing other pathologies
merged_data <- merged_data[merged_data[,"Pathology.Result.Updated.10.1"]%in%names_data,]

colnames(merged_data) <- c(colnames(merged_data)[1:length(score_of_interest)],"Pathology")
melted_data <- melt(merged_data)
summarized_melted <- summarySE(melted_data, measurevar="value", groupvars=c("Pathology","variable"))

cols <- c("khaki3","red")

ggplot(summarized_melted, aes(x=factor(Pathology, level = names_data), y=value, fill=variable)) +
  geom_errorbar(data=summarized_melted, mapping=aes(x=factor(Pathology, level = names_data), ymin=value-se, ymax=value+se),
                size = 1.5, width = .6, position = "dodge", color = rep(cols,4)) +
  geom_point(size = 2, position = position_dodge(width = .6)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.ontop = TRUE
  ) +
  ylim(-1.1,1)
  # theme_bw()
  
for (pathology in names_data){
  temp_data <- as.matrix(merged_data[merged_data[,"Pathology"] == pathology,score_of_interest])
  melted_data <- melt(temp_data)[,2:3]
  melted_data$X2 <- ordered(melted_data$X2, levels = score_of_interest)
  
  x <- ggboxplot(melted_data, x = "X2", y = "value", 
                 color = "X2", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 order = score_of_interest,
                 ylab = "value", xlab = "Treatment", main = pathology)
  # plot(x)
  
  res.aov <- aov(value ~ X2, data = melted_data)
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  print("")
  print("")
}

#########################
# figure 3D

merged_data <- read.csv("/inputs/merged_data_luad.csv", row.names = 1)
# stages 3 and above are listed as 3 in merged_data file

merged_data[,"tumor_stage_num"] <- as.integer(merged_data[,"tumor_stage_num"])
# selecting only primary tumor
merged_data <- merged_data[merged_data[,"sample_type_id"]%in%c("1"),]

score_of_interest <- c("ADH1B_CAF","FAP_CAF")
for (score in score_of_interest){
  merged_data[,score] <- scale(merged_data[,c(score)])
}
merged_data <- merged_data[,c(score_of_interest,"tumor_stage_num")]
merged_data <- merged_data[!is.na(merged_data[,"tumor_stage_num"]),]

merged_data[merged_data[,"tumor_stage_num"] == 3,"tumor_stage_num"]  <- 2

stages <- c(1,2)

# changing column name from "tumor_stage_num" to "stages"
colnames(merged_data) <- c(colnames(merged_data)[1:length(score_of_interest)],"stages")

merged_data[,"stages"] <- as.character(merged_data[,"stages"])

melted_data <- melt(merged_data, id="stages")
summarized_melted <- summarySE(melted_data, measurevar="value", groupvars=c("stages","variable"))

cols <- c("khaki3","red")

# if adjusting number of stanges, "color = rep(cols,2))" must be changed
ggplot(summarized_melted, aes(x=stages, y=value, fill=variable)) +
  geom_errorbar(data=summarized_melted, mapping=aes(x=stages, ymin=value-se, ymax=value+se),
                size = 1.5, width = .6, position = "dodge", color = rep(cols,2)) +
  geom_point(size = 2, position = position_dodge(width = .6)) +
  ylim(-0.5,0.5) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.ontop = TRUE
  )

t.test(merged_data[merged_data[,"stages"] == 2,"ADH1B_CAF"],merged_data[merged_data[,"stages"] == 2,"FAP_CAF"])

for (stage in stages){
  temp_data <- as.matrix(merged_data[merged_data[,"stages"] == stage,score_of_interest])
  melted_data <- melt(temp_data)[,2:3]
  melted_data$X2 <- ordered(melted_data$X2, levels = score_of_interest)
  x <- ggboxplot(melted_data, x = "X2", y = "value", 
                 color = "X2", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 order = score_of_interest,
                 ylab = "value", xlab = "Treatment", main = stage)
  # plot(x)
  
  res.aov <- aov(value ~ X2, data = melted_data)
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  print("")
  print("")
}

#########################
# figure 3G

merged_data <- read.csv("/inputs/merged_data_luad.csv", row.names = 1, stringsAsFactors = F)
sigs <- c("Fib","ADH1B_CAF","FAP_CAF","LCAM","delta")

# selecting only primary tumor
merged_data <- merged_data[merged_data[,"sample_type_id"]%in%c("1"),]

# various residuals, etc
delta_caf <- merged_data[,"FAP_CAF"] - merged_data[,"ADH1B_CAF"]

delta_to_lcam_residuals <- lm(merged_data[,"LCAM"] ~ delta_caf)

merged_data <- cbind(merged_data,NA,NA)
merged_data[,"delta"] <- delta_caf
merged_data[,sigs] <- scale(merged_data[,sigs])

summary(delta_to_lcam_residuals)
#########################

#########################
#########################
#########################
#########################
sqrt(0.4292)
#########################
#########################
#########################
#########################
ggplot(merged_data, aes(x=delta, y=LCAM)) +
  geom_point(size=3) +
  theme_bw() +
  ggtitle("R = 0.66")



