`%notin%` <- Negate(`%in%`)

library(gplots)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(RColorBrewer)



###################
# adh1b versus fap staining

plus_stains <- read.csv("/inputs/S5.csv",stringsAsFactors = F,fileEncoding="UTF-8-BOM")

ggplot(plus_stains, aes(x = ADH1B/100, y = Adj.FAP/100, color = stroma)) +
  geom_point(aes(shape = barrier), size = 3) +
  scale_color_manual(values = c("a" = "forestgreen", "f" = "red")) +
  scale_shape_manual(values=c(16,1)) +
  theme_bw() +
  # geom_hline(yintercept = 0, linetype = "solid") +
  # geom_vline(xintercept = 0, linetype = "solid") +
  scale_x_continuous(trans='log10', limits = c(0.001, 1)) +
  scale_y_continuous(trans='log10', limits = c(0.001, 1))
  # scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  # geom_hline(yintercept = 0.2, linetype = "dashed", color="gray") +
  # geom_vline(xintercept = 0.2, linetype = "dashed", color="gray")
  # xlim(0,1) +
  # ylim(0,1)

# power test showing that low probability of finding samples above 50 & 50
thresholds <- seq(5,75,5)

for (thresh in thresholds){
  print(thresh)
  print(mean(plus_stains$Adj.FAP > thresh)*mean(plus_stains$ADH1B > thresh)*nrow(plus_stains))
}

thresholds <- seq(5,50,5)

# testing by hypergeometric
obsv_p_vals <- list()
for (thresh in thresholds){
  # thresh <- runif(1)*100
  # upper right quadrant [2,2]
  x = sum(plus_stains$Adj.FAP > thresh & plus_stains$ADH1B > thresh)
  # positive y axis samples over threshold, if adj.fap is y axis
  m = sum(plus_stains$Adj.FAP > thresh)
  # samples below threshold on y axis, if adj.fap is y axis
  n = sum(plus_stains$Adj.FAP <= thresh)
  # total sample numbers
  k = sum(plus_stains$ADH1B > thresh)
  
  obsv_p_vals[as.character(thresh)] <- phyper(x, m, n, k, lower.tail = TRUE, log.p = FALSE)
}

# p values
obsv_p_vals
# p values with correction
for (num in 1:length(obsv_p_vals)){
  print(obsv_p_vals[num])
  print(obsv_p_vals[[num]]*length(thresholds))
}





###################
# cluster adjacent versus tumor enrichment

sum_tissue <- read.csv("/inputs/S2.csv")

fibs <- c(9,8,4,20,24,26,22,19,11)
for (fib in colnames(sum_tissue)){
  if (fib == "X" | fib == "tissue"){
    next
  }
  tumor_fib <- sum_tissue[sum_tissue[,"tissue"]=="tumor",fib]
  adj_fib <- sum_tissue[sum_tissue[,"tissue"]=="adjacent",fib]
  
  print(fib)
  print(wilcox.test(adj_fib,tumor_fib, exact=F))
}





###################
# infiltration graphing and analysis
brewer.pal(n = 8, name = "Dark2")

quant_table <- "/inputs/S7.csv"
infil_table <- read.csv(quant_table)
colnames(infil_table) <- c("patient",c(colnames(infil_table)[2:ncol(infil_table)]))

quant_type <- c("tum_strom_ratio_cd8","tumor_cd8","stroma_cd8","tum_strom_ratio_cd3","tumor_cd3","stroma_cd3")
quant_type <- c("tum_strom_ratio_cd8","tumor_cd8","tum_strom_ratio_cd3","tumor_cd3")

for (quant in quant_type){
  
  print(quant)
  
  if (quant == "tum_strom_ratio_cd8"){
    binwidth=0.07
    ylim = c(0.05,1.5)
  }else if(quant == "tumor_cd8"){
    binwidth=0.1
    ylim = c(10,1000)
  }else if(quant == "tum_strom_ratio_cd3"){
    binwidth=0.06
    ylim = c(0.05,1)
  }else{
    binwidth=0.1
    ylim = c(30,3000)
  }
  
  # plotting aSMA versus quant
  asma_infil <-  infil_table[!is.na(infil_table[,"asma"]),]
  g <- ggplot(asma_infil, aes(x=asma, y=get(quant))) +
    theme_bw() +
    geom_point(size = 3) +
    labs(title = quant, x = "", y = "") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 40)) +
    scale_y_continuous(trans='log10')
  
  print(quant)
  print(g)
  print(cor.test(asma_infil[,"asma"],asma_infil[,quant]))
  
  # plotting adh1b versus fap
  g <- ggplot(infil_table,aes(x=CAF,y=get(quant),fill=CAF))
  g <- g + geom_dotplot(binaxis="y",stackdir="center", binwidth=binwidth) +
    coord_cartesian(ylim = ylim, expand = TRUE) +
    theme_bw() +
    labs(title = quant, x = "", y = "") +
    scale_fill_brewer(palette="Dark2") +
    scale_y_continuous(trans='log10')
  
  # print(g)
  
  # plotting adh1b barrier + versus -
  new_infil_table <- melt(infil_table[,c("CAF","barrier",quant)])
  new_infil_table <- new_infil_table[new_infil_table[,"CAF"]=="a",]
  
  temp <- ggplot(new_infil_table, aes(x=barrier, y=value)) +
    coord_cartesian(ylim = ylim, expand = TRUE) +
    theme_bw() +
    labs(title = quant) +
    geom_dotplot(binaxis='y', stackdir='center', 
                 position=position_dodge(0.8), binwidth=binwidth) +
    scale_fill_brewer(palette="Dark2") +
    scale_y_continuous(trans='log10')
  
  # print(temp)
  
  # tumor cd8 density t test for adh1b barrier
  adh1b_infil <- infil_table[infil_table[,"CAF"] == "a",]
  temp_test <- t.test(adh1b_infil[adh1b_infil[,"barrier"] == "n",quant],
                      adh1b_infil[adh1b_infil[,"barrier"] == "y",quant])
  
  # print(temp_test$p.value)
}





###################
# infiltration graphing and analysis
tiling_table <- "/inputs/S10.csv"

# or

tiling_table <- "/inputs/S11.csv"

tiling_table <- read.csv(tiling_table)

vars <- c("TS_density","Stromal_CD3_density","Tumor_CD3_density","TS_CD3_ratio")

g <- ggplot(tiling_table, aes(x=Barrier_scores, y=Stromal_CD8_density)) +
  theme_bw() +
  geom_point() +
  ylab("Stromal_CD8_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000 ))
print(cor.test(tiling_table[,"Barrier_scores"], tiling_table[,"Stromal_CD8_density"]))
print(g)

g <- ggplot(tiling_table, aes(x=Barrier_scores, y=Tumor_CD8_density)) +
  theme_bw() +
  geom_point() +
  ylab("Tumor_CD8_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2000 ))
print(cor.test(tiling_table[,"Barrier_scores"], tiling_table[,"Tumor_CD8_density"]))
print(g)

######

g <- ggplot(tiling_table, aes(x=Barrier_scores, y=Stromal_CD3_density)) +
  theme_bw() +
  geom_point() +
  ylab("Stromal_CD3_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10000 ))
print(cor.test(tiling_table[,"Barrier_scores"], tiling_table[,"Stromal_CD3_density"]))
print(g)

g <- ggplot(tiling_table, aes(x=Barrier_scores, y=Tumor_CD3_density)) +
  theme_bw() +
  geom_point() +
  ylab("Tumor_CD3_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000 ))
print(cor.test(tiling_table[,"Barrier_scores"], tiling_table[,"Tumor_CD3_density"]))
print(g)





###################
# infiltration graphing and analysis
tiling_table <- "/inputs/S8.csv"

# or

tiling_table <- "/inputs/S9.csv"

tiling_table <- read.csv(tiling_table)

tiling_table <- tiling_table[!is.na(tiling_table[,"Stromal_CD3_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"Tumor_CD3_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"Stromal_CD8_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"Tumor_CD8_density"]),]

tiling_table <- tiling_table[!is.na(tiling_table[,"ASMA_stroma_area"]),]

g <- ggplot(tiling_table, aes(x=ASMA_stroma_area, y=get("Stromal_CD8_density"))) +
  theme_bw() +
  geom_point() +
  ylab("Stromal_CD8_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000))
cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"Stromal_CD8_density"])
print(g)

g <- ggplot(tiling_table, aes(x=ASMA_stroma_area, y=get("Stromal_CD3_density"))) +
  theme_bw() +
  geom_point() +
  ylab("Stromal_CD3_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6000))
cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"Stromal_CD3_density"])
print(g)

tiling_table <- "D:/Dropbox/Grout et al/figures/supplemental info/table_8_asma_caf_tiling_24.csv"
tiling_table <- read.csv(tiling_table)
tiling_table <- tiling_table[!is.na(tiling_table[,"Tumor_CD8_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"Tumor_CD3_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"ASMA_stroma_area"]),]

cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"Tumor_CD3_density"])
cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"Tumor_CD8_density"])

tiling_table <- "D:/Dropbox/Grout et al/figures/supplemental info/table_9_asma_caf_tiling_54.csv"
tiling_table <- read.csv(tiling_table)
tiling_table <- tiling_table[!is.na(tiling_table[,"Tumor_CD8_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"Tumor_CD3_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"ASMA_stroma_area"]),]

cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"Tumor_CD3_density"])
cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"Tumor_CD8_density"])






