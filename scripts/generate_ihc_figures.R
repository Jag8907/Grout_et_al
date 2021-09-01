`%notin%` <- Negate(`%in%`)

library(gplots)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(RColorBrewer)

###################
# adh1b versus fap staining

plus_stains <- read.csv("intermediate_files/S5.csv",stringsAsFactors = F,fileEncoding="UTF-8-BOM")

ggplot(plus_stains, aes(x = ADH1B_coverage_stroma/100, y = adjusted_FAP/100, color = predominant_CAF)) +
  geom_point(aes(shape = MYH11_aSMA_CAF), size = 3) +
  scale_color_manual(values = c("ADH1B" = "forestgreen", "FAP" = "red")) +
  scale_shape_manual(values=c(16,1)) +
  theme_bw() +
  scale_x_continuous(trans='log10', limits = c(0.001, 1)) +
  scale_y_continuous(trans='log10', limits = c(0.001, 1))

# power test showing that low probability of finding samples above 50 & 50
thresholds <- seq(5,75,5)

for (thresh in thresholds){
  print(thresh)
  print(mean(plus_stains$adjusted_FAP > thresh)*mean(plus_stains$ADH1B_coverage_stroma > thresh)*nrow(plus_stains))
}

thresholds <- seq(5,50,5)

# testing by hypergeometric
obsv_p_vals <- list()
for (thresh in thresholds){
  # thresh <- runif(1)*100
  # upper right quadrant [2,2]
  x = sum(plus_stains$adjusted_FAP > thresh & plus_stains$ADH1B_coverage_stroma > thresh)
  # positive y axis samples over threshold, if adjusted_FAP is y axis
  m = sum(plus_stains$adjusted_FAP > thresh)
  # samples below threshold on y axis, if adjusted_FAP is y axis
  n = sum(plus_stains$adjusted_FAP <= thresh)
  # total sample numbers
  k = sum(plus_stains$ADH1B_coverage_stroma > thresh)
  
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

sum_tissue <- read.csv("intermediate_files/S3.csv", row.names = 1)

for (fib in colnames(sum_tissue)){
  if (fib == "tissue"){
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

quant_table <- "intermediate_files/S7.csv"
infil_table <- read.csv(quant_table)
colnames(infil_table) <- c("patient",c(colnames(infil_table)[2:ncol(infil_table)]))

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
  asma_infil <-  infil_table[!is.na(infil_table[,"asma_stroma_coverage"]),]
  g <- ggplot(asma_infil, aes(x=asma_stroma_coverage, y=get(quant))) +
    theme_bw() +
    geom_point(size = 3) +
    labs(title = quant, x = "", y = "") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 40)) +
    scale_y_continuous(trans='log10')
  
  print(g)
  print(cor.test(asma_infil[,"asma_stroma_coverage"],asma_infil[,quant]))
  
  # plotting adh1b versus fap
  g <- ggplot(infil_table,aes(x=CAF,y=get(quant),fill=CAF))
  g <- g + geom_dotplot(binaxis="y",stackdir="center", binwidth=binwidth) +
    coord_cartesian(ylim = ylim, expand = TRUE) +
    theme_bw() +
    labs(title = quant, x = "", y = "") +
    scale_fill_brewer(palette="Dark2") +
    scale_y_continuous(trans='log10')
  
  # print(g)
  
  # plotting adh1b MYH11_aSMA_CAF + versus -
  new_infil_table <- melt(infil_table[,c("CAF","MYH11_aSMA_CAF",quant)])
  new_infil_table <- new_infil_table[new_infil_table[,"CAF"]=="ADH1B",]
  
  temp <- ggplot(new_infil_table, aes(x=MYH11_aSMA_CAF, y=value)) +
    coord_cartesian(ylim = ylim, expand = TRUE) +
    theme_bw() +
    labs(title = quant) +
    geom_dotplot(binaxis='y', stackdir='center', 
                 position=position_dodge(0.8), binwidth=binwidth) +
    scale_fill_brewer(palette="Dark2") +
    scale_y_continuous(trans='log10')
  
  # print(temp)
  
  # tumor cd8 density t test for adh1b MYH11_aSMA_CAF
  adh1b_infil <- infil_table[infil_table[,"CAF"] == "ADH1B",]
  temp_test <- t.test(adh1b_infil[adh1b_infil[,"MYH11_aSMA_CAF"] == "ABSENT",quant],
                      adh1b_infil[adh1b_infil[,"MYH11_aSMA_CAF"] == "PRESENT",quant])
  
  # print(temp_test$p.value)
}





###################
# infiltration graphing and analysis
tiling_table <- "intermediate_files/S10.csv"

# or

tiling_table <- "intermediate_files/S11.csv"

tiling_table <- read.csv(tiling_table)

vars <- c("tumor_stroma_log_cd8_ratio","stromal_CD3_density","tumor_CD3_density","tumor_stroma_log_cd3_ratio")

g <- ggplot(tiling_table, aes(x=MYH11.aSMA..scores, y=stromal_CD8_density)) +
  theme_bw() +
  geom_point() +
  ylab("stromal_CD8_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000 ))
print(cor.test(tiling_table[,"MYH11.aSMA..scores"], tiling_table[,"stromal_CD8_density"]))
print(g)

g <- ggplot(tiling_table, aes(x=MYH11.aSMA..scores, y=tumor_CD8_density)) +
  theme_bw() +
  geom_point() +
  ylab("tumor_CD8_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2000 ))
print(cor.test(tiling_table[,"MYH11.aSMA..scores"], tiling_table[,"tumor_CD8_density"]))
print(g)

######

g <- ggplot(tiling_table, aes(x=MYH11.aSMA..scores, y=stromal_CD3_density)) +
  theme_bw() +
  geom_point() +
  ylab("stromal_CD3_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10000 ))
print(cor.test(tiling_table[,"MYH11.aSMA..scores"], tiling_table[,"stromal_CD3_density"]))
print(g)

g <- ggplot(tiling_table, aes(x=MYH11.aSMA..scores, y=tumor_CD3_density)) +
  theme_bw() +
  geom_point() +
  ylab("tumor_CD3_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000 ))
print(cor.test(tiling_table[,"MYH11.aSMA..scores"], tiling_table[,"tumor_CD3_density"]))
print(g)





###################
# infiltration graphing and analysis
tiling_table <- "intermediate_files/S8.csv"

# or

tiling_table <- "intermediate_files/S9.csv"

tiling_table <- read.csv(tiling_table)

tiling_table <- tiling_table[!is.na(tiling_table[,"stromal_CD3_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"tumor_CD3_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"stromal_CD8_density"]),]
tiling_table <- tiling_table[!is.na(tiling_table[,"tumor_CD8_density"]),]

tiling_table <- tiling_table[!is.na(tiling_table[,"aSMA_stroma_area"]),]

g <- ggplot(tiling_table, aes(x=aSMA_stroma_area, y=get("stromal_CD8_density"))) +
  theme_bw() +
  geom_point() +
  ylab("stromal_CD8_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000))
cor.test(tiling_table[,"aSMA_stroma_area"], tiling_table[,"stromal_CD8_density"])
print(g)

g <- ggplot(tiling_table, aes(x=aSMA_stroma_area, y=get("stromal_CD3_density"))) +
  theme_bw() +
  geom_point() +
  ylab("stromal_CD3_density") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6000))
cor.test(tiling_table[,"aSMA_stroma_area"], tiling_table[,"stromal_CD3_density"])
print(g)

cor.test(tiling_table[,"aSMA_stroma_area"], tiling_table[,"tumor_CD3_density"])
cor.test(tiling_table[,"aSMA_stroma_area"], tiling_table[,"tumor_CD8_density"])








