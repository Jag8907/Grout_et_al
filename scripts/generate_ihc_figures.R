`%notin%` <- Negate(`%in%`)

library(gplots)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(RColorBrewer)

###################

# aSMA coverage bar plot

main_path <- "D:/Dropbox/Grout et al/0. REVISIONS for CANCER DISCOVERY/Response to reviewers/Supplementary files"
bar_fap_asma <- read.csv(paste(main_path, "/Table S5.csv", sep= ""), stringsAsFactors = F, fileEncoding="UTF-8-BOM")

temp_mat <- matrix(NA, nrow = nrow(bar_fap_asma), ncol = 4)
colnames(temp_mat) <- c("patient_ID","aSMA_pos","neg","adjusted_FAP")

temp_mat[,"aSMA_pos"] <- bar_fap_asma[,"FAP_pos_aSMA_pos"]
temp_mat[,"neg"] <- 100-temp_mat[,"aSMA_pos"]
temp_mat[,"patient_ID"] <- bar_fap_asma[,"patient_ID"]
temp_mat[,"adjusted_FAP"] <- bar_fap_asma[,"adjusted_FAP"]

temp_mat <- temp_mat[!is.na(temp_mat[,"aSMA_pos"]),]

temp_mat <- temp_mat[order(temp_mat[,"patient_ID"]),]

temp_mat <- melt(temp_mat[,c("aSMA_pos","neg")])
colnames(temp_mat) <- c("patient_ID","stain","value")

temp_mat[,"value"] <- temp_mat[,"value"]/100

ggplot(temp_mat, aes(fill=stain, y=value, x=patient_ID)) + 
  geom_bar(position="fill", stat="identity", width = 0.75) +
  theme_bw() +
  scale_fill_manual("legend", values = c("asmapos" = "red", "neg" = "black"))

###################
# adh1b versus fap staining

main_path <- "D:/Dropbox/Grout et al/0. REVISIONS for CANCER DISCOVERY/Response to reviewers/Supplementary files"
quant_table <- paste(main_path, "/Table S5.csv", sep= "")

plus_stains <- read.csv(quant_table,stringsAsFactors = F,fileEncoding="UTF-8-BOM")

ggplot(plus_stains, aes(x = ADH1B_coverage_stroma/100, y = adjusted_FAP/100, color = predominant_CAF)) +
  geom_point(aes(shape = MYH11_aSMA_CAF), size = 3) +
  scale_color_manual(values = c("ADH1B" = "forestgreen", "FAP" = "red")) +
  scale_shape_manual(values=c(16,1)) +
  theme_bw() +
  scale_x_continuous(trans='log10', limits = c(0.001, 1)) +
  scale_y_continuous(trans='log10', limits = c(0.001, 1)) +
  theme(panel.grid.minor = element_blank())

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

main_path <- "D:/Dropbox/Grout et al/0. REVISIONS for CANCER DISCOVERY/Response to reviewers/Supplementary files"

# quant_table <- paste(main_path, "/Table S7 foxp3 update.csv", sep= "")
quant_table <- paste(main_path, "/Table S7.csv", sep= "")

infil_table <- read.csv(quant_table)
colnames(infil_table) <- c("patient",c(colnames(infil_table)[2:ncol(infil_table)]))

quant_type <- c("tum_strom_ratio_cd8","tumor_cd8","stroma_cd8",
                "tum_strom_ratio_cd3","tumor_cd3","stroma_cd3",
                "tum_strom_ratio_foxp3","tumor_foxp3","stroma_foxp3")

quant_type <- c("tum_strom_ratio_cd3","tumor_cd3","tum_strom_ratio_cd8","tumor_cd8")


quant <- "tum_strom_ratio_cd3"
# quant <- "tumor_foxp3"

log_graphs_e <- -4:3
breaks <- 10^(log_graphs_e)
minor_breaks <- rep(1:9, length(breaks)-1)*(10^rep(log_graphs_e[1:length(breaks)-1], each=9))

for (quant in quant_type){
  
  print(quant)
  
  if (quant == "tum_strom_ratio_cd8"){
    binwidth=0.08
    ylim = c(0.03,1)
    
  }else if(quant == "tumor_cd8"){
    binwidth=0.1
    ylim = c(10,1000)
    
  }else if(quant == "stroma_cd8"){
    binwidth=0.08
    ylim = c(100,3000)
  
  }else if(quant == "tum_strom_ratio_cd3"){
    binwidth=0.08
    ylim = c(0.03,1)
    
  }else if(quant == "tumor_cd3"){
    binwidth=0.09
    ylim = c(30,3000)
    
  }else if(quant == "stroma_cd3"){
    binwidth=0.1
    ylim = c(100,10000)
    
  }else if(quant == "tum_strom_ratio_foxp3"){
    binwidth=0.08
    ylim = c(0.03,1)
    
  }else if(quant == "tumor_foxp3"){
    binwidth=0.08
    ylim = c(10,300)
    
  }else if(quant == "stroma_foxp3"){
    binwidth=0.12
    ylim = c(10,1000)
  }
  
  # plotting aSMA versus quant
  asma_infil <-  infil_table[!is.na(infil_table[,"asma_stroma_coverage"]),]
  g <- ggplot(asma_infil, aes(x=asma_stroma_coverage, y=get(quant))) +
    coord_cartesian(ylim = ylim, expand = TRUE) +
    theme_bw() +
    geom_point(size = 3) +
    labs(title = quant, x = "", y = "") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 40)) +
    scale_y_continuous(trans='log10') +
    theme(panel.grid.minor = element_blank())
  
  print(g)
  temp <- cor.test(asma_infil[,"asma_stroma_coverage"],asma_infil[,quant])
  print(temp$p.value)
  print(temp$estimate)
  
  # plotting adh1b versus fap
  g <- ggplot(infil_table,aes(x=CAF,y=get(quant),fill=CAF))
  g <- g + geom_dotplot(binaxis="y",stackdir="center",binwidth=binwidth,stackratio=0.8,dotsize=1) +
    scale_fill_manual(values=c("ADH1B"="grey", "FAP"="black")) +
    coord_cartesian(ylim = ylim, expand = TRUE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(title = quant, x = "", y = "") +
    # scale_fill_brewer(palette="Dark2") +
    theme(legend.position="none") +
    scale_y_continuous(trans='log10')
  
  print(g)
  
  temp_adh1b <- infil_table[,"CAF"]=="ADH1B"
  temp_fap <- infil_table[,"CAF"]=="FAP"
  
  temp_test <- t.test(infil_table[temp_adh1b,quant],infil_table[temp_fap,quant])
  print(paste(quant, "adh1b v fap p value", round(temp_test$p.value,4)))
  
  # plotting adh1b MYH11_aSMA_CAF + versus -
  new_infil_table <- melt(infil_table[,c("stage","CAF","MYH11_aSMA_CAF",quant)])
  # new_infil_table <- new_infil_table[new_infil_table[,"CAF"]=="ADH1B",]
  new_infil_table <- new_infil_table[new_infil_table[,"stage"]=="low",]
  
  temp <- ggplot(new_infil_table, aes(x=MYH11_aSMA_CAF, y=value)) +
    coord_cartesian(ylim = ylim, expand = TRUE) +
    theme_bw() +
    labs(title = quant) +
    geom_dotplot(binaxis='y', stackdir='center', 
                 position=position_dodge(0.8), binwidth=binwidth) +
    # scale_fill_brewer(palette="Dark2") +
    scale_y_continuous(trans='log10') +
    theme(panel.grid.minor = element_blank())
  
  print(temp)
  
  # tumor cd8 density t test for adh1b MYH11_aSMA_CAF
  temp_test <- t.test(new_infil_table[new_infil_table[,"MYH11_aSMA_CAF"] == "ABSENT","value"],
                      new_infil_table[new_infil_table[,"MYH11_aSMA_CAF"] == "PRESENT","value"])
  
  print(paste(quant, "myh11 p value", round(temp_test$p.value,4)))
}






###################
# collagen 11 and 12
# infiltration graphing and analysis

pat_38_col11 <- "G:/My Drive/Paper website/Scripts/intermediate_files/S13.csv"
# stromal_aSMA_coverage
# stromal_Coll11_coverage
# stromal_Coll4_coverage
# stromal_FAP_coverage

pat_38_col12 <- "G:/My Drive/Paper website/Scripts/intermediate_files/S14.csv"
# stromal_aSMA_coverage
# stromal_Coll12A_coverage
# stromal_FAP_coverage

pat_58_col11 <- "G:/My Drive/Paper website/Scripts/intermediate_files/S15.csv"
# stromal_aSMA_coverage
# stromal_Coll12A_coverage
# stromal_Coll4_coverage
# stromal_FAP_coverage

table_names <- c("pat_38_col11","pat_38_col12","pat_58_col11")
tables <- c(pat_38_col11,pat_38_col12,pat_58_col11)

collagens <- c("stromal_Coll12A_coverage","stromal_Coll4_coverage","stromal_Coll11_coverage")

# tiling_table[tiling_table < 1e-2] <- 1e-2

n <- 1

for (table in tables){
  
  table_name <- table_names[n]
  n <- n + 1
  
  tiling_table <- table
  tiling_table <- read.csv(tiling_table, row.names = 1)
  vars <- colnames(tiling_table)
  temp_cols <- collagens[collagens%in%vars]
  
  for (col in temp_cols){
    
    if (col == "stromal_Coll12A_coverage"){
      
      ylimit <- c(0,60)
      yint <- 60
      
    }else if(col == "stromal_Coll4_coverage"){
      
      ylimit <- c(0,50)
      yint <- 50
      
    }else{
      
      ylimit <- c(0,12)
      yint <- 12.5
      
    }
    
    xlimit_fap <- c(0,50)
    
    g <- ggplot(tiling_table, aes(x=stromal_aSMA_coverage, y=get(col))) +
      theme_bw() +
      geom_point(size=2) +
      geom_vline(xintercept = 50, color = "black") +
      geom_hline(yintercept = yint, color = "black") +
      ylim(ylimit) +
      xlim(xlimit_fap) +
      ylab(col) +
      ggtitle(paste(table_name, max(tiling_table[,"stromal_aSMA_coverage"]), max(tiling_table[,col]))) +
      theme(panel.grid.minor = element_blank())
    print(col)
    print(cor.test(tiling_table[,"stromal_aSMA_coverage"], tiling_table[,col], method = "spearman"))
    print(g)
    
  }
}







for (col in collagens){
  
  if (col == "stromal_Coll12A_coverage"){
    xlimit <- c(1e-2,50)
    ylimit <- c(1e-2,10)
    
    g <- ggplot(tiling_table, aes(x=stromal_aSMA_coverage, y=get(col))) +
      theme_bw() +
      geom_point(size=2) +
      ylab(col) +
      theme(panel.grid.minor = element_blank())
      # scale_x_log10(limits = xlimit) +
      # scale_y_log10(limits = ylimit)
    print(cor.test(tiling_table[,"stromal_aSMA_coverage"], tiling_table[,col]))
    print(g)
    
  }else if(col == "stromal_Coll4_coverage"){
    xlimit <- c(1e-2,50)
    ylimit <- c(1e-2,10)
    
    g <- ggplot(tiling_table, aes(x=stromal_aSMA_coverage, y=get(col))) +
      theme_bw() +
      geom_point(size=2) +
      ylab(col) +
      theme(panel.grid.minor = element_blank())
      # scale_x_log10(limits = xlimit) +
      # scale_y_log10(limits = ylimit)
    print(cor.test(tiling_table[,"stromal_aSMA_coverage"], tiling_table[,col]))
    print(g)
    
  }else{
    xlimit <- c(1e-2,50)
    ylimit <- c(1e-2,50)
    
    g <- ggplot(tiling_table, aes(x=stromal_aSMA_coverage, y=get(col))) +
      theme_bw() +
      geom_point(size=2) +
      ylab(col) +
      theme(panel.grid.minor = element_blank())
      # scale_x_log10(limits = xlimit) +
      # scale_y_log10(limits = ylimit)
    print(cor.test(tiling_table[,"stromal_aSMA_coverage"], tiling_table[,col]))
    print(g)
    
  }
}

g <- ggplot(tiling_table, aes(x=ASMA_stroma_area, y=stromal_Coll11_coverage)) +
  theme_bw() +
  geom_point(size=2) +
  ylab("stromal_Coll11_coverage") +
  theme(panel.grid.minor = element_blank()) +
  scale_x_log10(limits = c(1e-2,10)) +
  scale_y_log10(limits = c(1e-2,10))
  # scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 4000 ))
print(cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"stromal_Coll11_coverage"]))
print(g)

g <- ggplot(tiling_table, aes(x=ASMA_stroma_area, y=stromal_Coll4_coverage)) +
  theme_bw() +
  geom_point(size=2) +
  ylab("stromal_Coll4_coverage") +
  theme(panel.grid.minor = element_blank()) +
  scale_x_log10(limits = c(1e-2,10)) +
  scale_y_log10(limits = c(1e-2,100))
print(cor.test(tiling_table[,"ASMA_stroma_area"], tiling_table[,"stromal_Coll4_coverage"]))
print(g)


tiling_table <- "G:/My Drive/Paper website/Scripts/intermediate_files/S13.csv"

tiling_table <- read.csv(tiling_table, row.names = 1)

vars <- c("stromal_aSMA_coverage","stromal_Coll12A_coverage")

g <- ggplot(tiling_table, aes(x=stromal_aSMA_coverage, y=stromal_Coll12A_coverage)) +
  theme_bw() +
  geom_point(size=2) +
  ylab("stromal_Coll12A_coverage") +
  theme(panel.grid.minor = element_blank()) +
  scale_x_log10(limits = c(1e-1,100)) +
  scale_y_log10(limits = c(1e-2,100))
# scale_x_continuous(expand = c(0, 0), limits = c(-1, 4)) +
# scale_y_continuous(expand = c(0, 0), limits = c(0, 4000 ))
print(cor.test(tiling_table[,"stromal_aSMA_coverage"], tiling_table[,"stromal_Coll12A_coverage"]))
print(g)

###################
# infiltration graphing and analysis
tiling_table <- "intermediate_files/S10.csv"

# or

tiling_table <- "intermediate_files/S11.csv"
tiling_table <- "G:/My Drive/Paper website/Scripts/intermediate_files/S11.csv"

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








