`%notin%` <- Negate(`%in%`)

library(Matrix)
library(rhdf5)
library(plotly)
library(ggplot2)
library(data.table)
library(scDissector)
library(pheatmap)

main_path <- "D:/Dropbox/Grout et al/0. REVISIONS for CANCER DISCOVERY/Final editor edits/Paper website/Scripts"

setwd(main_path)
gene_lists <- paste(main_path,"/inputs/gene_lists.txt",sep="")
data_path <- paste(main_path,"/inputs",sep="")

# updated lists for external dataset analysis
fap_asma <- "DERL3,CYP26A1,SULF1,C1QTNF3,CXCL14,SUGCT,PLAU,HOPX,HTRA3,EPYC,MMP1,SCUBE3"
temp_meso <- "MT1G,WT1,CCL7,BDKRB1,CALB2,SERPINB2,HP,ITLN1,TM4SF1,CPA4,KRT8,KRT18"
temp_adh1b <- "CH25H,FRZB,TGM2,A2M,PPP1R14A,PTGDS,PLEKHH2,FMO2,PRELP,SLIT2,SCN7A,DST,SEPP1,LTBP4,ALDH1A1,ADH1B,SOD3,ABCA8,PLAC9,FBLN1,GSN,APOD,C7,C3,APOE"
temp_fap <- "COL11A1,CTSK,CHN1,COL6A3,VCAN,SPARC,COL12A1,THBS2,COL1A2,CTHRC1,POSTN,DIO2,COL1A1,RIN2,MMP11,COL6A1,COL5A2,CRABP2,LRRC17,ADAM12,TWIST1,LRRC15,CILP,GJB2,TNFAIP6,INHBA,GREM1,CST1,IGFL2,MMP7,P4HA3,PLPP4"
myh11_sig <- "IGFBP2,ITGBL1,COL4A1,COL4A2,LTBP2,LTBP1,ELN,MYH11,WIF1,PCSK1N,WFDC1,MYLK,F2R,CITED2,ROBO2,LUZP2,TSLP,FAM150A,TYRP1,CST2,ANO4,ENTPD1-AS1,KRT17,WNT11,ATRNL1,COL9A1,BRF2"
alv <- "G0S2,SRGN,RGCC,FIGF,GPM6B,MYH10,LIMCH1,INMT,FHL1,TIMP3,MAMDC2,MACF1,TCF21,AOC3"

temp_meso <- strsplit(temp_meso,",")
temp_adh1b <- strsplit(temp_adh1b,",")
temp_fap <- strsplit(temp_fap,",")
fap_asma <- strsplit(fap_asma,",")
temp_myh11 <- strsplit(myh11_sig,",")
alv <- strsplit(alv,",")

# graphing parameters
log_graphs_e <- -4:1
breaks <- 10^(log_graphs_e)
minor_breaks <- rep(1:9, length(breaks)-1)*(10^rep(log_graphs_e[1:length(breaks)-1], each=9))
reg <- 1e-4

#####################################

# run after running kim code

fib_cells <- rownames(temp_u)

downsample=function(u,min_umis){
  all_op=1:nrow(u)
  base_tab=rep(0,nrow(u))
  names(base_tab)=all_op
  
  downsamp_one=function(x,n){
    tab=base_tab
    tab2=table(sample(rep(all_op,x),size = n,replace=F))
    tab[names(tab2)]=tab2
    return(tab)
  }
  
  cell_mask=colnames(u)[colSums(u,na.rm=T)>min_umis]
  message("Downsampling ", length(cell_mask), " cells to ",min_umis," UMIs")
  
  chunk_size=5000
  breaks=unique(c(seq(from=1,to = length(cell_mask),by = chunk_size),length(cell_mask)+1))
  for (i in 1:(length(breaks)-1)){
    if (i==1){
      ds=Matrix(apply(u[,cell_mask[breaks[i]:(breaks[i+1]-1)],drop=F],2,downsamp_one,min_umis))
    }else{
      ds=cBind(ds,Matrix(apply(u[,cell_mask[breaks[i]:(breaks[i+1]-1)],drop=F],2,downsamp_one,min_umis)))
    }
    
  }
  rownames(ds)=rownames(u)
  colnames(ds)=cell_mask
  return(ds)
}

min_umis <- 2000

original_u <- u

# adjust to only include cells after filterin
fib_cells <- fib_cells[fib_cells%in%colnames(u)]
temp_u <- u[,fib_cells]
u <- temp_u

u <- downsample(u, min_umis)

##
# for Kim et al only for observing MYH11 cells
u <- u[,myh11_cells[myh11_cells%in%colnames(u)]]

u <- u[,other_cells[other_cells%in%colnames(u)]]
##

write_path <- "C:/Users/johna/Desktop/"

setwd("G:/My Drive/Paper website/Scripts")

# source files
gene_lists <- "inputs/gene_lists.txt"
stroma_ldm <- "inputs/ldm_2000_v2.RData"
data_dir <- "inputs"
sample_annots <- read.csv(file.path(data_dir,"sample_annots.csv"),r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv(file.path(data_dir,"model_lung_190606_annot_lists.csv"),r=1,h=1,stringsAsFactors = F)
annots_list$norm_group[annots_list$lineage=="MNP"] <- "MNP"
lcam_metadata <- "inputs/cell_metadata.rd"
lcam_id_conversion <- "inputs/conversion.txt"

# genes to remove from gene module analysis
t1=read.table(gene_lists,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
junk <- c(gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)

gene_list1_adj <- temp_adh1b[[1]]
gene_list2_adj <- temp_fap[[1]]

ingenes <- c(gene_list1_adj,gene_list2_adj)
ingenes <- ingenes[ingenes%in%rownames(u)]

myColor <- colorRampPalette(c("blue", "red"))(4)

# MYH11 focus
ingenes <- temp_myh11[[1]]
gene_list2_adj <- temp_myh11[[1]]
gene_list1_adj <- temp_myh11[[1]]

myh11_cells <- myh11_cells[myh11_cells%in%colnames(u)]
u <- u[,myh11_cells]

other_cells <- other_cells[other_cells%in%colnames(u)]
u <- u[,other_cells]
u <- u[,sample(other_cells, 50)]

colgrad_abs_file=system.file("extdata", "colors_paul.txt", package="scDissector")
colgrad_abs<<-read.table(colgrad_abs_file,stringsAsFactors=F)[,1]

breaks <- seq(0,8,length = length(colgrad_abs))

sort_genes <- gene_list2_adj[gene_list2_adj%in%rownames(u)]
# sort_genes <- sort_genes[sort_genes %notin% c("ELN")]

sort_temp <- colSums(log2(u[sort_genes,]+1))
sorted_cells <- names(sort(sort_temp,decreasing  = T))

pheatmap(u[ingenes,sorted_cells],cluster_rows = F,
         cluster_cols = F, show_colnames = F, breaks = breaks, color = colgrad_abs)

pheatmap(u[ingenes,sorted_cells],cluster_rows = F,
         cluster_cols = F, show_colnames = F, breaks = c(0,1,2,4,8), color = myColor)

#####################################
# kim analysis

t1=read.table(gene_lists,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
junk <- c(gl$HLA, gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)

gl$myh11 <- temp_myh11[[1]]
gl$ADH1B_CAF <- temp_adh1b[[1]]
gl$FAP_CAF <- temp_fap[[1]]
gl$fap_asma <- fap_asma[[1]]
gl$alv <- alv[[1]]
gl$meso <- temp_meso[[1]]

u <- fread(file=paste(data_path,"/Kim_cells_of_interest_v2.csv",sep=""))
gene_names <- unlist(u$Index)
cell_names <- colnames(u)[2:ncol(u)]

u <- u[u$Index%notin%junk,]

####

cell_annots <- read.csv(paste(data_path, "/kim_cell_annots.csv", sep=""), row.names = 1)
staging <- read.csv(paste(data_path, "/kim_staging.csv", sep=""), row.names = 1)

cell_types <- as.character(unique(cell_annots[,"Cell_type"]))
cell_subtypes <- as.character(unique(cell_annots[,"Cell_subtype"]))

# columns are cell types
# rowns are patients
patient_count_matrix <- matrix(NA,nrow = nrow(staging),ncol = length(cell_subtypes)+2)
rownames(patient_count_matrix) <- staging[,"Samples"]
colnames(patient_count_matrix) <- c(cell_subtypes, "ADH1B", "FAP")

samp <- as.character(staging[,"Samples"][1])
type <- colnames(patient_count_matrix)[1]
type <- "ADH1B"
for (samp in staging[,"Samples"]){
  for (type in colnames(patient_count_matrix)){
    if (is.na(type)){
      next
    }
    # matchning sample
    type_count <- cell_annots[cell_annots[,"Sample"] == samp,]
    # matching cell type
    type_count <- type_count[type_count[,"Cell_subtype"] == type,]
    
    # removing introduced NA values - result from trying to match FAP
    type_count <- type_count[!is.na(type_count[,1]),]

    patient_count_matrix[samp,type] <- nrow(type_count)
  }
}
patient_count_matrix

####

for (num in 1:length(gl)){
  gl[[num]] <- gl[[num]][gl[[num]]%in%u$Index]
}

gl$extended_fib <- gl$extended_fib[gl$extended_fib%notin%c(gl$FAP_CAF,gl$ADH1B_CAF)]

temp_fib_mat <- u[u$Index%in%gl$extended_fib,]
temp_fib_mat$Index <- NULL
temp_fib_mat <- as.matrix(temp_fib_mat)
dataFib <- c()

for (cell in cell_names){
  temp_sum <- sum(u[,..cell])
  if (temp_sum <= 1000){
    temp_calc <- 0
    dataFib <- c(dataFib, temp_calc)
  }else{
    temp_calc <- sum(temp_fib_mat[,cell])/temp_sum
    dataFib <- c(dataFib, temp_calc)
  }
}

temp_pvc_mat <- u[u$Index%in%gl$pan_PvC,]
temp_pvc_mat$Index <- NULL
temp_pvc_mat <- as.matrix(temp_pvc_mat)
dataPvC <- c()

for (cell in cell_names){
  temp_sum <- sum(u[,..cell])
  if (temp_sum <= 1000){
    temp_calc <- 0
    dataPvC <- c(dataPvC, temp_calc)
  }else{
    temp_calc <- sum(temp_pvc_mat[,cell])/temp_sum
    dataPvC <- c(dataPvC, temp_calc)
  }
}

# graphing parameters
log_graphs_e <- -4:1
breaks <- 10^(log_graphs_e)
minor_breaks <- rep(1:9, length(breaks)-1)*(10^rep(log_graphs_e[1:length(breaks)-1], each=9))

temp <- data.frame(dataFib,dataPvC)

to_regularize <- list("dataFib","dataPvC")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

ggplot(temp, aes(x=dataFib, y=dataPvC)) +
  geom_point(size=1,alpha=0.2) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,5e-2),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,3e-1),breaks = breaks, minor_breaks = minor_breaks) +
  geom_vline(xintercept = 2e-2) +
  geom_hline(yintercept = 1.5e-3) +
  geom_vline(xintercept = 3e-1, color = "black") +
  geom_hline(yintercept = 5e-2, color = "black") +
  theme(panel.grid.minor = element_blank())

silico_gated <- c(dataFib > 2e-2 & dataPvC < 1.5e-3)

print(length(silico_gated[silico_gated])/(nrow(temp)-nrow(temp[rowSums(temp) == 0,])))

u <- fread(file=paste(data_path,"/Kim_cells_of_interest_v2.csv",sep=""), select = cell_names[silico_gated])
u <- as.matrix(u)
rownames(u) <- gene_names
u <- u[rownames(u)[rownames(u)%notin%junk],]

dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)
dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataADH1B <- colSums(u[gl$ADH1B_CAF,])/colSums(u)
dataFAP <- colSums(u[gl$FAP_CAF,])/colSums(u)
dataMYH11 <- colSums(u[gl$myh11,])/colSums(u)
dataMESO <- colSums(u[gl$meso,])/colSums(u)
dataASMA <- colSums(u[gl$fap_asma,])/colSums(u)
dataALV <- colSums(u[gl$alv,])/colSums(u)

# splitting by stage
rm(temp_pvc_mat)
rm(temp_fib_mat)

staging <- read.csv(paste(data_path, "/kim_staging.csv", sep=""), row.names = 1)
cell_stages <- c()
cell_samples <- c()

head(staging[,"Stages"])
early_stages <- c("IA","IIA","IA3","IIIA")

temp_u <- u

cell <- colnames(u)[2]
for (cell in colnames(u)){
  prefix <- strsplit(cell,"_")[[1]][3]
  cell_stage <- staging[prefix,"Stages"]
  
  if (cell_stage%in%early_stages){
    cell_stages <- c(cell_stages, "early")
  }else{
    cell_stages <- c(cell_stages, "late")
  }
  cell_samples <- c(cell_samples, prefix)
}
new_rownames <- c(rownames(u),"stage","sample")
u <- rbind(u, NA, NA)
rownames(u) <- new_rownames
u["stage",] <- cell_stages
u["sample",] <- cell_samples
head(u["stage",])

group.colors <- c(early = "#203FE0", late = "#B00E0E")

temp <- data.frame(dataFib,dataADH1B,dataFAP,dataASMA,dataMYH11,dataALV,dataMESO,u["stage",])
colnames(temp) <- c("dataFib","dataADH1B","dataFAP","dataASMA","dataMYH11","dataALV","dataMESO","stage")

to_regularize <- list("dataFib","dataADH1B","dataFAP","dataASMA","dataMYH11","dataALV","dataMESO")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

ggplot(temp, aes(x=dataADH1B, y=dataALV)) +
  geom_point(size=1, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1.5e-2) +
  ggtitle("Kim") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataALV[dataALV < 1.5e-2]),]

ggplot(temp, aes(x=dataADH1B, y=dataMESO)) +
  geom_point(size=1, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 2e-3) +
  ggtitle("Kim") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMESO[dataMESO < 2e-3]),]

ggplot(temp, aes(x=dataADH1B, y=dataMYH11)) +
  geom_point(size=1, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1e-2) +
  ggtitle("Kim") +
  theme(panel.grid.minor = element_blank())

##
# for making myh11 truthplot
myh11_cells <- names(dataMYH11[dataMYH11 > 1e-2])
dataMYH11 <- dataMYH11[dataMYH11 > 1e-2]

other_cells <- names(dataMYH11[dataMYH11 < 1e-2])
dataMYH11 <- dataMYH11[dataMYH11 < 1e-2]

myh11_u <- temp
##

temp <- temp[names(dataMYH11[dataMYH11 < 1e-2]),]

ggplot(temp, aes(x=dataADH1B, y=dataFAP, color=stage)) +
  geom_point(size=1, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 5e-1, color = "black") +
  ggtitle("Kim") +
  theme(panel.grid.minor = element_blank())

ggplot(temp, aes(x=dataASMA, y=dataFAP)) +
  geom_point(size=1, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  ggtitle("Kim") +
  theme(panel.grid.minor = element_blank())

#####################################
# laugh analysis

gene_lists <- "G:/My Drive/Paper website/Scripts/inputs/gene_lists.txt"
t1=read.table(gene_lists,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
junk <- c(gl$HLA, gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)

gl$myh11 <- temp_myh11[[1]]
gl$ADH1B_CAF <- temp_adh1b[[1]]
gl$FAP_CAF <- temp_fap[[1]]
gl$fap_asma <- fap_asma[[1]]
gl$alv <- alv[[1]]
gl$meso <- temp_meso[[1]]

load(paste(data_path,"/laughney_data_v2.RData", sep=""))

u <- u[rownames(u)[rownames(u)%notin%junk],]

for (num in 1:length(gl)){
  gl[[num]] <- gl[[num]][gl[[num]]%in%rownames(u)]
}

gl$extended_fib <- gl$extended_fib[gl$extended_fib%notin%c(gl$FAP_CAF,gl$ADH1B_CAF)]

dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)

temp <- data.frame(dataFib,dataPvC)

to_regularize <- list("dataFib","dataPvC")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

ggplot(temp, aes(x=dataFib, y=dataPvC)) +
  geom_point(size=2,alpha=0.2) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,5e-2),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,3e-1),breaks = breaks, minor_breaks = minor_breaks) +
  geom_vline(xintercept = 1.5e-2) +
  geom_hline(yintercept = 1.5e-3) +
  geom_hline(yintercept = 5e-2, color = "black") +
  geom_vline(xintercept = 3e-1, color = "black") +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("Laughney")

silico_gated <- c(dataFib > 1.5e-2 & dataPvC < 1.5e-3)

print(length(silico_gated[silico_gated])/nrow(temp))

u <- u[,silico_gated]

##
# splitting by stage
staging <- read.csv(paste(data_path, "/laughney_staging.csv", sep=""), row.names = 1)
cell_stages <- c()
for (cell in colnames(u)){
  prefix <- strsplit(cell,"_")[[1]][1]
  staging[prefix,"t_num"]
  cell_stages <- c(cell_stages,staging[prefix,"t_num"])
}
new_rownames <- c(rownames(u),"stages")
u <- rbind(u, NA)
rownames(u) <- new_rownames
u["stages",] <- cell_stages
##

dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)
dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataADH1B <- colSums(u[gl$ADH1B_CAF,])/colSums(u)
dataFAP <- colSums(u[gl$FAP_CAF,])/colSums(u)
dataMYH11 <- colSums(u[gl$myh11,])/colSums(u)
dataMESO <- colSums(u[gl$meso,])/colSums(u)
dataASMA <- colSums(u[gl$fap_asma,])/colSums(u)
dataALV <- colSums(u[gl$alv,])/colSums(u)

temp <- data.frame(dataFib,dataADH1B,dataFAP,dataASMA,dataMYH11,dataALV,dataMESO)

to_regularize <- list("dataFib","dataADH1B","dataFAP","dataASMA","dataMYH11","dataALV","dataMESO")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

temp[colnames(u)[u["stages",]==1],"stage"] <- "early"
temp[colnames(u)[u["stages",]==2],"stage"] <- "late"

group.colors <- c(early = "#203FE0", late = "#B00E0E")

temp <- temp[order(temp[,"stage"],decreasing=T),]

ggplot(temp, aes(x=dataADH1B, y=dataALV)) +
  geom_point(size= 2, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1.5e-2) +
  ggtitle("Laughney") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataALV[dataALV < 1.5e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMESO)) +
  geom_point(size= 2, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 2e-3) +
  ggtitle("Laughney") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMESO[dataMESO < 2e-3]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMYH11)) +
  geom_point(size= 2, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1e-2) +
  ggtitle("Laughney") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMYH11[dataMYH11 < 1e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataFAP, color=stage)) +
  geom_point(size=2, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 5e-1, color = "black") +
  ggtitle("Laughney") +
  theme(panel.grid.minor = element_blank())

ggplot(temp, aes(x=dataASMA, y=dataFAP)) +
  geom_point(size=2, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  ggtitle("Laughney") +
  theme(panel.grid.minor = element_blank())

#####################################
# wu analysis

gene_lists <- "G:/My Drive/Paper website/Scripts/inputs/gene_lists.txt"
t1=read.table(gene_lists,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
junk <- c(gl$HLA, gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)

gl$myh11 <- temp_myh11[[1]]
gl$ADH1B_CAF <- temp_adh1b[[1]]
gl$FAP_CAF <- temp_fap[[1]]
gl$fap_asma <- fap_asma[[1]]
gl$alv <- alv[[1]]
gl$meso <- temp_meso[[1]]

# wu saved as u
load(paste(data_path,"/wu_data_v3.RData", sep=""))

u <- u[rownames(u)[rownames(u)%notin%junk],]

for (num in 1:length(gl)){
  gl[[num]] <- gl[[num]][gl[[num]]%in%rownames(u)]
}

gl$extended_fib <- gl$extended_fib[gl$extended_fib%notin%c(gl$FAP_CAF,gl$ADH1B_CAF)]

dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)

temp <- data.frame(dataFib,dataPvC)

to_regularize <- list("dataFib","dataPvC")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

ggplot(temp, aes(x=dataFib, y=dataPvC)) +
  geom_point(size=1,alpha=0.2) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,5e-2),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,3e-1),breaks = breaks, minor_breaks = minor_breaks) +
  geom_vline(xintercept = 2e-2) +
  geom_hline(yintercept = 1.5e-3) +
  geom_hline(yintercept = 5e-2, color = "black") +
  geom_vline(xintercept = 3e-1, color = "black") +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("Wu")

silico_gated <- c(dataFib > 2e-2 & dataPvC < 1.5e-3)

print(length(silico_gated[silico_gated])/nrow(temp))

u <- u[,silico_gated]

dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)
dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataADH1B <- colSums(u[gl$ADH1B_CAF,])/colSums(u)
dataFAP <- colSums(u[gl$FAP_CAF,])/colSums(u)
dataMYH11 <- colSums(u[gl$myh11,])/colSums(u)
dataMESO <- colSums(u[gl$meso,])/colSums(u)
dataASMA <- colSums(u[gl$fap_asma,])/colSums(u)
dataALV <- colSums(u[gl$alv,])/colSums(u)

temp <- data.frame(dataFib,dataADH1B,dataFAP,dataASMA,dataMYH11,dataALV,dataMESO)

to_regularize <- list("dataFib","dataADH1B","dataFAP","dataASMA","dataMYH11","dataALV","dataMESO")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

ggplot(temp, aes(x=dataADH1B, y=dataALV)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1.5e-2) +
  ggtitle("Wu") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataALV[dataALV < 1.5e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMESO)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 2e-3) +
  ggtitle("Wu") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMESO[dataMESO < 2e-3]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMYH11)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1e-2) +
  ggtitle("Wu") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMYH11[dataMYH11 < 1e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataFAP)) +
  geom_point(size=1, alpha = 0.5, color="#B00E0E") +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 5e-1, color = "black") +
  ggtitle("Wu") +
  theme(panel.grid.minor = element_blank())

ggplot(temp, aes(x=dataASMA, y=dataFAP)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  ggtitle("Wu") +
  theme(panel.grid.minor = element_blank())

#####################################
# lambs analysis

gene_lists <- "G:/My Drive/Paper website/Scripts/inputs/gene_lists.txt"
t1=read.table(gene_lists,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
junk <- c(gl$HLA, gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)

gl$myh11 <- temp_myh11[[1]]
gl$ADH1B_CAF <- temp_adh1b[[1]]
gl$FAP_CAF <- temp_fap[[1]]
gl$fap_asma <- fap_asma[[1]]
gl$alv <- alv[[1]]
gl$meso <- temp_meso[[1]]

# lambs saved as u_lambs
load(paste(data_path,"/lambs_tum_only.RData", sep=""))

u <- u_lambs
rm(u_lambs)

u <- u[rownames(u)[rownames(u)%notin%junk],]

# removing patients 5 and 8, large cell carcinoma and pleiomorphic
cell <- colnames(u)[1]
cells_to_remove <- c()
for (cell in colnames(u)){
  prefix <- as.character(strsplit(cell,"_")[[1]][1])
  if (prefix%in%c("79","80","81","92","93","94","NA")){
    cells_to_remove <- c(cells_to_remove, cell)
  }
}
u <- u[,colnames(u)[!colnames(u)%in%cells_to_remove]]

for (num in 1:length(gl)){
  gl[[num]] <- gl[[num]][gl[[num]]%in%rownames(u)]
}

gl$extended_fib <- gl$extended_fib[gl$extended_fib%notin%c(gl$FAP_CAF,gl$ADH1B_CAF)]

dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)

temp <- data.frame(dataFib,dataPvC)

to_regularize <- list("dataFib","dataPvC")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

ggplot(temp, aes(x=dataFib, y=dataPvC)) +
  geom_point(size=1,alpha=0.2) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,5e-2),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,3e-1),breaks = breaks, minor_breaks = minor_breaks) +
  geom_vline(xintercept = 2.5e-2) +
  geom_hline(yintercept = 1.5e-3) +
  geom_hline(yintercept = 5e-2, color = "black") +
  geom_vline(xintercept = 3e-1, color = "black") +
  theme(panel.grid.minor = element_blank())

silico_gated <- c(dataFib > 2.5e-2 & dataPvC < 1.5e-3)

print(length(silico_gated[silico_gated])/nrow(temp))

u <- u[,silico_gated]

# reading cell stage
staging <- read.csv(paste(data_path, "/lambs_staging.csv", sep=""), row.names = 1)
cell_stages <- c()
cell_tumor_type <- c()
cell <- colnames(u)[1]

for (cell in colnames(u)){
  prefix <- as.character(strsplit(cell,"_")[[1]][1])
  cell_stages <- c(cell_stages,staging[prefix,"t_num"])
  cell_tumor_type <- c(cell_tumor_type,as.character(staging[prefix,"Type"]))
}

names(cell_stages) <- colnames(u)
names(cell_tumor_type) <- colnames(u)

dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)
dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataADH1B <- colSums(u[gl$ADH1B_CAF,])/colSums(u)
dataFAP <- colSums(u[gl$FAP_CAF,])/colSums(u)
dataMYH11 <- colSums(u[gl$myh11,])/colSums(u)
dataMESO <- colSums(u[gl$meso,])/colSums(u)
dataASMA <- colSums(u[gl$fap_asma,])/colSums(u)
dataALV <- colSums(u[gl$alv,])/colSums(u)

library(RColorBrewer)
cols <- brewer.pal(4,"Spectral")

tum_to_col_num <- c(AD = 1, Mix = 2, LCC = 3, SCC = 4)

# color ID based on tumor type
color_id <- rev(cols)[tum_to_col_num[cell_tumor_type]]

# color id based on stage
cell_stages <- as.numeric(gsub(3,4,cell_stages))
cell_stages <- as.numeric(gsub(2,4,cell_stages))
names(cell_stages) <- colnames(u)
# red is high stage, blue is low stage
color_id <- rev(cols)[cell_stages]
names(color_id) <- colnames(u)

temp <- data.frame(dataFib,dataADH1B,dataFAP,dataASMA,dataALV,dataMYH11,dataMESO,"stage"=NA)

to_regularize <- list("dataFib","dataADH1B","dataFAP","dataASMA","dataMYH11","dataALV","dataMESO")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

# temp[,"dataFAP"] <- pmax(temp[,"dataFAP"],1e-3)

temp[names(cell_stages)[cell_stages==1],"stage"] <- "early"
temp[names(cell_stages)[cell_stages==4],"stage"] <- "late"

group.colors <- c(early = "#203FE0", late = "#B00E0E")

ggplot(temp, aes(x=dataADH1B, y=dataALV)) +
  geom_point(size=1.8, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1e-2) +
  ggtitle("Lambrechts") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataALV[dataALV < 1e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMESO)) +
  geom_point(size=1.8, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 5e-3) +
  ggtitle("Lambrechts") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMESO[dataMESO < 5e-3]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMYH11)) +
  geom_point(size=1.8, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1e-2) +
  ggtitle("Lambrechts") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMYH11[dataMYH11 < 1e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataFAP)) +
  geom_point(size=1.8, color="#B00E0E") +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 5e-1, color = "black") +
  ggtitle("Lambrechts") +
  theme(panel.grid.minor = element_blank())

ggplot(temp, aes(x=dataASMA, y=dataFAP)) +
  geom_point(size=1.8, alpha = 1) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  ggtitle("Lambrechts") +
  theme(panel.grid.minor = element_blank())

#####################################
# in house analysis

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
junk <- c(gl$HLA, gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)

gl$myh11 <- temp_myh11[[1]]
gl$ADH1B_CAF <- temp_adh1b[[1]]
gl$FAP_CAF <- temp_fap[[1]]
gl$fap_asma <- fap_asma[[1]]
gl$alv <- alv[[1]]
gl$meso <- temp_meso[[1]]

load(paste(data_path,"/in_house_v3.RData", sep=""))
load(paste(data_path,"/ldm_1000.RData", sep=""))
cluster_id <- ldm$dataset$cell_to_cluster

fap_cells <- names(cluster_id[cluster_id%in%c("22","24","26")])
adh1b_cells <- names(cluster_id[cluster_id%in%c("20")])
alv_cells <- names(cluster_id[cluster_id%in%c("8","4")])
msc_cells <- names(cluster_id[cluster_id%in%c("9")])
myh11_cells <- names(cluster_id[cluster_id%in%c("11")])
clu_cells <- names(cluster_id[cluster_id%in%c("19")])
meso_cells <- names(cluster_id[cluster_id%in%c("1")])
pvc_cells <- names(cluster_id[cluster_id%in%c("7","25","5","12","17")])

u <- u[rownames(u)[rownames(u)%notin%junk],]

for (num in 1:length(gl)){
  gl[[num]] <- gl[[num]][gl[[num]]%in%rownames(u)]
}

gl$extended_fib <- gl$extended_fib[gl$extended_fib%notin%c(gl$FAP_CAF,gl$ADH1B_CAF)]

dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)
dataMYH11 <- colSums(u[gl$myh11,])/colSums(u)
dataMESO <- colSums(u[gl$meso,])/colSums(u)

temp <- data.frame(dataFib,dataPvC)

cluster_id <- cluster_id[names(cluster_id)[names(cluster_id)%in%rownames(temp)]]
temp <- temp[names(cluster_id),]

to_regularize <- list("dataFib","dataPvC")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

ggplot(temp, aes(x=dataFib, y=dataPvC)) +
  geom_point(size=1,alpha=0.2) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,5e-2),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,3e-1),breaks = breaks, minor_breaks = minor_breaks) +
  geom_vline(xintercept = 2e-2) +
  geom_hline(yintercept = 1.5e-3) +
  geom_hline(yintercept = 5e-2, color = "black") +
  geom_vline(xintercept = 3e-1, color = "black") +
  theme(panel.grid.minor = element_blank())

silico_gated <- c(dataFib > 2e-2 & dataPvC < 1.5e-3)
silico_gated[is.na(silico_gated)] <- FALSE

print(length(silico_gated[silico_gated])/nrow(temp))

u <- u[,silico_gated]

cluster_id <- cluster_id[names(cluster_id)[names(cluster_id)%in%colnames(u)]]
u <- u[,names(cluster_id)]

##
# splitting by stage
s1=read.table("G:/My Drive/scRNAseq_analysis/compiled/clustering_data_lung4/sample_sets.txt",stringsAsFactors = F,row.names = 1)
il=strsplit(s1[,1],",")
names(il)=rownames(s1)

samp_names <- il$JG_stroma_consented_tumor

stage_2_plus_samps <- as.character(c(246,248,260,400,412))
stage_1_samps <- samp_names[samp_names%notin%stage_2_plus_samps]

stage_1_cells <- c()
stage_2_plus_cells <- c()
for (cell in colnames(u)){
  prefix <- strsplit(cell,"_")[[1]][1]
  if (prefix%in%stage_1_samps){
    stage_1_cells <- c(stage_1_cells, cell)
  }else{
    stage_2_plus_cells <- c(stage_2_plus_cells, cell)
  }
}
##

dataPvC <- colSums(u[gl$pan_PvC,])/colSums(u)
dataFib <- colSums(u[gl$extended_fib,])/colSums(u)
dataADH1B <- colSums(u[gl$ADH1B_CAF,])/colSums(u)
dataFAP <- colSums(u[gl$FAP_CAF,])/colSums(u)
dataMYH11 <- colSums(u[gl$myh11,])/colSums(u)
dataMESO <- colSums(u[gl$meso,])/colSums(u)
dataASMA <- colSums(u[gl$fap_asma,])/colSums(u)
dataALV <- colSums(u[gl$alv,])/colSums(u)

temp <- data.frame(dataADH1B,dataFAP,dataFib,dataASMA)
temp_filter <- temp[!(rowSums(temp)<=reg*2),]

temp <- data.frame(dataFib,dataADH1B,dataMESO,dataFAP,dataASMA,dataALV,dataMYH11,"stage"=NA)

to_regularize <- list("dataFib","dataADH1B","dataMESO","dataFAP","dataASMA","dataALV","dataMYH11")
for (set in to_regularize){
  temp[,set] <- pmax(temp[,set],reg)
}

# temp[,"dataFAP"] <- pmax(temp[,"dataFAP"],1e-3)

alv_cells_cluster <- alv_cells[alv_cells%in%names(dataASMA)]
msc_cells_cluster <- msc_cells[msc_cells%in%names(dataASMA)]
fap_cells_cluster <- fap_cells[fap_cells%in%names(dataASMA)]
adh1b_cells_cluster <- adh1b_cells[adh1b_cells%in%names(dataASMA)]
myh11_cells_cluster <- myh11_cells[myh11_cells%in%names(dataASMA)]
clu_cells_cluster <- clu_cells[clu_cells%in%names(dataASMA)]
meso_cells_cluster <- meso_cells[meso_cells%in%names(dataASMA)]
pvc_cells_cluster <- pvc_cells[pvc_cells%in%names(dataASMA)]

temp[fap_cells_cluster,"stage"] <- "red"
temp[adh1b_cells_cluster,"stage"] <- "beige"
temp[alv_cells_cluster,"stage"] <- "blue"
temp[msc_cells_cluster,"stage"] <- "pink"
temp[myh11_cells_cluster,"stage"] <- "gray"
temp[clu_cells_cluster,"stage"] <- "gray"
temp[meso_cells_cluster,"stage"] <- "brown"
temp[pvc_cells_cluster,"stage"] <- "green"
temp[is.na(temp[,"stage"]),"stage"] <- "black"

temp <- temp[order(temp[,"stage"],decreasing=TRUE),]

group.colors <- c(blue = "#203FE0", red = "#B00E0E", green = "#00FF00", beige = "#A19A1D",
                  pink = "#B989C9", gray = "#808080", brown = "#694000", black = "#000000")

ggplot(temp, aes(x=dataADH1B, y=dataALV)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1.5e-2) +
  ggtitle("In house") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataALV[dataALV < 1.5e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMESO)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 2e-3) +
  ggtitle("In house") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMESO[dataMESO < 2e-3]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataMYH11)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 1e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  geom_hline(yintercept = 1e-2) +
  ggtitle("In house") +
  theme(panel.grid.minor = element_blank())

temp <- temp[names(dataMYH11[dataMYH11 < 1e-2]),]
temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataFAP, color=stage)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 5e-1, color = "black") +
  theme(panel.grid.minor = element_blank())

temp[stage_1_cells,"stage"] <- "early"
temp[stage_2_plus_cells,"stage"] <- "late"

temp <- temp[order(temp[,"stage"],decreasing=TRUE),]

group.colors <- c(early = "#203FE0", late = "#B00E0E")

temp <- temp[!is.na(temp[,1]),]

ggplot(temp, aes(x=dataADH1B, y=dataFAP, color=stage)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 5e-1, color = "black") +
  theme(panel.grid.minor = element_blank())

ggplot(temp, aes(x=dataASMA, y=dataFAP)) +
  geom_point(size=1, alpha = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(1e-3,5e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(limits = c(1e-4,1e-1),breaks = breaks, minor_breaks = minor_breaks) +
  scale_color_manual(values=group.colors) +
  theme(legend.position="none") +
  geom_hline(yintercept = 5e-1, color = "black") +
  geom_vline(xintercept = 1e-1, color = "black") +
  ggtitle("In house") +
  theme(panel.grid.minor = element_blank())




