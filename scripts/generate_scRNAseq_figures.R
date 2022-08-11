# if (download_data){
#   download_files(pipeline_path)
# }
# 
# ldm_filename=paste(pipeline_path,"inputs","ldm_2000_v2.RData",sep="/")
# lcam_metadata_filename=paste(pipeline_path,"inputs","cell_metadata.rd",sep="/")
# 
# load(ldm_filename,envir = env)
# load(lcam_metadata_filename,envir = env)

setwd("G:/My Drive/Paper website/Scripts")

# libraries
library(TCGA2STAT)
library(reshape2)
library(ggplot2)
library(gplots)
library(skmeans)
library(Matrix.utils)
library(mixtools)
library(httpuv)
library(shiny)
library(robustbase)
library(devtools)
library(scDissector)
library(matrixStats)
library(Matrix)
library(scales)
library(RColorBrewer)
library(pheatmap)
library(seriation)

write_path <- "outputs/"

# source files
gene_lists <- "inputs/gene_lists_temp.txt"
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

########################
# functions

`%notin%` <- Negate(`%in%`)
fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
fisher.z2r <- function(z) {(exp(2 * z) - 1)/(1 + exp(2 * z))}

open_plot=function(path="figures/",fn,plot_type="png",width=NA,height=NA,bg="white"){
  if (plot_type=="png"){
    if (is.na(width)){
      width=400
    }
    if (is.na(height)){
      height=400
    }
    png(paste(path,"/",fn,".png",sep=""),width,height,bg=bg)
  }else if(plot_type=="pdf"){
    if (is.na(width)){
      width=4
    }
    if (is.na(height)){
      height=4
    }
    pdf(paste(path,"/",fn,".pdf",sep=""),width,height,bg=bg)
  }
  par(cex.axis=1.5)
}

colgrad_abs_file=system.file("extdata", "colors_paul.txt", package="scDissector")

if (colgrad_abs_file==""){
  colgrad_abs_file="extdata/colors_viridis.txt"
}

colgrad_abs<<-read.table(colgrad_abs_file,stringsAsFactors=F)[,1]

close_plot=function(){
  dev.off()
}

Relative_or_Absolute="Relative"
colgrad<<-c(colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100))
sample_cols<<-rep(paste("#",read.table("inputs/sample_colors.txt",
                                       stringsAsFactors = F)[,1],sep=""),10)

make_truth_plots=function(ldm,clusters_to_exclude,gene_list1,gene_list2,samples=0,
                          zlim=c(0,2.5),main_title="stroma",clust_num=200,samp_num=100){
  
  clusters=setdiff(ldm$cluster_order,clusters_to_exclude)
  ds=ldm$dataset$ds[[match("2000",ldm$dataset$ds_numis)]]
  ds=ds[,ldm$dataset$cell_to_cluster[colnames(ds)]%in%clusters]
  l=list()
  
  cluster_id <- ldm$dataset$cell_to_cluster
  sample_id <- ldm$dataset$cell_to_sample
  
  if (samples == 0){
    samples <- unique(sample_id)
  }

  # sample and cluster normalization
  for (cluster in clusters){
    
    temp_mask <- c()
    for (samp in samples){
      
      cluster_mask=colnames(ds)[cluster_id[colnames(ds)]==cluster]
      sample_mask=colnames(ds)[sample_id[colnames(ds)]==samp]
      
      mask <- cluster_mask[cluster_mask%in%sample_mask]
      
      if (length(mask) == 0){
        next
      }
      if (samp_num < length(mask)){
        mask <- sample(mask,size = samp_num,replace = F)
      }
      temp_mask <- c(temp_mask, mask)
    }
    
    if (length(temp_mask) > clust_num){
      l[[paste(cluster)]]=sample(temp_mask,clust_num)
    }else{
      l[[paste(cluster)]]=sample(temp_mask,length(temp_mask))
    }
  }
  
  cells=unlist(l)
  ds=ds[,cells]
  
  gene_list1_adj=gene_list1
  gene_list2_adj=gene_list2
  
  open_plot(path=write_path,fn = "truth_color_legend",plot_type = "png",
            width = 500,height = 100)
  
  par(mar=c(1,1,1,1))
  image(matrix(1:100,100,1),col=colgrad_abs,pch=20,cex=3,axes=F)
  close_plot()
  
  open_plot(path=write_path,fn = main_title,plot_type = "png",width = 4000,height = 2000)
  
  # sorting cells based on input genes
  if (sort_cells){
    if (length(sort_genes) == 1){
      sorted_cells <- names(sort(log2(ldm$dataset$umitab[sort_genes,]+1),decreasing  = T))
    }else{
      sorted_cells <- names(sort(colSums(log2(ldm$dataset$umitab[sort_genes,]+1)),decreasing  = T))
    }
    ldm$dataset$umitab <- ldm$dataset$umitab[,sorted_cells]
    ldm$dataset$cell_to_cluster <- ldm$dataset$cell_to_cluster[sorted_cells]
    ldm$dataset$cell_to_sample <- ldm$dataset$cell_to_sample[sorted_cells]
    ds <- ds[,sorted_cells[sorted_cells%in%colnames(ds)]]
  }
  
  plot_truth_heatmap(ds,cell_to_sample =ldm$dataset$cell_to_sample[colnames(ds)],
                     cell_to_cluster = ldm$dataset$cell_to_cluster[colnames(ds)],
                     ingenes=gene_list1_adj,inclusts=clusters,
                     zlim=zlim,cols=colgrad_abs,plot_batch_bar=F,gene_text_cex=2,cluster_text_cex=2)
  close_plot()
}

sort_cells <- FALSE

########################
# loading ldm

load(stroma_ldm)

ldm$dataset$umitab <- ldm$dataset$umitab[!(rownames(ldm$dataset$umitab)%in%junk),]

cluster_id <- ldm$dataset$cell_to_cluster
sample_id <- ldm$dataset$cell_to_sample

samples <- unique(sample_id)



# convert_clust <- read.table("G:/My Drive/sc_scripts/Gene sets/clust_to_nam v2.txt")
# rownames(convert_clust) <- as.character(convert_clust[,1])
# 
# new_mat <- matrix(cluster_id)
# rownames(new_mat) <- names(cluster_id)
# new_mat <- cbind(new_mat, NA)
# new_mat[,2] <- as.character(convert_clust[new_mat[,1],2])
# write.csv(new_mat, "C:/Users/johna/Desktop/cluster_annotations.csv")

########################
# module analysis fibroblast

get_avg_gene_to_gene_cor=function(ds,cell_to_sample,samples=NULL,weighted=F,
                                  min_number_of_cell_per_sample=5,min_umi_counts_per_samples=5,
                                  showShinyProgressBar=F,session=NULL){
  
  zmat=matrix(0,nrow(ds),nrow(ds),dimnames=list(rownames(ds),rownames(ds)))
  samples_inp=unique(cell_to_sample)
  if (is.null(samples)){
    samples=samples_inp
  }else{
    samples=intersect(samples,samples_inp)
  }
  samples=samples[table(cell_to_sample)[samples]>=min_number_of_cell_per_sample]
  if (weighted){
    w=table(cell_to_sample)[samples]
    w=w/sum(w)
  }else{
    w=rep(1/length(samples),length(samples))
    names(w)=samples
  }
  i=0
  get_cor_per_sample=function(){
    for (samp in samples){
      i=i+1
      if (showShinyProgressBar){
        setProgress(i,session = session)
      }else{
        print(samp)
      }
      dsi=ds[Matrix::rowSums(ds[,cell_to_sample==samp,drop=F],na.rm=T)>min_umi_counts_per_samples,cell_to_sample==samp,drop=F]
      z=fisher.r2z(.99*cor(as.matrix(Matrix::t(dsi)),use="complete.obs"))
      rm(dsi)
      z[is.na(z)]=0
      zmat[rownames(z),colnames(z)]= zmat[rownames(z),colnames(z)] + w[samp]*z
      rm(z)
      gc()
    }
    return(zmat)
  }
  if (showShinyProgressBar){
    withProgress( {zmat=get_cor_per_sample()},1,length(samples),1,message = "Computing Correlations",session = session)
  }else{
    zmat=get_cor_per_sample()
  }
  return(fisher.z2r(zmat/sum(w)))
}

get_avg_module_to_gene_cor=function(ds,genes,modules_list,cell_to_sample,samples=NULL,
                                    weighted=F,min_number_of_cell_per_sample=5,min_umi_counts_per_samples=5){
  if (is.null(names(modules_list))){
    names(modules_list)=1:length(modules_list)
  }
  zmat=matrix(0,length(modules_list),length(genes),dimnames=list(names(modules_list),genes))
  samples_inp=unique(cell_to_sample)
  if (is.null(samples)){
    samples=samples_inp
  }else{
    samples=intersect(samples,samples_inp)
  }
  if (weighted){
    w=table(cell_to_sample)[samples]
    w=w/sum(w)
  }else{
    rep(1/length(samples),length(samples))
  }
  samples=samples[table(cell_to_sample)[samples]>=min_number_of_cell_per_sample]
  if (length(sample)==0){
    return()
  }
  for (samp in samples){
    print(samp)
    genes=genes[Matrix::rowSums(ds[genes,cell_to_sample==samp,drop=F],na.rm=T)>=min_umi_counts_per_samples]
    dsi=log2(1+ds[genes,cell_to_sample==samp,drop=F])
    get_one_mod_sum=function(x,ds){colSums(log2(1+as.matrix(ds[x,,drop=F])))}
    ds_mods_i=t(sapply(modules_list,get_one_mod_sum,ds[,cell_to_sample==samp,drop=F]))
    z=fisher.r2z(.99*cor(as.matrix(Matrix::t(ds_mods_i)),as.matrix(Matrix::t(dsi))))
    z[is.na(z)]=0
    zmat[rownames(z),colnames(z)]=z+ w[samp]*zmat[rownames(z),colnames(z),drop=F]
    rm(z)
    rm(dsi)
    gc()
  }
  zmat[,setdiff(colnames(zmat),genes)]=NA
  return(fisher.z2r(zmat/sum(w)))
}

get_gene_cormap=function(ldm,ds_version="2000",cells=NA){
  ds=ldm$dataset$ds[[match(ds_version,ldm$dataset$ds_numis)]]
  if (!all(is.na(cells))){
    ds=ds[,intersect(colnames(ds),cells)]
  }
  cormat=cor_analysis(ds,ldm)
  return(cormat)
}

save_gene_cor_map=function(cormat,modules_version,zbreaks=c(-1,seq(-.5,.5,l=99),1),
                           ser_method="GW",cor_cols=colorRampPalette(c("blue","white","red"))(100)){
  par(mar=c(3,3,3,3))
  ord=get_order(seriate(as.dist(1-cormat),ser_method))
  image(cormat[ord,ord],col=cor_cols,breaks=zbreaks,axes=F)
  mtext(text = colnames(cormat)[ord],side = 1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord],side = 3,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord],side = 4,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  box()
  dev.off()
}

gene_cor_analysis=function(ldm,ds_version,min_varmean_per_gene=0.15,
                           min_number_of_UMIs=50,genes_to_exclude=c(),
                           clusters=NULL,samples=NULL,weighted=F,modules_list=NULL,selected_genes=NULL){
  
  if (is.null(clusters)){
    clusters=colnames(ldm$model$models)
  }
  ds=ldm$dataset$ds[[match(ds_version,ldm$dataset$ds_numis)]]
  ds=ds[,ldm$dataset$cell_to_cluster[colnames(ds)]%in%clusters]
  s1=Matrix::rowSums(ds,na.rm=T) 
  s2=Matrix::rowSums(ds^2,na.rm=T)
  mask1=s1>=min_number_of_UMIs&(!rownames(ds)%in%genes_to_exclude)
  message(sum(mask1)," genes passed expression threshold")
  m1=s1[mask1]/ncol(ds)
  m2=s2[mask1]/ncol(ds)
  v1=m2-m1^2
  x=log10(m1)
  breaks=seq(min(x,na.rm=T),max(x,na.rm=T),.2)
  lv=log2(v1/m1)
  llv=split(lv,cut(x,breaks))
  mask_llv=sapply(llv,length)>0   
  z=sapply(llv[mask_llv],min,na.rm=T)
  b=breaks[-length(breaks)]
  b=b[mask_llv]
  lo=loess(z~b)
  high_var_genes=names(which((lv-predict(lo,newdata =x))>min_varmean_per_gene))
  if (length(selected_genes)>0){
    tempList <- names(lv-predict(lo,newdata =x))
    high_var_genes=tempList[tempList%in%selected_genes]
  }
  message(length(high_var_genes)," High var genes")
  if (is.null(modules_list)){
    cormat=get_avg_gene_to_gene_cor(log2(1+ds[high_var_genes,]),ldm$dataset$cell_to_sample[colnames(ds)],samples=samples,weighted = weighted) 
  }else{
    cormat=get_avg_module_to_gene_cor(ds,genes=high_var_genes,modules_list = modules_list,ldm$dataset$cell_to_sample[colnames(ds)],samples=samples,weighted = weighted) 
  }
}

parse_modules=function(ldm,cormat,ds_version,modules_version="",
                       nmods=50,reg=1e-6,mod_size=4,min_mod_cor=0.1,zlim=c(-.9,.9),
                       ord_viz_method="OLO_complete",figures_path="",tables_path=""){
  
  gene_mask2=names(which(apply(cormat,1,quantile,1-mod_size/ncol(cormat),na.rm=T)>=min_mod_cor))
  ds=ldm$dataset$ds[[match(ds_version,ldm$dataset$ds_numis)]]
  mods=cutree(hclust(as.dist(1-cormat[gene_mask2,gene_mask2])),nmods)
  modsl=split(names(mods),mods)
  modsums=aggregate(ds[gene_mask2,],groupings=mods,fun="sum")
  cormat=cor(t(as.matrix(modsums)))
  zbreaks=c(-1,seq(zlim[1],zlim[2],l=99),1)
  cor_cols=colorRampPalette(c("blue","white","red"))(100)
  par(mar=c(3,3,3,3))
  ord=get_order(seriate(as.dist(1-cormat),method="OLO_complete"))
  cormat=cormat[ord,ord]
  colnames(cormat)=1:nmods
  rownames(cormat)=1:nmods
  modsl=modsl[ord]
  names(modsl)=1:nmods
  ord_viz=get_order(seriate(as.dist(1-cormat),method=ord_viz_method))
  par(mar=c(1,1,1,1))
  image(cormat[ord_viz,ord_viz],col=cor_cols,breaks=zbreaks,axes=F)
  mtext(text = colnames(cormat)[ord_viz],side = 1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord_viz],side = 3,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord_viz],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord_viz],side = 4,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  box()
  dev.off()
  write.table(file=paste(tables_path,modules_version,"_modules.txt",sep=""),sapply(modsl,paste,collapse=","),row.names = T,col.names = F,quote=F)
  
  return(modsl)
}

run_gene_modules_analysis=function(){
  
  genes_to_exclude=c(grep("RP",rownames(ileum_ldm$dataset$umitab),val=T),grep("MT-",rownames(ileum_ldm$dataset$umitab),val=T))
  
  cormat=gene_cor_analysis(ileum_ldm,"2000",.1,30,genes_to_exclude=genes_to_exclude,
                           clusters=unlist(ileum_ldm$cluster_sets$T),samples = setdiff(c(inflamed_samples,uninflamed_samples),c("68","69")))
  
  parse_modules(ileum_ldm,cormat,"2000","ileum_T",nmods =50,figures_path=supp_figures_path,tables_path = tables_path,zlim = c(-.5,.5))
  
  genes_to_exclude=c(grep("RP",rownames(ileum_ldm$dataset$umitab),val=T),grep("MT-",rownames(ileum_ldm$dataset$umitab),val=T))
  
  cormat=gene_cor_analysis(ileum_ldm,"2000",.1,30,genes_to_exclude=genes_to_exclude,
                           clusters=unlist(ileum_ldm$cluster_sets$MNP),samples = setdiff(c(inflamed_samples,uninflamed_samples),c("68","69")))
  
  parse_modules(ileum_ldm,cormat,"2000","ileum_MNP",nmods =50,figures_path=supp_figures_path,tables_path = tables_path,zlim = c(-.5,.5))
}

# removing all genes not in the transcription factor list
genes_to_exclude <- rownames(ldm$dataset$umitab)[rownames(ldm$dataset$umitab)%notin%gl$transcription]

test_run=function(){
  genes_to_exclude=c(genes_to_exclude,junk)
  genes_to_exclude <- unique(genes_to_exclude)
  clusters = c(9,8,4,20,24,26,22,19,11)
  cormat=gene_cor_analysis(ldm,"2000",0.9,10,genes_to_exclude=genes_to_exclude,clusters=clusters)
  save_gene_cor_map(cormat,"all_cells2",ser_method = "OLO_complete")
  modules_version="all_cells2"
  ser_method = "OLO_complete"
  temp <<- parse_modules(ldm,cormat,"2000","transcription",nmods = 15,mod_size=4,min_mod_cor=0.15,figures_path=write_path,tables_path=write_path)
}

test_run()

########################
# expression by patient x gene

# plot the fold expression increase over all non-selected clusters of each patient_cluster of cluster of interest
tumor_samps <- as.character(c(240,242,244,246,248,250,252,253,254,258,260,400,412,417,920))

myh11_ecm <- "VWA1,TINAGL1,PODN,GAS6,THSD4,PXDN,LTBP1,CRELD1,COL9A1,COL4A2,COL4A1,COL27A1,TGFBI,IGFBP3,LTBP2,ELN,COL4A5"
fap_ecm <- "NID2,COL5A3,COL18A1,OGN,VCAN,COL6A3,THBS2,SPARC,POSTN,COL5A2,COL5A1,COL3A1,COL1A1,COL12A1,COL11A1,CILP,MFAP5,LUM,TNFAIP6,IGFBP5,FNDC1"
# adh1b_ecm <- "PAPLN,IGFBP4,FMOD,FBLN2,VWA5A,SRPX,MGP,MATN2,LAMA2,COL14A1,PRELP"
# alv_ecm <- "WISP2,PRG4,NPNT,COL21A1,CTGF,FBLN1,DCN,TNXB,SVEP1,SPARCL1,LTBP4,COL6A5,COL13A1,ABI3BP,TNC"

myh11_fap_merge <- "VWA1,TINAGL1,PODN,GAS6,THSD4,PXDN,LTBP1,CRELD1,COL9A1,COL4A2,COL4A1,COL27A1,TGFBI,IGFBP3,LTBP2,ELN,COL4A5,VCAN,COL6A3,THBS2,SPARC,POSTN,COL5A2,COL5A1,COL3A1,COL1A1,COL12A1,COL11A1,CILP,MFAP5,LUM,TNFAIP6,IGFBP5,FNDC1"

ecm_lists <- c(myh11_ecm,fap_ecm)

num <- 1
for (num in 1:length(ecm_lists)){
  
  ecm_genes <- ecm_lists[num]
  ecm_genes <- myh11_fap_merge
  ecm_genes <- strsplit(ecm_genes,",")[[1]]
  
  temp_dataset <- ldm$dataset
  
  noise_adjusted <- pmax(temp_dataset$counts-temp_dataset$noise_counts,0)
  tumor_samps <- as.character(c(240,242,244,246,248,250,252,253,254,258,260,400,412,417,920))
  noise_adjusted <- noise_adjusted[tumor_samps,,]
  noise_adjusted <- noise_adjusted + 1
  
  # merging FAP+ and FAP+aSMA+ CAFs
  temp_dataset$cell_to_cluster <- gsub("^26$","22",temp_dataset$cell_to_cluster)
  temp_dataset$cell_to_cluster <- gsub("^4$","8",temp_dataset$cell_to_cluster)
  
  # clusters to select for
  # i = as.character(c(8,20,22,24,11))
  i = as.character(c(20,22,24,11))
  # i = as.character(c(20,22,11))
  
  merged_cluster_id <- temp_dataset$cell_to_cluster
  merged_sample_id <- temp_dataset$cell_to_sample
  samples <- unique(merged_sample_id)
  
  u=temp_dataset$umitab
  u=u[,merged_cluster_id[colnames(u)]%in%i]
  l=list()
  
  # sample normalization
  for (cluster in i){
    for (samp in samples){
      cluster_mask=colnames(u)[merged_cluster_id[colnames(u)]==cluster]
      sample_mask=colnames(u)[merged_sample_id[colnames(u)]==samp]
      
      mask <- cluster_mask[cluster_mask%in%sample_mask]
      if (length(mask) == 0){
        next
      }
      if (500 < length(mask)){
        l[[paste(cluster,samp)]]=sample(mask,size = 500,replace = F)
      }else{
        l[[paste(cluster,samp)]]=sample(mask,size = length(mask),replace = F)
      }
    }
  }
  cells=unlist(l)
  u=u[,cells]
  # u <- t(t(u)/rowSums(t(u)))
  
  # myh11_excluded <- as.character(c(20,22,24,8))
  # fap_asma_excluded <- as.character(c(20,24,11,8))
  
  myh11_excluded <- as.character(c(20,22,24))
  fap_asma_excluded <- as.character(c(20,24,11))
  
  myh11_excluded <- as.character(c(20,24))
  fap_asma_excluded <- as.character(c(20,24))
  
  # ##
  # myh11_excluded <- as.character(c(20,22))
  # fap_asma_excluded <- as.character(c(20,11))
  # ##
  
  exclusion_lists <- list(myh11_excluded, fap_asma_excluded)
  excluded_clust <- exclusion_lists[[num]]
  
  clust_of_interest <- c("11","22")[num]
  
  # ##
  # excluded_clust <- "8"
  # ##
  
  # cells of interest selected
  temp_cell_names <- colnames(u)[merged_cluster_id[colnames(u)]%in%excluded_clust]
  # data of genes and cells of interest selected
  temp_u <- rowSums(u[ecm_genes,temp_cell_names])/sum(u[,temp_cell_names])
  # taking mean gene value across cells in cluster
  temp_means <- temp_u + 1e-5
  
  ##
  excluded_clust <- exclusion_lists[[num]]
  ##
  
  # clust_of_interest <- i[i%notin%excluded_clust]
  
  # selecting only cells in clusters of interest
  clust_cells <- names(merged_cluster_id[merged_cluster_id%in%clust_of_interest])
  # counts of clusters of interest in each sample
  clust_counts <- table(merged_sample_id[clust_cells])
  # selecting only samples with more than 30 cells in it
  samps <- names(clust_counts[clust_counts > 30])
  
  plot_matrix <- matrix(NA, nrow = length(unique(samps)), ncol = length(ecm_genes))
  rownames(plot_matrix) = unique(samps)
  colnames(plot_matrix) = ecm_genes
  
  ####
  new_u <- noise_adjusted
  ####
  
  # new_u <- t(t(new_u)/rowSums(t(new_u)))
  
  # samp <- samps[1]
  # for (samp in samps){
  #   cell_samps <- names(merged_sample_id[merged_sample_id==samp])
  #   temp_u <- rowSums(new_u[ecm_genes,cell_samps]/sum(new_u[,cell_samps]))
  #   
  #   plot_matrix[samp,ecm_genes] <- temp_u
  # }
  
  ####
  ecm_proportion <- new_u[samps,ecm_genes,clust_of_interest]/rowSums(new_u[samps,,clust_of_interest])
  
  fold_change_over_other <- t(t(ecm_proportion)/temp_means)
  fold_change_over_other <- log2(fold_change_over_other)
  ####
  
  # fold_change_over_other <- t(t(plot_matrix)/temp_means)
  ecm_proportion[,"POSTN"]
  temp_means["POSTN"]
  fold_change_over_other[,"POSTN"]
  ecm_proportion[,"POSTN"]/temp_means["POSTN"]
  
  # fold_change_over_other <- log2(fold_change_over_other)
  
  breaks = seq(ceiling(min(fold_change_over_other)),floor(max(fold_change_over_other)),length.out=8)
  
  if (num == 1){
    breaks = seq(-3,3,length.out=7)
  }else{
    breaks = seq(-3,3,length.out=7)
  }
  
  pheatmap(fold_change_over_other,
           cluster_rows = F,
           cluster_cols = F,
           main = paste("Cluster of interest:", clust_of_interest),
           scale = "none",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
           breaks = breaks)
  
}

########################
# Figure 1B

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list2=gl$lineage_genes
gene_list1=gl$lineage_genes

endo <- c(28,23,3,27,14,21,2,15,10)
peri <- c(5,12)
sm <- c(17,7,25)
fibs <- c(9,8,4,20,24,26,22,19,11)
meso <- c(1)

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- c(fibs, meso, sm, peri, endo)
ldm$cluster_order <- as.character(ldm$cluster_order)
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=150,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_1b",zlim=c(0,2.5))

########################
# Figure S1C

u <- ldm$dataset$umitab
par(mar=c(2,7,1,2))
clusters <- c(1:28)
boxplot(split(log10(colSums(u)),cluster_id),las=2,ylab="log10(#UMIs)",ylim=c(2.7,4.4))

########################
# Figure S1F

figure_S1F=function(fn="figure_S1F"){
  
  gene_mask=apply(ldm$model$models,1,max)>5e-5&rownames(ldm$model$models)%notin%junk
  
  m=log2(1e-5+ldm$model$models[gene_mask,unlist(ldm$cluster_sets[names(ldm$cluster_sets)%notin%"Not good"])])
  
  m=m-rowMeans(m)
  cormat=cor(m)
  d=as.dist(1-cormat)
  order=seriate(d,method = "GW")
  ord=get_order(order)
  
  large_margin=4
  small_margin=.5
  open_plot(path = write_path,fn,plot_type = "pdf",width = 5,height = 5)
  
  par(mar=c(large_margin,large_margin,small_margin,small_margin))
  image(cormat[ord,ord],col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  
  box(lwd=2)
  mtext(rownames(cormat)[ord],1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  mtext(rownames(cormat)[ord],2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  close_plot()
  
  open_plot(path = write_path,fn="figure_colorscale",plot_type = "pdf",width = 4,height = 1)
  par(mar=c(2,1,0,1))
  image(matrix(seq(-1,1,l=100),,1),col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  axis(1)
  close_plot()
}

ldm$cluster_sets <- list()

ldm$cluster_sets$Fib <- c(9,8,4,20,24,26,22,19,11)

allClust <- as.integer(unlist(unique(ldm$cluster_sets)))
cluster_to_cluster_set_with_pdc <- unique(as.vector(unlist(ldm$cluster_sets)))
ldm$cluster_sets$"Not good" <- allClust[allClust%notin%cluster_to_cluster_set_with_pdc]
ldm$cluster_sets$unannotated <- c()

figure_S1F("fib_no_meso_S1F")

ldm$cluster_sets$Meso <- c(1)

allClust <- as.integer(unlist(unique(ldm$cluster_sets)))
cluster_to_cluster_set_with_pdc <- unique(as.vector(unlist(ldm$cluster_sets)))
ldm$cluster_sets$"Not good" <- allClust[allClust%notin%cluster_to_cluster_set_with_pdc]
ldm$cluster_sets$unannotated <- c()

figure_S1F("fib_S1F")

ldm$cluster_sets$PvC <- c(17,5,12,7,25)

allClust <- as.integer(unlist(unique(ldm$cluster_sets)))
cluster_to_cluster_set_with_pdc <- unique(as.vector(unlist(ldm$cluster_sets)))
ldm$cluster_sets$"Not good" <- allClust[allClust%notin%cluster_to_cluster_set_with_pdc]
ldm$cluster_sets$unannotated <- c()

figure_S1F("pvc_fib_S1F")

ldm$cluster_sets$LEC <- c(10)
ldm$cluster_sets$BEC <- c(28,23,3,27,14,21,2,15)

allClust <- as.integer(unlist(unique(ldm$cluster_sets)))
cluster_to_cluster_set_with_pdc <- unique(as.vector(unlist(ldm$cluster_sets)))
ldm$cluster_sets$"Not good" <- allClust[allClust%notin%cluster_to_cluster_set_with_pdc]
ldm$cluster_sets$unannotated <- c()

figure_S1F("all_stroma_S1F")

########################
# Figure 1C

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list2=gl$fib_genes
gene_list1=gl$fib_genes

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- as.character(c(9,8,4,20,24,26,22,11,19,1))
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=250,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_1C",zlim=c(0,2.5))

########################
# Figure 1F

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

mhcii <- gl$mhcii
  
ccl19 <- c("CCL19", "CCL21", "RBP5","TNFSF10","VCAM1")
adh1b <- rev(c("CFD","RGCC","GPC3","INMT","MYH10","CES1","CAV1","LIMCH1","ITGBL1","ROBO2","FIGF"))

# with mhcii
gene_list2= c(ccl19, adh1b, mhcii)
gene_list1= c(ccl19, adh1b, mhcii)
sort_genes <- c("CCL19")

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- c("20")
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

sort_cells <- TRUE
make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=1000,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_1f_mhcii",zlim=c(0,2.5))
sort_cells <- FALSE

# # without mhcii
# gene_list2= c(ccl19, adh1b)
# gene_list1= c(ccl19, adh1b)
# sort_genes <- c("CCL19")

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- c("20")
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

sort_cells <- TRUE
make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=1000,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_1f",zlim=c(0,2.5))
sort_cells <- FALSE

########################
# figure 2A top panel

tmp <- table(ldm$dataset$cell_to_cluster,ldm$dataset$cell_to_sample)

normal_samps <- as.character(c(239,241,243,245,247,249,251,257,259,401,411,416))
tumor_samps <- as.character(c(240,242,244,246,248,250,252,253,254,258,260,400,412,417,920))

cluster_order <- as.character(c(1,9,4,8,20,19,11,24,26,22))
cluster_order <- as.character(c(1,9,4,8,20,19,11,24,26,22))
tmp <- tmp[cluster_order,]
tmp <- t(t(tmp)/rowSums(t(tmp)))

tmp_tumor <- tmp[,tumor_samps]
patient_order <- tmp_tumor["20",]-colSums(tmp_tumor[c("24","26","22"),])

tmp_tumor <- tmp_tumor[,names(sort(patient_order))]
barplot(tmp_tumor,col=c("saddlebrown","pink","royalblue3","royalblue1","khaki3","gray70",
                        "gray50","red1","red3","red3"))

# figure 2A bottom panel

tmp_normal <- tmp[,normal_samps]
normal_order <- as.character(c(401,247,259,245,411,249,257,251,239,416,243,241))
tmp_normal <- tmp_normal[,normal_order]
barplot(tmp_normal,col=c("saddlebrown","pink","royalblue3","royalblue1","khaki3","gray70",
                         "gray50","red1","red3","red3"))

# figure 1D
  
merged_tmp_tumor <- rowSums(tmp_tumor)/sum(tmp_tumor)
merged_tmp_normal <- rowSums(tmp_normal)/sum(tmp_normal)

merged_tmp <- rbind(merged_tmp_normal,merged_tmp_tumor)
barplot(t(merged_tmp),col=c("saddlebrown","pink","royalblue4","royalblue1","khaki3","gray70",
                            "gray50","red1","firebrick","firebrick"))

# figure 1D revision

stage_2_plus <- as.character(c(400,248,260,246,412))
stage_1 <- tumor_samps[tumor_samps%notin%stage_2_plus]

for (temp in list(stage_2_plus,stage_1)){
  
  tumor_samps <- temp
  tmp_tumor <- tmp[,tumor_samps]
  
  merged_tmp_tumor <- rowSums(tmp_tumor)/sum(tmp_tumor)
  merged_tmp_normal <- rowSums(tmp_normal)/sum(tmp_normal)
  
  merged_tmp <- rbind(merged_tmp_normal,merged_tmp_tumor)
  barplot(t(merged_tmp),col=c("saddlebrown","pink","royalblue4","royalblue1","khaki3","gray70",
                              "gray50","red1","firebrick","firebrick"))
  
}

tumor_samps <- as.character(c(240,242,244,246,248,250,252,253,254,258,260,400,412,417,920))
tmp_tumor <- tmp[,tumor_samps]

# figure 1C adjacent and tumor enrichment

ratio_clusts <- (1/colSums(merged_tmp))*t(merged_tmp)
ratio_clusts <- ratio_clusts[as.character(c(9,8,4,20,24,26,22,11,19,1)),]
barplot(t(ratio_clusts),col=c("gray85","gray40"))

# figure S2B adjacent and tumor enrichment
tmp <- table(ldm$dataset$cell_to_cluster,ldm$dataset$cell_to_sample)
cluster_order <- as.character(c(5,12,17,25,7))
tmp <- tmp[cluster_order,]
tmp <- t(t(tmp)/rowSums(t(tmp)))
tmp_tumor <- tmp[,tumor_samps]
tmp_normal <- tmp[,normal_samps]
merged_tmp_tumor <- rowSums(tmp_tumor)/sum(tmp_tumor)
merged_tmp_normal <- rowSums(tmp_normal)/sum(tmp_normal)
merged_tmp <- rbind(merged_tmp_normal,merged_tmp_tumor)

tmp_normal <- tmp_normal[,normal_order]
tmp_tumor <- tmp_tumor[,names(sort(patient_order))]

ratio_clusts <- (1/colSums(merged_tmp))*t(merged_tmp)
ratio_clusts <- ratio_clusts[rev(as.character(c(7,25,17,12,5))),]
barplot(t(ratio_clusts),col=c("gray85","gray40"))

# figure S2A adjacent and tumor enrichment
tmp <- table(ldm$dataset$cell_to_cluster,ldm$dataset$cell_to_sample)
cluster_order <- as.character(c(14,2,15,27,21,3,23,28,10))
tmp <- tmp[cluster_order,]
tmp <- t(t(tmp)/rowSums(t(tmp)))
tmp_tumor <- tmp[,tumor_samps]
tmp_normal <- tmp[,normal_samps]
merged_tmp_tumor <- rowSums(tmp_tumor)/sum(tmp_tumor)
merged_tmp_normal <- rowSums(tmp_normal)/sum(tmp_normal)
merged_tmp <- rbind(merged_tmp_normal,merged_tmp_tumor)

tmp_normal <- tmp_normal[,normal_order]
tmp_tumor <- tmp_tumor[,names(sort(patient_order))]

ratio_clusts <- (1/colSums(merged_tmp))*t(merged_tmp)
ratio_clusts <- ratio_clusts[as.character(c(15,2,14,27,21,3,23,28)),]
barplot(t(ratio_clusts),col=c("gray85","gray40"))

########################
# Figure NA
# histograms

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

# names(gl) <- c(1:length(names(gl)))

msc <- c(gl$'11')
alv <- c(gl$'3')
fap <- c(gl$'14',gl$'15')

gene_sets <- list(msc, alv, fap)
set_names <- c("MSC score", "Alveolar score", "FAP score")

i = rev(c(4,8,9,20,24,26,22))

u=ldm$dataset$umitab
u=u[,cluster_id[colnames(u)]%in%i]

for (num in 1:length(gene_sets)){
  gene_set <- gene_sets[[num]]
  all_genes <- gene_set[gene_set%in%row.names(u)]
  l=list()
  
  # sample normalization
  for (cluster in i){
    for (samp in samples){
      
      cluster_mask=colnames(u)[cluster_id[colnames(u)]==cluster]
      sample_mask=colnames(u)[sample_id[colnames(u)]==samp]
      mask <- cluster_mask[cluster_mask%in%sample_mask]
      
      if (length(mask) == 0){
        next
      }
      if (100 < length(mask)){
        l[[paste(cluster,samp)]]=sample(mask,size = 100,replace = F)
      }else{
        l[[paste(cluster,samp)]]=mask
      }
    }
  }
  
  cells=unlist(l)
  u=u[,cells]
  
  # normalizing by percent total expression
  u <- t(t(u)/rowSums(t(u)))
  # u <- u + 1e-04
  # u <- log10(u)
  
  mods <- list(c(all_genes))
  names(mods) <- set_names[num]
  mod_scores <- colSums(u[unlist(mods[1]),])
  
  layout(matrix(seq(length(mods)),nrow=1))
  par(mgp=c(2,1,0),mar=c(3,1,2,1))
  
  new_i <- rev(i)
  annot_ord <- new_i
  
  cols <- rev(c("red4","red4","red1","khaki3","purple","royalblue1","royalblue3"))
  
  subtype_2_col <- rgb(t(col2rgb(cols)*0.7/255))
  names(subtype_2_col) <- annot_ord
  
  for(mod_iter in names(mods)){
    # xlim <- range(log10(10^-3.5+mod_scores[cluster_id[colnames(u)]%in%new_i]),na.rm=T)
    xlim <- c(-3.5, -0.5)
    plot(1:10,col="white",xlim = xlim, ylim=c(0,8),ylab="",yaxt="n",bty="L",xlab="")
    
    for(iter in rev(seq(length(annot_ord)))){
      cell_mask <- cluster_id[colnames(u)]==annot_ord[iter]
      d <- density(log10(10^-4+mod_scores[cell_mask]), adjust=1)
      polygon(d$x,iter-1+d$y,col=subtype_2_col[iter])
      abline(h=iter-1)
      mtext(paste("MSC_score",": ",mod_iter,sep=""))
    }
  }
}

########################
# figure S1D

# read sample sets
t1=read.table(gene_lists,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

cell_umi_tot <- c()
for (pat in ldm$dataset$numis_before_filtering){
  cell_umi_tot <- c(cell_umi_tot, pat)
}
ldm$dataset$insilico_gating_scores$MITO
mito_umi_count <- ldm$dataset$insilico_gating_scores$MITO * cell_umi_tot

normal_samps <- as.character(c(239,241,243,245,247,249,251,257,259,401,411,416))
tumor_samps <- as.character(c(240,242,244,246,248,250,252,253,254,258,260,400,412,417,920))
samps <- as.character(c(normal_samps,tumor_samps))
sample_summary <- matrix(NA, nrow = 1, ncol = length(samps), dimnames = list(rownames("Mito"), samps))

for (samp in samps){
  umi_tot <- sum(cell_umi_tot[names(ldm$dataset$numis_before_filtering[[samp]])])
  mito_tot <- sum(mito_umi_count[names(ldm$dataset$numis_before_filtering[[samp]])])
  sample_summary[,samp] <- mito_tot/umi_tot
}

melted_sample_summary <- melt(sample_summary)
melted_sample_summary[,"Var2"] <- as.character(melted_sample_summary[,"Var2"])

melted_sample_summary$Var2 <- factor(melted_sample_summary$Var2, levels = c(normal_samps,tumor_samps))

ggplot(melted_sample_summary, aes(y=value, x=Var2)) +
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  theme_bw()

########################
# figure S1E

t1=read.table(gene_lists,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list1 <- gl$modules_all_stroma
gene_list2 <- gl$modules_all_stroma

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- rev(as.character(c(9,8,4,20,24,26,22,11,19,1,17,25,7,5,12,14,2,15,27,21,3,23,28,10)))
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=250,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_s1e",zlim=c(0,4))

########################
# Figure S2A

gene_list1 <- c("EDNRB","IL1RL1","FENDRR","TBX3","CA4","GJA5","HEY1","DKK2","CCL2","IL6","CSF3","SELE",
                "ACKR1","FUT7","INSR","ESM1")
gene_list2 <- gene_list1

# clusters to select for
i = as.character(c(15,14,2,27,21,3,23,28))

gene_sets <- list("gene_list1")
for (set in gene_sets){
  
  all_genes <- gene_list1
  
  gene_list2=all_genes
  gene_list1=all_genes
  
  clusters <- unique(ldm$dataset$cell_to_cluster)
  ldm$cluster_order <- as.character(i)
  clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]
  
  u=ldm$dataset$umitab
  u=u[,cluster_id[colnames(u)]%in%i]
  l=list()
  
  for (cluster in i){
    temp_mask <- c()
    for (samp in samples){
      cluster_mask=colnames(u)[cluster_id[colnames(u)]==cluster]
      sample_mask=colnames(u)[sample_id[colnames(u)]==samp]
      
      mask <- cluster_mask[cluster_mask%in%sample_mask]
      if (length(mask) == 0){
        next
      }
      if (150 < length(mask)){
        mask <- sample(mask,size = 150,replace = F)
      }
      temp_mask <- c(temp_mask, mask)
    }
    l[[paste(cluster,samp)]]=temp_mask
  }
  
  cells=unlist(l)
  u=u[,cells]
  u <- t(t(u)/rowSums(t(u)))
  
  all_genes <- all_genes[all_genes%in%rownames(u)]
  
  bars_matrix <- matrix(nrow=(length(i)), ncol=length(all_genes))
  colnames(bars_matrix) <- all_genes
  rownames(bars_matrix) <- i
  
  for (num in 1:length(i)){
    temp_num <- i[num]
    temp_cell_names <- colnames(u)[cluster_id[colnames(u)]==temp_num]
    temp_u <- u[all_genes,temp_cell_names]
    temp_means <- rowMeans(temp_u)
    bars_matrix[num,] <- temp_means
  }
  
  min_val <- sort(unique(as.vector(bars_matrix)))[2]
  fold_bars_matrix <- t(t(bars_matrix)/colMeans(bars_matrix))
  breaks = seq(0,2.5,length.out=8)
  
  temp <- pheatmap(fold_bars_matrix,
                   cluster_rows = F,
                   cluster_cols = F,
                   scale = "none",
                   color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaks)),
                   breaks = breaks)
  temp
}

########################
# Figure S2B

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list1 <- gl$pvc_modules
gene_list2 <- gene_list1

# clusters to select for
i = as.character(c(17,25,7,5,12))

gene_sets <- list("gene_list1")

for (set in gene_sets){
  
  all_genes <- gene_list1
  
  gene_list2=all_genes
  gene_list1=all_genes
  
  clusters <- unique(ldm$dataset$cell_to_cluster)
  ldm$cluster_order <- as.character(i)
  clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]
  
  u=ldm$dataset$umitab
  u=u[,cluster_id[colnames(u)]%in%i]
  l=list()
  
  for (cluster in i){
    temp_mask <- c()
    for (samp in samples){
      cluster_mask=colnames(u)[cluster_id[colnames(u)]==cluster]
      sample_mask=colnames(u)[sample_id[colnames(u)]==samp]
      
      mask <- cluster_mask[cluster_mask%in%sample_mask]
      if (length(mask) == 0){
        next
      }
      if (150 < length(mask)){
        mask <- sample(mask,size = 150,replace = F)
      }
      temp_mask <- c(temp_mask, mask)
    }
    l[[paste(cluster,samp)]]=temp_mask
  }
  
  cells=unlist(l)
  u=u[,cells]
  u <- t(t(u)/rowSums(t(u)))
  
  all_genes <- all_genes[all_genes%in%rownames(u)]
  
  bars_matrix <- matrix(nrow=(length(i)), ncol=length(all_genes))
  colnames(bars_matrix) <- all_genes
  rownames(bars_matrix) <- i
  
  for (num in 1:length(i)){
    temp_num <- i[num]
    temp_cell_names <- colnames(u)[cluster_id[colnames(u)]==temp_num]
    temp_u <- u[all_genes,temp_cell_names]
    temp_means <- rowMeans(temp_u)
    bars_matrix[num,] <- temp_means
  }
  
  min_val <- sort(unique(as.vector(bars_matrix)))[2]
  fold_bars_matrix <- t(t(bars_matrix)/colMeans(bars_matrix))
  breaks = seq(0,2.5,length.out=8)
  
  temp <- pheatmap(fold_bars_matrix,
                   cluster_rows = T,
                   cluster_cols = T,
                   scale = "none",
                   color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaks)),
                   breaks = breaks)
  temp
  
  # write(temp$tree_col$labels[temp$tree_col$order], paste(write_path, set, "_from_graph.txt", sep=""))
}

########################
# Figure S3E

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list2=gl$iCAF_myCAF_apCAF
gene_list1=gl$iCAF_myCAF_apCAF

i = c(8,20,24,26,22,11,1)
set_name <- "fig_s3e"

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- as.character(i)
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=500,samp_num=200,
                 gene_list1,gene_list2,main_title=set_name,zlim=c(0,2.5))

########################
# Figure S4A

# # cell type to cell type clustering
# 
# u <- ldm$dataset$umitab
# 
# LEC <- c(10)
# BEC <- c(28,23,3,27,14,21,2,15)
# PvC <- c(5,12,7,25,17)
# # Meso <- c(1)
# Fib <- c(9,8,4,20,24,26,22,11,19)
# 
# # selecting only fibs
# clusts <- c(Fib)
# 
# normal_samps <- as.character(c(239,241,243,245,247,249,251,257,259,401,411,416))
# tumor_samps <- as.character(c(240,242,244,246,248,250,252,253,254,258,260,400,412,417,920))
# samps <- as.character(c(normal_samps,tumor_samps))
# 
# #selecting only tumor
# samps <- as.character(c(tumor_samps))
# samp_by_clust <- matrix(NA, nrow = length(clusts), ncol = length(samps), dimnames = list(clusts, samps))
# 
# # troubleshooting
# clust <- clusts[1]
# samp <- samps[1]
# 
# for (clust in clusts){
#   for (samp in samps){
#     samp_cells <- names(sample_id)[sample_id == samp]
#     clust_cells <- names(cluster_id)[cluster_id == clust]
#     cells <- clust_cells[clust_cells%in%samp_cells]
#     samp_by_clust[as.character(clust),as.character(samp)] <- length(cells)
#   }
# }
# # patient IDs are now rows and clusters columns
# samp_by_clust <- t(samp_by_clust)
# 
# adj_samp_by_clust <- samp_by_clust
# 
# for (clust in clusts){
#   if (clust%in%c(LEC,BEC)){
#     print('endo')
#     print(clust)
#     adj_samp_by_clust[,as.character(clust)] <- samp_by_clust[,as.character(clust)]/rowSums(samp_by_clust[,as.character(c(LEC,BEC))])
#   }else if(clust%in%Fib){
#     print('fib')
#     print(clust)
#     adj_samp_by_clust[,as.character(clust)] <- samp_by_clust[,as.character(clust)]/rowSums(samp_by_clust[,as.character(Fib)])
#   }else{
#     print('pvc')
#     print(clust)
#     adj_samp_by_clust[,as.character(clust)] <- samp_by_clust[,as.character(clust)]/rowSums(samp_by_clust[,as.character(PvC)])
#   }
# }
# 
# # regularizing
# adj_samp_by_clust <- adj_samp_by_clust + 0.001
# adj_samp_by_clust <- log10(adj_samp_by_clust)
# # renaming clusts
# cluster_names <- "G:/My Drive/scRNAseq_analysis/compiled/clustering_data_lung4/model_lung_stroma_210412_28_cluster_sets.txt"
# temp <- read.delim(cluster_names,header=TRUE)
# rownames(temp) <- as.character(temp[,"node"])
# 
# colnames(adj_samp_by_clust) <- temp[colnames(adj_samp_by_clust),2]
# 
# pheatmap(adj_samp_by_clust,cluster_rows = F,cluster_cols = T)
# 
# library(corrplot)
# # colnames(adj_samp_by_clust) <- c("MSC","Inf. alv. fib.","Alv. fib.","ADH1B+ CAF","FAP+ CAF","FAP+ CAF","FAP+ CAF","MYH11+ CAF","CLU+ CAF")
# cor_clusts <- cor(adj_samp_by_clust)
# 
# # # cor tests
# # new_clust_by_samp <- adj_samp_by_clust
# # new_clust_by_samp <- cbind(new_clust_by_samp, "new_caf1"=NA)
# # new_clust_by_samp[,"new_caf1"] <- rowMeans(new_clust_by_samp[,colnames(new_clust_by_samp)=="CAF1"])
# # new_clust_by_samp <- t(new_clust_by_samp)
# # 
# # cor.test(new_clust_by_samp["ADH1B",],new_clust_by_samp['CAF2',],)$p.value*3
# # cor.test(new_clust_by_samp["ADH1B",],new_clust_by_samp['new_caf1',],)$p.value*3
# # cor.test(new_clust_by_samp["new_caf1",],new_clust_by_samp['CAF2',])$p.value*3
# new_order <- c("MSC","Inf. alv. fib.","MYH11+aSMA+ CAF","ADH1B","Alv. fib.","CLU+ fib.","FAP+ CAF","FAP+aSMA+ CAF_26","FAP+aSMA+ CAF_22")
# cor_clusts <- cor_clusts[new_order,new_order]
# 
# corrplot(cor_clusts, method = 'color', type = "upper", number.font = 1,
#          tl.col = "black", tl.srt = 45, col=rev(brewer.pal(n=8, name="RdBu")))

########################
# Figure S4C

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list1 <- c()

for (module in 1:15){
  gene_list1 <- c(gene_list1, gl[as.character(module)])
}

gene_list1=unlist(gene_list1, use.names = FALSE)
gene_list2=unlist(gene_list1, use.names = FALSE)

gene_list1 <- rev(gene_list1)
gene_list2 <- rev(gene_list2)

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- as.character(c(9,8,4,20,24,26,22))
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=250,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_s5a",zlim=c(0,4))

########################
# figure 3A and 5E
# immunomod and ecm graphing

gene_sets <- c("ecm_genes","imm_genes","chemo_cyto_only","growth_only","contractile")

for (set in gene_sets){
  if (set == "ecm_genes"){
    all_genes <- gl$ecm
    all_genes <- rev(all_genes)
    i = c(8,4,20,24,26,22,11)
  }else if (set == "contractile"){
    all_genes <- gl$contractile
    all_genes <- rev(all_genes)
    i = c(8,4,20,24,26,22,11)
  }else if (set == "chemo_cyto_only"){
    all_genes <- gl$chemo_cyto_only
    all_genes <- rev(all_genes)
    i = c(8,4,20,24,26,22,11)
  }else if (set == "growth_only"){
    all_genes <- gl$growth_only
    all_genes <- rev(all_genes)
    i = c(8,4,20,24,26,22,11)
  }else{
    all_genes <- gl$ligands
    i = c(20,24,26,22,11)
  }
  
  cluster_id <- ldm$dataset$cell_to_cluster
  
  gene_list2=all_genes
  gene_list1=all_genes
  
  clusters <- unique(ldm$dataset$cell_to_cluster)
  ldm$cluster_order <- as.character(i)
  clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]
  
  u=ldm$dataset$umitab
  u=u[,cluster_id[colnames(u)]%in%i]
  l=list()
  
  for (cluster in i){
    temp_mask <- c()
    for (samp in samples){
      cluster_mask=colnames(u)[cluster_id[colnames(u)]==cluster]
      sample_mask=colnames(u)[sample_id[colnames(u)]==samp]
      
      mask <- cluster_mask[cluster_mask%in%sample_mask]
      if (length(mask) == 0){
        next
      }
      if (150 < length(mask)){
        mask <- sample(mask,size = 150,replace = F)
      }
      temp_mask <- c(temp_mask, mask)
    }
    l[[paste(cluster,samp)]]=temp_mask
  }
  
  cells=unlist(l)
  u=u[,cells]
  u <- t(t(u)/rowSums(t(u)))
  
  all_genes <- all_genes[all_genes%in%rownames(u)]
  
  bars_matrix <- matrix(nrow=(length(i)), ncol=length(all_genes))
  colnames(bars_matrix) <- all_genes
  rownames(bars_matrix) <- i
  
  for (num in 1:length(i)){
    temp_num <- i[num]
    temp_cell_names <- colnames(u)[cluster_id[colnames(u)]==temp_num]
    temp_u <- u[all_genes,temp_cell_names]
    temp_means <- rowMeans(temp_u)
    bars_matrix[num,] <- temp_means
  }
  
  min_val <- sort(unique(as.vector(bars_matrix)))[2]
  fold_bars_matrix <- t(t(bars_matrix)/colMeans(bars_matrix))
  breaks = seq(0,2.5,length.out=8)
  
  if (set%in%c("contractile","chemo_cyto_only","growth_only")){
    temp <- pheatmap(fold_bars_matrix,
                     cluster_rows = F,
                     cluster_cols = T,
                     scale = "none",
                     color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaks)),
                     breaks = breaks)
    temp
  }else{
    temp <- pheatmap(fold_bars_matrix,
                     cluster_rows = F,
                     cluster_cols = F,
                     scale = "none",
                     color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaks)),
                     breaks = breaks)
    temp
  }
  print(colnames(fold_bars_matrix))
}

########################
# Figure mhcii with endo

# immune data saved as idm
load("G:/My Drive/sc_scripts/tumor_immune.RData")

# DC1 = cluster 45
# mature DC = cluster 41

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

temp_dataset <- ldm$dataset

# merging FAP+ and FAP+aSMA+ CAFs
temp_dataset$cell_to_cluster <- gsub("^26$","24",temp_dataset$cell_to_cluster)

temp_dataset$cell_to_cluster <- gsub("^4$","8",temp_dataset$cell_to_cluster)

temp_dataset$cell_to_cluster <- gsub("^2$","15",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^14$","15",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^27$","15",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^21$","15",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^3$","15",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^23$","15",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^28$","15",temp_dataset$cell_to_cluster)

gene_list1 <- rev(c("CD74",gl$'mhcii'))
gene_list2 <- rev(c("CD74",gl$'mhcii'))

# # clusters to select for
# i = as.character(c(9,8,20,22,24,11,15,10))
# cols <- c("pink","royalblue3","khaki3","red1","red3","gray","green","purple","black")
# i = as.character(c(20,22,24,11,15,10))
# cols <- c("khaki3","red1","red3","gray","green","purple","black")
# i = as.character(c(9,8,20,24,22,11))
# cols <- c("pink","royalblue3","khaki3","red1","red3","gray","black")
# i = as.character(c(20,22,24,11))
# cols <- c("khaki3","red1","red3","gray","black")

ADH1B_caf_cells <- names(temp_dataset$cell_to_cluster[temp_dataset$cell_to_cluster=="20"])
ccl19_pos_cells <- temp_dataset$umitab[,temp_dataset$umitab["CCL19",] > 0]
ccl19_pos_cells <- ccl19_pos_cells[,ADH1B_caf_cells[ADH1B_caf_cells%in%colnames(ccl19_pos_cells)]]
ccl19_cell_names <- colnames(ccl19_pos_cells)
temp_dataset$cell_to_cluster[ccl19_cell_names] <- gsub("^20$","99",temp_dataset$cell_to_cluster[ccl19_cell_names])

# all fib with endo
i = as.character(c(99,9,8,4,20,22,24,11,15,10,1))
cols <- c("black","pink","royalblue3","royalblue1","khaki3","red1","red3","darkgray","green","purple","brown")
# cafs with endo
# i = as.character(c(99,20,22,24,11,15,10))
# cols <- c("black","khaki3","red1","red3","darkgray","green","purple")
# # all fib
# i = as.character(c(99,9,8,20,24,22,11))
# cols <- c("black","pink","royalblue3","khaki3","red1","red3","darkgray")
# cafs
# i = as.character(c(99,20,22,24,11,1))
# cols <- c("black","khaki3","red1","red3","darkgray","brown")
# # adh1b only
# i = as.character(c(99,20))
# cols <- c("black","khaki3")

merged_cluster_id <- temp_dataset$cell_to_cluster
merged_sample_id <- temp_dataset$cell_to_sample
samples <- unique(merged_sample_id)

u=temp_dataset$umitab
u=u[,merged_cluster_id[colnames(u)]%in%i]
l=list()

# sample normalization
for (cluster in i){
  for (samp in samples){
    cluster_mask=colnames(u)[merged_cluster_id[colnames(u)]==cluster]
    sample_mask=colnames(u)[merged_sample_id[colnames(u)]==samp]
    
    mask <- cluster_mask[cluster_mask%in%sample_mask]
    if (length(mask) == 0){
      next
    }
    if (100 < length(mask)){
      l[[paste(cluster,samp)]]=sample(mask,size = 100,replace = F)
    }else{
      l[[paste(cluster,samp)]]=sample(mask,size = length(mask),replace = F)
    }
  }
}
cells=unlist(l)
u=u[,cells]
u <- t(t(u)/rowSums(t(u)))

genesAll <- gene_list1
bars_matrix <- matrix(nrow=(length(i)), ncol=length(genesAll))
colnames(bars_matrix) <- genesAll
rownames(bars_matrix) <- i

for (num in 1:length(i)){
  temp_num <- i[num]
  temp_cell_names <- colnames(u)[merged_cluster_id[colnames(u)]==temp_num]
  temp_u <- u[genesAll,temp_cell_names]
  temp_means <- rowMeans(temp_u)
  bars_matrix[num,] <- temp_means
}

min_val <- sort(unique(as.vector(bars_matrix)))[2]

genes_to_keep <- rev(gene_list1)
bars_matrix <- bars_matrix[,genes_to_keep]
bars_matrix_reg <- bars_matrix + 1e-4

# adding immune and selecting only DC1 cells
immune_cells <- idm$dataset$umitab
immune_cells <- immune_cells[,names(idm$dataset$cell_to_cluster)[idm$dataset$cell_to_cluster == "45"]]
immune_cells <- t(t(immune_cells)/rowSums(t(immune_cells)))
mhc_immune <- rowSums(immune_cells[colnames(bars_matrix_reg),])
mhc_immune <- mhc_immune + 1e-4
bars_matrix_reg <- rbind(bars_matrix_reg,mhc_immune)

melted_bars <- melt(bars_matrix_reg)
melted_bars[,1] <- as.character(melted_bars[,1])
# changing format to counts per 100000 UMI
melted_bars[,3] <- melted_bars[,3]*10000

# i = as.character(c(99,20,22,24,11,15,10))
# cols <- c("black","khaki3","red1","red3","darkgray","green","purple")

# i = as.character(c(99,20,22,24,11))
# cols <- c("black","khaki3","red1","red3","darkgray")

# i = as.character(c(99,20,22,24,11,1))
# cols <- c("black","khaki3","red1","red3","darkgray","tan4")

i = as.character(c(99,1,15,10,"mhc_immune"))
cols <- c("black","tan4","darkgreen","purple","lightgray")

melted_bars <- melted_bars[melted_bars[,1]%in%i,]

clust_ord <- as.character(i)
melted_bars$Var1 <- factor(melted_bars$Var1, levels = clust_ord)

# graphing parameters
log_graphs_e <- 0:5
breaks <- 10^(log_graphs_e)
minor_breaks <- rep(1:9, length(breaks)-1)*(10^rep(log_graphs_e[1:length(breaks)-1], each=9))

ggplot(melted_bars, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  scale_fill_manual("legend",values = cols) +
  theme_bw() +
  scale_y_continuous(trans='log10', breaks = c(1,10,100,1000,10000,100000)) +
  # scale_y_log10(limits = c(1,1e5)) +
  # scale_y_log10(limits = c(1,1e1)) +
  geom_hline(yintercept = 1e5) +
  theme(panel.grid.minor = element_blank())

########################
# Figure 2F

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

temp_dataset <- ldm$dataset

# merging FAP+ and FAP+aSMA+ CAFs
temp_dataset$cell_to_cluster <- gsub("^26$","22",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^24$","22",temp_dataset$cell_to_cluster)

temp_dataset$cell_to_cluster <- gsub("^4$","8",temp_dataset$cell_to_cluster)

gene_list1 <- gl$'fig_2f_genes'
gene_list2 <- gl$'fig_2f_genes'

# clusters to select for
i = as.character(c(9,8,20,22))

merged_cluster_id <- temp_dataset$cell_to_cluster
merged_sample_id <- temp_dataset$cell_to_sample
samples <- unique(merged_sample_id)

u=temp_dataset$umitab
u=u[,merged_cluster_id[colnames(u)]%in%i]
l=list()

# sample normalization
for (cluster in i){
  for (samp in samples){
    cluster_mask=colnames(u)[merged_cluster_id[colnames(u)]==cluster]
    sample_mask=colnames(u)[merged_sample_id[colnames(u)]==samp]
    
    mask <- cluster_mask[cluster_mask%in%sample_mask]
    if (length(mask) == 0){
      next
    }
    if (100 < length(mask)){
      l[[paste(cluster,samp)]]=sample(mask,size = 100,replace = F)
    }else{
      l[[paste(cluster,samp)]]=sample(mask,size = length(mask),replace = F)
    }
  }
}
cells=unlist(l)
u=u[,cells]
u <- t(t(u)/rowSums(t(u)))

genesAll <- gene_list1
bars_matrix <- matrix(nrow=(length(i)), ncol=length(genesAll))
colnames(bars_matrix) <- genesAll
rownames(bars_matrix) <- i

for (num in 1:length(i)){
  temp_num <- i[num]
  temp_cell_names <- colnames(u)[merged_cluster_id[colnames(u)]==temp_num]
  temp_u <- u[genesAll,temp_cell_names]
  temp_means <- rowMeans(temp_u)
  bars_matrix[num,] <- temp_means
}
min_val <- sort(unique(as.vector(bars_matrix)))[2]

genes_to_keep <- rev(gene_list1)
bars_matrix <- bars_matrix[,genes_to_keep]
bars_matrix_reg <- bars_matrix + 1e-4

melted_bars <- melt(bars_matrix_reg)
melted_bars[,1] <- as.character(melted_bars[,1])
# changing format to counts per 100000 UMI
melted_bars[,3] <- melted_bars[,3]*10000

cols <- c("purple","royalblue3","khaki3","red1","red2","red3")

clust_ord <- as.character(c(9,8,20,22))
melted_bars$Var1 <- factor(melted_bars$Var1, levels = clust_ord)

ggplot(melted_bars, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  scale_fill_manual("legend",values = cols) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  theme(panel.grid.minor = element_blank())



########################
# figure 3A and 5E
# immunomod and ecm gene selection
# figure 2f must be run first

# clusters
i = c(9,8,22,11)
# minimum expression
x <- 1/2000*0.05
# fold change
fold_change <- 3

########################
# first pass gene selection
# immune ligands

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

all_genes <- gl$all_ligands

u=temp_dataset$umitab
u=u[,merged_cluster_id[colnames(u)]%in%i]
l=list()

for (cluster in i){
  temp_mask <- c()
  for (samp in samples){
    
    cluster_mask=colnames(u)[merged_cluster_id[colnames(u)]==cluster]
    sample_mask=colnames(u)[merged_sample_id[colnames(u)]==samp]
    
    mask <- cluster_mask[cluster_mask%in%sample_mask]
    
    if (length(mask) == 0){
      next
    }
    if (150 < length(mask)){
      mask <- sample(mask,size = 150,replace = F)
    }
    temp_mask <- c(temp_mask, mask)
  }
  l[[paste(cluster,samp)]]=temp_mask
}

cells=unlist(l)
u=u[,cells]

# normalizing by percent total expression
u <- t(t(u)/rowSums(t(u)))
all_genes <- all_genes[all_genes%in%row.names(u)]
all_de_genes <- c()

# generating cluster mean expression values for genes
gene_mat <- matrix(NA, nrow = length(i), ncol = length(all_genes),
                   dimnames = list(c(i), all_genes))
i <- as.character(i)

for (num in i){
  cluster1 <- colnames(u)[colnames(u)%in%names(merged_cluster_id[merged_cluster_id == num])]
  gene_mat[num,all_genes] <- rowMeans(u[all_genes,cluster1])
}

# checking for minimum expression
expr_thresh <- colMaxs(gene_mat)>x
gene_mat <- gene_mat[,expr_thresh]

# isolating genes with expression fold change x over at least one other cluster
all_de_genes <- c()
for (num in i){
  for (var in i){
    if (num == var){
      next
    }
    comp1genes <- c()
    comp2genes <- c()
    comp1 <- (gene_mat[num,] > gene_mat[var,]*fold_change)
    comp1genes <- names(comp1)[comp1]
    
    all_de_genes <- c(all_de_genes,comp1genes)
  }
}
all_de_genes <- unique(all_de_genes)
genesAll <- all_de_genes

mmp_exclude <- c("MMP1","MMP12","MMP13","MMP2","MMP7","MMP9")
genesAll <- genesAll[genesAll%notin%mmp_exclude]

# write(genesAll,paste(write_path, "imm_genes.txt", sep=""))
# write(genesAll,paste("C:/Users/johna/Desktop/imm_genes.txt", sep=""))

########################
# first pass gene selection
# transcription

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

all_genes <- gl$transcription

u=temp_dataset$umitab
u=u[,merged_cluster_id[colnames(u)]%in%i]
l=list()

for (cluster in i){
  temp_mask <- c()
  for (samp in samples){
    
    cluster_mask=colnames(u)[merged_cluster_id[colnames(u)]==cluster]
    sample_mask=colnames(u)[merged_sample_id[colnames(u)]==samp]
    
    mask <- cluster_mask[cluster_mask%in%sample_mask]
    
    if (length(mask) == 0){
      next
    }
    if (150 < length(mask)){
      mask <- sample(mask,size = 150,replace = F)
    }
    temp_mask <- c(temp_mask, mask)
  }
  l[[paste(cluster,samp)]]=temp_mask
}

cells=unlist(l)
u=u[,cells]

# normalizing by percent total expression
u <- t(t(u)/rowSums(t(u)))
all_genes <- all_genes[all_genes%in%row.names(u)]
all_de_genes <- c()

# generating cluster mean expression values for genes
gene_mat <- matrix(NA, nrow = length(i), ncol = length(all_genes),
                   dimnames = list(c(i), all_genes))
i <- as.character(i)

for (num in i){
  cluster1 <- colnames(u)[colnames(u)%in%names(merged_cluster_id[merged_cluster_id == num])]
  gene_mat[num,all_genes] <- rowMeans(u[all_genes,cluster1])
}

# checking for minimum expression
expr_thresh <- colMaxs(gene_mat)>x
gene_mat <- gene_mat[,expr_thresh]

# isolating genes with expression fold change x over at least one other cluster
all_de_genes <- c()
for (num in i){
  for (var in i){
    if (num == var){
      next
    }
    comp1genes <- c()
    comp2genes <- c()
    comp1 <- (gene_mat[num,] > gene_mat[var,]*fold_change)
    comp1genes <- names(comp1)[comp1]
    
    all_de_genes <- c(all_de_genes,comp1genes)
  }
}
all_de_genes <- unique(all_de_genes)
genesAll <- all_de_genes

write.csv(genesAll,paste(write_path, "transcription_genes.csv", sep=""))

########################
# first pass gene selection
# ecm ligands

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

all_genes <- gl$all_ecm

u=temp_dataset$umitab
u=u[,merged_cluster_id[colnames(u)]%in%i]
l=list()

for (cluster in i){
  temp_mask <- c()
  for (samp in samples){
    
    cluster_mask=colnames(u)[merged_cluster_id[colnames(u)]==cluster]
    sample_mask=colnames(u)[merged_sample_id[colnames(u)]==samp]
    
    mask <- cluster_mask[cluster_mask%in%sample_mask]
    
    if (length(mask) == 0){
      next
    }
    if (150 < length(mask)){
      mask <- sample(mask,size = 150,replace = F)
    }
    temp_mask <- c(temp_mask, mask)
  }
  l[[paste(cluster,samp)]]=temp_mask
}

cells=unlist(l)
u=u[,cells]

# normalizing by percent total expression
u <- t(t(u)/rowSums(t(u)))
all_genes <- all_genes[all_genes%in%row.names(u)]
all_de_genes <- c()

# generating cluster mean expression values for genes
gene_mat <- matrix(NA, nrow = length(i), ncol = length(all_genes),
                   dimnames = list(c(i), all_genes))
i <- as.character(i)

for (num in i){
  cluster1 <- colnames(u)[colnames(u)%in%names(merged_cluster_id[merged_cluster_id == num])]
  gene_mat[num,all_genes] <- rowMeans(u[all_genes,cluster1])
}

# checking for minimum expression
expr_thresh <- colMaxs(gene_mat)>x
gene_mat <- gene_mat[,expr_thresh]

# isolating genes with expression fold change x over at least one other cluster
all_de_genes <- c()
for (num in i){
  for (var in i){
    if (num == var){
      next
    }
    comp1genes <- c()
    comp2genes <- c()
    comp1 <- (gene_mat[num,] > gene_mat[var,]*fold_change)
    comp1genes <- names(comp1)[comp1]
    
    all_de_genes <- c(all_de_genes,comp1genes)
  }
}
all_de_genes <- unique(all_de_genes)
genesAll <- all_de_genes

# write(genesAll,paste(write_path, "ecm_genes.txt", sep=""))

########################
# first pass gene selection
# contractile

t1=read.table(gene_lists,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

all_genes <- gl$all_contractile

u=temp_dataset$umitab
u=u[,merged_cluster_id[colnames(u)]%in%i]
l=list()

for (cluster in i){
  temp_mask <- c()
  for (samp in samples){
    
    cluster_mask=colnames(u)[merged_cluster_id[colnames(u)]==cluster]
    sample_mask=colnames(u)[merged_sample_id[colnames(u)]==samp]
    
    mask <- cluster_mask[cluster_mask%in%sample_mask]
    
    if (length(mask) == 0){
      next
    }
    if (150 < length(mask)){
      mask <- sample(mask,size = 150,replace = F)
    }
    temp_mask <- c(temp_mask, mask)
  }
  l[[paste(cluster,samp)]]=temp_mask
}

cells=unlist(l)
u=u[,cells]

# normalizing by percent total expression
u <- t(t(u)/rowSums(t(u)))
all_genes <- all_genes[all_genes%in%row.names(u)]
all_de_genes <- c()

# generating cluster mean expression values for genes
gene_mat <- matrix(NA, nrow = length(i), ncol = length(all_genes),
                   dimnames = list(c(i), all_genes))
i <- as.character(i)

for (num in i){
  cluster1 <- colnames(u)[colnames(u)%in%names(merged_cluster_id[merged_cluster_id == num])]
  gene_mat[num,all_genes] <- rowMeans(u[all_genes,cluster1])
}

# checking for minimum expression
expr_thresh <- colMaxs(gene_mat)>x
gene_mat <- gene_mat[,expr_thresh]

# isolating genes with expression fold change x over at least one other cluster
all_de_genes <- c()
for (num in i){
  for (var in i){
    if (num == var){
      next
    }
    comp1genes <- c()
    comp2genes <- c()
    comp1 <- (gene_mat[num,] > gene_mat[var,]*fold_change)
    comp1genes <- names(comp1)[comp1]
    
    all_de_genes <- c(all_de_genes,comp1genes)
  }
}
all_de_genes <- unique(all_de_genes)
genesAll <- all_de_genes

write(genesAll,paste(write_path, "ecm_genes.txt", sep=""))

########################
# Figure 3B

load(lcam_metadata)

temp_1 <- sample_annots

cell_metadata <- cbind(cell_metadata,temp_1[as.character(cell_metadata$sample_ID),c("patient_ID","tissue","library_chemistry","prep","prime")])
cell_metadata <- cbind(cell_metadata,annots_list[as.character(cell_metadata$cluster_assignment),])
prep_mask <- cell_metadata$prime=="3"
cell_metadata <- cell_metadata[prep_mask,]
patient_tissue <- apply(temp_1[,c("patient_ID","tissue")],1,paste,collapse="_")
names(patient_tissue) <- rownames(temp_1)
patient_tissue <- patient_tissue[as.character(cell_metadata$sample_ID)]
tab <- table(patient_tissue,cell_metadata$sub_lineage)

# removing 'junk' cells
tab <- tab[,colnames(tab)%notin%"epi_endo_fibro_doublet"]
tab <- tab/rowSums(tab)

for(norm_group in c("T","B&plasma","MNP","lin_neg")){
  if (norm_group == "lin_neg"){
    tab_tmp <- tab[,annots_list$norm_group[match(colnames(tab),annots_list$sub_lineage)]==norm_group]
    tab_tmp <- tab_tmp[,-ncol(tab_tmp)]
    tab_tmp <- tab_tmp / rowSums(tab[,-ncol(tab)])
    tab[,colnames(tab_tmp)] <- tab_tmp
  }else{
    tab_tmp <- tab[,annots_list$norm_group[match(colnames(tab),annots_list$sub_lineage)]==norm_group]
    tab_tmp <- tab_tmp[,-ncol(tab_tmp)]
    tab_tmp <- tab_tmp / rowSums(tab[,-ncol(tab)])
    tab[,colnames(tab_tmp)] <- tab_tmp
  }
}

LCAM <- c("T_activated","IgG","MoMac-II")
LCAM_score <- rowSums(log(tab[,LCAM]+1e-2))

resting_clusts <- c("B","cDC2","AZU1_mac","Tcm/naive_II","cDC1","AM")
resting_score <- rowSums(log(tab[,resting_clusts]+1e-2))

pat_ord_tumor <- order((LCAM_score-resting_score)[grep("Tumor",rownames(tab),v=T)])
pat_ord_normal <- order((LCAM_score-resting_score)[grep("Normal",rownames(tab),v=T)])

tab_tumor <- tab[grep("Tumor",rownames(tab)),]
tab_normal <- tab[grep("Normal",rownames(tab)),]

mat_normal <- t(tab_normal[pat_ord_normal,rev(c(resting_clusts,LCAM))])
mat_tumor <- t(tab_tumor[pat_ord_tumor,rev(c(resting_clusts,LCAM))])
clust.means <- rowMeans(cbind(mat_normal,mat_tumor))
mat_normal <- log2((1e-2+mat_normal)/(1e-2+clust.means))
mat_tumor <- log2((1e-2+mat_tumor)/(1e-2+clust.means))

thresh <- 2
mat_tumor[mat_tumor < -thresh] <- -thresh
mat_tumor[mat_tumor > thresh] <- thresh
mat_tumor <- mat_tumor/thresh
mat_tumor <- round((mat_tumor+1)/2*49)+1
col_tumor <- bluered(50)[min(mat_tumor):max(mat_tumor)]

h <- seq(-0.5,8.5)[4]/8

colnames(mat_tumor) <- gsub("_Tumor","",colnames(mat_tumor))
saved_lcam_order <- colnames(mat_tumor)

t1=read.table(lcam_id_conversion,row.names = 1,
              stringsAsFactors = F)

mat_tumor <- mat_tumor[,colnames(mat_tumor)%in%rownames(t1)]
colnames(mat_tumor) <- t1[colnames(mat_tumor),]

t1=read.table(gene_lists,row.names = 1,
              stringsAsFactors = F)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

adh1b <- gl$'adh1b'
adh1b_myh11 <- gl$'adh1b_myh11'
fap <- gl$'fap'

l_adh1b <- length(adh1b)
l_adh1b_myh11 <- length(adh1b_myh11)
l_fap <- length(fap)

# aligning adh1b and fap values for heatmap
fib <- c(rep(1,l_adh1b),
         rep(1,l_adh1b_myh11),
         rep(50,l_fap))
pats <- c(adh1b,adh1b_myh11,fap)
names(fib) <- pats

mat_tumor <- mat_tumor[,colnames(mat_tumor)%in%pats]
fib <- fib[names(fib)%in%colnames(mat_tumor)]
pats <- pats[pats%in%colnames(mat_tumor)]

mat_tumor <- rbind(mat_tumor,Fib=NA)
mat_tumor["Fib", pats] <- fib

fib_sorted <- sort(mat_tumor["Fib",])
mat_tumor <- mat_tumor[,names(fib_sorted)]

for (row in rownames(mat_tumor)){
  mat_tumor[row,] <- scale(mat_tumor[row,])
}

image(t(mat_tumor),col=col_tumor,xaxt="n",yaxt="n")
mtext(colnames(mat_tumor),side=1,at=seq(0,1,1/(ncol(mat_tumor)-1)),las=2,line=0.25,cex=0.7)
mtext("Tumor samples",cex=2,line=0.5)
box()

# correlation value calculation
names(LCAM_score) <- gsub("_Tumor","",names(LCAM_score))
saved_lcam_order <- as.character(saved_lcam_order)

t1=read.table(lcam_id_conversion,row.names = 1,
              stringsAsFactors = F)
saved_lcam_order <- t1[saved_lcam_order,]
saved_lcam_order <- saved_lcam_order[!is.na(saved_lcam_order)]
names(LCAM_score) <- t1[names(LCAM_score),]
LCAM_score <- LCAM_score[!is.na(names(LCAM_score))]

saved_lcam_order <- as.character(saved_lcam_order[saved_lcam_order%in%colnames(mat_tumor)])
cor.test(mat_tumor["Fib",saved_lcam_order], LCAM_score[saved_lcam_order])

########################

########################
# Figure NA

cd34_mat <- c(1,0,0,0,0,1)
cd10_mat <- c(0,1,0,0,0,0)
adh1b_mat <- c(1,1,1,0,0,0)
fap_mat <- c(0,0,0,1,1,0)
asma_mat <- c(0,0,0,0,1,1)
myh11_mat <- c(0,0,0,0,0,1)
ihc_stain_matrix <-cbind(cd34_mat, cd10_mat, adh1b_mat, fap_mat, asma_mat, myh11_mat)
colnames(ihc_stain_matrix) <- c("CD34","CD10","ADH1B","FAP","aSMA","MYH11")

pheatmap(ihc_stain_matrix,
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c("#FFFFFF","#000000"))



