# if (download_data){
#   download_files(pipeline_path)
# }
# 
# ldm_filename=paste(pipeline_path,"inputs","ldm_2000_v2.RData",sep="/")
# lcam_metadata_filename=paste(pipeline_path,"inputs","cell_metadata.rd",sep="/")
# 
# load(ldm_filename,envir = env)
# load(lcam_metadata_filename,envir = env)


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

write_path <- "D:/Google Drive/Paper website/Scripts/outputs/"

# source files
ecm_immunomod <- "D:/Google Drive/Paper website/Scripts/inputs/immunomod/"
genes_to_exclude <- "D:/Google Drive/Paper website/Scripts/inputs/excluded_genes.txt"
image_genes <- "D:/Google Drive/Paper website/Scripts/inputs/image_gene_lists.txt"
lig_heatmap_genes <- "D:/Google Drive/Paper website/Scripts/inputs/lig_rec_ligands.txt"
ecm_heatmap_genes <- "D:/Google Drive/Paper website/Scripts/inputs/ecm.txt"
modules_3a_genes <- "D:/Google Drive/Paper website/Scripts/inputs/modules_fig_3a.txt"
modules_all_stroma <- "D:/Google Drive/Paper website/Scripts/inputs/modules_all_stroma.txt"
sample_sets <- "D:/Google Drive/Paper website/Scripts/inputs/sample_sets.txt"
lcam_cafs <- "D:/Google Drive/Paper website/Scripts/inputs/lcam_cafs.txt"
stroma_ldm <- "D:/Google Drive/Paper website/Scripts/inputs/ldm_2000_v2.RData"
# LCAM source codes
data_dir <- "D:/Google Drive/Paper website/Scripts/inputs"
sample_annots <- read.csv(file.path(data_dir,"sample_annots.csv"),r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv(file.path(data_dir,"model_lung_190606_annot_lists.csv"),r=1,h=1,stringsAsFactors = F)
annots_list$norm_group[annots_list$lineage=="MNP"] <- "MNP"
lcam_metadata <- "D:/Google Drive/Paper website/Scripts/inputs/cell_metadata.rd"
lcam_id_conversion <- "D:/Google Drive/Paper website/Scripts/inputs/conversion.txt"

# genes to remove from gene module analysis
t1=read.table(genes_to_exclude,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
junk <- c(gl$RPs, gl$Junk, gl$Ribosome, gl$IgVar, gl$Mito)

########################
# functions

`%notin%` <- Negate(`%in%`)
fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
fisher.z2r <- function(z) {(exp(2 * z) - 1)/(1 + exp(2 * z))}








# get_gene_symbol_convetors=function(path="D:/Google Drive/scDissector-master/inst/extdata/"){
#   hgnc<-read.delim(paste(path,"hgnc_complete_set.txt",sep=""),header = T,stringsAsFactors = F)
#   old_symbol=ifelse(hgnc[,"prev_symbol"]=="",hgnc[,"symbol"],hgnc[,"prev_symbol"])
#   l_old_symbol=strsplit(old_symbol,"\\|")
#   old_symbol2=unlist(l_old_symbol)
#   
#   new_symbol=hgnc[,"symbol"]
#   new_symbol2=rep(new_symbol,sapply(l_old_symbol,length))
#   
#   gene_symbol_old2new<-new_symbol2
#   gene_symbol_new2old<-old_symbol
# }
# get_gene_symbol_convetors()












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
  colgrad_abs_file="../extdata/colors_viridis.txt"
}

colgrad_abs<<-read.table(colgrad_abs_file,stringsAsFactors=F)[,1]

close_plot=function(){
  dev.off()
}

Relative_or_Absolute="Relative"
colgrad<<-c(colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100))
sample_cols<<-rep(paste("#",read.table("D:/Google Drive/Paper website/Scripts/inputs/sample_colors.txt",
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
    if (length(sort_genes == 1)){
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

test_run=function(){
  genes_to_exclude=junk
  clusters = c(9,8,4,20,24,26,22)
  cormat=gene_cor_analysis(ldm,"2000",0.9,10,genes_to_exclude=genes_to_exclude,clusters=clusters)
  save_gene_cor_map(cormat,"all_cells2",ser_method = "OLO_complete")
  modules_version="all_cells2"
  ser_method = "OLO_complete"
  temp <<- parse_modules(ldm,cormat,"2000","all_fibs",nmods = 15,mod_size=4,min_mod_cor=0.15,figures_path=write_path,tables_path=write_path)
}

test_run()

########################
# Figure 1B

t1=read.table(image_genes,
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
# Figure 1C

ldm$cluster_sets <- list()

ldm$cluster_sets$LEC <- c(10)
ldm$cluster_sets$BEC <- c(28,23,3,27,14,21,2,15)
ldm$cluster_sets$PvC <- c(5,12,7,25)
ldm$cluster_sets$Meso <- c(1)
ldm$cluster_sets$Fib <- c(9,8,4,20,24,26,22,19,11,17)

allClust <- as.integer(unlist(unique(ldm$cluster_sets)))
cluster_to_cluster_set_with_pdc <- unique(as.vector(unlist(ldm$cluster_sets)))
ldm$cluster_sets$"Not good" <- allClust[allClust%notin%cluster_to_cluster_set_with_pdc]
ldm$cluster_sets$unannotated <- c()

figure_1c=function(){
  
  gene_mask=apply(ldm$model$models,1,max)>5e-5&rownames(ldm$model$models)%notin%junk
  
  m=log2(1e-5+ldm$model$models[gene_mask,unlist(ldm$cluster_sets[names(ldm$cluster_sets)%notin%"Not good"])])
  
  m=m-rowMeans(m)
  cormat=cor(m)
  d=as.dist(1-cormat)
  order=seriate(d,method = "GW")
  ord=get_order(order)
  
  large_margin=4
  small_margin=.5
  open_plot(path = write_path,fn="figure_1c",plot_type = "pdf",width = 5,height = 5)
  
  par(mar=c(large_margin,large_margin,small_margin,small_margin))
  image(cormat[ord,ord],col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  
  box(lwd=2)
  mtext(rownames(cormat)[ord],1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  mtext(rownames(cormat)[ord],2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  close_plot()
  
  open_plot(path = write_path,fn="figure_1c_colorscale",plot_type = "pdf",width = 4,height = 1)
  par(mar=c(2,1,0,1))
  image(matrix(seq(-1,1,l=100),,1),col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  axis(1)
  close_plot()
}

figure_1c()

########################
# Figure 1D

t1=read.table(image_genes,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list2=gl$fib_genes
gene_list1=gl$fib_genes

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- as.character(c(9,8,4,20,24,26,22,11,19,1))
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=250,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_1d",zlim=c(0,2.5))

########################
# Figure 1G

ccl19 <- c("CCL19", "CCL21", "RBP5","TNFSF10","VCAM1")
adh1b <- rev(c("CFD","RGCC","GPC3","INMT","MYH10","CES1","CAV1","LIMCH1","ITGBL1","ROBO2","FIGF"))

gene_list2= c(ccl19, adh1b)
gene_list1= c(ccl19, adh1b)
sort_genes <- c("CCL19")

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- c("20")
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

sort_cells <- TRUE
make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=1000,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_1g",zlim=c(0,2.5))
sort_cells <- FALSE

########################
# figure 2A left panel

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

# figure 2A right panel

tmp_normal <- tmp[,normal_samps]
normal_order <- as.character(c(401,247,259,245,411,249,257,251,239,416,243,241))
tmp_normal <- tmp_normal[,normal_order]
barplot(tmp_normal,col=c("saddlebrown","pink","royalblue3","royalblue1","khaki3","gray70",
                         "gray50","red1","red3","red3"))

# figure 1E

merged_tmp_tumor <- rowSums(tmp_tumor)/sum(tmp_tumor)
merged_tmp_normal <- rowSums(tmp_normal)/sum(tmp_normal)

merged_tmp <- rbind(merged_tmp_normal,merged_tmp_tumor)
barplot(t(merged_tmp),col=c("saddlebrown","pink","royalblue4","royalblue1","khaki3","gray70",
                            "gray50","red1","firebrick","firebrick"))

# figure 1D adjacent and tumor enrichment

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

tmp_normal <- tmp_normal[,normal_order]
tmp_tumor <- tmp_tumor[,names(sort(patient_order))]

ratio_clusts <- (1/colSums(merged_tmp))*t(merged_tmp)
ratio_clusts <- ratio_clusts[as.character(c(17,25,7,5,12)),]
barplot(t(ratio_clusts),col=c("gray85","gray40"))

########################
# Figure 3A

t1=read.table(modules_3a_genes,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

names(gl) <- c(1:length(names(gl)))

msc <- c(gl$'11')
alv <- c(gl$'3')
fap <- c(gl$'14',gl$'15')

gene_sets <- list(msc, alv, fap)
set_names <- c("MSC score", "Alveolar score", "FAP score")

i = rev(c(9,4,8,20,24,26,22))

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
  
  cols <- rev(c("red4","red4","red1","khaki3","royalblue1","royalblue3","purple"))
  
  subtype_2_col <- rgb(t(col2rgb(cols)*0.7/255))
  names(subtype_2_col) <- annot_ord
  
  for(mod_iter in names(mods)){
    xlim <- range(log10(10^-3.5+mod_scores[cluster_id[colnames(u)]%in%new_i]),na.rm=T)
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

# read all saved sample sets
s1=read.table(sample_sets,stringsAsFactors = F,row.names = 1)
il=strsplit(s1[,1],",")
names(il)=rownames(s1)
# saving the normal and tumor stroma sets
samps <- paste("data_",il$stroma,".rd",sep="")
# compiled samples
samps <- paste("D:/Google Drive/Paper website/Scripts/inputs/stroma_lung_compiled",samps,sep="/")
samples <- samps
# table with mito genes
t1=read.table(genes_to_exclude,stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)
#samples
normal_samps <- as.character(c(239,241,243,245,247,249,251,257,259,401,411,416))
tumor_samps <- as.character(c(240,242,244,246,248,250,252,253,254,258,260,400,412,417,920))
samps <- as.character(c(normal_samps,tumor_samps))
sample_summary <- matrix(NA, nrow = nrow(u1), ncol = length(samps), dimnames = list(rownames(u1), samps))
# loading the samples
u1=NA
samp_env=new.env()
samp <- samples[1]
load(samp,envir=samp_env)
u1=samp_env$umitab
# loading
for (num in 1:length(samples)){
  samp_env=new.env()
  
  samp <- samples[num]
  samp_num <- samps[num]
  
  print(samp)
  load(samp,envir=samp_env)
  sample_summary[,samp_num] <- rowMeans(samp_env$umitab)
}

u <- u1
rm(u1)

mito_sums <- colSums(sample_summary[gl$Mito,])/colSums(sample_summary)
mito_sums <- melt(mito_sums)
mito_sums <- cbind(mito_sums, rownames(mito_sums))
colnames(mito_sums) <- c("value","patients")
mito_sums$patients <- factor(mito_sums$patients, levels = mito_sums$patients)

ggplot(mito_sums, aes(y=value, x=patients)) +
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  theme_bw()
rm(u)

########################
# figure S1E

t1=read.table(modules_all_stroma,stringsAsFactors = F)
gl=c(t1[,1])

gene_list1 <- gl
gene_list2 <- gl

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

t1=read.table(image_genes,
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

t1=read.table(image_genes,
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
# Figure S5A

t1=read.table(modules_3a_genes,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list2=unlist(gl, use.names = FALSE)
gene_list1=unlist(gl, use.names = FALSE)

gene_list1 <- rev(gene_list1)
gene_list2 <- rev(gene_list2)

clusters <- unique(ldm$dataset$cell_to_cluster)
ldm$cluster_order <- as.character(c(9,8,4,20,24,26,22))
clusters_to_exclude <- clusters[clusters%notin%ldm$cluster_order]

make_truth_plots(ldm,clusters_to_exclude=clusters_to_exclude,clust_num=250,samp_num=100,
                 gene_list2,gene_list1,main_title="fig_s5a",zlim=c(0,4))

########################
# figure 3E and 5D
# immunomod and ecm graphing

t1=read.table(image_genes,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_sets <- c("ecm_genes","imm_genes")

for (set in gene_sets){
  all_genes <- gl[set][[1]]
  
  if (set == "ecm_genes"){
    all_genes <- read.table(paste(ecm_immunomod,"ecm_genes.txt",sep=""))
    all_genes <- all_genes[[1]]
    all_genes <- rev(all_genes)
    i = c(8,4,20,24,26,22,11)
  }else{
    all_genes <- read.table(paste(ecm_immunomod,"imm_genes.txt",sep=""))
    all_genes <- all_genes[[1]]
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
  
  temp <- pheatmap(fold_bars_matrix,
                   cluster_rows = F,
                   cluster_cols = F,
                   scale = "none",
                   color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaks)),
                   breaks = breaks)
  temp
}

########################
# Figure 3B
# Gene expression table

temp_dataset <- ldm$dataset

# merging FAP+ and FAP+aSMA+ CAFs
temp_dataset$cell_to_cluster <- gsub("^26$","22",temp_dataset$cell_to_cluster)
temp_dataset$cell_to_cluster <- gsub("^24$","22",temp_dataset$cell_to_cluster)

# modules to look at
t1=read.table(modules_3a_genes,
              stringsAsFactors = F,row.names = 1)
gl=strsplit(t1[,1],",")
names(gl)=rownames(t1)

gene_list2=unlist(gl, use.names = FALSE)
gene_list1=unlist(gl, use.names = FALSE)

gene_list1 <- rev(c(gene_list1,"APOD","PLA2G2A"))
gene_list2 <- rev(c(gene_list2,"APOD","PLA2G2A"))

# clusters to select for
i = as.character(c(9,8,4,20,22))
set_name <- "modules_all"
all_genes <- gene_list1

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

genes_to_keep <- rev(c("ACTA2","COL1A1","COL3A1","BGN","ASPN",
                       "MACF1","TCF21","FIGF","C3","APOD","PLA2G2A"))
bars_matrix <- bars_matrix[,genes_to_keep]
bars_matrix_reg <- bars_matrix + 1e-4

melted_bars <- melt(bars_matrix_reg)
melted_bars[,1] <- as.character(melted_bars[,1])
# changing format to counts per 100000 UMI
melted_bars[,3] <- melted_bars[,3]*10000

cols <- c("purple","royalblue3","royalblue1","khaki3","red1","red2","red3")

clust_ord <- as.character(c(9,8,4,20,22))
melted_bars$Var1 <- factor(melted_bars$Var1, levels = clust_ord)

ggplot(melted_bars, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  scale_fill_manual("legend",values = cols) +
  theme_bw() +
  scale_y_continuous(trans='log10')

########################
# figure 3E and 5D
# ecm and immunomodulatory genes

# clusters
i = c(20,22,11)

# minimum expression
x <- 1/2000*0.05
# fold change
fold_change <- 3

########################
# figure 3E first pass gene selection

# immune
t1=read.table(lig_heatmap_genes,
              stringsAsFactors = F)
all_genes <- unlist(strsplit(t1[,2],","))

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

########################
# figure 5D first pass gene selection

# ecm
t1=read.table(ecm_heatmap_genes,
              stringsAsFactors = F)
all_genes <- unlist(strsplit(t1[,2],","))

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
# Figure 3F

load(lcam_metadata)

cell_metadata <- cbind(cell_metadata,sample_annots[as.character(cell_metadata$merged_sample_id),
                                                   c("patient_ID","tissue","library_chemistry","prep","Project","prime")])
cell_metadata <- cbind(cell_metadata,annots_list[as.character(cell_metadata$cluster_assignment),])
prep_mask <- cell_metadata$library_chemistry=="V2" & cell_metadata$prep=="beads" & cell_metadata$Project=="lung" & cell_metadata$prime=="3"
patient_tissue <- apply(sample_annots[,c("patient_ID","tissue")],1,paste,collapse="_")
names(patient_tissue) <- rownames(sample_annots)
patient_tissue <- patient_tissue[as.character(cell_metadata$merged_cluster_id)]
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

# without regard to prep first:
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

# removing zilionis and lambrechts data sets and selecting only tumor
exclude_tums <- grep("zilionis",colnames(mat_tumor),)
exclude_tums <- c(exclude_tums,grep("Lambrechts",colnames(mat_tumor),))
include_tums <- c(grep("Tumor",colnames(mat_tumor),))
include_tums <- include_tums[-exclude_tums]
colnames(mat_tumor) <- gsub("_Tumor","",colnames(mat_tumor))
mat_tumor <- mat_tumor[,include_tums]
saved_lcam_order <- colnames(mat_tumor)

t1=read.table(lcam_id_conversion,row.names = 1,
              stringsAsFactors = F)
mat_tumor <- mat_tumor[,colnames(mat_tumor)%in%rownames(t1)]
colnames(mat_tumor) <- t1[colnames(mat_tumor),]

t1=read.table(lcam_cafs,row.names = 1,
              stringsAsFactors = F)
lcam_caf_ids <- unlist(strsplit(t1[,1],","))

adh1b <- strsplit(t1["adh1b",],",")[[1]]
adh1b_myh11 <- strsplit(t1["adh1b_myh11",],",")[[1]]
fap <- strsplit(t1["fap",],",")[[1]]

l_adh1b <- length(adh1b)
l_adh1b_myh11 <- length(adh1b_myh11)
l_fap <- length(fap)

# aligning adh1b and fap values for heatmap
fib <- c(rep(1,l_adh1b),
         rep(1,l_adh1b_myh11),
         rep(50,l_fap))
pats <- as.character(lcam_caf_ids)
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
saved_lcam_order <- saved_lcam_order[saved_lcam_order%in%colnames(mat_tumor)]
cor.test(mat_tumor["Fib",saved_lcam_order], LCAM_score[saved_lcam_order])

########################







