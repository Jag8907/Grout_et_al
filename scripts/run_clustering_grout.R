library(methods)
scClustering_dir="martin_et_al_cell_2019/clustering/"
source(paste(scClustering_dir,"/clustering3.r",sep=""))


rps=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/ribsomial_genes.csv",sep=""),stringsAsFactors = F)[,1])
igk=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGK_genes.csv",sep=""),stringsAsFactors = F)[,1])
igh=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGH_genes.csv",sep=""),stringsAsFactors = F)[,1])
igl=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGL_genes.csv",sep=""),stringsAsFactors = F)[,1])
ep=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/EP_markers.csv",sep=""),stringsAsFactors = F)[,1])
variable_seg_light_chain=c(igk,igl,igh)
malat=c("MALAT1")
xist="XIST"
jchain="JCHAIN"
#highly_expressed=c("IGKC","IGHA1","JCHAIN")
cell_cycle=strsplit("STMN1,CDC20,ANP32E,CKS1B,NUF2,ASPM,UBE2T,CENPF,RRM2,SGOL2,SGOL1,SMC4,CENPE,CCNA2,CCNB1,PTTG1,HMMR,MXD3,HIST1H4C,CENPW,H2AFV,CKS2,SMC2,ZWINT,CDK1,MKI67,C12or\
                    f75,CDKN3,NUSAP1,CCNB2,KIAA0101,PRC1,ARL6IP1,PLK1,CENPN,AURKB,TOP2A,KPNA2,HN1,TK1,BIRC5,TYMS,TPX2,UBE2C,CALM3,UBE2S,GTSE1,CLSPN,CDCA7,CENPU,GMNN,MCM7,CCDC34,CDCA5,HELLS,TMEM106C\
                    ,IDH2,SNRNP25,CDT1,PCNA,UHRF1,ASF1B",",")[[1]]
hla=strsplit("HLA-F,HLA-G,HLA-A,HLA-E,HLA-C,HLA-B,HLA-DRA,HLA-DRB5,HLA-DRB1,HLA-DQA1,HLA-DQB1,HLA-DQB1-AS1,HLA-DQA2,HLA-DQB2,HLA-DOB,HLA-DMB,HLA-DMA,HLA-DOA,HLA-DPA1,HLA-DPB1","\
             ,")[[1]]
mts=c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB')
metallothionein=c('MT1HL1','MT2A','MT1E','MT1M','MT1A','MT1B','MT1F','MT1G','MT1H','MT1X')
stress=strsplit("DUSP10,CSRNP1,IRF1,HSP90AB1,DNAJA1,TUBB4B,DDIT4,RGCC,SOCS1,UBB,PMAIP1,SERTAD1,MYADM,HSPA6,HSPE1,HSPA1A,HSPA1B,BAG3,HSPH1,HSP90AA1,DNAJB1,PIK3R3",",")[[1]]
mtrnr=strsplit("MTRNR2L11,MTRNR2L12,MTRNR2L13,MTRNR2L6,MTRNR2L10,MTRNR2L8,MTRNR2L7,MTRNR2L5,MTRNR2L4,MTRNR2L1,MTRNR2L3",",")[[1]]
am_genes=strsplit("S100A4,FN1,PPBP,SDC2,CTSL,CD163,RNASE1,CCL13,CCL18,GNLY,CXCL8,HBEGF,SCGB3A2,CD74,SCGB3A1,HLA-DQB1,HLA-DQA2,HLA-DQB2,XIST,SFTPC,LSP1,SCGB1A1,EGR2,SFTPA1,PHLDA1\
                  ,IFI27,CCL3,TGM2",",")[[1]]

gene_lists=c()
gene_lists["Epi"]="AGR3,SFTA3,TMC5,ELF3,CLDN7,EPCAM,FOLR1,NKX2-1,TSTD1,SFTA2,CPM,CLDN4,ATF7IP2,PARD6B,TNNC1,KRT8,CEACAM6,GGTLC1,KRT19,FXYD3,CLIC3,KRT7,MSLN,SUSD2,AQP4,CLIC5,SCEL,GPRC5A,TACSTD2,KRT18,RAB11FIP1,CLDN18,COL4A3,ALCAM,EFNA1,ANKRD1,ARHGEF26,SLC39A8,AGER,CST6,UPK3B,RTKN2,GOLGA8B,RP11-27M24.1,VEGFA,ANOS1,DENND3,ADRB2,ANKRD29,SPOCK2,C4BPA,SFTPD,CLDN3,CTSE,SLC34A2,MPZL2,MALL,S100A14,PGC,SERPINA1,SFTPA1,CD24,NAPSA,AGR2,SFTPA2,CXCL17,MUC1,SFTPB,WFDC2,PPP6R1,SCGB3A1,SCGB3A2,SPINK5,SPINK1,FGGY,PIGR,CRYM,SMIM22,AQP3"
gene_lists["Mito"]="MTRNR2L1,MT-CO2,MTRNR2L8,MTRNR2L12,MT-CO3,MT-ND4,MT-CYB,MT-ND3,MT-ATP6,MT-CO1,NEAT1,MT-ND2,MT-ND1,MT-ATP8,MT-ND4L,MT-ND5,MT-ND6"
gene_lists["Tcell"]="CD69,CD2,GZMB,CTSW,XCL1,GZMA,GNLY,CD7,CCL5,PRF1,KLRD1,TRDC,AC092580.4,KLRB1,NKG7,KLRC1,GZMK,CD160,FCGR3A,CD247,XCL2,ID2,B3GNT7,CMC1,SPON2,FGFBP2,HOPX,PLAC8,MYOM2,CLIC3,IFITM2,KLRF1,KLRC2,KRT86,KRT81,GZMH,CST7,HCST,LAIR2"
gene_lists["RBC"]="CSTA,HBD,S100A7,KRT17,NPPC,PAGE2,HBA1,HBB,HBA2"
gene_lists["Plasma"]="ANKRD28,HERPUD1,C15orf48,IGHA1,JCHAIN,IGKV4-1,IGLC7,AL928768.3,IGHA2,IGHM,IGLV6-57,ISG20,SEC11C,SSR4,IGKV1-12,TNFRSF17,IGLC6,IGLL1,FKBP11,IGLL5,IGLC5,MZB1,DERL3,IGLC2,IGLC3,IGHD,IGHE,IGHEP1,IGKV1-5,IGKC,IGLV3-1,IGHG2,IGHGP,IGHG4,IGHG3,IGHG1"

gene_lists=strsplit(gene_lists,",")

insilico_gating=list()

insilico_gating$MITO=list()
insilico_gating$MITO$genes=gene_lists$Mito
insilico_gating$MITO$interval=c(0,0.25)

insilico_gating$Plasma=list()
insilico_gating$Plasma$genes=gene_lists$Plasma
insilico_gating$Plasma$interval=c(0,0.05)

insilico_gating$RBC=list()
insilico_gating$RBC$genes=gene_lists$RBC
insilico_gating$RBC$interval=c(0,0.2)

insilico_gating$Tcell=list()
insilico_gating$Tcell$genes=gene_lists$Tcell
insilico_gating$Tcell$interval=c(0,0.015)

insilico_gating$Epi=list()
insilico_gating$Epi$genes=gene_lists$Epi
insilico_gating$Epi$interval=c(0,0.01)

samps=c("239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","257","258","259","260","400","401","411","412")
libs=annots[samps,"amp_batch_ID"]
print(libs)
data_l=read_multiple_mtx("data/",libs,cell_interval = c(800,1e6),noise_interval = c(100,1e6))
save(data_l,file="tmp/stroma_data_l.rd")

model_name="lung_stroma_210412_28"
pdf(paste(model_name,".pdf",sep=""))
cluster(data_l_path="tmp/stroma_data_l.rd",model_name = model_name,k=28,load_seed=F,
	running_mode="LSF_seeding",
#	running_mode="local",
        params=list(
          train_set_size  = 7500,
          test_set_size = 2500,
          insilico_gating = insilico_gating,
          genes_excluded = c(mts,mtrnr,malat,xist,metallothionein,variable_seg_light_chain,hla,jchain),
          seeding_varmean_quantile=.92,
          min_umis_per_seeding_gene=50,
          genes_excluded_from_seeding = c(rps),
          init_min_umis_quantile_range=c(0,.5),
          fixed_downsampling_min_umis=NA,
          reg=5e-6,n_init_seeds=1000,
          max_n_cores=1,
          init_method="TGL_kmeans",
          km_reg=.2,
          init_alpha_noise=0.04,
          thresh_fraction_cells_moved=0.005))
dev.off()
