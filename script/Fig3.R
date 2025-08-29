## library
library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)
library(patchwork)
library(cowplot)
library(DoubletFinder)

#### process ###################################################################################
## celldex ref data
ref<-readRDS("../ref.rds")

## function to remove doublet cell by DoubletFinder
removeDoublets <- function(data) {
  sweep <- paramSweep_v3(data,  sct = F)
  sweep <- summarizeSweep(sweep, GT = FALSE)
  bcmvn <- find.pK(sweep)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  # find nExp
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  # DoubletRate = ncol(data)*8*1e-6 #https://www.jianshu.com/p/6770c6a05287
  DoubletRate = ncol(data)*25*1e-6 # manmade
  nExp_poi <- round(DoubletRate * ncol(data@assays$RNA@data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  # find Doublet
  data <-
    doubletFinder_v3(
      data,
      PCs = 1:10,
      pN = 0.25,
      pK = mpK,
      nExp = nExp_poi,
      reuse.pANN = FALSE,
      sct = F
    )
  return(data)
}

## function to process raw data
process<-function(samplename){
  data<-Read10X(paste0(path,samplename))
  data<-CreateSeuratObject(data,project = samplename)
  
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
  data[["percent.rp"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")
  data[["gene.per.umi"]] <- log10(data$nFeature_RNA)/log10(data$nCount_RNA)
  
  data <- subset(data,subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 10 & gene.per.umi>0.8)
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- ScaleData(data,features=VariableFeatures(data)) 
  data <- RunPCA(data, features = VariableFeatures(object = data))
  p<-ElbowPlot(data)
  ggsave(p,filename = paste0(samplename,"_elbowplot.pdf"),width = 5,height = 5)
  
  data<-removeDoublets(data)
  return(data)
}

## procee sample data
F12<-process("F12")
F17<-process("F17")
F19<-process("F19")
F15<-process("F15")

## save merged data
saveRDS(list(F12=F12,F17=F17,F19=F19,F15=F15),file = "NC-CS.F12-17-19-15.rds")

data<-readRDS("NC-CS.F12-17-19-15.rds")

## filter doublet cell 
{
  F12<-data$F12
  F12$doubletFinder<-F12$DF.classifications_0.25_0.19_4371
  F12$pANN_0.25_0.19_4371<-NULL
  F12$DF.classifications_0.25_0.19_4371<-NULL
  F12<-subset(F12,doubletFinder=="Singlet")
  
  F17<-data$F17
  F17$doubletFinder<-F17$DF.classifications_0.25_0.005_2864
  F17$pANN_0.25_0.005_2864<-NULL
  F17$DF.classifications_0.25_0.005_2864<-NULL
  F17<-subset(F17,doubletFinder=="Singlet")
  
  F15<-data$F15
  F15$doubletFinder<-F15$DF.classifications_0.25_0.13_1517
  F15$pANN_0.25_0.13_1517<-NULL
  F15$DF.classifications_0.25_0.13_1517<-NULL
  F15<-subset(F15,doubletFinder=="Singlet")
  
  F19<-data$F19
  F19$doubletFinder<-F19$DF.classifications_0.25_0.27_1458
  F19$pANN_0.25_0.27_1458<-NULL
  F19$DF.classifications_0.25_0.27_1458<-NULL
  F19<-subset(F19,doubletFinder=="Singlet")
}

## sample integrate
{
  sct_func<-function(x){
    x<-SCTransform(x,vars.to.regress = c("percent.mt"),verbose = TRUE,conserve.memory = TRUE)
    return(x)
  }
  
  data.list <- list(F12=F12,F17=F17,F19=F19,F15=F15)
  
  remove(data,F12,F15,F17,F19)
  
  data.list <- lapply(X = data.list, FUN = sct_func)
  features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features)
  data <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  remove(data.list)
}

## add sample annotation
{
  sample<-data$orig.ident
  sample[sample=="F12"]<-"L-NC-F12"
  sample[sample=="F17"]<-"L-NC-F17"
  sample[sample=="F15"]<-"L-CS-F15"
  sample[sample=="F19"]<-"L-CS-F19"
  
  sample<-factor(sample,levels = c("L-NC-F12", "L-NC-F17", "L-CS-F15", "L-CS-F19"))
  data$sample<-sample
  
  PYMT<-data@assays$SCT@data["transgene-PyMT",]
  data$PYMT<-PYMT
  tissue<-as.character(lapply(as.character(data$sample), function(x){return(strsplit(x,"-")[[1]][1])}))
  data$tissue<-tissue
  group<-as.character(lapply(as.character(data$sample), function(x){return(strsplit(x,"-")[[1]][2])}))
  data$group<-factor(group,levels = c(  'NC',  'CS'  ))
}

## standard pipline
data <- RunPCA(data, verbose = FALSE)
ElbowPlot(data,ndims = 50)
# JackStrawPlot(data)
data <- RunUMAP(data, dims = 1:30)
data <- FindNeighbors(data, dims = 1:30, verbose = FALSE,reduction = "umap")
data <- FindClusters(data, verbose = FALSE,resolution=0.8)

## function of Singler annotation
annoBySingler<-function(data,ref.se=ref$ImmGen,fineOrMain="fine",addClusterToAnno=TRUE,check.cluster=TRUE){
  # ref.se<-ImmGenData()
  matrix <- GetAssayData(data, slot = 'data')
  if(fineOrMain=="fine"){
    singler.cluster <-SingleR(matrix,ref = ref.se,
                              labels = ref.se$label.fine,
                              clusters = data@meta.data$seurat_clusters)
  }else if(fineOrMain=="main"){
    singler.cluster <-SingleR(matrix,ref = ref.se,
                              labels = ref.se$label.main,
                              clusters = data@meta.data$seurat_clusters)
  }
  celltype = data.frame(ClusterID = rownames(singler.cluster),celltype = singler.cluster$labels,stringsAsFactors = F)
  ## load to meta.data
  data@meta.data$celltype = 'NA'
  if(addClusterToAnno){
    for (i in 1:nrow(celltype)) {
      if(check.cluster){
        data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- paste0("c",i,"_",celltype$celltype[i])
      }else{
        data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- paste0("c",i-1,"_",celltype$celltype[i])
      }
      
    }
  }else{
    for (i in 1:nrow(celltype)) {
      data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]
    }
  }
  
  return(data)
}

## Singler annotation
data<-annoBySingler(data,ref.se=ref$MRG,fineOrMain="main",addClusterToAnno=FALSE)
data$celltype_MRG<-data$celltype
data<-annoBySingler(data,ref.se=ref$ImmGen,fineOrMain="main",addClusterToAnno=FALSE)
data$celltype_IGD<-data$celltype

## save data
saveRDS(data,"Lung.cs-nc.integrated.rds")

#### cell assignment ###################################################################################
## CHANGE LABEL
data$celltype<-data$celltype_IGD
data@meta.data[data@meta.data$seurat_clusters==13,]$celltype<-"Epithelial cells"
data@meta.data[data@meta.data$seurat_clusters==30,]$celltype<-"Epithelial cells"
data@meta.data[data@meta.data$seurat_clusters==15,]$celltype<-"MoMacDc"
data@meta.data[data@meta.data$celltype=="Macrophages",]$celltype<-"MoMacDc"
data@meta.data[data@meta.data$celltype=="DC",]$celltype<-"MoMacDc"
data@meta.data[data@meta.data$celltype=="Monocytes",]$celltype<-"MoMacDc"
data@meta.data[data@meta.data$celltype=="NKT",]$celltype<-"T cells"
data@meta.data[data@meta.data$celltype=="Tgd",]$celltype<-"T cells"
data@meta.data[data@meta.data$celltype=="ILC",]$celltype<-"T cells"
data@meta.data[data@meta.data$celltype=="Stromal cells",]$celltype<-"Fibroblasts"


p1<-DimPlot(data,reduction = 'umap',group.by = "celltype_IGD",repel=T, label=T, label.size=4,cols = c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")) )
p2<-DimPlot(data,reduction = 'umap',group.by = "celltype_MRG",repel=T, label=T, label.size=4,cols = c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")) )
p3<-DimPlot(data,reduction = 'umap',group.by = "celltype_map",repel=T, label=T, label.size=4,cols = c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")) )

data$celltype<-factor(data$celltype,levels = c("B cells","NK cells","T cells","Neutrophils","Basophils","MoMacDc","Fibroblasts","Endothelial cells","Epithelial cells"))
DimPlot(data,reduction = 'umap',group.by = "celltype",repel=T, label=T, label.size=4,cols = c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")) )


## all markers
{
  allfeat <- c(
    #celltype
    "Ptprc", #CD45
    "Cd53",  #CD53\
    #B
    "Cd79a", #CD79 b
    "Ms4a1", #CD20 B
    #ILC
    #"Itgae",#CD103
    #"Ncam1",#cd56
    # "Il4",
    "Il5",
    # "Il9",
    # "Il17",
    #"Cxcr3",#cs183
    # "Ccl3",
    #NK
    "Ncr1" ,  #CD335 NK
    #NK &ILC
    "Klrd1", #CD94 NK
    "Klre1", #CD94 NK
    
    "Il2rb",#CD122
    #T NKT
    "Cd226",#CD226
    #T
    "Cd3d",#CD3
    "Cd3e",#CD3
    "Cd3g",#CD3
    "Cd4",
    "Cd8a",
    #Gran
    "Cd24a",#CD24
    
    "Csf3r",#CD114
    
    "Retnlg",
    "S100a8",
    "S100a9",
    "Cxcr2",  #cd128
    "Il3ra",#CD123 Baso
    #Bas
    "Il13",
    "Il4",
    "Mcpt8",
    "Ccl3",
    "Gata2",
    #DC
    
    
    #"Lamp3",#CD208
    "Cd207",
    "Cd209a",
    "Itgax",#CD11c
    "Cd86",#CD86h
    #M
    "Adgre1",#F4/80
    "Cd68",  #CD68 M
    # "Bst2",#CD317
    # "Cxcr3",
    
    # "Entpd1",#CD39 M
    "Mrc1",   #CD206 M
    "Csf1r",  #CD115 M
    "Fcgr1" ,  #CD64 M
    "Fcgr3" ,  #CD16 M
    # "Lilra5",#Cd85
    
    #Mono
    "Itgal",#CD11A
    # "Cd14",
    #Fibr
    # "Il1r1",#CD121a
    #"Lrp1",#CD91
    "Pdgfra","Pdgfrb",
    'Col1a1',
    'Col1a2',
    'Igf2',
    #stro
    # "Cd248",
    # "Thy1",#CD90
    # "Jag1",#CD339
    # "Bst1",#cs157
    # "Stag1",
    # "Stag2",
    # "Stag3",
    # "Stim1",
    # "Stim2",
    #endo
    # "Cd34",#
    "Pecam1",#Cd31
    "Vegfa",# Vegf
    "Cdh5",#CD144
    "Eng",#CD105
    "Erg","Cldn5",
    #ALVEOLAR
    "Sftpa1", "Sftpb", "Ager",
    #"Aqp4",
    #EPI
    "Epcam",#CD326
    #"Krt14",
    #"Krt19",
    # "Krt7","Krt8", #cytokeratin
    #CC
    "PYMT"
  )
} 

DotPlot(data,features = allfeat,group.by = "celltype",cols = c("grey","#E41A1C"))+
  coord_flip()+theme(axis.text.x =  element_text(angle = 90, hjust = 1,vjust = 0.5))
DotPlot(data,features = allfeat,assay = "SCT",group.by = "celltype",cols = c("grey","#E41A1C"))+
  coord_flip()+theme(axis.text.x =  element_text(angle = 90, hjust = 1,vjust = 0.5))
DotPlot(data,features = allfeat,group.by = "seurat_clusters",cols = c("grey","#E41A1C"))+
  coord_flip()+theme(axis.text.x =  element_text(angle = 90, hjust = 1,vjust = 0.5))
DotPlot(data,features = allfeat,assay = "RNA",group.by = "seurat_clusters",cols = c("grey","#E41A1C"))+
  coord_flip()+theme(axis.text.x =  element_text(angle = 90, hjust = 1,vjust = 0.5))

saveRDS(data,"Lung.cs-nc.integrated.rds")
data.bk<-data

## function of cell percent plot
plotBarLine<-function(data,name,plot2sam.width=3.5,plot4sam.width=4,plotline.width=4){
  data.bk<-data
  data<-data.bk@meta.data
  sample<-levels(data$sample)[1:4]
  # ct<-unique(data$celltype)
  ct<-levels(data$celltype)
  df<-c()
  for (c in ct) {
    out<-c()
    for (s in sample) {
      num<-filter(data,celltype==c&sample==s)%>% nrow()
      out<-c(out,num)
    }
    #out<-out/sum(out)
    df<-rbind(df,out)
  }
  colnames(df)<-sample
  rownames(df)<-ct
  
  df2<-apply(df,2,function(x){return(x/sum(x))})
  
  df<-data.frame(value=c(as.numeric(df2[,1]),as.numeric(df2[,2]),as.numeric(df2[,3]),as.numeric(df2[,4])),
                 type=factor(c(rownames(df2),rownames(df2),rownames(df2),rownames(df2)),levels = ct),
                 sample=factor(
                   c(rep('L-NC-F12',nrow(df2)),
                     rep('L-NC-F17',nrow(df2)),
                     rep('L-CS-F15',nrow(df2)),
                     rep('L-CS-F19',nrow(df2))),
                   levels = c('L-NC-F12', 'L-NC-F17', 'L-CS-F15', 'L-CS-F19')))
  
  p1<-ggplot(df,aes(x=sample,y=value,fill=type))+
    geom_bar(stat = 'identity')+
    theme_bw() +
    scale_fill_manual(values=c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2"))) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10,angle = 90, hjust = 1,vjust = 0.5),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black")
    )
  ggsave(p1,filename = paste0(name,".cellprop-bar.4.pdf"),width = plot4sam.width,height = 5)
  
  L_NC<-(df2[,1]+df2[,2])/2
  L_CS<-(df2[,3]+df2[,4])/2
  
  df<-data.frame(value=c(as.numeric(L_NC),as.numeric(L_CS)),
                 type=factor(c(names(L_NC),names(L_CS)),levels = ct),
                 sample=factor(c(rep('NC',length(L_NC)),rep('CS',length(L_CS))),
                               levels = c("NC","CS")))
  library(RColorBrewer)
  
  p1<-ggplot(df,aes(x=sample,y=value,fill=type))+
    geom_bar(stat = 'identity')+
    theme_bw() +
    scale_fill_manual(values=c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2"))) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black")
    )
  ggsave(p1,filename = paste0(name,".cellprop-bar.pdf"),width = plot2sam.width,height = 5)
  
  
  p2<-ggplot(df,aes(x=sample,y=value,group=type,shape=type))+
    geom_line(aes(color=type),size=0.5)+
    geom_point(size=3,aes(colour=type))+
    scale_shape_manual(values = c(0:20))+
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.line = element_line(colour = "black"),
      panel.border=element_rect(fill=NA,color="black", size=0.2, linetype="solid")
    )+
    scale_color_manual(values =c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")))
  ggsave(p2,filename = paste0(name,".cellprop-line.pdf"),width = plotline.width,height = 5)
}
## function of plot fisher bar plot
plotFisherBar<-function(data,name,plot.width=8){
  data<-data@meta.data
  # ct<-unique(data$celltype)
  ct<-levels(data$celltype)
  p_df<-c()
  for (cell in ct) {
    print(cell)
    
    wt<-filter(data,group=="NC")%>%nrow()
    wt_1<-filter(data,celltype==cell&group=="NC")%>%nrow()
    wt_2<-filter(data,celltype!=cell&group=="NC")%>%nrow()
    
    nc<-filter(data,group=="CS")%>%nrow()
    nc_1<-filter(data,celltype==cell&group=="CS")%>%nrow()
    nc_2<-filter(data,celltype!=cell&group=="CS")%>%nrow()
    
    wt_nc<-fisher.test(matrix(c(wt_1,wt_2,nc_1,nc_2),nrow=2))$p.value
    p_df<-rbind(p_df,wt_nc)
  }
  colnames(p_df)<-c("NC-CS")
  rownames(p_df)<-ct
  write.csv(p_df,file = paste0(name,".celltype.fisher.p.csv"))
  
  value<-c()
  sample<-c()
  celltype<-c()
  for (cell in ct) {
    print(cell)
    
    wt<-filter(data,group=="NC")%>%nrow()
    wt_1<-filter(data,celltype==cell&group=="NC")%>%nrow()/wt
    
    nc<-filter(data,group=="CS")%>%nrow()
    nc_1<-filter(data,celltype==cell&group=="CS")%>%nrow()/nc
    
    value<-c(value,wt_1,nc_1)
    sample<-c(sample,"NC","CS")
    celltype<-c(celltype,rep(cell,2))
  }
  
  plotdata<-data.frame(value=value,
                       sample=factor(sample,levels = c("NC","CS")),celltype=factor(celltype,levels = ct))
  
  p<-ggplot(plotdata,aes(x=celltype,y=value,fill=sample))+
    geom_bar(stat = 'identity',position="dodge",width = 0.7)+
    scale_fill_manual(values=brewer.pal(12,"Set3"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = -45,hjust = 0))+
    ylab("Fraction")
  ggsave(p,filename = paste0(name,".celltype.cell-sample-fisher.pdf"),width = plot.width,height = 4)
  
  ## 4 sample
  value<-c()
  sample<-c()
  celltype<-c()
  for (cell in ct) {
    print(cell)
    
    WT_F13<-filter(data,sample=="L-NC-F12")%>%nrow()
    WT_F13<-filter(data,celltype==cell&sample=="L-NC-F12")%>%nrow()/WT_F13
    WT_F14<-filter(data,sample=="L-NC-F17")%>%nrow()
    WT_F14<-filter(data,celltype==cell&sample=="L-NC-F17")%>%nrow()/WT_F14
    
    NC_F12<-filter(data,sample=="L-CS-F15")%>%nrow()
    NC_F12<-filter(data,celltype==cell&sample=="L-CS-F15")%>%nrow()/NC_F12
    NC_F17<-filter(data,sample=="L-CS-F19")%>%nrow()
    NC_F17<-filter(data,celltype==cell&sample=="L-CS-F19")%>%nrow()/NC_F17
    
    value<-c(value,WT_F13,WT_F14,NC_F12,NC_F17)
    sample<-c(sample,"NC-F12","NC-F17","CS-F15","CS-F19")
    celltype<-c(celltype,rep(cell,4))
  }
  
  plotdata<-data.frame(value=value,
                       sample=factor(sample,levels = c("NC-F12","NC-F17","CS-F15","CS-F19")),celltype=factor(celltype,levels = ct))
  
  p<-ggplot(plotdata,aes(x=celltype,y=value,fill=sample))+
    geom_bar(stat = 'identity',position="dodge",width = 0.7,color ="black")+
    scale_fill_manual(values=c( "#CCFF99","#99CC66","#FF9999","#FF6666"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = -45,hjust = 0))+
    ylab("Fraction")
  
  ggsave(p,filename = paste0(name,".celltype.cell-4sample-fisher.pdf"),width = plot.width,height = 4)
  
}


#### Fig3b ####################################
data.bk<-readRDS('Lung.cs-nc.integrated.rds')

data.bk$celltype_use<-as.character(data.bk$celltype)
data.bk@meta.data[data.bk@meta.data$celltype=="T cells",]$celltype_use<-"T_ILC_NKT"
data.bk$celltype_use<-factor(data.bk$celltype_use,levels = c("B cells","NK cells","T_ILC_NKT","Neutrophils","Basophils","MoMacDc","Fibroblasts","Endothelial cells","Epithelial cells"))

# celltype
p<-DimPlot(data.bk,reduction = 'umap',group.by = "celltype_use",repel=T, label=T, label.size=4 ,cols = cols1,)+
  theme_blank()+ggtitle("Main")
ggsave(p,filename = "Fig3b.pdf",width = 6.3,height = 5)

#### Fig3c ####################################
# plot percent barplot
data.bk$celltype<-data.bk$celltype_use
plotBarLine(data.bk,"Fig.3c",plot2sam.width=4,plot4sam.width=5,plotline.width=5)
plotFisherBar(data.bk,"Fig.3c",plot.width = 8)

#### Fig3d ####################################
# marker
{
  allfeat <- c(
    #celltype
    "Ptprc", #CD45
    "Cd53",  #CD53\
    #B
    "Cd79a", #CD79 b
    "Ms4a1", #CD20 B
    #ILC
    #"Itgae",#CD103
    #"Ncam1",#cd56
    # "Il4",
    "Il5",
    # "Il9",
    # "Il17",
    #"Cxcr3",#cs183
    # "Ccl3",
    #NK
    "Ncr1" ,  #CD335 NK
    #NK &ILC
    "Klrd1", #CD94 NK
    "Klre1", #CD94 NK
    
    "Il2rb",#CD122
    #T NKT
    "Cd226",#CD226
    #T
    "Cd3d",#CD3
    "Cd3e",#CD3
    "Cd3g",#CD3
    "Cd4",
    "Cd8a",
    #Gran
    "Cd24a",#CD24
    
    "Csf3r",#CD114
    
    "Retnlg",
    "S100a8",
    # "S100a9",
    "Cxcr2",  #cd128
    "Il3ra",#CD123 Baso
    #Bas
    "Il13",
    "Il4",
    "Mcpt8",
    "Ccl3",
    "Gata2",
    #DC
    
    
    #"Lamp3",#CD208
    "Cd207",
    "Cd209a",
    "Itgax",#CD11c
    "Cd86",#CD86h
    #M
    "Adgre1",#F4/80
    "Cd68",  #CD68 M
    # "Bst2",#CD317
    # "Cxcr3",
    
    # "Entpd1",#CD39 M
    "Mrc1",   #CD206 M
    "Csf1r",  #CD115 M
    "Fcgr1" ,  #CD64 M
    "Fcgr3" ,  #CD16 M
    # "Lilra5",#Cd85
    
    #Mono
    "Itgal",#CD11A
    # "Cd14",
    #Fibr
    # "Il1r1",#CD121a
    #"Lrp1",#CD91
    "Pdgfra","Pdgfrb",
    'Col1a1',
    'Col1a2',
    'Igf2',
    #stro
    # "Cd248",
    # "Thy1",#CD90
    # "Jag1",#CD339
    # "Bst1",#cs157
    # "Stag1",
    # "Stag2",
    # "Stag3",
    # "Stim1",
    # "Stim2",
    #endo
    # "Cd34",#
    "Pecam1",#Cd31
    "Vegfa",# Vegf
    "Cdh5",#CD144
    "Eng",#CD105
    # "Erg",
    "Cldn5",
    #ALVEOLAR
    "Sftpa1", "Sftpb", "Ager",
    #"Aqp4",
    #EPI
    "Epcam",#CD326
    #"Krt14",
    #"Krt19",
    # "Krt7","Krt8", #cytokeratin
    #CC
    "PYMT"
  )
} 
# plot
p<-DotPlot(data.bk,features = allfeat,group.by = "celltype_use",cols = c("grey","#E41A1C"),assay = "RNA")+
  theme(axis.text.x =  element_text(angle = 90, hjust = 1,vjust = 0.5))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave(p,filename = "Fig.3d.pdf",width = 9,height = 4)


