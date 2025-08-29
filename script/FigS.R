#### Fig S6 #################################################################

### S6A
## process T cell cluster
name<-"T cells"
data<-subset(data.bk,celltype==name&celltype_modify!="B cells")
DimPlot(data,group.by = "seurat_clusters")
p1<-DimPlot(data,group.by = "group")
p2<-DimPlot(data,group.by = "sample")
ggsave(plot_grid(p1,p2),filename = paste0(name,".sample.pdf"),width = 14,height = 6)


data<-FindVariableFeatures(data,nfeatures = 2000)
data <- ScaleData(data,vars.to.regress = "percent.mt")
data <- RunPCA(data, verbose = FALSE)
ElbowPlot(data,ndims = 50)+ylim(1,10)+xlim(1,50)

data <- RunUMAP(data, dims = 1:20)
data <- FindNeighbors(data, dims = 1:20, verbose = FALSE)
data <- FindClusters(data, verbose = FALSE,resolution=0.5)
DimPlot(data,label = T)

{library(clustree)
  data2 <- FindClusters(
    object = data,
    resolution = c(seq(.01,0.8,.02))
  )
  clustree(data2@meta.data, prefix = "integrated_snn_res.")
}


data<-annoBySingler(data,ref.se=ref$ImmGen,fineOrMain="fine",check.cluster = FALSE,addClusterToAnno = TRUE)

data$celltype<-factor(data$celltype,levels = c("c0_T cells (T.4NVE44-49D-11A-)",  "c1_T cells (T.4NVE44-49D-11A-)",            "c2_T cells (T.Tregs)"   ,        
                                               "c3_T cells (T.8Nve)"    ,         "c4_NKT (NKT.4-)"  ,               "c5_T cells (T.4NVE44-49D-11A-)" ,
                                               "c6_T cells (T.8EFF.OT1LISO)" ,    "c7_Tgd (Tgd.mat.VG2+)" ,          "c8_NKT (NKT.4-)"     ,           
                                               "c9_ILC (ILC2)","c10_T cells (T.4NVE44-49D-11A-)",
                                               "c11_T cells (T.4NVE44-49D-11A-)", "c12_T cells (T.8EFF.OT1LISO)" ,   "c13_T cells (T.CD4CONTROL)"     ,
                                               "c14_T cells (T.4NVE44-49D-11A-)" ,"c15_T cells (T.4Mem)"  ))
DimPlot(data,reduction = 'umap',group.by = "celltype",pt.size = 1,repel=T, label=T, label.size=4,
        cols = c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")))
VlnPlot(data,features = "Cd3d",group.by = "celltype")+NoLegend()


## modify a B cell cluster
# { 
#   feat<-c("Ptprc","Cd52","Cd79a","Cd19","Ms4a1","Cd3d","Cd3g","Cd3e","Cd4","Cd8a")
#   DotPlot(data,group.by = "celltype",features = feat,assay = "RNA")+coord_flip()+NoLegend()+theme(axis.text.x =  element_text(angle = -45, hjust = 0,vjust = 0.5))
#   
#   data.bk$celltype_modify<-"NA"
#   data.bk$celltype_fine<-"NA"
#   
#   data.bk@meta.data[data@meta.data[data@meta.data$celltype=="B cells (B.GC)",]%>%rownames(),]$celltype_modify<-"B cells"
#   data.bk@meta.data[data@meta.data[data@meta.data$celltype=="B cells (B.GC)",]%>%rownames(),]$celltype_fine<-"B cells (B.GC)"
#   DimPlot(data.bk,group.by = "celltype_modify")
#   data<-subset(data,celltype!="B cells (B.GC)")
# }

## rename
data$celltype_fine<-data$celltype
data$celltype<-""
data@meta.data[data@meta.data$celltype_fine=="ILC (ILC2)",]$celltype<-"ILC"
data@meta.data[data@meta.data$celltype_fine=="NKT (NKT.4-)",]$celltype<-"NKT"
data@meta.data[data@meta.data$celltype_fine=="T cells (T.4Mem)",]$celltype<-"Tmem"
data@meta.data[data@meta.data$celltype_fine=="T cells (T.Tregs)",]$celltype<-"Treg"
data@meta.data[data@meta.data$celltype_fine=="T cells (T.8Nve)",]$celltype<-"T.8Nve"
data@meta.data[data@meta.data$celltype_fine=="Tgd (Tgd.mat.VG2+)",]$celltype<-"Tgd"
data@meta.data[data@meta.data$celltype_fine=="T cells (T.8EFF.OT1LISO)",]$celltype<-"T.8Eff"
data@meta.data[data@meta.data$celltype_fine=="T cells (T.CD4CONTROL)",]$celltype<-"T.4Nve"
data@meta.data[data@meta.data$celltype_fine=="T cells (T.4NVE44-49D-11A-)",]$celltype<-"T.4Nve"
data$celltype<-factor(data$celltype,levels = c('T.4Nve',"T.8Nve","T.8Eff","Tmem","Treg","Tgd","NKT","ILC"))

### S6A 
p1<-DimPlot(data,group.by = "celltype",label = T,repel = T,cols = cols2,pt.size = 0.8)+theme_blank()
ggsave(p1,filename = "S6A.pdf",width = 4.8,height = 4)

### S6B 
## function for plot
plotBarLine<-function(data,name,plot2sam.width=3.5,plot4sam.width=4,plotline.width=4){
  cols<-c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5","#BC80BD", "#CCEBC5", "#FFED6F")
  
  
  data<-data@meta.data
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
  
  df_4sample<-data.frame(value=c(as.numeric(df2[,1]),as.numeric(df2[,2]),as.numeric(df2[,3]),as.numeric(df2[,4])),
                         type=factor(c(rownames(df2),rownames(df2),rownames(df2),rownames(df2)),levels = ct),
                         sample=factor(
                           c(rep('L-NC-F12',nrow(df2)),
                             rep('L-NC-F17',nrow(df2)),
                             rep('L-CS-F15',nrow(df2)),
                             rep('L-CS-F19',nrow(df2))),
                           levels = c('L-NC-F12', 'L-NC-F17', 'L-CS-F15', 'L-CS-F19')))
  
  p1<-ggplot(df_4sample,aes(x=sample,y=value,fill=type))+
    geom_bar(stat = 'identity')+
    theme_bw() +
    scale_fill_manual(values=cols) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10,angle = 90, hjust = 1,vjust = 0.5),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black")
    )+
    ylab("Fraction")+
    xlab("Sample")
  ggsave(p1,filename = paste0(name,".cellprop-bar.4.pdf"),width = plot4sam.width,height = 5)
 ######################################################################
  group<-levels(data$group)
  # ct<-unique(data$celltype)
  ct<-levels(data$celltype)
  df<-c()
  for (c in ct) {
    out<-c()
    for (s in group) {
      num<-filter(data,celltype==c&group==s)%>% nrow()
      out<-c(out,num)
    }
    #out<-out/sum(out)
    df<-rbind(df,out)
  }
  colnames(df)<-group
  rownames(df)<-ct
  
  df2<-apply(df,2,function(x){return(x/sum(x))})
  
  df_2sample<-data.frame(value=c(as.numeric(df2[,1]),as.numeric(df2[,2])),
                         type=factor(c(rownames(df2),rownames(df2)),levels = ct),
                         group=factor(
                           c(rep('NC',nrow(df2)),
                             rep('CS',nrow(df2))
                           ),
                           levels = c('NC', 'CS')))
 #########################33
  # L_NC<-(df2[,1]+df2[,2])/2
  # L_CS<-(df2[,3]+df2[,4])/2
  # 
  # df_2sample<-data.frame(value=c(as.numeric(L_NC),as.numeric(L_CS)),
  #                type=factor(c(names(L_NC),names(L_CS)),levels = ct),
  #                group=factor(c(rep('NC',length(L_NC)),rep('CS',length(L_CS))),
  #                              levels = c("NC","CS")))
  library(RColorBrewer)
  
  p1<-ggplot(df_2sample,aes(x=group,y=value,fill=type))+
    geom_bar(stat = 'identity')+
    theme_bw() +
    scale_fill_manual(values=cols) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black")
    )+
    ylab("Fraction")+
    xlab("Group")
  ggsave(p1,filename = paste0(name,".cellprop-bar.pdf"),width = plot2sam.width,height = 5)
  
  
  p2<-ggplot(df,aes(x=group,y=value,group=type,shape=type))+
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
    scale_color_manual(values =cols)+
    ylab("Fraction")+
    xlab("Group")
  ggsave(p2,filename = paste0(name,".cellprop-line.pdf"),width = plotline.width,height = 5)
}
plotFisherBar<-function(data,name,plot.width=8){
  cols<-c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5","#BC80BD", "#CCEBC5", "#FFED6F")
  
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
    scale_fill_manual(values=c("#FBB4AE", "#B3CDE3"))+
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
    scale_fill_manual(values=c( "#FBB4AE", "#FBB4AE","#B3CDE3","#B3CDE3"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = -45,hjust = 0))+
    ylab("Fraction")
  
  ggsave(p,filename = paste0(name,".celltype.cell-4sample-fisher.pdf"),width = plot.width,height = 4)
  
}
plotFisherBarInAll<-function(data,data.bk,name,plot.width=4,plot.height=4){
  cols<-c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5","#BC80BD", "#CCEBC5", "#FFED6F")
  
  data<-data@meta.data
  ct<-levels(data$celltype)
  
  value<-c()
  sample<-c()
  celltype<-c()
  for (cell in ct) {
    print(cell)
    
    # wt<-filter(data,group=="NC")%>%nrow()
    wt<-filter(data.bk@meta.data,group=="NC")%>%nrow()
    wt_1<-filter(data,celltype==cell&group=="NC")%>%nrow()/wt
    
    # nc<-filter(data,group=="CS")%>%nrow()
    nc<-filter(data.bk@meta.data,group=="CS")%>%nrow()
    nc_1<-filter(data,celltype==cell&group=="CS")%>%nrow()/nc
    
    value<-c(value,wt_1,nc_1)
    sample<-c(sample,"NC","CS")
    celltype<-c(celltype,rep(cell,2))
  }
  
  plotdata<-data.frame(value=value,
                       sample=factor(sample,levels = c("NC","CS")),celltype=factor(celltype,levels = ct))
  
  p<-ggplot(plotdata,aes(x=celltype,y=value,fill=sample))+
    geom_bar(stat = 'identity',position="dodge",width = 0.7)+
    scale_fill_manual(values=c("#FBB4AE", "#B3CDE3"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = -45,hjust = 0))+
    ylab("Fraction")
  ggsave(p,filename = paste0(name,".celltype.cell-group.pdf"),width = plot.width,height = plot.height)
}

# bar
plotBarLine(data,"S6B",plot2sam.width=3,plot4sam.width=3.5,plotline.width=3.5,cols = cols2)
plotFisherBar(data,"S6B",plot.width=8,cols = c("#1F78B4", "#B2DF8A"))
plotFisherBarInAll(data,data.bk,"S6B",plot.width=4,plot.height=4,cols=c("#1F78B4", "#B2DF8A"))


### S6C 
# marker
feat.t<-c("Cd3d","Cd3g","Cd4","Cd8a","Lef1","Sell",
          "Gzmb","Gzmk",
          'Calcrl',"Ccr7",
          "Ctla4","Foxp3",
          'Trgv2','Trdv4',
          'Klrk1',"Klrd1",
          "Gata3","Il2ra")
p<-DotPlot(data,
           # scale = F,
           features = feat.t,cols = c("white","red"),group.by = "celltype",assay = "RNA")+theme(axis.text.x =  element_text(angle = -45, hjust = 0,vjust = 0.5))
ggsave(p,filename = "S6C.pdf",width = 7,height = 3)

### S6FGH 待补充 




#### Fig S7 ##################################################################
### S7A
name<-"MoMacDc"
data<-subset(data.bk,celltype=="MoMacDc")
DimPlot(data,group.by = "seurat_clusters")

feat<-c("Lamp3",#CD208
        "Cd207",
        "Cd209a",
        "Itgax",#CD11c
        "Cd86",#CD86h
        "Irf8",#TF
        #M
        "Adgre1",#F4/80
        "Cd68",  #CD68 M
        "Bst2",#CD317
        "Cxcr3",
        
        "Entpd1",#CD39 M
        "Mrc1",   #CD206 M
        "Csf1r",  #CD115 M
        "Fcgr1" ,  #CD64 M
        "Fcgr3" ,  #CD16 M
        "Lilra5",#Cd85
        
        #Mono
        "Itgal"#CD11A
)
# DotPlot(data,features = feat)
VlnPlot(data,features = feat,flip = T,stack = T,group.by = "celltype")

## re-pipline
data<-FindVariableFeatures(data,nfeatures = 2000)
data <- ScaleData(data,vars.to.regress = "percent.mt")
data <- RunPCA(data, verbose = FALSE)
# ElbowPlot(data,ndims = 50)+ylim(1,10)+xlim(1,50)

data <- RunUMAP(data, dims = 1:20)
data <- FindNeighbors(data, dims = 1:20, verbose = FALSE)
data <- FindClusters(data, verbose = FALSE,resolution=0.1,n.iter = 20,n.start = 6)
DimPlot(data,label = T)

data<-annoBySingler(data,ref.se=ImmGenData(),fineOrMain="fine",check.cluster = FALSE,addClusterToAnno = TRUE)

DimPlot(data,group.by = "celltype", cols = c(brewer.pal(12,"Set3")))
data<-subset(data,celltype!="c4_Neutrophils (GN.ARTH)")

# remove Neutrophils re-run
data<-FindVariableFeatures(data,nfeatures = 2000)
data <- ScaleData(data,vars.to.regress = "percent.mt")
data <- RunPCA(data, verbose = FALSE)
# ElbowPlot(data,ndims = 50)+ylim(1,10)+xlim(1,50)

data <- RunUMAP(data, dims = 1:20)
data <- FindNeighbors(data, dims = 1:20, verbose = FALSE)
data <- FindClusters(data, verbose = FALSE,resolution=0.1,n.iter = 20,n.start = 6)
DimPlot(data,label = T)

## singleR anno
data<-annoBySingler(data,ref.se=ImmGenData(),fineOrMain="fine",check.cluster = FALSE,addClusterToAnno = TRUE)
DimPlot(data,group.by = "celltype", cols = c(brewer.pal(12,"Set3")))

## celltype assignment
data$celltype<-""
data@meta.data[data@meta.data$seurat_clusters==0,]$celltype<-'Mac-1'
data@meta.data[data@meta.data$seurat_clusters==1,]$celltype<-'Mono'
data@meta.data[data@meta.data$seurat_clusters==2,]$celltype<-'DC-1'
data@meta.data[data@meta.data$seurat_clusters==3,]$celltype<-'Mac-2'
data@meta.data[data@meta.data$seurat_clusters==4,]$celltype<-'DC-2'
data@meta.data[data@meta.data$seurat_clusters==5,]$celltype<-'DC-3'
data$celltype<-factor(data$celltype)

## plot
p1<-DimPlot(data,group.by = "celltype",label = T,repel = T,cols = cols2,pt.size = 0.9)+theme_blank()
ggsave(p1,filename = "S7A.pdf",width = 4.8,height = 4)

### S7B
plotBarLine(data,"S7B",plot2sam.width=3,plot4sam.width=3.5,plotline.width=3.5,cols = cols2)
plotFisherBar(data,"S7B",plot.width=4,cols = c("#1F78B4", "#B2DF8A"))

### S7C
## DEG 
deg<-FindAllMarkers(data,only.pos = T,test.use = "wilcox",min.pct = 0.5)
deg<-filter(deg,p_val_adj<0.05&avg_log2FC>1)
## DEG top15
DEG<-deg%>%group_by(cluster)%>%top_n(20)

# downsample for expr
sub<-subset(data,downsample=50)
mat<-GetAssayData(sub,slot = "scale.data")
cluster_info<-factor(sort(sub$celltype,decreasing = F))
mat<-as.matrix(mat[intersect(DEG$gene,rownames(mat)), names(cluster_info)])%>%unique.data.frame

# anno
top_anno <- HeatmapAnnotation(cluster = anno_block( gp = gpar(fill = cols),labels_gp = gpar(cex = 0.1, col = "white")))

mark_gene <-intersect(c(unique(top5$gene),highlight_gene),rownames(mat))
gene_pos <-which(rownames(mat) %in% mark_gene)
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = mark_gene))

# scale
mat2<-apply(mat,1,function(x) return(((x-mean(x)))/sd(x)))
# plot
Heatmap(t(mat2),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = cluster_info,
        top_annotation = top_anno,
        # right_annotation  = row_anno,
        # jitter=TRUE,
        show_heatmap_legend = TRUE
        #column_title = NULL
)
### S7EFGH 待补充 


#### Fig S8###################################################################

### S8A
data<-subset(data.bk,celltype=="Endothelial cells")
data@active.ident<-factor(data$group)

p<-DotPlot(data,features = c('Tjp1','Ocln','Cldn5'),assay = "RNA",scale = T,cols = c("grey","salmon"))+coord_flip()+xlab("Gene of \nvascular permeability") +ylab("group")
ggsave(p,filename = "S8A.pdf",width = 4,height =2.5)

### S8B
## GSEA enrichment
library(singleseqgset)
library(heatmap3)
library(msigdbr)

path <- msigdbr(species="Mus musculus",category="C2", subcategory = "KEGG")
path.names <- unique(path$gs_name)
path.sets <- vector("list",length=length(path.names))
names(path.sets) <- path.names
for (i in names(path.sets)) {
  path.sets[[i]] <- pull(path[path$gs_name==i,"gene_symbol"])
}
# remove_prefix
n<-6
name<-names(path.sets)
names(path.sets)<-lapply(name, function(x){substr(x,n,nchar(x))})%>%as.character()

logfc.data <- logFC(cluster.ids=as.character(data@active.ident),expr.mat=GetAssayData(data,assay = "RNA"))
gse.res <- wmw_gsea(expr.mat=GetAssayData(data,assay = "RNA"),cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=path.sets)


# filter
res.stats <- gse.res[["GSEA_statistics"]]%>%na.omit()
res.pvals <- gse.res[["GSEA_p_values"]]

res.pvals <- apply(res.pvals,2,p.adjust,method="fdr")
res.pvals <- res.pvals[rownames(res.stats ),]

p<-pheatmap::pheatmap(res.stats[c("TIGHT_JUNCTION","GAP_JUNCTION","ADHERENS_JUNCTION"),],scale = "none",
                      display_numbers = (ifelse(res.pvals[c("TIGHT_JUNCTION","GAP_JUNCTION","ADHERENS_JUNCTION"),] <0.05, "*", "")),border_color = "white",
                      cluster_cols = F,cluster_rows = F,fontsize_number = 20,fontsize_row = 10)
ggsave(p,filename = "S8B.pdf",width = 2.8,height =1.4)

### S8B 待补充
#### Fig S11 待补充 #######################################################################
#### Fig S12 待补充 ###############################################################
#### Fig S13 待补充 #########