library(SeuratDisk)
library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)
library(patchwork)
library(cowplot)
library(DoubletFinder)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(RColorBrewer)


#### process ########################################################################
# load data
ref<-readRDS("../ref.rds")
nc_cs<-readRDS("NC-CS.F12-17-19-15.rds")

# remove doublet cell
F12<-nc_cs$F12
F12$doubletFinder<-F12$DF.classifications_0.25_0.19_4371
F12$pANN_0.25_0.19_4371<-NULL
F12$DF.classifications_0.25_0.19_4371<-NULL
F12<-subset(F12,doubletFinder=="Singlet")

F17<-nc_cs$F17
F17$doubletFinder<-F17$DF.classifications_0.25_0.005_2864
F17$pANN_0.25_0.005_2864<-NULL
F17$DF.classifications_0.25_0.005_2864<-NULL
F17<-subset(F17,doubletFinder=="Singlet")

F15<-nc_cs$F15
F15$doubletFinder<-F15$DF.classifications_0.25_0.13_1517
F15$pANN_0.25_0.13_1517<-NULL
F15$DF.classifications_0.25_0.13_1517<-NULL
F15<-subset(F15,doubletFinder=="Singlet")

F19<-nc_cs$F19
F19$doubletFinder<-F19$DF.classifications_0.25_0.27_1458
F19$pANN_0.25_0.27_1458<-NULL
F19$DF.classifications_0.25_0.27_1458<-NULL
F19<-subset(F19,doubletFinder=="Singlet")

## merge
data<-merge(F12,y=c(F17,F15,F19))
remove(nc_cs,F12,F17,F15,F19)

## SCTransform but no intergrate
data.sct<-SCTransform(data,vars.to.regress = c("percent.mt"),verbose = TRUE,conserve.memory = TRUE)

data.sct <- RunPCA(data.sct, verbose = FALSE)
ElbowPlot(data.sct,ndims = 50)
data.sct <- RunUMAP(data.sct, features = VariableFeatures(data.sct))
data.sct <- FindNeighbors(data.sct, dims = 1:30, verbose = FALSE)
data.sct <- FindClusters(data.sct, verbose = FALSE,resolution=0.5)

sample<-data.sct$orig.ident
sample[sample=="F12"]<-"L-NC-F12"
sample[sample=="F17"]<-"L-NC-F17"
sample[sample=="F15"]<-"L-CS-F15"
sample[sample=="F19"]<-"L-CS-F19"
sample<-factor(sample,levels = c("L-NC-F12", "L-NC-F17", "L-CS-F15", "L-CS-F19"))
data.sct$sample<-sample


p2<-DimPlot(data.sct,reduction = 'umap',group.by = 'sample')
p1<-DimPlot(data.sct,reduction = 'umap',group.by = 'seurat_clusters')
ggsave(plot_grid(p2,p1,ncol =2),filename = 'cs-nc.main.pdf',width = 13,height = 6)

data<-data.sct
remove(data.sct)

## anno by SingleR ImmGen
ref.se<-ref$ImmGen
matrix <- GetAssayData(data, slot = 'data')
singler.cluster <-SingleR(matrix,ref = ref.se,
                          labels = ref.se$label.main,
                          clusters = data@meta.data$seurat_clusters)
celltype = data.frame(ClusterID = rownames(singler.cluster),celltype = singler.cluster$labels,stringsAsFactors = F)
## load to meta.data
data@meta.data$celltype_ImmGen = 'NA'
for (i in 1:nrow(celltype)) {
  data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype_ImmGen'] <- celltype$celltype[i]
}
## plot
p<-DimPlot(data,reduction = 'umap',group.by = "celltype_ImmGen",repel=T, label=T, label.size=4)
ggsave(p,filename = "cs-nc.celltype_IMMGEN.pdf",width = 9,height = 8)

## anno by SingleR MRG
ref.se<-ref$MRG
matrix <- GetAssayData(data, slot = 'data')
singler.cluster <-SingleR(matrix,ref = ref.se,
                          labels = ref.se$label.main,
                          clusters = data@meta.data$seurat_clusters)
celltype = data.frame(ClusterID = rownames(singler.cluster),celltype = singler.cluster$labels,stringsAsFactors = F)
## load to meta.data
data@meta.data$celltype_MRD = 'NA'
for (i in 1:nrow(celltype)) {
  data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype_MRD'] <- celltype$celltype[i]
}
# plot
p<-DimPlot(data,reduction = 'umap',group.by = "celltype_MRD",repel=T, label=T, label.size=4)
# ggsave(p,filename = "cs-nc.celltype_MRD.pdf",width = 9,height = 8)

## recheck doublet cell
library(scran)
mat<-GetAssayData(data,slot = 'counts')
out<-doubletCells(mat)
data$doubletCellScore<-log10(out+1)

## Alveolar marker
feat<-c("Sftpa1", "Sftpb", "Ager", "Aqp4")
VlnPlot(data,features = feat,pt.size = 0,ncol = 2)
FeaturePlot(data,features = feat,ncol=2)
DimPlot(data,reduction = 'umap',group.by = 'seurat_clusters',label = T,label.size = 5,repel = T)

## celltype assignment
data$celltype<-data$celltype_ImmGen
data@meta.data[data@meta.data$seurat_clusters==6,]$celltype<-"Alveolar"
data@meta.data[data@meta.data$seurat_clusters==10,]$celltype<-"Alveolar"

# plot
p<-DimPlot(data,reduction = 'umap',group.by = 'celltype')
ggsave(p,filename = 'cs-nc.celltype.pdf',width = 7,height = 6)

## check pymt expression
FeaturePlot(data,"transgene-PyMT")
VlnPlot(data,"transgene-PyMT")
PYMT<-data@assays$SCT@data["transgene-PyMT",]
data$PYMT<-PYMT

# assignment cancer cells
for (i in c(1:nrow(data@meta.data))) {
  if(data@meta.data$PYMT[i]!=0){
    data@meta.data$celltype[i]<-"Cancer cells"
  }
}
saveRDS(data,'cs-nc.rds')

#### neu ##############################################################################
data.Neu<-subset(data,celltype=="Neutrophils")
DimPlot(data.Neu)

# re SCT
data.Neu<-SCTransform(data.Neu,vars.to.regress = c("percent.mt"),verbose = TRUE,conserve.memory = TRUE)
data.Neu <- RunPCA(data.Neu, verbose = FALSE)
ElbowPlot(data.Neu,ndims = 50)
data.Neu <- RunUMAP(data.Neu, features = VariableFeatures(data.Neu))
data.Neu <- FindNeighbors(data.Neu, dims = 1:20, verbose = FALSE)
data.Neu<- FindClusters(data.Neu, verbose = FALSE,resolution=0.1)
DimPlot(data.Neu,label = T)


data<-data.Neu

{
  library(clustree)
  data2 <- FindClusters(
    object = data,
    resolution = c(seq(.01,0.5,.02))
  )
  clustree(data2@meta.data, prefix = "SCT_snn_res.")
}

ref.se<-ref$ImmGen
matrix <- GetAssayData(data, slot = 'data')
singler.cluster <-SingleR(matrix,ref = ref.se,
                          labels = ref.se$label.fine,
                          clusters = data@meta.data$seurat_clusters)
celltype = data.frame(ClusterID = rownames(singler.cluster),celltype = singler.cluster$labels,stringsAsFactors = F)

## load to meta.data
data@meta.data$Neu_IMMGEN = 'NA'
for (i in 1:nrow(celltype)) {
  data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'Neu_IMMGEN'] <- paste0("c",i-1,"-",celltype$celltype[i])
}
## plot
DimPlot(data,group.by = "Neu_IMMGEN")

#### new cluster by clustree ###########################################################

data$new_cluster<-""
data@meta.data[data@meta.data$seurat_clusters==0,]$new_cluster<-0
data@meta.data[data@meta.data$seurat_clusters==1,]$new_cluster<-1
data@meta.data[data@meta.data$seurat_clusters==2,]$new_cluster<-1
data@meta.data[data@meta.data$seurat_clusters==3,]$new_cluster<-2
data@meta.data[data@meta.data$seurat_clusters==4,]$new_cluster<-3
data$seurat_clusters<-factor(data$new_cluster)
data@active.ident<-data$seurat_clusters

markers<-FindAllMarkers(data,only.pos = TRUE,min.pct = 0.5,test.use = "wilcox",downsample=100)
deg <- markers %>% filter(p_val_adj<0.05&avg_log2FC>0.8) %>% group_by(cluster)
# DoHeatmap(subset(data,downsample=50), features = deg$gene)

DEG<-deg%>%filter(!grepl('^Rp[sl]',gene)) %>%filter(!grepl('^Gm[1-9]',gene)) %>%filter(!grepl('^Hb[ab]-',gene))%>%filter(!grepl('^mt-',gene))%>%group_by(cluster)
top10 <- DEG %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
# DoHeatmap(subset(data,downsample=100), features = top10$gene)

sub<-subset(data,downsample=100)
mat<-GetAssayData(sub)
cluster_info<-factor(sort(sub$seurat_clusters,decreasing = F))
mat<-as.matrix(mat[top10$gene, names(cluster_info)])
top_anno <- HeatmapAnnotation(cluster = anno_block( gp = gpar(fill = RColorBrewer::brewer.pal(6,"Set3")),# 设置填充色
                                                    labels_gp = gpar(cex = 0.1, col = "white")))

mat2<-apply(mat,1,function(x) return(((x-mean(x)))/sd(x)))
Heatmap(t(mat2),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        
        column_split = cluster_info,
        top_annotation = top_anno)

data$celltype<-""
data@meta.data[data@meta.data$seurat_clusters==0,]$celltype<-"c0: Ccl4/Il1r2"
data@meta.data[data@meta.data$seurat_clusters==1,]$celltype<-"c1: Fos/Btg2"
data@meta.data[data@meta.data$seurat_clusters==2,]$celltype<-"c2: Tmsb10/Rack1"
data@meta.data[data@meta.data$seurat_clusters==3,]$celltype<-"c3: NGP/Camp"

p<-DimPlot(data,group.by ="celltype",label = F,label.size = 5,repel = T ,pt.size = 0.8)
ggsave(p,filename = "Neu_newcluster.celltype.pdf",width = 6.8,height = 5)

data.Neu<-data
saveRDS(data.Neu,"cs-nc.Neu_newcluster.rds")


## new name 
data$celltype_new<-data$celltype
data$celltype<-""
data@meta.data[data@meta.data$celltype_new=="c0: Ccl4/Il1r2",]$celltype<-"N1"
data@meta.data[data@meta.data$celltype_new=="c1: Fos/Btg2",]$celltype<-"N2"
data@meta.data[data@meta.data$celltype_new=="c2: Tmsb10/Rack1",]$celltype<-"N3"
data@meta.data[data@meta.data$celltype_new=="c3: NGP/Camp",]$celltype<-"N4"

data$celltype<-factor(data$celltype)

data$group=factor(data$sam,levels = c("NC","CS"))
data$sample<-factor(data$sample,levels = c("L-NC-F12", "L-NC-F17", "L-CS-F15" ,"L-CS-F19"))

### for scvelo  cellrank
write.csv(data@meta.data,file = "Neu.metadata.csv")
write.csv(data@reductions$umap@cell.embeddings,file = "Neu.umap.csv")
write.csv(data@reductions$pca@cell.embeddings,file = "Neu.pca.csv")
data@active.assay<-"RNA"
SaveLoom(data,filename = "Neu.loom")


# barcode switch mapping data.bk
data@meta.data$names<-rownames(data@meta.data)

lis<-data@meta.data[data@meta.data$orig.ident=='F15',]$names
lis.change<-gsub("_3","_4",lis)
data@meta.data[data@meta.data$orig.ident=='F15',]$names<-lis.change

lis<-data@meta.data[data@meta.data$orig.ident=='F19',]$names
lis.change<-gsub("_4","_3",lis)
data@meta.data[data@meta.data$orig.ident=='F19',]$names<-lis.change

saveRDS(data,"../Neu.rds")

#### Fig5a ###########################################
# Neu NC CS

data.Neu<-readRDS("Neu.rds")

p<-DimPlot(data.Neu,reduction = 'umap',group.by = "celltype",
           pt.size = 0.8,repel=T, label=T, label.size=4 ,
           cols = RColorBrewer::brewer.pal(4,"Set2"),split.by = "group")+
  theme_blank()+ggtitle("Celltype/group")
ggsave(p,filename = "Fig5a.pdf",width = 7.5,height = 3.7)

#### Fig5b ###########################################
plotBarLine(data.Neu,"Fig5b",plot2sam.width=4,plot4sam.width=5,plotline.width=5)
plotFisherBar(data.Neu,"Fig5b",plot.width = 8)

#### Fig5c ###########################################
data<-data.Neu
## deg
markers<-FindAllMarkers(data,only.pos = TRUE,min.pct = 0.5,test.use = "wilcox",downsample=100)
deg <- markers %>% filter(p_val_adj<0.05&avg_log2FC>0.8) %>% group_by(cluster)
write.csv(deg,file = "PLOT/Neu.deg.csv")

deg<-read.csv("Neu.deg.csv")
## filter gene
DEG<-deg%>%filter(!grepl('^Rp[sl]',gene)) %>%filter(!grepl('^Gm[1-9]',gene)) %>%filter(!grepl('^Hb[ab]-',gene))%>%filter(!grepl('^mt-',gene))%>%group_by(cluster)
top10 <- DEG %>% group_by(cluster) %>% top_n(n = 10, wt = p_val_adj)

## get expr
sub<-subset(data,downsample=2000)
mat<-GetAssayData(sub)
cluster_info<-factor(sort(sub$celltype,decreasing = F))
mat<-as.matrix(mat[top10$gene, names(cluster_info)])

## anno of heatmap
top_anno <- HeatmapAnnotation(cluster = anno_block( gp = gpar(fill = RColorBrewer::brewer.pal(4,"Set2")),# 设置填充色
                                                    labels_gp = gpar(cex = 0.1, col = "white")))
## add interested gene marker for plot
N1<-top10[top10$cluster==0,]$gene
N1<-c("Ccl3","Ccl4","Cxcl2","Il1r2",N1)%>%unique()
N2<-top10[top10$cluster==1,]$gene
add<-c("Fcgr3","Cd244a","Lyn","Cxcr2","Grk2","Igf1r")
N2<-c(add,N2)%>%unique()
N1N2<-c('Ccl3','Ccl4','Cxcl2','Il1r2',"Nfkbia",'Slfn2','Cebpb','Srgn','Hcar2','Fth1','Tnfaip3','Gadd45b','S100a9','Mpp7','Cd14','Arg1','Arg2','Cxcr2','Fos','Btg2','Csf3r','Zfp36','Fcgr3','Il1b','Txnip','dusp1','Sftpc','Junb','Ier2','Lst1','Tnfaip2','Cebpd')
N3<-top10[top10$cluster==2,]$gene[1:10]
N4<-top10[top10$cluster==3,]$gene[1:10]
gene<-c(N1N2,N3,N4)

# remove cd3
gene<-gene[gene!="Cd3g"]
gene<-gene[gene!="Cd3d"]

## downsample for expression
sub<-subset(data,downsample=100)
mat<-GetAssayData(sub)
cluster_info<-factor(sort(sub$celltype,decreasing = F))
intersect(gene,rownames(mat))

## expr mat
mat<-as.matrix(mat[intersect(gene,rownames(mat)), names(cluster_info)])
top_anno <- HeatmapAnnotation(cluster = anno_block( gp = gpar(fill = RColorBrewer::brewer.pal(4,"Set2")),# 设置填充色
                                                    labels_gp = gpar(cex = 0.1, col = "white")))

## scale
mat2<-apply(mat,1,function(x) return(((x-mean(x)))/sd(x)))
pdf("Fig5c.pdf",width = 7,height = 12)
Heatmap(t(mat2),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = cluster_info,
        top_annotation = top_anno)
dev.off()

#### Fig5d ###########################################
p<-VlnPlot(data.Neu,features = c("Ccl3","Ccl4","Il1r2","Cxcl2","Cebpb"),ncol = 5,
           group.by = "group",assay = "RNA",pt.size = 0)+ylab("Group")+xlab("")
ggsave(p,filename ="Fig5d.pdf",width = 7,height = 4 )


#### Fig5f ###########################################

# Neu marker expression

p<-DotPlot(data.Neu,features = c("Ccl3","Ccl4","Il1r2","Cxcl2","Cebpb"),
           group.by = "group",scale = T,
           cols = c("grey","#E41A1C"),assay = "RNA")+
  ggtitle("Neutrophils")+ylab("Group")+xlab("")
ggsave(p,filename ="Fig5f.pdf",width = 5.5,height = 2 )

#### Fig5g ###########################################
# table prepare
exp<-read.table('data/GSE109467/gene-TPM-matrix.txt')
colnames(exp)<-c("GMP","preNeutrophil","Immature Neutrophil","Mature Neutrophil","Blood Neutrophil 1","Blood Neutrophil 2")


sub<-subset(data,downsample=100)
mat<-GetAssayData(sub)
cluster_info<-factor(sort(sub$celltype,decreasing = F))
mat<-as.matrix(mat[top10$gene, names(cluster_info)])
top_anno <- HeatmapAnnotation(cluster = anno_block( gp = gpar(fill = RColorBrewer::brewer.pal(4,"Set2")),# 设置填充色
                                                    labels_gp = gpar(cex = 0.1, col = "white")))
mat2<-apply(mat,1,function(x) return(((x-mean(x)))/sd(x)))
Heatmap(t(mat2),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = cluster_info,
        top_annotation = top_anno)
gene<-top10%>%select(cluster,gene)


ref<-t(mat2)
ref.y<-gene$cluster%>%as.numeric()

X<-log2(exp[gene$gene,]+1)
X<-apply(X,1,function(x) return(((x-mean(x)))/sd(x)))%>%t()
X<-X%>%na.omit()
ref<-ref[rownames(X),]

Mat<-log2(exp[gene$gene,]+1)
Mat<-apply(Mat,1,function(x) return(((x-mean(x)))/sd(x)))%>%t()
Mat<-Mat%>%na.omit()
write.table(Mat,"Mat.txt",sep = "\t")

ref<-ref[rownames(Mat),]
write.table(ref,"ref.txt",sep = "\t")
write.table(ref.y,"ref.y.txt",sep = "\t")
Mat<-Mat%>%na.omit()
write.table(Mat,"Mat.txt",sep = "\t",quote = F)
ref<-ref[rownames(Mat),]
write.table(ref,"ref.txt",sep = "\t",quote = F)
write.table(ref.y,"ref.y.txt",sep = "\t",quote = F)

## for CIBERSORTx https://cibersortx.stanford.edu/
## get CIBERSORTx_Job7_Results.txt
## plot by Excel in CIBERSORTx_Job7_Results.noblood.xlsx



#### Fig5f ################################################
## for scvelo  cellrank
write.csv(data@meta.data,file = "Neu.metadata.csv")
write.csv(data@reductions$umap@cell.embeddings,file = "Neu.umap.csv")
write.csv(data@reductions$pca@cell.embeddings,file = "Neu.pca.csv")
data@active.assay<-"RNA"
SaveLoom(data,filename = "Neu.loom")

## shall code for velocyto
# velocyto run10x /home/DATA/zhengjie/DATA/yuzuoren/F12_PYMT2 /home/old_home/public_data/GENOME_REF/mm10-ranger-pymt/genes/genes.gtf -@ 20
# velocyto run10x /home/DATA/zhengjie/DATA/yuzuoren/F17_PYMT2 /home/old_home/public_data/GENOME_REF/mm10-ranger-pymt/genes/genes.gtf -@ 20
# velocyto run10x /home/DATA/zhengjie/DATA/yuzuoren/F15_PYMT2 /home/old_home/public_data/GENOME_REF/mm10-ranger-pymt/genes/genes.gtf -@ 20
# velocyto run10x /home/DATA/zhengjie/DATA/yuzuoren/F19_PYMT2 /home/old_home/public_data/GENOME_REF/mm10-ranger-pymt/genes/genes.gtf -@ 20

## analysis by python notebook: script/scVelo.ipynb
## analysis by python notebook: script/cellrank.ipynb
# note: the name of class N1 and N2 maybe changed by hand


