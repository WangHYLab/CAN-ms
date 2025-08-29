#### Fig7b A ################################################

FeaturePlot(data.Neu,features = "Cebpb")
ggsave(filename ="Fig7b-a.pdf",width = 5.5,height = 2 )



#### Fig7b B ################################################

data<-readRDS("Neu.rds")

exprMat <- GetAssayData(data,slot = "counts")%>%as.matrix()
cellInfo <- data.frame(data@meta.data)

cellInfo <- cellInfo[,c("sample","group","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")

## Initialize settings
scenicOptions <- initializeScenic(org="mgi", dbDir="/home/DATA/zhengjie/ref/SCENIC", nCores=10) 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#scenicOptions <-readRDS("scenicOptions.Rds")

## co-expr
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

## Build and score the GRN 
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- initializeScenic(org="mgi", dbDir="/home/DATA/zhengjie/ref/SCENIC", nCores=1) 
runSCENIC_3_scoreCells(scenicOptions, exprMat_log,skipTsne = T) 
runSCENIC_4_aucell_binarize(scenicOptions)

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Fos"]
viewMotifs(tableSubset)

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Cebpb"]
viewMotifs(tableSubset)

tsneAUC(scenicOptions, aucType="AUC")
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

## AUC
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

anno_colors = list(
  group = c(CS="#E41A1C",NC= "#377EB8"),
  celltype = c( 'N1-Ccl4/Il1r2'="#66C2A5" ,'N2-Fos/Btg2'="#FC8D62",'N3-Tmsb10/Rack1'= "#8DA0CB", 'N4-NGP/Camp'="#E78AC3")
) ## change to N2 N1 N3 N4 by hand later

anno_col= cellInfo[,2:3]
pheatmap::pheatmap(as.matrix(auc),
                   scale = "row",
                   annotation_col = anno_col,
                   clustering_method = "ward.D2",
                   annotation_colors = anno_colors,
                   breaks=seq(-2, 2, length.out = 100),
                   treeheight_row = 10,
                   treeheight_col = 16,
                   cutree_cols = 3,
                   show_colnames = F)

# group
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,
                   treeheight_row=10,
                   treeheight_col=10, 
                   cluster_cols = F,
                   border_color=NA)


topRegulators <- reshape2::melt(apply(regulonActivity_byCellType_Scaled, 2, function(x) cbind(sort(x[x>0], decreasing=TRUE))))[c("L1","Var1", "value")]; colnames(topRegulators) <- c("CellType","Regulon", "RelativeActivity")

viewTable(topRegulators)

## binary
minPerc <- .0
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))

binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

pheatmap::pheatmap(binaryActPerc_subset, color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA)


## activity of TF
library(Seurat)

dr_coords <- Embeddings(data, reduction="umap")

tfs <- c("Fos","Jund","Atf3","Cebpb","Egr1","Klf2","Ets1")

par(mfrow=c(3,3))

AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")

