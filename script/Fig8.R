#### Fig8a ################################################
# CELLCHAT 
library(CellChat)

cellchat<-createCellChat(data.bk,assay = "integrated",group.by = "celltype_use")
{
  bk<-data.bk@meta.data[data.bk$celltype_use=="Neutrophils",]%>%rownames()
  n<-data.Neu@meta.data$names
  
  p<-intersect(n,bk)
  # p<-n
  meta<-data.Neu@meta.data
  rownames(meta)<-meta$names
  meta<-meta[p,]
  
  data.bk$celltype_merge_Neu<-"na"
  for (i in c(1:nrow(meta))) {
    data.bk@meta.data[rownames(meta[i,]),"celltype_merge_Neu"]<-as.character(meta[i,"celltype"])
  }
  
  ## have PYMT expression, so assignment to cancer cells
  data.bk@meta.data['GCCAGGTCAAGCGAGT-1_2',]$celltype_merge_Neu<-"Cancer cells"
  
  data<-subset(data.bk,celltype_merge_Neu=="Cancer cells"|celltype_merge_Neu=="N1"|celltype_merge_Neu=="N2"|celltype_merge_Neu=="N3"|celltype_merge_Neu=="N4")
  cellchat.all<-cellchat
  
  cellchat<-createCellChat(data,assay = "integrated",group.by = "celltype_merge_Neu")
  
  
  ## sub Neu celltype with all other cells
  bk<-data.bk@meta.data[data.bk$celltype_cancer=="Neutrophils",]%>%rownames()
  n<-data.Neu@meta.data$names
  
  p<-intersect(n,bk)
  # p<-n
  meta<-data.Neu@meta.data
  rownames(meta)<-meta$names
  meta<-meta[p,]
  
  data.bk$celltype_merge_Neu<-as.character(data.bk$celltype_cancer)
  for (i in c(1:nrow(meta))) {
    data.bk@meta.data[rownames(meta[i,]),"celltype_merge_Neu"]<-as.character(meta[i,"celltype"])
  }
  # data.bk@meta.data['GCCAGGTCAAGCGAGT-1_2',]$celltype_merge_Neu<-"Cancer cells"
  
  data<-subset(data.bk,celltype_merge_Neu!="Neutrophils")
  cellchat<-createCellChat(data,assay = "integrated",group.by = "celltype_merge_Neu")
  
}

CellChatDB <- CellChatDB.mouse

{ #add CCL4-CCR1
  CCL4_CCR1<-c('CCL4_CCR1','CCL','Ccl4','Ccr1','','','' ,'' ,'PMID:23233369','Secreted Signaling','Ccl4  - Ccr1')
  aa<-CellChatDB$interaction
  aa<-rbind(aa,CCL4_CCR1)
  rownames(aa)<-aa$interaction_name
  CellChatDB$interaction<-aa
  
}

CellChatDB.use <- subsetDB(CellChatDB, search =c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor") )    # use Secreted Signaling for cell-cell communication analysis

cellchat@DB <- CellChatDB.use 

showDatabaseCategory(CellChatDB)

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

# cellchat <- computeCommunProb(cellchat) # the pri code have bug, so fix
computeCommunProb_modify<-function (object, type = c("triMean", "truncatedMean", "median"), 
                                    trim = NULL, LR.use = NULL, raw.use = TRUE, population.size = FALSE, 
                                    do.fast = TRUE, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1) {
  # type <- match.arg(type)
  # FunMean <- switch(type, triMean = triMean, truncatedMean = function(x) mean(x, 
  #                                                                             trim = trim, na.rm = TRUE), median = function(x) median(x, 
  #                                                                                                                                     na.rm = TRUE))
  FunMean <- triMean
  if (raw.use) {
    data <- as.matrix(object@data.signaling)
  }
  else {
    data <- object@data.project
  }
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  else {
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  # my.sapply <- ifelse(test = future::nbrOfWorkers() == 1, yes = sapply,
  #                     no = future.apply::future_sapply)
  my.sapply <-sapply
  
  ptm = Sys.time()
  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!\n         You may need to drop unused levels using 'droplevels' function. e.g.,\n         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }
  if (all(data[1:5, ] == floor(data[1:5, ]))) {
    stop("Please check your input data matrix and ensure that you use the normalized data instead of count data!")
  }
  data.use <- data/max(data)
  nC <- ncol(data.use)
  if (do.fast) {
    data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)
    data.use.avg <- t(data.use.avg[, -1])
    colnames(data.use.avg) <- levels(group)
    
    dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
    dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
    
    dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                     data.use.avg, pairLRsig, type = "A")
    dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                     data.use.avg, pairLRsig, type = "I")
    
    dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor
    dataLavg2 <- t(replicate(nrow(dataLavg), as.numeric(table(group))/nC))
    dataRavg2 <- dataLavg2
    index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != 
                             "")
    index.antagonist <- which(!is.na(pairLRsig$antagonist) & 
                                pairLRsig$antagonist != "")
    
    Prob <- array(0, dim = c(numCluster, numCluster, nLR))
    Pval <- array(0, dim = c(numCluster, numCluster, nLR))
    set.seed(seed.use)
    permutation <- replicate(nboot, sample.int(nC, size = nC))
    data.use.avg.boot <- my.sapply(X = 1:nboot, FUN = function(nE) {
      groupboot <- group[permutation[, nE]]
      data.use.avgB <- aggregate(t(data.use), list(groupboot), 
                                 FUN = FunMean)
      data.use.avgB <- t(data.use.avgB[, -1])
      return(data.use.avgB)
    }, simplify = FALSE)
    pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())
    for (i in 1:nLR) {
      dataLR <- Matrix::crossprod(matrix(dataLavg[i, ], 
                                         nrow = 1), matrix(dataRavg[i, ], nrow = 1))
      P1 <- dataLR^n/(Kh^n + dataLR^n)
      if(any(is.na(P1))){
        Pnull = P1
        Prob[, , i] <- Pnull
        p = 1
        Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, 
                              byrow = FALSE)
      }else{
        if (sum(P1) == 0) {
          Pnull = P1
          Prob[, , i] <- Pnull
          p = 1
          Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, 
                                byrow = FALSE)
        }
        else {
          if (is.element(i, index.agonist)) {
            data.agonist <- computeExpr_agonist(data.use = data.use.avg, 
                                                pairLRsig, cofactor_input, index.agonist = i, 
                                                Kh = Kh, n = n)
            P2 <- Matrix::crossprod(matrix(data.agonist, 
                                           nrow = 1))
          }
          else {
            P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          if (is.element(i, index.antagonist)) {
            data.antagonist <- computeExpr_antagonist(data.use = data.use.avg, 
                                                      pairLRsig, cofactor_input, index.antagonist = i, 
                                                      Kh = Kh, n = n)
            P3 <- Matrix::crossprod(matrix(data.antagonist, 
                                           nrow = 1))
          }
          else {
            P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          if (population.size) {
            P4 <- Matrix::crossprod(matrix(dataLavg2[i, 
            ], nrow = 1), matrix(dataRavg2[i, ], nrow = 1))
          }
          else {
            P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          Pnull = P1 * P2 * P3 * P4
          Prob[, , i] <- Pnull
          Pnull <- as.vector(Pnull)
          Pboot <- sapply(X = 1:nboot, FUN = function(nE) {
            data.use.avgB <- data.use.avg.boot[[nE]]
            dataLavgB <- computeExpr_LR(geneL[i], data.use.avgB, 
                                        complex_input)
            dataRavgB <- computeExpr_LR(geneR[i], data.use.avgB, 
                                        complex_input)
            dataRavgB.co.A.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                              data.use.avgB, pairLRsig[i, , drop = FALSE], 
                                                              type = "A")
            dataRavgB.co.I.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                              data.use.avgB, pairLRsig[i, , drop = FALSE], 
                                                              type = "I")
            dataRavgB <- dataRavgB * dataRavgB.co.A.receptor/dataRavgB.co.I.receptor
            dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
            P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
            if (is.element(i, index.agonist)) {
              data.agonist <- computeExpr_agonist(data.use = data.use.avgB, 
                                                  pairLRsig, cofactor_input, index.agonist = i, 
                                                  Kh = Kh, n = n)
              P2.boot <- Matrix::crossprod(matrix(data.agonist, 
                                                  nrow = 1))
            }
            else {
              P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            if (is.element(i, index.antagonist)) {
              data.antagonist <- computeExpr_antagonist(data.use = data.use.avgB, 
                                                        pairLRsig, cofactor_input, index.antagonist = i, 
                                                        Kh = Kh, n = n)
              P3.boot <- Matrix::crossprod(matrix(data.antagonist, 
                                                  nrow = 1))
            }
            else {
              P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            if (population.size) {
              groupboot <- group[permutation[, nE]]
              dataLavg2B <- as.numeric(table(groupboot))/nC
              dataLavg2B <- matrix(dataLavg2B, nrow = 1)
              dataRavg2B <- dataLavg2B
              P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
            }
            else {
              P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
            }
            Pboot = P1.boot * P2.boot * P3.boot * P4.boot
            return(as.vector(Pboot))
          })
          Pboot <- matrix(unlist(Pboot), nrow = length(Pnull), 
                          ncol = nboot, byrow = FALSE)
          nReject <- rowSums(Pboot - Pnull >= 0)
          p = nReject/nboot
          Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, 
                                byrow = FALSE)
        }
        
      }
      
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
  }
  else {
    dataL <- computeExpr_LR(geneL, data.use, complex_input)
    dataR <- computeExpr_LR(geneR, data.use, complex_input)
    dataR.co.A.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                  data.use, pairLRsig, type = "A")
    dataR.co.I.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                  data.use, pairLRsig, type = "I")
    dataR <- dataR * dataR.co.A.receptor/dataR.co.I.receptor
    dataLavg <- aggregate(t(dataL), list(group), FUN = FunMean)
    dataLavg <- t(dataLavg[, -1])
    rownames(dataLavg) <- geneL
    dataRavg <- aggregate(t(dataR), list(group), FUN = FunMean)
    dataRavg <- t(dataRavg[, -1])
    rownames(dataRavg) <- geneR
    dataL.binary = (dataL > 0) * 1
    dataR.binary = (dataR > 0) * 1
    dataLavg2 <- aggregate(t(dataL.binary), list(group), 
                           FUN = sum)
    dataLavg2 <- t(dataLavg2[, -1])/nC
    dataRavg2 <- aggregate(t(dataR.binary), list(group), 
                           FUN = sum)
    dataRavg2 <- t(dataRavg2[, -1])/nC
    index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != 
                             "")
    index.antagonist <- which(!is.na(pairLRsig$antagonist) & 
                                pairLRsig$antagonist != "")
    set.seed(seed.use)
    permutation <- replicate(nboot, sample.int(nC, size = nC))
    Prob <- array(0, dim = c(numCluster, numCluster, nLR))
    Pval <- array(0, dim = c(numCluster, numCluster, nLR))
    pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())
    for (i in 1:nLR) {
      dataLR <- Matrix::crossprod(matrix(dataLavg[i, ], 
                                         nrow = 1), matrix(dataRavg[i, ], nrow = 1))
      P1 <- dataLR^n/(Kh^n + dataLR^n)
      if (sum(P1) == 0) {
        Pnull = P1
        Prob[, , i] <- Pnull
        p = 1
        Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, 
                              byrow = FALSE)
      }
      else {
        if (is.element(i, index.agonist)) {
          data.agonist <- computeExprGroup_agonist(data.use = data.use, 
                                                   pairLRsig, cofactor_input, group = group, 
                                                   index.agonist = i, Kh = Kh, FunMean = FunMean, 
                                                   n = n)
          P2 <- Matrix::crossprod(matrix(data.agonist, 
                                         nrow = 1))
        }
        else {
          P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (is.element(i, index.antagonist)) {
          data.antagonist <- computeExprGroup_antagonist(data.use = data.use, 
                                                         pairLRsig, cofactor_input, group = group, 
                                                         index.antagonist = i, Kh = Kh, FunMean = FunMean, 
                                                         n = n)
          P3 <- Matrix::crossprod(matrix(data.antagonist, 
                                         nrow = 1))
        }
        else {
          P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (population.size) {
          P4 <- Matrix::crossprod(matrix(dataLavg2[i, 
          ], nrow = 1), matrix(dataRavg2[i, ], nrow = 1))
        }
        else {
          P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        Pnull = P1 * P2 * P3 * P4
        Prob[, , i] <- Pnull
        Pnull <- as.vector(Pnull)
        dataL.i <- dataL[i, ]
        dataR.i <- dataR[i, ]
        dataL2.i <- dataL.binary[i, ]
        dataR2.i <- dataR.binary[i, ]
        Pboot <- my.sapply(X = 1:nboot, FUN = function(nE) {
          groupboot <- group[permutation[, nE]]
          dataLavgB <- aggregate(matrix(dataL.i, ncol = 1), 
                                 list(groupboot), FUN = FunMean)
          dataLavgB <- t(dataLavgB[, -1])
          dataLavgB <- matrix(dataLavgB, nrow = 1)
          dataRavgB <- aggregate(matrix(dataR.i, ncol = 1), 
                                 list(groupboot), FUN = FunMean)
          dataRavgB <- t(dataRavgB[, -1])
          dataRavgB <- matrix(dataRavgB, nrow = 1)
          dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
          P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
          if (is.element(i, index.agonist)) {
            data.agonist <- computeExprGroup_agonist(data.use = data.use, 
                                                     pairLRsig, cofactor_input, group = groupboot, 
                                                     index.agonist = i, Kh = Kh, FunMean = FunMean, 
                                                     n = n)
            P2.boot <- Matrix::crossprod(matrix(data.agonist, 
                                                nrow = 1))
          }
          else {
            P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          if (is.element(i, index.antagonist)) {
            data.antagonist <- computeExprGroup_antagonist(data.use = data.use, 
                                                           pairLRsig, cofactor_input, group = groupboot, 
                                                           index.antagonist = i, Kh = Kh, FunMean = FunMean, 
                                                           n = n)
            P3.boot <- Matrix::crossprod(matrix(data.antagonist, 
                                                nrow = 1))
          }
          else {
            P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          dataLavg2B <- by(matrix(dataL2.i, ncol = 1), 
                           groupboot, sum)/nC
          dataLavg2B <- matrix(dataLavg2B, nrow = 1)
          dataRavg2B <- by(matrix(dataR2.i, ncol = 1), 
                           groupboot, sum)/nC
          dataRavg2B <- matrix(dataRavg2B, nrow = 1)
          if (population.size) {
            P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
          }
          else {
            P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
          }
          Pboot = P1.boot * P2.boot * P3.boot * P4.boot
          return(as.vector(Pboot))
        })
        Pboot <- matrix(unlist(Pboot), nrow = length(Pnull), 
                        ncol = nboot, byrow = FALSE)
        nReject <- rowSums(Pboot - Pnull >= 0)
        p = nReject/nboot
        Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, 
                              byrow = FALSE)
      }
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
  }
  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list(prob = Prob, pval = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")
  object@options$parameter <- list(type.mean = type, trim = trim, 
                                   raw.use = raw.use, population.size = population.size, 
                                   nboot = nboot, seed.use = seed.use, Kh = Kh, n = n)
  object@net <- net
  return(object)
} 


cellchat <- computeCommunProb_modify(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


# CCI circle 
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,targets.use = 10, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength (Cancer cell)")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,targets.use = 4, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength (Neutrophils)")

# CXCL plot
signal_path<-"CXCL"
pairLR.use <- extractEnrichedLR(cellchat, signaling = signal_path)
netVisual_bubble(cellchat, pairLR.use = pairLR.use)
netVisual_individual(cellchat, signaling = signal_path, pairLR.use = "CXCL2_CXCR2", layout = "circle")
netVisual_individual(cellchat, signaling = signal_path, pairLR.use = "CXCL2_CXCR1", layout = "circle")
netVisual_individual(cellchat, signaling = signal_path, pairLR.use = "CXCL12_CXCR4", layout = "circle")
netVisual_individual(cellchat, signaling = signal_path, pairLR.use = "CXCL1_CXCR2", layout = "circle")
netVisual_heatmap(cellchat, signaling = signal_path, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = signal_path)

signal_path<-c("CXCL","CCL")
pairLR.use <- extractEnrichedLR(cellchat, signaling = signal_path)
netVisual_bubble(cellchat, pairLR.use = pairLR.use,targets.use = 10,title.name = "CXCL/CCL interaction (Cancercell target)")
netVisual_bubble(cellchat, pairLR.use = pairLR.use,targets.use = 4,title.name = "CXCL/CCL interaction (Neutrophils target)")
netAnalysis_contribution(cellchat, signaling = signal_path)


## cellcaht  prob.original>0.005
cellchat<-cellchat.all
signal_path<-c("CXCL","CCL")
pairLR <- extractEnrichedLR(cellchat, signaling = signal_path)
pairLR.use<-c("CXCL2_CXCR2","CCL3_CCR1","CCL4_CCR1")%>%as.data.frame()
colnames(pairLR.use)<-colnames(pairLR)
df<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,targets.use = 3,sources.use = c(8,9),
                     n.colors = 5,
                     min.quantile =0.2,
                     title.name = "CXCL/CCL interaction (Neutrophils target)",remove.isolate=T,return.data = T)
df<-df$communication%>%as.data.frame()%>%dplyr::select(source.target,prob.original,interaction_name_2)
df<-df%>%filter(interaction_name!='CCL3_CCR1'|source!="N2")  ### filter

## ggplot 
{  
  PLOT<-function(df,color.heatmap = "Spectral", 
                 n.colors = 10, direction = -1, thresh = 0.05, comparison = NULL, 
                 group = NULL, remove.isolate = FALSE, max.dataset = NULL, 
                 min.dataset = NULL, min.quantile = 0, max.quantile = 1, line.on = TRUE, 
                 line.size = 0.2, color.text.use = TRUE, color.text = NULL, 
                 title.name = NULL, font.size = 10, font.size.title = 10, 
                 show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
                 angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE){
    
    if (is.null(vjust.x) | is.null(hjust.x)) {
      angle = c(0, 45, 90)
      hjust = c(0, 1, 1)
      vjust = c(0, 1, 0.5)
      vjust.x = vjust[angle == angle.x]
      hjust.x = hjust[angle == angle.x]
    }
    
    g<-ggplot(df, aes(x = source.target, y = interaction_name_2, 
                      color = prob, size = pval)) + geom_point(pch = 16) + 
      theme_linedraw() + theme(panel.grid.major = element_blank()) + 
      theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
                                       vjust = vjust.x), axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
    values <- c(1, 2, 3)
    names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
    
    
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)})
    color.use <-color.use[length(color.use):1]
    
    g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), 
                          breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
                                                                                   sort(unique(df$pval))], name = "p-value")
    if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
      g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                      na.value = "white", limits = c(quantile(df$prob, 
                                                                              0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
                                      breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
                                                                                           1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.8, barheight=3,
                                                                                                                                                                      title = "Commun. Prob."))
    }
    else {
      g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                      na.value = "white") + guides(color = guide_colourbar(barwidth = 1, barheight=2,
                                                                                           title = "Commun. Prob."))
    }
    g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
      theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
    if (grid.on) {
      if (length(unique(df$source.target)) > 1) {
        g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                                               0.5, 1), lwd = 0.1, colour = color.grid)
      }
      if (length(unique(df$interaction_name_2)) > 1) {
        g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                               0.5, 1), lwd = 0.1, colour = color.grid)
      }
    }
    if (!is.null(title.name)) {
      g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
    }
    if (!is.null(comparison)) {
      if (line.on) {
        xintercept = seq(0.5 + length(dataset.name[comparison]), 
                         length(group.names0) * length(dataset.name[comparison]), 
                         by = length(dataset.name[comparison]))
        g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                            color = "grey60", size = line.size)
      }
      if (color.text.use) {
        if (is.null(group)) {
          group <- 1:length(comparison)
          names(group) <- dataset.name[comparison]
        }
        if (is.null(color.text)) {
          color <- ggPalette(length(unique(group)))
        }
        else {
          color <- color.text
        }
        names(color) <- names(group[!duplicated(group)])
        color <- color[group]
        dataset.name.order <- levels(df$source.target)
        dataset.name.order <- stringr::str_match(dataset.name.order, 
                                                 "\\(.*\\)")
        dataset.name.order <- stringr::str_sub(dataset.name.order, 
                                               2, stringr::str_length(dataset.name.order) - 
                                                 1)
        xtick.color <- color[dataset.name.order]
        g <- g + theme(axis.text.x = element_text(colour = xtick.color))
      }
    }
    return(g)
  }
  
  p<-PLOT(df)
  ggsave(p,filename = "Fig8a.pdf",width = 3,height = 3)
}


#### Fig8b #######################################################################
## remove baso effect
library(CellChat)
data<-subset(data.bk,celltype_cancer!='Basophils')
data$celltype_cancer<-factor(data$celltype_cancer,levels = c("B cells","NK cells","T_ILC_NKT","Neutrophils","MoMacDc","Fibroblasts","Endothelial cells","Epithelial cells","Cancer cells"))

cellchat<-createCellChat(data,assay = "integrated",group.by = "celltype_cancer")
CellChatDB <- CellChatDB.mouse

# add CCL4-CCR1
CCL4_CCR1<-c('CCL4_CCR1','CCL','Ccl4','Ccr1','','','' ,'' ,'PMID:23233369','Secreted Signaling','Ccl4  - Ccr1')
aa<-CellChatDB$interaction
aa<-rbind(aa,CCL4_CCR1)
rownames(aa)<-aa$interaction_name
CellChatDB$interaction<-aa

## select db
CellChatDB.use <- subsetDB(CellChatDB, search =c("Secreted Signaling","Cell-Cell Contact","ECM-Receptor") )    # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use 

## run cellchat
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb_modify(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

signal_path<-"CCL"
pairLR.use <- extractEnrichedLR(cellchat, signaling = signal_path)

## plot 
pdf("Fig8b.pdf",width = 10,height = 10)
netVisual_individual(cellchat, signaling = signal_path, pairLR.use = "CCL3_CCR1", layout = "circle")
netVisual_individual(cellchat, signaling = signal_path, pairLR.use = "CCL4_CCR1", layout = "circle")
dev.off()
