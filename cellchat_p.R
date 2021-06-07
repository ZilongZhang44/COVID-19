Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")

library(Seurat)
library(tidyverse)  
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
rm(list = ls())
options(stringsAsFactors = FALSE)
  

load("nCoV.integrated2.Rdata")
  scRNA <- subset(nCoV.integrated2,disease=="Y")
  Idents(scRNA) <- 'celltype_new'

  table(scRNA$orig.ident)



cov_p <- createCellChat(scRNA@assays$RNA@data, meta = scRNA@meta.data, group.by = "celltype_new")


if(F){

  if(F){


      cellchat <- cov_p

      cellchat <- setIdent(cellchat, ident.use = "celltype_new")
      groupSize <- as.numeric(table(cellchat@idents))     # number of cells in each cell group
      

      CellChatDB <- CellChatDB.human 
      showDatabaseCategory(CellChatDB)

      CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

      # set the used database in the object
      cellchat@DB <- CellChatDB.use
      

      cellchat <- subsetData(cellchat) # subset the expression data of signaling genes
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- projectData(cellchat, PPI.human)
      
    }
    

    if(F){

      cellchat <- computeCommunProb(cellchat, raw.use = FALSE,   population.size = TRUE)
      # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      df.net <- subsetCommunication(cellchat)
      write.csv(df.net, "net_lr.csv")
      

      cellchat <- computeCommunProbPathway(cellchat)
      df.netp <- subsetCommunication(cellchat, slot.name = "netP")
      write.csv(df.netp, "net_pathway.csv")
      
      saveRDS(cellchat,"cov_p.rds")
    }
 
    
  
   
    

    if(F){
      cellchat <- read_rds("cov_p.rds")
      cellchat <- aggregateNet(cellchat)
      groupSize <- as.numeric(table(cellchat@idents))
      par(mfrow = c(1,2), xpd=TRUE)
      netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                       label.edge= F, title.name = "Number of interactions")
      netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                       label.edge= F, title.name = "Interaction weights/strength")
      if (F){
      par(mfrow = c(1,1), xpd=TRUE)
      interaction_weight <-  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                                              label.edge= F, title.name = "Interaction weights/strength")
      ggsave("interaction_weight.pdf", interaction_weight, width = 5, height = 6)
      }

      par(mfrow = c(1,1), xpd=TRUE)
      p_net <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                                label.edge= F, title.name = "Interaction weights/strength")


      mat <- cellchat@net$count
      par(mfrow = c(3,4), xpd=TRUE)
      for (i in 1:nrow(mat)) {
        # i = 1
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                         arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
      }

      
      mat <- cellchat@net$weight
      par(mfrow = c(3,4), xpd=T)
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                         arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
      }

    }

    if(F){
      cellchat@netP$pathways      # show all signaling pathways
      pathways.show <- c("COMPLEMENT")  # select one signaling pathway to visualising

      # Hierarchy plot
      levels(cellchat@idents)    # show all celltype
      vertex.receiver = c(1,2,3,6) # define a numeric vector giving the index of the celltype as targets
      netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

      
      # Circle plot
      par(mfrow=c(1,1))
      netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

      
      # Chord diagram
      par(mfrow=c(1,1))
      netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

      
      # Heatmap
      par(mfrow=c(1,1))
      netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

      
      netAnalysis_contribution(cellchat, signaling = pathways.show)
      pairLR.COMPLEMENT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)

      
      # Hierarchy plot
      LR.show <- pairLR.COMPLEMENT[1,] # show one ligand-receptor pair
      vertex.receiver = c(1,2,3,6) # a numeric vector
      netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, 
                           vertex.receiver = vertex.receiver)

      
      # Circle plot
      netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
      # save as TIL/CXCL_circle2.pdf
      
      # Chord diagram
      netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


      # Access all the signaling pathways showing significant communications
      pathways.show.all <- cellchat@netP$pathways

      levels(cellchat@idents)
      vertex.receiver = c(1,2,3,6)
      dir.create("all_pathways_com")

      for (i in 1:length(pathways.show.all)) {
        # Visualize communication network associated with both signaling pathway and individual L-R pairs
        netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
                  vertex.receiver = vertex.receiver, layout = "hierarchy")
        # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
        gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
        ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
               plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
      }
    
    }

    if(F){
      ### Bublle plot
      levels(cellchat@idents)
      # show all the significant interactions (L-R pairs)
      p = netVisual_bubble(cellchat, sources.use = c(4,5), 
                           targets.use = c(1,2,3,6), remove.isolate = FALSE)
      ggsave("cov_p_bubble.pdf", p, width = 8, height = 14)

      if(F){
        netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), 
                         signaling = c("CCL","CXCL"), remove.isolate = FALSE)
        # show all the significant interactions (L-R pairs) based on user's input
        pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
        netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), 
                         pairLR.use = pairLR.use, remove.isolate = TRUE)
      }
      
      ## Plot the signaling gene expression distribution
      p = plotGeneExpression(cellchat, signaling = "CXCL")
      ggsave("CXCL_GeneExpression_vln.pdf", p, width = 8, height = 8)
      p = plotGeneExpression(cellchat, signaling = "CXCL", type = "dot")
      ggsave("CXCL_GeneExpression_dot.pdf", p, width = 8, height = 6)
    }

    if(F){

      cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

      netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                        width = 15, height = 6, font.size = 10)

      ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5)
      ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5)
      ht1 + ht2

    }
    
    if(F){

      selectK(cellchat, pattern = "outgoing")

      nPatterns = 3 
      cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, 
                                                width = 8, height = 12, font.size = 6)

      
      # river plot
      netAnalysis_river(cellchat, pattern = "outgoing")

      
      # dot plot
      netAnalysis_dot(cellchat, pattern = "outgoing")


      selectK(cellchat, pattern = "incoming")

      nPatterns = 3
      cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, 
                                                width = 5, height = 12, font.size = 10)

      netAnalysis_river(cellchat, pattern = "incoming")


      netAnalysis_dot(cellchat, pattern = "incoming")

    }
    

    if(F){
      cellchat <- computeNetSimilarity(cellchat, type = "functional")
      cellchat <- netEmbedding(cellchat, type = "functional")
      cellchat <- netClustering(cellchat, type = "functional")
      p = netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
      ggsave("Manifold_functional_cluster.pdf", p, width = 8, height = 6)


      cellchat <- computeNetSimilarity(cellchat, type = "structural")
      cellchat <- netEmbedding(cellchat, type = "structural")
      cellchat <- netClustering(cellchat, type = "structural")
      p = netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
      ggsave("Manifold_structural_cluster.pdf", p, width = 8, height = 6)
    }
    

    
    
####role####
    gg1 <- netAnalysis_signalingRole_scatter(cellchat)
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    # Signaling role analysis on the cell-cell communication networks of interest
    
    
    g<- netAnalysis_signalingRole_scatter(cellchat,label.size = 5)
  }
  ggsave("figure3_B",g,width = 12.5,height = 11)
  
