
if(F){
  library(Seurat)
  library(monocle)
  library(tidyverse)
  library(patchwork)
  library(future)
  rm(list=ls())
  
  setwd("~/covid")
  dir.create("monocle/Monocle4")
  
  
  ####analysis####

  load("macro_p.Rdata")
  sco <- macro_p   
  sco$celltype_new <- as.factor(as.character(sco$celltype_new))
  

  data <- GetAssayData(sco, assay = "RNA", slot = "counts")

  pd <- new('AnnotatedDataFrame', data = sco@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  mycds <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
  

  mycds <- estimateSizeFactors(mycds)
  mycds <- estimateDispersions(mycds, cores=8)
  

  order.genes <- VariableFeatures(sco)
  saveRDS(order.genes, "monocle/Monocle/order.genes.seurat.rds")

  if(F){
    tmp <- order.genes

    ex1 <- c(cc.genes$s.genes, cc.genes$g2m.genes)
    ex2 <- grep('^MT-', rownames(mycds), value = T)
    ex3 <- grep('^RP[LS]', rownames(mycds), value = T)
    ex <- unique(c(ex1, ex2, ex3))
    tmp <- tmp[!tmp %in% ex]

    include <- c('CD3D','CD4','CD8A')
    tmp <- union(tmp, include)
  }

  order.genes <- readRDS("~/project/Monocle/order.genes.seurat.rds")
  mycds <- setOrderingFilter(mycds, order.genes)
  p <- plot_ordering_genes(mycds)
  ggsave("Monocle/order.genes.seurat.pdf", p, width = 8, height = 6)

  disp_table <- dispersionTable(mycds)
  order.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit) %>% 
    pull(gene_id) %>% as.character()
  mycds <- setOrderingFilter(mycds, order.genes)
  p <- plot_ordering_genes(mycds)
  ggsave("monocle/Monocle/order.genes.monocle.pdf", p, width = 8, height = 6)
  

  mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
  
  mycds <- orderCells(mycds)
  mycds <- orderCells(mycds, root_state = 4)
 
  
  ####plot#### 
  


  load("monocle_p_4.Rdata")
  mycds <- mycds_p_4
  

  plot1 <- plot_cell_trajectory(mycds, color_by = "State")


  plot2 <- plot_cell_trajectory(mycds, color_by = "celltype_new")

  

  plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")

  

  plot4 <- plot_cell_trajectory(mycds, color_by = "group")


  plotc <- plot1|plot3|plot2|plot4

  ggsave("monocle_trajectory.pdf", plotc, width = 8, height = 6)
  


  
  
  
  
}
