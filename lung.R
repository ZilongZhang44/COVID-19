####读取数据，构建对象####
library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(AnnotationHub)  
library(org.Hs.eg.db)   

setwd("D:/准备文章/covid19/code/covid/GSE145926_RAW")

C141<-Read10X_h5("GSM4339769_C141_filtered_feature_bc_matrix.h5")
C142<-Read10X_h5("GSM4339770_C142_filtered_feature_bc_matrix.h5")
C143<-Read10X_h5("GSM4339771_C143_filtered_feature_bc_matrix.h5")
C144<-Read10X_h5("GSM4339772_C144_filtered_feature_bc_matrix.h5")
C145<-Read10X_h5("GSM4339773_C145_filtered_feature_bc_matrix.h5")
C146<-Read10X_h5("GSM4339774_C146_filtered_feature_bc_matrix.h5")
C148<-Read10X_h5("GSM4475051_C148_filtered_feature_bc_matrix.h5")
C149<-Read10X_h5("GSM4475052_C149_filtered_feature_bc_matrix.h5")
C152<-Read10X_h5("GSM4475053_C152_filtered_feature_bc_matrix.h5")
C51<- Read10X_h5("GSM4475048_C51_filtered_feature_bc_matrix.h5")
C52<- Read10X_h5("GSM4475049_C52_filtered_feature_bc_matrix.h5")
C100<- Read10X_h5("GSM4475050_C100_filtered_feature_bc_matrix.h5")

setwd("D:/covid19/code/covid")

GSM3660650 <- Read10X("GSM3660650")

C141<-CreateSeuratObject(counts = C141, project = "C141",min.cells = 3, min.features = 200)
C142<-CreateSeuratObject(counts = C142, project = "C142",min.cells = 3, min.features = 200)
C143<-CreateSeuratObject(counts = C143, project = "C143",min.cells = 3, min.features = 200)
C144<-CreateSeuratObject(counts = C144, project = "C144",min.cells = 3, min.features = 200)
C145<-CreateSeuratObject(counts = C145, project = "C145",min.cells = 3, min.features = 200)
C146<-CreateSeuratObject(counts = C146, project = "C146",min.cells = 3, min.features = 200)
C148<-CreateSeuratObject(counts = C148, project = "C148",min.cells = 3, min.features = 200)
C149<-CreateSeuratObject(counts = C149, project = "C149",min.cells = 3, min.features = 200)
C152<-CreateSeuratObject(counts = C152, project = "C152",min.cells = 3, min.features = 200)
C51<-CreateSeuratObject(counts = C51, project = "C51",min.cells = 3, min.features = 200)
C52<-CreateSeuratObject(counts = C52, project = "C52",min.cells = 3, min.features = 200)
C100<-CreateSeuratObject(counts = C100, project = "C100",min.cells = 3, min.features = 200)
GSM3660650<-CreateSeuratObject(counts = GSM3660650, project = "GSM3660650",min.cells = 3, min.features = 200)

C141$group<-"mild"
C142$group<-"mild"
C143$group<-"severe"
C144$group<-"mild"
C145$group<-"severe"
C146$group<-"severe"
C148$group<-"severe"
C149$group<-"severe"
C152$group<-"severe"
C51$group <- "healthy"
C52$group <- "healthy"
C100$group <- "healthy"
GSM3660650$group <- "healthy"

####QC####
C141[["percent.mt"]]<-PercentageFeatureSet(C141,pattern = "^MT")
p141<-VlnPlot(C141, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C142[["percent.mt"]]<-PercentageFeatureSet(C142,pattern = "^MT")
p142<-VlnPlot(C142, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C143[["percent.mt"]]<-PercentageFeatureSet(C143,pattern = "^MT")
p143<-VlnPlot(C143, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C144[["percent.mt"]]<-PercentageFeatureSet(C144,pattern = "^MT")
p144<-VlnPlot(C144, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C145[["percent.mt"]]<-PercentageFeatureSet(C145,pattern = "^MT")
p145<-VlnPlot(C145, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C146[["percent.mt"]]<-PercentageFeatureSet(C146,pattern = "^MT")
p146<-VlnPlot(C146, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C148[["percent.mt"]]<-PercentageFeatureSet(C148,pattern = "^MT")
p148<-VlnPlot(C148, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C149[["percent.mt"]]<-PercentageFeatureSet(C149,pattern = "^MT")
p149<-VlnPlot(C149, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C152[["percent.mt"]]<-PercentageFeatureSet(C152,pattern = "^MT")
p152<-VlnPlot(C152, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C51[["percent.mt"]]<-PercentageFeatureSet(C51,pattern = "^MT")
p51<-VlnPlot(C51, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C52[["percent.mt"]]<-PercentageFeatureSet(C52,pattern = "^MT")
p52<-VlnPlot(C52, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C100[["percent.mt"]]<-PercentageFeatureSet(C100,pattern = "^MT")
p100<-VlnPlot(C100, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GSM3660650[["percent.mt"]]<-PercentageFeatureSet(GSM3660650,pattern = "^MT")
pGSM3660650<-VlnPlot(GSM3660650, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

C141 <- subset(C141, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C142 <- subset(C142, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C143 <- subset(C143, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C144 <- subset(C144, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C145 <- subset(C145, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C146 <- subset(C146, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C148 <- subset(C148, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C149 <- subset(C149, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C152 <- subset(C152, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C51 <- subset(C51, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C52 <- subset(C52, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C100 <- subset(C100, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
GSM3660650 <- subset(GSM3660650, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)

C141 <- NormalizeData(C141)
C142 <- NormalizeData(C142)
C143 <- NormalizeData(C143)
C144 <- NormalizeData(C144)
C145 <- NormalizeData(C145)
C146 <- NormalizeData(C146)
C148 <- NormalizeData(C148)
C149 <- NormalizeData(C149)
C152 <- NormalizeData(C152)
C51 <- NormalizeData(C51)
C52 <- NormalizeData(C52)
C100 <- NormalizeData(C100)
GSM3660650<- NormalizeData(GSM3660650)


C141 <- FindVariableFeatures(C141, selection.method = "vst", nfeatures = 2000)
C142 <- FindVariableFeatures(C142, selection.method = "vst", nfeatures = 2000)
C143 <- FindVariableFeatures(C143, selection.method = "vst", nfeatures = 2000)
C144 <- FindVariableFeatures(C144, selection.method = "vst", nfeatures = 2000)
C145 <- FindVariableFeatures(C145, selection.method = "vst", nfeatures = 2000)
C146 <- FindVariableFeatures(C146, selection.method = "vst", nfeatures = 2000)
C148 <- FindVariableFeatures(C148, selection.method = "vst", nfeatures = 2000)
C149 <- FindVariableFeatures(C149, selection.method = "vst", nfeatures = 2000)
C152 <- FindVariableFeatures(C152, selection.method = "vst", nfeatures = 2000)
C51 <- FindVariableFeatures(C51, selection.method = "vst", nfeatures = 2000)
C52 <- FindVariableFeatures(C52, selection.method = "vst", nfeatures = 2000)
C100 <- FindVariableFeatures(C100, selection.method = "vst", nfeatures = 2000)
GSM3660650 <- FindVariableFeatures(GSM3660650, selection.method = "vst", nfeatures = 2000)

nCoV.list <- list(C141,C142, C143,C144,C145,C146,C148,C149,C152,C51,C52,C100,GSM3660650)

#save(nCov.list,C141,C142, C143,C144,C145,C146,C148,C149,C152,C51,C52,C100,GSM3660650,file="./data_processed.Rdata")
####数据整合####
#load("./data_processed.Rdata")
nCoV <- FindIntegrationAnchors(object.list = nCoV.list, dims = 1:50)
#save(nCoV,file="nCov.Rdata")
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))
save(nCoV.integrated,file="interated.Rdata")
####metadata####
#load("interated.Rdata")
samples = read.delim2("meta.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
sample_info = as.data.frame(colnames(nCoV.integrated))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = nCoV.integrated@meta.data$orig.ident
sample_info = dplyr::left_join(sample_info,samples)
rownames(sample_info) = sample_info$ID
nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)

####QC integrated####
DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
ElbowPlot(object = nCoV.integrated,ndims = 100)

###cluster
nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2)

nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)

nCoV.integrated <- RunUMAP(object = nCoV.integrated, dims = 1:50)
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)

nCoV.integrated@misc$markers <- FindAllMarkers(object = nCoV.integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
write.table(nCoV.integrated@misc$markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')




find_markers<- read.table(file = "marker_MAST.txt",sep = '\t')

DefaultAssay(nCoV.integrated) <- "RNA"
# find markers for every cluster compared to all remaining cells, report only the positive ones

VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
supplementary_1 <- VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
ggsave(filename ="./pdf/supplmentary_1.pdf", plot =  supplementary_1)


markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF','DCN',
            'FCGR3A','TREM2','KRT18')
hc.markers=nCoV.integrated@misc$markers
hc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) -> top30
var.genes = c(nCoV.integrated@assays$RNA@var.features,top30$gene,markers)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"),features = var.genes)

#save(nCoV.integrated,file="before_downstream.Rdata")

####annotation####
setwd("../code/covid/")
load("before_downstream.Rdata")


DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)

DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'sample_new')

DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'sample_new', ncol = 4)

DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'group')

DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'group')

DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'disease')

DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'disease', ncol = 2)

markers = c('TPPP3','KRT18','CD68','FCGR3B','CD1C','CLEC9A','LILRA4','TPSB2','CD3D','KLRD1','MS4A1','IGHG4')

pp_temp = FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
supplement_lung_2  <- CombinePlots(plots = plots,ncol = 4,legend = 'right')
supplement_lung_2



pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + theme(axis.line = element_line(size = 0.6))
fig_1_2<- pp
fig_1_2

library(ggplot2)
nCoV.integrated[["cluster"]] <- Idents(object = nCoV.integrated)
big.cluster = nCoV.integrated@meta.data
organ.summary = table(big.cluster$sample_new,big.cluster$cluster)
#write.table(organ.summary,file = '1-nCoV-percentage-sample.txt',quote = FALSE,sep = '\t')
library(ggplot2)
library(dplyr)
library(reshape2)
organ.summary = read.delim2("1-nCoV-percentage-sample.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
organ.summary$group = rownames(organ.summary)
organ.summary.dataframe = melt(organ.summary)
colnames(organ.summary.dataframe) = c('group','cluster','cell')
samples_name_new = c('HC1','HC2','HC3','HC4','M1','M2','M3','S1','S2','S3','S4','S5','S6')
organ.summary.dataframe$group = factor(organ.summary.dataframe$group,labels = samples_name_new,levels = samples_name_new)
organ.summary.dataframe$cell = as.numeric(organ.summary.dataframe$cell)

if (FALSE){
  

  new_order = c('0','1','2','3','4','5','7','8','9','11','12','16','17','20','21','23','28','19','24','27','13','25','15','6','10','30','18','22','26','29','31','14'
                )
  organ.summary.dataframe$cluster = factor(organ.summary.dataframe$cluster,levels = new_order,labels = new_order)
  cols = c('#32b8ec','#60c3f0','#8ccdf1','#cae5f7','#92519c','#b878b0','#d7b1d2','#e7262a','#e94746','#eb666d','#ee838f','#f4abac','#fad9d9')
  pp = ggplot(data=organ.summary.dataframe, aes(x=cluster, y=cell, fill=group)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Fraction of sample per cluster (%)') +
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
          legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) +
    scale_fill_manual(values = cols)
  pp
}

#mast
FeaturePlot(nCoV.integrated,features = 'TPSB2') 
#pdc
FeaturePlot(nCoV.integrated,features = 'LILRA4')
#plasma
FeaturePlot(nCoV.integrated,features = 'IGHG4')
#B cells
FeaturePlot(nCoV.integrated,features = 'MS4A1')

nCoV.integrated1 <- RenameIdents(object = nCoV.integrated,
                                 '0'='Macrophages',
                                 '1'='Macrophages',
                                 '2'='Macrophages',
                                 '3'='Macrophages',
                                 '4'='Macrophages',
                                 '5'='Macrophages',
                                 '6'='T cells',
                                 '7'='Macrophages',
                                 '8'='Macrophages',
                                 '9'='T cells',
                                 '10'='Macrophages',
                                 '11'='Macrophages',
                                 '12'='Macrophages',
                                 '13'='Ciliated cells',
                                 '14'='T cells',
                                 '15'='Secretory cells',
                                 '16'='Macrophages',
                                 '17'='Macrophages',
                                 '18'='NK cells',
                                 '19'='Neutrophil',
                                 '20'='Macrophages',
                                 '21'='mDC',
                                 '22'='Macrophages',
                                 '23'='Plasma cells',
                                 '24'='Epithelial',
                                 '25'='B cells', 
                                 '26'='Plasma cells',
                                 '27'='Macrophages',
                                 '28'='pDC',
                                 '29'='T cells',
                                 '30'='Mast cells'
                                 )

DimPlot(object = nCoV.integrated1, reduction = 'umap',label = T, label.size = 6, ncol = 3,
        repel = TRUE,combine = TRUE)

#save(nCoV.integrated1,file="final_downstream.Rdata")


####downstream####
load("final_downstream.Rdata")
#


if (F){
  nCoV.integrated1$celltype = Idents(nCoV.integrated1)
  
  
    
    nCoV.integrated1 = subset(nCoV.integrated1,subset = celltype != 'Doublets')
  
  
  
  nCoV_groups = c('Epithelial','Macrophages','Neutrophil','mDC','pDC','Mast','T cells','NK cells ',
                  'B cells','Plasma cells','Ciliated cells','Secretory cells')
  nCoV.integrated1$celltype = factor(nCoV.integrated1$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(nCoV.integrated1) = nCoV.integrated1$celltype
}

fig_lung_1 <-  DimPlot(object = nCoV.integrated, reduction = 'umap',label = T, label.size = 4, ncol = 3,
                       combine = TRUE) 
fig_lung_1

fig_lung_2 <- DimPlot(object = nCoV.integrated1, reduction = 'umap',label = T, label.size = 4, ncol = 3,
       combine = TRUE)
fig_lung_2


markers = c('TPPP3','KRT18','CD68','FCGR3B','CD1C','CLEC9A','LILRA4','TPSB2','CD3D','KLRD1','MS4A1','IGHG4')

if (FALSE){
  

  pp_temp = FeaturePlot(object = nCoV.integrated1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + 
                    theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
  pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
  fig_lung_3  <- CombinePlots(plots = plots,ncol = 4,legend = 'right')
  fig_lung_3
}
pp_temp = FeaturePlot(object = nCoV.integrated1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12),
                        plot.title = element_text(family = 'sans',face='italic',size=16),
                        legend.text = element_text(size = 16),legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
fig_lung_3  <- CombinePlots(plots = plots,ncol = 4,legend = 'right')
fig_lung_3

#ggsave("./pdf/sup_lung_2.pdf", plot = sup_lung_2)

#"SR-B1"
nCoV.integrated1$celltype <- Idents(nCoV.integrated1)
VlnPlot(nCoV.integrated1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')

DotPlot(nCoV.integrated1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')

fig_lung_4 <- DotPlot(nCoV.integrated1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'disease')
fig_lung_4


fig_lung_5 <- DotPlot(nCoV.integrated1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
fig_lung_5


fig_lung <- CombinePlots(plots = list(fig_lung_1,fig_lung_2,fig_lung_3,fig_lung_4,fig_lung_5))
fig_lung
ggsave("./pdf/lung.pdf", plot = fig_lung,width = 15, height = 18)



markers = c('TPPP3','KRT18','CD68','FCGR3B','CD1C','CLEC9A','LILRA4','TPSB2','CD3D','KLRD1','MS4A1','IGHG4')

supplement_lung_1 <- DimPlot(nCoV.integrated,group.by = 'orig.ident',reduction = 'pca')
supplement_lung_1
ggsave("./pdf/suplement_lung_1.pdf", plot = supplement_lung_1)

pp_temp = FeaturePlot(object = nCoV.integrated1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
pp
supplement_lung_2  <- CombinePlots(plots = plots,ncol = 4,legend = 'right')
supplement_lung_2


pp = DotPlot(nCoV.integrated1, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + theme(axis.line = element_line(size = 0.6))
supplement_lung_3<- pp
supplement_lung_3

if (F){
  pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',label = FALSE, label.size = 6,split.by = 'group', ncol = 3,repel = TRUE,combine = TRUE)
  pp_temp = pp_temp + theme(axis.title = element_text(size = 17),axis.text = element_text(size = 17),strip.text = element_text(family = 'arial',face='plain',size=17), legend.text = element_text(size = 17),axis.line = element_line(size = 1),axis.ticks = element_line(size = 0.8),legend.key.height = unit(1.4,"line"))
  pp_temp
  
  pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',label = FALSE, split.by = 'sample_new', label.size = 7,ncol = 5,repel = TRUE, combine = TRUE,pt.size = 1.5)
  pp_temp = pp_temp + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 22),strip.text = element_text(family = 'sans',face='plain',size=22),legend.text = element_text(size = 22),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line"),axis.line = element_line(size = 1.2),axis.ticks = element_line(size = 1.2))
  pp_temp
}








####macrophage####
#mac <- subset(nCoV.integrated1,celltype=='Macrophages')
#mac_iden<- c("0","1","2","3","4","5","7","8","10",'11','12','16','17',
        '20','22','27')
lung <- nCoV.integrated1
lung$celltype <- Idents(lung)
mac_dis <- FindMarkers(lung, ident.1 = "Y", 
                       group.by = 'disease', 
                       subset.ident = "Macrophages")
top100 <- rownames(mac_dis[1:100,])

top100

write.table(top100,file = "./enrichment/lung_top100.txt",quote = F,row.names = F)

####GO, KEGG####
go_lung <- read.table("./enrichment/lung_top100.txt",sep=" ")

go_lung <- t(go_lung)
keytypes(org.Hs.eg.db)

go_lung_id_trance <- bitr(go_lung, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)
write.table(go_lung_id_trance$ENSEMBL,"./enrichment/go_lung_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./enrichment/go_lung_id_trance.txt")
f <- f[c(1)]

EG2Ensembl = toTable(org.Hs.egENSEMBL)
f = f$V1
geneLists = data.frame(ensembl_id = f)
results = merge(geneLists,EG2Ensembl,by='ensembl_id',all.x  =T)
id =na.omit(results$gene_id)

All <- enrichGO(OrgDb="org.Hs.eg.db",gene=id,ont="ALL",readable=T)

dotplot(All,showCategory=10,title="Enrichment Go")+
  barplot(All, showCategory=15,title="EnrichmentGO") 

dotplot(All,showCategory=10,title="Enrichment Go")
barplot(All, showCategory=15,title="EnrichmentGO") 


barplot(All, showCategory=5,title="EnrichmentGO",split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free")

fig_lung_6 <- dotplot(All,split="ONTOLOGY",title ="Enrichment GO Top5",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
fig_lung_6


####KEGG
KEGG <- enrichKEGG(gene= id, organism  = 'hsa', pvalueCutoff = 0.05)
dotplot(KEGG,font.size=12)
barplot(KEGG,font.size=8) 
barplot(KEGG, font.size=12, showCategory=10, title="Enrichment KEGG Top10")+
  scale_size(rang=c(5.20))

dotplot(KEGG,showCategory=5,title="Enrichment KEGG Top5")
fig_lung_7<-dotplot(KEGG, font.size=12, showCategory=10, title="Enrichment KEGG Top10")+
  scale_size(rang=c(5.20))
fig_lung_7


lung_all <- CombinePlots(plots = list(fig_lung_1,fig_lung_2,fig_lung_3,
                                fig_lung_4,fig_lung_5,fig_lung_6,fig_lung_7))
lung_all

ggsave("./fig_2/A.pdf",plot=fig_lung_1)
ggsave("./fig_2/B.pdf",plot=fig_lung_2)
ggsave("./fig_2/C.pdf",plot=fig_lung_3)
ggsave("./fig_2/D.pdf",plot=fig_lung_4)
ggsave("./fig_2/E.pdf",plot=fig_lung_5)
ggsave("./fig_2/F.pdf",plot=fig_lung_6)
ggsave("./fig_2/G.pdf",plot=fig_lung_7)


ggsave("./pdf/fig_lung_all.pdf", plot = lung_all, width = 15, height = 18)
#ggsave("./pdf/fig_lung.pdf", plot = lung_all, width = 30, height = 36)


#save(All,KEGG,file="KEGG.Rdata")
load("KEGG.Rdata")

####
c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL") %in% rownames(mac_dis)


DotPlot(mac,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'disease')
