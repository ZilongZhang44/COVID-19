####clustering####
library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(AnnotationHub)  
library(org.Hs.eg.db)  

setwd("D:/准备文章/covid19/data/")

bladder_1<- Read10X(data.dir = "./GSE129845_bladder/human/1/")
bladder_2<- Read10X(data.dir = "./GSE129845_bladder/human/2/")
bladder_3<- Read10X(data.dir = "./GSE129845_bladder/human/3/")

bladder_1<-CreateSeuratObject(counts = bladder_1, project = "bladder_1",min.cells = 3, min.features = 200)
bladder_2<-CreateSeuratObject(counts = bladder_2, project = "bladder_2",min.cells = 3, min.features = 200)
bladder_3<-CreateSeuratObject(counts = bladder_3, project = "bladder_3",min.cells = 3, min.features = 200)

bladder_1[["percent.mt"]]<-PercentageFeatureSet(bladder_1,pattern = "^MT")
pbladder_1 <-VlnPlot(bladder_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
bladder_2[["percent.mt"]]<-PercentageFeatureSet(bladder_2,pattern = "^MT")
pbladder_2 <-VlnPlot(bladder_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
bladder_3[["percent.mt"]]<-PercentageFeatureSet(bladder_3,pattern = "^MT")
pbladder_3 <-VlnPlot(bladder_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbladder_1
pbladder_2
pbladder_3


bladder_1 <- subset(bladder_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 )
bladder_2 <- subset(bladder_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 )
bladder_3 <- subset(bladder_3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 )


bladder_1 <- NormalizeData(bladder_1)
bladder_2 <- NormalizeData(bladder_2)
bladder_3 <- NormalizeData(bladder_3)

bladder_1 <- FindVariableFeatures(bladder_1, selection.method = "vst",   mean.cutoff = c(0.0125, 3),
            dispersion.cutoff = c(0.5, Inf))
bladder_2 <- FindVariableFeatures(bladder_2, selection.method = "vst",mean.cutoff = c(0.0125, 3),
                                  dispersion.cutoff = c(0.5, Inf))
bladder_3 <- FindVariableFeatures(bladder_3, selection.method = "vst", mean.cutoff = c(0.0125, 3),
                                  dispersion.cutoff = c(0.5, Inf))

bladder_list <- list(bladder_1,bladder_2,bladder_3)

bladder <- FindIntegrationAnchors(object.list = bladder_list, dims = 1:50)
bladder <- IntegrateData(anchorset = bladder, dims = 1:50)



#save(bladder,file="bladder.Rdata")

#load("./GSE129845_bladder/bladder.Rdata")

DefaultAssay(bladder) <- "RNA"
bladder[['percent.mito']] <- PercentageFeatureSet(bladder, pattern = "^MT-")
bladder <- NormalizeData(object = bladder, normalization.method = "LogNormalize", scale.factor = 1e4)
bladder <- FindVariableFeatures(object = bladder, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
bladder <- ScaleData(bladder, verbose = FALSE, vars.to.regress = c("nCount_RNA","percent.mito"))
VlnPlot(object = bladder, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(object = bladder, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

DefaultAssay(bladder) <- "integrated"
bladder <- ScaleData(bladder, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
bladder <- RunPCA(bladder, verbose = FALSE,npcs = 100)
bladder <- ProjectDim(object = bladder)
ElbowPlot(object = bladder,ndims = 100)

###cluster
bladder <- FindNeighbors(object = bladder, dims = 1:15)
bladder_integrated <- FindClusters(object = bladder_integrated, resolution = 0.8)

bladder_integrated <- RunTSNE(object = bladder_integrated, dims = 1:15)
DimPlot(object = bladder_integrated, reduction = 'tsne',label = TRUE)

bladder_integrated <- RunUMAP(object = bladder_integrated, dims = 1:15)
DimPlot(object = bladder_integrated, reduction = 'umap',label = TRUE)
p_umap<- DimPlot(object = bladder_integrated, reduction = 'umap',label = TRUE)
ggsave("./GSE129845_bladder/pdf/bladder_umap_2.pdf", plot = p_umap)

save(bladder_integrated,file="./GSE129845_bladder/before_downstream.Rdata")
#load("before_downstream.Rdata")
DefaultAssay(bladder_integrated) <- "RNA"
# find markers for every cluster compared to all remaining cells, report only the positive ones
#bladder_integrated@misc$markers <- FindAllMarkers(object = bladder_integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')


####annotation####
FeaturePlot(bladder_integrated,features = "KRT18")
FeaturePlot(bladder_integrated,features = "KRT19")

#basal
FeaturePlot(bladder_integrated,features = "KRT5")
FeaturePlot(bladder_integrated,features = "KRT17")

#intermediate
FeaturePlot(bladder_integrated,features = "KRT13")

#Umbrella
FeaturePlot(bladder_integrated,features = "KRT20")
FeaturePlot(bladder_integrated,features = "UPK1A")
FeaturePlot(bladder_integrated,features = "UPK1B")
FeaturePlot(bladder_integrated,features = "UPK3B")
FeaturePlot(bladder_integrated,features = "UPK3A")
FeaturePlot(bladder_integrated,features = "UPK2")

#interstitial
FeaturePlot(bladder_integrated,features = "VIM")

#fibroblast
FeaturePlot(bladder_integrated,features = "S100A4")
FeaturePlot(bladder_integrated,features = "COL3A1")
FeaturePlot(bladder_integrated,features = "COL1A1")
FeaturePlot(bladder_integrated,features = "COL1A2")

#myofibroblast
FeaturePlot(bladder_integrated,features = "ACTA2")

#smooth muscle cell
FeaturePlot(bladder_integrated,features = "TAGLN")
FeaturePlot(bladder_integrated,features = "DES")
FeaturePlot(bladder_integrated,features = "CNN1")
FeaturePlot(bladder_integrated,features = "ACTG2")
FeaturePlot(bladder_integrated,features = "TPM2")

#endothelial
FeaturePlot(bladder_integrated,features = "SELE")
FeaturePlot(bladder_integrated,features = "PECAM1")
FeaturePlot(bladder_integrated,features = "VCAM1")
FeaturePlot(bladder_integrated,features = "CDH5")

#monocyte
FeaturePlot(bladder_integrated,features = "LYZ")
FeaturePlot(bladder_integrated,features = "MS4A7")
FeaturePlot(bladder_integrated,features = "CD14")

FeaturePlot(bladder_integrated,features = "CD209")

#T cell
FeaturePlot(bladder_integrated,features = "CD3D")
FeaturePlot(bladder_integrated,features = "CD3E")

#B cell
FeaturePlot(bladder_integrated,features = "MZB1")
FeaturePlot(bladder_integrated,features = "CD79A")

FeaturePlot(bladder_integrated,features = "GPM6A")


bladder_integrated_1 <- RenameIdents(object = bladder_integrated,
                                 '0'='Fibroblast',
                                 '1'='Basel cells',
                                 '2'='Basel cells',
                                 '3'='Basel cells',
                                 '4'='Intermediate cells',
                                 '5'='Basel cells',
                                 '6'='Fibroblast',
                                 '7'='Intermediate cells',
                                 '8'='Smooth muscle cells',
                                 '9'='Intermediate cells',
                                 '10'='Intermediate cells',
                                 '11'='Fibroblast',
                                 '12'='Endothelial cells',
                                 '13'='Monocytes',
                                 '14'='T cells',
                                 '15'='Umbrella cells',
                                 '16'='Myofibroblast',
                                 '17'='Interstitial cells',
                                 '18'='B cells'
                               
                              
                                 )

save(bladder_integrated_1,file="./GSE129845_bladder/bladder_final.Rdata")

####downstream####
load("./GSE129845_bladder/before_downstream.Rdata")
load("./GSE129845_bladder/bladder_final.Rdata")
DefaultAssay(bladder_integrated) <- "RNA"
DimPlot(object = bladder_integrated, reduction = 'umap',label = T)
fig_bladder_1 <- DimPlot(object = bladder_integrated, reduction = 'umap',label = T)
fig_bladder_1
ggsave(filename = "./GSE129845_bladder/pdf/bladder_clustering.pdf",plot=fig_bladder_1)

DimPlot(object = bladder_integrated_1, reduction = 'umap',label = T)
fig_bladder_2<-DimPlot(object =bladder_integrated_1, reduction = 'umap',label = T)
fig_bladder_2
ggsave(filename = "./GSE129845_bladder/pdf/bladder_annotated_final.pdf",plot=fig_bladder_2)

bladder_integrated_1$celltype = Idents(bladder_integrated_1)

markers = c("KRT5","UPK1A","KRT13","VIM","ACTA2","COL1A2","DES","SELE","LYZ","CD3D","MZB1")

pp_temp = FeaturePlot(object = bladder_integrated_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12),
                        plot.title = element_text(family = 'sans',face='italic',size=16),
                        legend.text = element_text(size = 16),legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,legend = 'right')
pp
fig_bladder_3 <- pp
fig_bladder_3

pp_temp = VlnPlot(object = bladder_integrated_1, features = markers,combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p )
pp = CombinePlots(plots = plots,legend = 'right')
pp
vln_marker <- pp
vln_marker


pp_dot = DotPlot(bladder_integrated_1, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp_dot = pp_dot + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + theme(axis.line = element_line(size = 0.6))

pp_dot
supplement_bladder_integrated_1 <- pp_dot

supplement_bladder_integrated_2 <- DoHeatmap(bladder_integrated_1, features = markers)
supplement_bladder_integrated_2


####


FeaturePlot(bladder_integrated,features = "NRP1")
FeaturePlot(bladder_integrated,features = "AXL")
FeaturePlot(bladder_integrated,features = "ACE2")
FeaturePlot(bladder_integrated,features = "TMPRSS2")
FeaturePlot(bladder_integrated,features = "FURIN")
FeaturePlot(bladder_integrated,features = "CTSL")
#FeaturePlot(bladder_integrated,features = "GRP78")
#FeaturePlot(bladder_integrated,features = "CD147")

VlnPlot(bladder_integrated_1,features = "NRP1",group.by = 'celltype')
VlnPlot(bladder_integrated_1,features = "AXL",group.by = 'celltype')
VlnPlot(bladder_integrated_1,features = "ACE2",group.by = 'celltype')
VlnPlot(bladder_integrated_1,features = "TMPRSS2",group.by = 'celltype')
VlnPlot(bladder_integrated_1,features = "FURIN",group.by = 'celltype')
VlnPlot(bladder_integrated_1,features = "CTSL",group.by = 'celltype')

VlnPlot(bladder_integrated_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
DotPlot(bladder_integrated_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
fig_bladder_4 <- DotPlot(bladder_integrated_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
fig_bladder_4

fig_bladder <- CombinePlots(plots = list(fig_bladder_1,fig_bladder_2,fig_bladder_3,fig_bladder_4))
fig_bladder

ggsave("./fig_bladder/A.pdf",plot=fig_bladder_1)
ggsave("./fig_bladder/B.pdf",plot=fig_bladder_2)
ggsave("./fig_bladder/C.pdf",plot=fig_bladder_3)
ggsave("./fig_bladder/D.pdf",plot=fig_bladder_4)

ggsave("./fig_bladder/fig_baldder.pdf",plot=fig_bladder)


####enrich_process####
sce = bladder_integrated_1
mat <- as.matrix(sce@assays$RNA@counts)



high_NRP1=colnames(subset(x = sce, subset = NRP1 > 0, slot = 'counts'))
high_AXL=colnames(subset(x = sce, subset = AXL > 0, slot = 'counts'))
high_FURIN=colnames(subset(x = sce, subset = FURIN > 0, slot = 'counts'))
high_CTSL=colnames(subset(x = sce, subset = CTSL > 0, slot = 'counts'))



highORlow_NRP1=ifelse(colnames(sce) %in% high_NRP1,'high','low')
highORlow_AXL=ifelse(colnames(sce) %in% high_AXL,'high','low')
highORlow_FURIN=ifelse(colnames(sce) %in% high_FURIN,'high','low')
highORlow_CTSL=ifelse(colnames(sce) %in% high_CTSL,'high','low')

table(highORlow_NRP1)
table(highORlow_AXL)
table(highORlow_FURIN)
table(highORlow_CTSL)

sce@meta.data$highORlow_NRP1=highORlow_NRP1
sce@meta.data$highORlow_AXL=highORlow_AXL
sce@meta.data$highORlow_FURIN=highORlow_FURIN
sce@meta.data$highORlow_CTSL=highORlow_CTSL



markers_NRP1 <- FindMarkers(sce, ident.1 = "high", 
                            group.by = 'highORlow_NRP1', 
                            subset.ident = c("Interstitial cells","Myofibroblast","Endothelial cells"))
top100_NRP1 <- rownames(markers_NRP1[1:100,])


markers_AXL <- FindMarkers(sce, ident.1 = "high", 
                           group.by = 'highORlow_AXL', 
                           subset.ident = c("Interstitial cells","Myofibroblast","T cells",
                                            "Monocytes","Smooth muscle cells","Fibroblast"))

top100_AXL <- rownames(markers_AXL[1:100,])


markers_FURIN <- FindMarkers(sce, ident.1 = "high", 
                             group.by = 'highORlow_FURIN', 
                             subset.ident = c("Endothelial cells"))

top100_FURIN <- rownames(markers_FURIN[1:100,])


markers_CTSL <- FindMarkers(sce, ident.1 = "high", 
                            group.by = 'highORlow_CTSL', 
                            subset.ident = c("Interstitial cells","Myofibroblast","T cells",
                                             "Monocytes","Smooth muscle cells","Fibroblast",
                                             "Umbrella cells", "Endothelial cells"
                                             ))

top100_CTSL <- rownames(markers_CTSL[1:100,])


write.table(top100_NRP1,file = "./GSE129845_bladder/top100_NRP1.txt",quote = F,row.names = F)
write.table(top100_AXL,file = "./GSE129845_bladder/top100_AXL.txt",quote = F,row.names = F)
write.table(top100_CTSL,file = "./GSE129845_bladder/top100_CTSL.txt",quote = F,row.names = F)
write.table(top100_FURIN,file = "./GSE129845_bladder/top100_FURIN.txt",quote = F,row.names = F)

####GO, NRP1####
go_NRP1 <- read.table("./GSE129845_bladder/top100_NRP1.txt",sep=" ")

go_NRP1 <- t(go_NRP1)
keytypes(org.Hs.eg.db)

go_NRP1_id_trance <- bitr(go_NRP1, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)
write.table(go_NRP1_id_trance$ENSEMBL,"./GSE129845_bladder/go_NRP1_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE129845_bladder/go_NRP1_id_trance.txt")
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

fig_NRP1_1 <- dotplot(All,split="ONTOLOGY",title ="Enrichment GO Top5",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
fig_NRP1_1

NRP1_go_dat <- All@result
NRP1_go <- NRP1_go_dat[NRP1_go_dat$p.adjust<0.05,"Description"]
NRP1_go
write.table(NRP1_go,"./GSE129845_bladder/NRP1_go.txt",row.names = F,col.names = F)


####KEGG_NRP1####
KEGG <- enrichKEGG(gene= id, organism  = 'hsa', pvalueCutoff = 0.05)
dotplot(KEGG,font.size=12)
barplot(KEGG,font.size=8) 
dotplot(KEGG,showCategory=5,title="Enrichment KEGG Top5")
fig_NRP2<-dotplot(KEGG, font.size=12, showCategory=10, title="Enrichment KEGG Top5")+
  scale_size(rang=c(5.20))
fig_NRP2

NRP1_kegg_dat <- KEGG@result
NRP1_kegg <- NRP1_kegg_dat[NRP1_kegg_dat$p.adjust<0.05,"Description"]

NRP1_kegg
write.table(NRP1_kegg,"./GSE129845_bladder/NRP1_kegg.txt",row.names = F,col.names = F)
####GO, AXL####
go_AXL <- read.table("./GSE129845_bladder/top100_AXL.txt",sep=" ")

go_AXL <- t(go_AXL)
keytypes(org.Hs.eg.db)

go_AXL_id_trance <- bitr(go_AXL, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)
write.table(go_AXL_id_trance$ENSEMBL,"./GSE129845_bladder/go_AXL_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE129845_bladder/go_AXL_id_trance.txt")
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

fig_AXL_1 <- dotplot(All,split="ONTOLOGY",title ="Enrichment GO Top5",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
fig_AXL_1

AXL_go_dat <- All@result
AXL_go <- AXL_go_dat[AXL_go_dat$p.adjust<0.05,"Description"]
AXL_go

####KEGG_AXL####
KEGG <- enrichKEGG(gene= id, organism  = 'hsa', pvalueCutoff = 0.05)
dotplot(KEGG,font.size=12)
barplot(KEGG,font.size=8) 
dotplot(KEGG,showCategory=5,title="Enrichment KEGG Top5")
fig_NRP1_2<-dotplot(KEGG, font.size=12, showCategory=10, title="Enrichment KEGG Top5")+
  scale_size(rang=c(5.20))
fig_NRP1_2

AXL_kegg_dat <- KEGG@result
AXL_kegg <- AXL_kegg_dat[AXL_kegg_dat$p.adjust<0.05,"Description"]

AXL_kegg


write.table(AXL_kegg,"./GSE129845_bladder/AXL_kegg.txt",row.names = F,col.names = F)

####GO, FURIN####
go_FURIN <- read.table("./GSE129845_bladder/top100_FURIN.txt",sep=" ")

go_FURIN <- t(go_FURIN)
keytypes(org.Hs.eg.db)

go_FURIN_id_trance <- bitr(go_FURIN, fromType = "SYMBOL", toType = "ENSEMBL",
                           OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_FURIN_id_trance$ENSEMBL,"./GSE129845_bladder/go_FURIN_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE129845_bladder/go_FURIN_id_trance.txt")
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

fig_FURIN_1 <- dotplot(All,split="ONTOLOGY",title ="Enrichment GO Top5",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
fig_FURIN_1

FURIN_go_dat <- All@result
FURIN_go <- FURIN_go_dat[FURIN_go_dat$p.adjust<0.05,"Description"]
FURIN_go

####KEGG_FURIN####
KEGG <- enrichKEGG(gene= id, organism  = 'hsa', pvalueCutoff = 0.05)
dotplot(KEGG,font.size=12)
barplot(KEGG,font.size=8) 
dotplot(KEGG,showCategory=5,title="Enrichment KEGG Top5")
fig_NRP1_2<-dotplot(KEGG, font.size=12, showCategory=10, title="Enrichment KEGG Top5")+
  scale_size(rang=c(5.20))
fig_NRP1_2

FURIN_kegg_dat <- KEGG@result
FURIN_kegg <- FURIN_kegg_dat[FURIN_kegg_dat$p.adjust<0.05,"Description"]

FURIN_kegg


####GO, CTSL####
go_CTSL <- read.table("./GSE129845_bladder/top100_CTSL.txt",sep=" ")

go_CTSL <- t(go_CTSL)
keytypes(org.Hs.eg.db)

go_CTSL_id_trance <- bitr(go_CTSL, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_CTSL_id_trance$ENSEMBL,"./GSE129845_bladder/go_CTSL_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE129845_bladder/go_CTSL_id_trance.txt")
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

fig_CTSL_1 <- dotplot(All,split="ONTOLOGY",title ="Enrichment GO Top5",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
fig_CTSL_1

CTSL_go_dat <- All@result
CTSL_go <- CTSL_go_dat[CTSL_go_dat$p.adjust<0.05,"Description"]
CTSL_go

####KEGG_CTSL####
KEGG <- enrichKEGG(gene= id, organism  = 'hsa', pvalueCutoff = 0.05)
dotplot(KEGG,font.size=12)
barplot(KEGG,font.size=8) 
dotplot(KEGG,showCategory=5,title="Enrichment KEGG Top5")
fig_NRP1_2<-dotplot(KEGG, font.size=12, showCategory=10, title="Enrichment KEGG Top5")+
  scale_size(rang=c(5.20))
fig_NRP1_2

CTSL_kegg_dat <- KEGG@result
CTSL_kegg <- CTSL_kegg_dat[CTSL_kegg_dat$p.adjust<0.05,"Description"]

CTSL_kegg

NRP1_bladder_go <- NRP1_go
AXL_bladder_go <- AXL_go
FURIN_bladder_go <- FURIN_go
CTSL_bladder_go <- CTSL_go 

NRP1_bladder_kegg <- NRP1_kegg
AXL_bladder_kegg <- AXL_kegg
FURIN_bladder_kegg <- FURIN_kegg
CTSL_bladder_kegg <- CTSL_kegg


save(NRP1_bladder_go,AXL_bladder_go,FURIN_bladder_go,CTSL_bladder_go,NRP1_bladder_kegg,
     AXL_bladder_kegg,FURIN_bladder_kegg,CTSL_bladder_kegg,file = "./GSE129845_bladder/enrichment_bladder.Rdata")


####enrichment result####
#load("./GSE109816_bladder/enrichment_bladder.Rdata")

#Reduce(intersect,list(NRP1_go,AXL_go,FURIN_go,CTSL_go))
#Reduce(intersect,list(NRP1_kegg,AXL_kegg,FURIN_kegg,CTSL_kegg))

#save(NRP1_go,AXL_go,FURIN_go,CTSL_go,NRP1_kegg,AXL_kegg,FURIN_kegg,CTSL_kegg,
#     file = "./GSE109816_bladder/enrichment_bladder.Rdata")

load("./GSE109816_heart/enrichment_heart.Rdata")
load("./GSE115469_liver/enrichment_liver.Rdata")

Reduce(intersect,list(NRP1_heart_go,NRP1_liver_go,NRP1_bladder_go))
Reduce(intersect,list(AXL_heart_go,AXL_liver_go,AXL_bladder_go))
Reduce(intersect,list(FURIN_heart_go,FURIN_bladder_go))
Reduce(intersect,list(CTSL_heart_go,CTSL_liver_go,CTSL_bladder_go))

Reduce(intersect,list(NRP1_heart_kegg,NRP1_liver_kegg,NRP1_bladder_kegg))


#ggsave("./GSE109816_heart/pdf/heart_fig.pdf", plot = fig_1, width = 15, height = 18)

if (FALSE){
  fig_1_1<- DimPlot(object = heart_1, reduction = 'umap',label = T)
  fig_1_1
  
  markers = c( "TNN",
               "RYR2",
               "XIRP2",
               "TECRL",
               "CMYA5",
               "ANKRD1",
               "MYH7",
               "NEBL",
               "MYH6",
               "TNNT2",
               "ACTN2",
               "NEXN",
               "CSRP3",
               "MYL3",
               "TRDN",
               "TPM1",
               
               "VWF",
               "IFI27",
               "PECAM1",
               "HLA-B",
               "CD74",
               
               "DCN",
               "C7",
               "ADH1B",
               "LUM",
               "PTN",
               "C1R",
               "SERPINF1",
               "CDH19",
               "ABCA8",
               
               "PTPRC",
               "LAPTM5",
               "MS4A6A",
               "AIF1",
               "ALOX5AP",
               "TYROBP",
               "FCER1G",
               "CD163",
               "F13A1",
               
               "RGS5",
               "PPP1R14A",
               "HIGD1B",
               "ACTA2")
  
  pp_temp = FeaturePlot(object = heart_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
  pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
  fig_heart_2  <- CombinePlots(plots = plots,ncol = 4,legend = 'right')
  fig_heart_2
}


