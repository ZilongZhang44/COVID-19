####clustering####
library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(AnnotationHub)  
library(org.Hs.eg.db)   

setwd("D:/准备文章/covid19/data/")

liver <- read.csv("./GSE115469_liver/GSE115469_Data.csv",header = T,row.names = 1)
liver<-CreateSeuratObject(counts = liver, project = "liver")

liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^MT-")
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#liver <- subset(liver, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 50)
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 10000)
#liver <- NormalizeData(liver)
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(liver), 10)
plot1 <- VariableFeaturePlot(liver)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
#s.genes <-cc.genes$s.genes
#g2m.genes<-cc.genes$g2m.genes
#liver <- CellCycleScoring(liver, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#all.genes <- rownames(liver)
#liver <- ScaleData(liver, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
liver <- ScaleData(liver,verbose = FALSE, vars.to.regress = c("nCount_RNA","percent.mito"))

#Eliminate batch effects with harmony and cell classification
liver <- RunPCA(liver, pc.genes = liver@var.genes, npcs = 20, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
if (FALSE){
liver <- liver %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(liver, 'harmony')
harmony_embeddings[1:5, 1:5]
liver <- liver %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.25) %>% 
  identity()
}

new.cluster.ids <- c(1, 2, 3, 4, 5)
names(new.cluster.ids) <- levels(liver)
liver <- RenameIdents(liver, new.cluster.ids)

#Calculating differentially expressed genes (DEGs) and Save rds file
liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(liver.markers,sep="\t",file="./GSE115469_liver/Seurat/0.25_20.xls")


liver <- ScaleData(liver, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
liver <- RunPCA(liver, verbose = FALSE,npcs = 100)
liver <- ProjectDim(object = liver)
ElbowPlot(object = liver,ndims = 100)

###cluster
liver <- FindNeighbors(object = liver, dims = 1:25)
liver <- FindClusters(object = liver, resolution = 0.8)

liver <- RunTSNE(object = liver, dims = 1:25)
DimPlot(object = liver, reduction = 'tsne',label = TRUE)


liver <- RunUMAP(object = liver, dims = 1:25)
DimPlot(object = liver, reduction = 'umap',label = TRUE)
p_umap<- DimPlot(object = liver, reduction = 'umap',label = TRUE)
ggsave("./GSE115469_liver/clustering.pdf", plot = p_umap)

#Some visual figure generation
#load("./GSE131685_liver/0.25_20.rds")
DimPlot(liver, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
#DimPlot(liver, reduction = "umap", group.by = "Phase", pt.size = .1)
DimPlot(liver, reduction = "umap", label = TRUE, pt.size = .1)
p_liver <- DimPlot(liver, reduction = "umap", label = TRUE, pt.size = .1)
#saveRDS(liver,file="./GSE115469_liver/0.25_20.rds")
#ggsave("./GSE131685_liver/liver_pdf/clusting.pdf",plot = p_liver)

#DoHeatmap(liver, features = c("SLC13A3","SLC34A1","GPX3","DCXR","SLC17A3","SLC22A8","SLC22A7","GNLY","NKG7","CD3D","CD3E","LYZ","CD14","KRT8","KRT18","CD24","VCAM1","UMOD","DEFB1","CLDN8","AQP2","CD79A","CD79B","ATP6V1G3","ATP6V0D2","TMEM213"))
#VlnPlot(liver, pt.size =0, idents= c(1,2,3), features = c("GPX3", "DCXR","SLC13A3","SLC34A1","SLC22A8","SLC22A7"))
#VlnPlot(liver, idents= c(8,10), features = c("AQP2", "ATP6V1B1","ATP6V0D2","ATP6V1G3"))

####annotation####
liver <- readRDS("./GSE115469_liver/0.25_20.rds")
#Hep tubule
FeaturePlot(liver,features = "CYP3A7")
FeaturePlot(liver,features = "CYP2A7")
FeaturePlot(liver,features = "CYP2A6")

FeaturePlot(liver,features = "SCD")
FeaturePlot(liver,features = "HMGCS1")
FeaturePlot(liver,features = "ACSS2")
FeaturePlot(liver,features = "TM7SF2")

FeaturePlot(liver,features = "SEC16B")
FeaturePlot(liver,features = "SLBP")
FeaturePlot(liver,features = "RND3")
FeaturePlot(liver,features = "PCK1")

FeaturePlot(liver,features = "BCHE")
FeaturePlot(liver,features = "G6PC")
FeaturePlot(liver,features = "GHR")
FeaturePlot(liver,features = "ALDH6A1")

FeaturePlot(liver,features = "RPP25L")
FeaturePlot(liver,features = "HSD11B1")
FeaturePlot(liver,features = "HAMP")
FeaturePlot(liver,features = "GHR")


FeaturePlot(liver,features = "HPR")
FeaturePlot(liver,features = "GSTA2")
FeaturePlot(liver,features = "AKR1C1")
FeaturePlot(liver,features = "MASP2")
#LSEC zone 1 (liver sinusoidal endothelial cells)
FeaturePlot(liver,features = "MGP")
FeaturePlot(liver,features = "SPARCL1")
FeaturePlot(liver,features = "TM4SF1")
FeaturePlot(liver,features = "CLEC14A")

#LSEC zone 2 and 3
FeaturePlot(liver,features = "CCL14")
FeaturePlot(liver,features = "CLEC1B")
FeaturePlot(liver,features = "FCN2")
FeaturePlot(liver,features = "S100A13")


#Portal endothelial
FeaturePlot(liver,features = "RAMP3")
FeaturePlot(liver,features = "INMT")
FeaturePlot(liver,features = "DNASE")
FeaturePlot(liver,features = "LIFR")

#Cholangiocytes
FeaturePlot(liver,features = "KRT7")
FeaturePlot(liver,features = "KRT19")
FeaturePlot(liver,features = "SOX9")
FeaturePlot(liver,features = "EPCAM")

#Stellate cells
FeaturePlot(liver,features = "ACTA2")
FeaturePlot(liver,features = "COL1A1")
FeaturePlot(liver,features = "RBP1")

#Inflammatory Macs
FeaturePlot(liver,features = "S100A8")
FeaturePlot(liver,features = "LYZ")
FeaturePlot(liver,features = "S100A9")
FeaturePlot(liver,features = "HLA-DPB1")

#Non-inflammatory Macs
FeaturePlot(liver,features = "CD5L")
FeaturePlot(liver,features = "MARCO")
FeaturePlot(liver,features = "VSIG4")

#CD3+ αβ T-cells
FeaturePlot(liver,features = "CD2")
FeaturePlot(liver,features = "CD3D")
FeaturePlot(liver,features = "TRAC")
FeaturePlot(liver,features = "GZMK")

#γδ T-cells 1
FeaturePlot(liver,features = "GNLY")
FeaturePlot(liver,features = "PTGDS")
FeaturePlot(liver,features = "GZMB")
FeaturePlot(liver,features = "TRDC")

#γδ T-cells 2
FeaturePlot(liver,features = "STMN1")
FeaturePlot(liver,features = "HMGB2")
FeaturePlot(liver,features = "TYMS")

#Mature B-cells
FeaturePlot(liver,features = "MS4A1")
FeaturePlot(liver,features = "LTB")
FeaturePlot(liver,features = "CD37")
FeaturePlot(liver,features = "CD79B")

#plasma cells
FeaturePlot(liver,features = "IGLC2")
FeaturePlot(liver,features = "IGHG1")
FeaturePlot(liver,features = "IGKC")

#NK-like cells
FeaturePlot(liver,features = "CD7")
FeaturePlot(liver,features = "KLRB1")
FeaturePlot(liver,features = "NKG7")

#erthyroid cells
FeaturePlot(liver,features = "HBB")
FeaturePlot(liver,features = "CA1")
FeaturePlot(liver,features = "ALAS2")

liver_1 <- RenameIdents(object = liver,
                      '0'='Hepatocytes',
                      '1'='Hepatocytes',
                      '2'='Hepatocytes',
                      '3'='CD3+ αβ T-cells',
                      '4'='Plasma cells',
                      '5'='CD3+ αβ T-cells',
                      '6'='γδ T-cells',
                      '7'='Inflammatory Macs',
                      '8'='CD3+ αβ T-cells',
                      '9'='Non-inflammatory Macs',
                      '10'='Inflammatory Macs',
                      '11'= 'LSEC cells',
                      '12'= 'LSEC cells',
                      '13'= 'LSEC cells',
                      '14'='B-cells',
                      '15'='Cholangiocytes',
                      '16'='Erthyroid cells',
                      '17'='Hepatocytes',
                      '18'='Hepatocytes',
                      '19' ='γδ T-cells',
                      '20' ='Hepatocytes',
                      '21'='Stellate cells'
                        )

fig_liver_1 <- DimPlot(object = liver, reduction = 'umap',label = T)
fig_liver_1

fig_liver_2<- DimPlot(object = liver_1, reduction = 'umap',label = T)
fig_liver_2
ggsave(filename = "./GSE115469_liver/pdf/liver_annotation.pdf",plot = fig_liver_2)


if (FALSE){

new_order = c('0','1','2','3','4','5','7','8','9','11','12','16','17','20','21')
organ.summary.dataframe$cluster = factor(organ.summary.dataframe$cluster,levels = new_order,labels = new_order)
cols = c('#32b8ec','#60c3f0','#8ccdf1','#cae5f7','#92519c','#b878b0','#d7b1d2','#e7262a','#e94746','#eb666d','#ee838f','#f4abac','#fad9d9')
pp = ggplot(data=organ.summary.dataframe, aes(x=cluster, y=cell, fill=group)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Fraction of sample per cluster (%)') +
  theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
        legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = cols)
pp
}
if (FALSE){
markers = c("CYP3A7","CYP2A7","CYP2A6","SCD","HMGCS1","ACSS2","TM7SF2","SEC16B","SLBP","RND3","PCK1","BCHE",
            "G6PC","GHR","ALDH6A1","RPP25L","HSD11B1","HAMP","GHR","HPR","GSTA2","AKR1C1","MASP2","MGP",
            "SPARCL1","TM4SF1","CLEC14A","CCL14","CLEC1B","FCN2","S100A13","RAMP3",
            "INMT","DNASE","LIFR","KRT7","KRT19","SOX9","EPCAM","ACTA2","COL1A1","RBP1","S100A8","LYZ","S100A9","HLA-DPB1","CD5L","MARCO","VSIG4",
            "CD2","CD3D","TRAC","GZMK","GNLY","PTGDS","GZMB","TRDC","STMN1","HMGB2","TYMS","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC",
            "CD7","KLRB1","NKG7","HBB","CA1","ALAS2")
}

markers = c("BCHE","TM4SF1","KRT7","ACTA2","LYZ","CD5L","GZMK","GNLY","MS4A1",
            "IGHG1","CA1")

pp_temp = FeaturePlot(object = liver_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,legend = 'right')
pp
fig_liver_3 <- pp


pp_dot = DotPlot(liver_1, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp_dot = pp_dot + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + theme(axis.line = element_line(size = 0.6))

pp_dot
supplement_liver_1 <- pp_dot

supplement_liver_2 <- DoHeatmap(liver_1, features = markers)
supplement_liver_2
#save(liver_1,file = "./GSE115469_liver/annotation.Rdata")



####downstream####
liver <- readRDS("./GSE115469_liver/0.25_20.rds")
fig_liver_1 <- DimPlot(object = liver, reduction = 'umap',label = T)
fig_liver_1




load("./GSE115469_liver/annotation.Rdata")
liver_1$celltype <- Idents(liver_1)

fig_liver_2<- DimPlot(object = liver_1, reduction = 'umap',label = T)
fig_liver_2


markers = c("BCHE","TM4SF1","KRT7","ACTA2","LYZ","CD5L","GZMK","GNLY","MS4A1",
            "IGHG1","CA1")

pp_temp = FeaturePlot(object = liver_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12),
                        plot.title = element_text(family = 'sans',face='italic',size=16),
                        legend.text = element_text(size = 16),legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
pp
fig_liver_3 <- pp




VlnPlot(liver_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
#fig_1_3 <- VlnPlot(liver_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"),group.by = 'celltype')
DotPlot(liver_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
fig_liver_4 <- DotPlot(liver_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
fig_liver_4
fig_4 <- CombinePlots(plots = list(fig_liver_1,fig_liver_2,fig_liver_3,fig_liver_4))
fig_4

ggsave("./fig_liver/A.pdf",plot=fig_liver_1)
ggsave("./fig_liver/B.pdf",plot=fig_liver_2)
ggsave("./fig_liver/C.pdf",plot=fig_liver_3)
ggsave("./fig_liver/D.pdf",plot=fig_liver_4)

ggsave("./fig_liver/fig_4.pdf",plot=fig_4)


ggsave("./GSE131685_liver/liver_pdf/liver_fig.pdf", plot = fig_1, width = 15, height = 18)

####enrich_process####
sce = liver_1
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
                            subset.ident = c('Stellate cells','Cholangiocytes',
                                             'LSEC cells','Non-inflammatory Macs'))
top100_NRP1 <- rownames(markers_NRP1[1:100,])


markers_AXL <- FindMarkers(sce, ident.1 = "high", 
                           group.by = 'highORlow_AXL', 
                           subset.ident = c('Stellate cells','Non-inflammatory Macs'))

top100_AXL <- rownames(markers_AXL[1:100,])


markers_FURIN <- FindMarkers(sce, ident.1 = "high", 
                             group.by = 'highORlow_FURIN', 
                             subset.ident = c("Hepatocytes"))

top100_FURIN <- rownames(markers_FURIN[1:100,])


markers_CTSL <- FindMarkers(sce, ident.1 = "high", 
                            group.by = 'highORlow_CTSL', 
                            subset.ident = c('Stellate cells','Erthyroid cells',
                                             'Cholangiocytes',
                                             'LSEC cells','Non-inflammatory Macs',
                                             'Inflammatory Macs','Hepatocytes'))

top100_CTSL <- rownames(markers_CTSL[1:100,])



#write.table(top100_NRP1,file = "./GSE115469_liver/top100_NRP1.txt",quote = F,row.names = F)
#write.table(top100_AXL,file = "./GSE115469_liver/top100_AXL.txt",quote = F,row.names = F)
#write.table(top100_CTSL,file = "./GSE115469_liver/top100_CTSL.txt",quote = F,row.names = F)
#write.table(top100_FURIN,file = "./GSE115469_liver/top100_FURIN.txt",quote = F,row.names = F)

####GO, NRP1####
go_NRP1 <- read.table("./GSE115469_liver/top100_NRP1.txt",sep=" ")

go_NRP1 <- t(go_NRP1)
keytypes(org.Hs.eg.db)

go_NRP1_id_trance <- bitr(go_NRP1, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)
write.table(go_NRP1_id_trance$ENSEMBL,"./GSE115469_liver/go_NRP1_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE115469_liver/go_NRP1_id_trance.txt")
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
#write.table(NRP1_go,"./GSE115469_liver/NRP1_go.txt",row.names = F,col.names = F)


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
write.table(NRP1_kegg,"./GSE115469_liver/NRP1_kegg.txt",row.names = F,col.names = F)



####GO, AXL####
go_AXL <- read.table("./GSE115469_liver/top100_AXL.txt",sep=" ")

go_AXL <- t(go_AXL)
keytypes(org.Hs.eg.db)

go_AXL_id_trance <- bitr(go_AXL, fromType = "SYMBOL", toType = "ENSEMBL",
                         OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_AXL_id_trance$ENSEMBL,"./GSE115469_liver/go_AXL_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE115469_liver/go_AXL_id_trance.txt")
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


#write.table(AXL_kegg,"./AXL_kegg.txt",row.names = F,col.names = F)


####GO, CTSL####
go_CTSL <- read.table("./GSE115469_liver/top100_CTSL.txt",sep=" ")

go_CTSL <- t(go_CTSL)
keytypes(org.Hs.eg.db)

go_CTSL_id_trance <- bitr(go_CTSL, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_CTSL_id_trance$ENSEMBL,"./GSE115469_liver/go_CTSL_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE115469_liver/go_CTSL_id_trance.txt")
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












####enrichment####
NRP1_liver_go <- NRP1_go
AXL_liver_go <- AXL_go
CTSL_liver_go <- CTSL_go

NRP1_liver_kegg <- NRP1_kegg
AXL_liver_kegg <- AXL_kegg
CTSL_liver_kegg <- CTSL_kegg

save(NRP1_liver_go,AXL_liver_go,CTSL_liver_go,NRP1_liver_kegg,AXL_liver_kegg,CTSL_liver_kegg,
     file = "./GSE115469_liver/enrichment_liver.Rdata")

intersect(NRP1_heart_go,
         NRP1_liver_go)

