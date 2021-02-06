####load####
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

heart.data <- read.csv("./GSE109816_heart/GSE109816_normal_heart_umi_matrix.csv",header = T,row.names = 1)
heart<-CreateSeuratObject(counts = heart.data, project = "heart",min.cells = 3, min.features = 500)

heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^MT-")
VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
heart <- subset(heart, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 72)
heart <- NormalizeData(heart, normalization.method = "LogNormalize", scale.factor = 10000)
heart <- NormalizeData(heart)
heart <- FindVariableFeatures(heart, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(heart), 10)
plot1 <- VariableFeaturePlot(heart)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


#s.genes <-cc.genes$s.genes
#g2m.genes<-cc.genes$g2m.genes
#heart <- CellCycleScoring(heart, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#all.genes <- rownames(heart)
#heart <- ScaleData(heart, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
heart <- ScaleData(heart,verbose = FALSE, vars.to.regress = c("nCount_RNA","percent.mito"))

#Eliminate batch effects with harmony and cell classification
heart <- RunPCA(heart, pc.genes = heart@var.genes, npcs = 20, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
if (FALSE){
  heart <- heart %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
  harmony_embeddings <- Embeddings(heart, 'harmony')
  harmony_embeddings[1:5, 1:5]
  heart <- heart %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.25) %>% 
    identity()


new.cluster.ids <- c(1, 2, 3, 4, 5)
names(new.cluster.ids) <- levels(heart)
heart <- RenameIdents(heart, new.cluster.ids)
}
#Calculating differentially expressed genes (DEGs) and Save rds file
heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(heart.markers,sep="\t",file="./GSE109816_heart/Seurat/0.25_20.xls")


heart <- ScaleData(heart, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
heart <- RunPCA(heart, verbose = FALSE,npcs = 100)
heart <- ProjectDim(object = heart)
ElbowPlot(object = heart,ndims = 100)

###cluster
heart <- FindNeighbors(object = heart, dims = 1:25)
heart <- FindClusters(object = heart, resolution = 0.8)

heart <- RunTSNE(object = heart, dims = 1:25)
DimPlot(object = heart, reduction = 'tsne',label = TRUE)


heart <- RunUMAP(object = heart, dims = 1:25)
DimPlot(object = heart, reduction = 'umap',label = TRUE)
p_umap<- DimPlot(object = heart, reduction = 'umap',label = TRUE)
ggsave("./GSE109816_heart/clustering.pdf", plot = p_umap)

#Some visual figure generation
#load("./GSE131685_heart/0.25_20.rds")
#DimPlot(heart, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
#DimPlot(heart, reduction = "umap", group.by = "Phase", pt.size = .1)
DimPlot(heart, reduction = "umap", label = TRUE, pt.size = .1)
p_heart <- DimPlot(heart, reduction = "umap", label = TRUE, pt.size = .1)
#saveRDS(heart,file="./GSE109816_heart/0.25_20.rds")
#ggsave("./GSE131685_heart/heart_pdf/clusting.pdf",plot = p_heart)

#DoHeatmap(heart, features = c("SLC13A3","SLC34A1","GPX3","DCXR","SLC17A3","SLC22A8","SLC22A7","GNLY","NKG7","CD3D","CD3E","LYZ","CD14","KRT8","KRT18","CD24","VCAM1","UMOD","DEFB1","CLDN8","AQP2","CD79A","CD79B","ATP6V1G3","ATP6V0D2","TMEM213"))
#VlnPlot(heart, pt.size =0, idents= c(1,2,3), features = c("GPX3", "DCXR","SLC13A3","SLC34A1","SLC22A8","SLC22A7"))
#VlnPlot(heart, idents= c(8,10), features = c("AQP2", "ATP6V1B1","ATP6V0D2","ATP6V1G3"))

####annotation####
heart <- readRDS("./GSE109816_heart/0.25_20.rds")
#CM
FeaturePlot(heart,features = "TNN")
FeaturePlot(heart,features = "RYR2")
FeaturePlot(heart,features = "XIRP2")
FeaturePlot(heart,features = "TECRL")
FeaturePlot(heart,features = "CMYA5")
FeaturePlot(heart,features = "ANKRD1")
FeaturePlot(heart,features = "MYH7")
FeaturePlot(heart,features = "NEBL")
FeaturePlot(heart,features = "MYH6")
FeaturePlot(heart,features = "TNNT2")
FeaturePlot(heart,features = "ACTN2")
FeaturePlot(heart,features = "NEXN")
FeaturePlot(heart,features = "CSRP3")
FeaturePlot(heart,features = "MYL3")
FeaturePlot(heart,features = "TRDN")
FeaturePlot(heart,features = "TPM1")



#EC
FeaturePlot(heart,features = "VWF")
FeaturePlot(heart,features = "IFI27")
FeaturePlot(heart,features = "PECAM1")
FeaturePlot(heart,features = "HLA-B")
FeaturePlot(heart,features = "CD74")


#FB
FeaturePlot(heart,features = "DCN")
FeaturePlot(heart,features = "C7")
FeaturePlot(heart,features = "ADH1B")
FeaturePlot(heart,features = "LUM")
FeaturePlot(heart,features = "PTN")
FeaturePlot(heart,features = "C1R")
FeaturePlot(heart,features = "SERPINF1")
FeaturePlot(heart,features = "CDH19")
FeaturePlot(heart,features = "ABCA8")

#MP
FeaturePlot(heart,features = "PTPRC")
FeaturePlot(heart,features = "LAPTM5")
FeaturePlot(heart,features = "MS4A6A")
FeaturePlot(heart,features = "AIF1")
FeaturePlot(heart,features = "ALOX5AP")
FeaturePlot(heart,features = "TYROBP")
FeaturePlot(heart,features = "FCER1G")
FeaturePlot(heart,features = "CD163")
FeaturePlot(heart,features = "F13A1")

#SMC
FeaturePlot(heart,features = "RGS5")
FeaturePlot(heart,features = "PPP1R14A")
FeaturePlot(heart,features = "HIGD1B")
FeaturePlot(heart,features = "ACTA2")





heart_1 <- RenameIdents(object = heart,
                        '0'='Endothelial cells',
                        '1'='Cardiomyocytes',
                        '2'='Cardiomyocytes',
                        '3'='Cardiomyocytes',
                        '4'='Smooth muscle cells',
                        '5'='Endothelial cells',
                        '6'='Endothelial cells',
                        '7'='Fibroblasts',
                        '8'='Smooth muscle cells',
                        '9'='Cardiomyocytes',
                        '10'='Macrophages',
                        '11'='Macrophages',
                        '12'='Endothelial cells',
                        '13'='Macrophages',
                        '14'='Endothelial cells'
                    
)




fig_heart_1<- DimPlot(object = heart_1, reduction = 'umap',label = T)
fig_heart_1

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
}
markers = c( 
  "TECRL",
  "VWF",
  "LUM",
  "LAPTM5",
  "RGS5")
pp_temp = FeaturePlot(object = heart_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
fig_heart_2  <- CombinePlots(plots = plots,legend = 'right')
fig_heart_2




pp_temp = FeaturePlot(object = heart_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
#pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
#pp

pp_dot = DotPlot(heart_1, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp_dot = pp_dot + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + theme(axis.line = element_line(size = 0.6))

pp_dot
supplement_heart_1 <- pp_dot
ggsave("./GSE109816_heart/pdf/supplement_heart_1.pdf", 
       plot = supplement_heart_1, width = 15, height = 18)

supplyment_heart_2 <- DoHeatmap(heart, features = markers)
supplyment_heart_2


#save(heart_1,file = "./GSE109816_heart/annotation.Rdata")

####downstream####
heart <- readRDS("./GSE109816_heart/0.25_20.rds")
fig_heart_1<- DimPlot(object = heart, reduction = 'umap',label = T)
fig_heart_1


load("./GSE109816_heart/annotation.Rdata")
heart_1$celltype <- Idents(heart_1)

fig_heart_2<- DimPlot(object = heart_1, reduction = 'umap',label = T)
fig_heart_2

markers = c( 
  "TECRL",
  "VWF",
  "LUM",
  "LAPTM5",
  "RGS5")
pp_temp = FeaturePlot(object = heart_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12),
                        plot.title = element_text(family = 'sans',face='italic',size=16),
                        legend.text = element_text(size = 16),legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
fig_heart_3  <- CombinePlots(plots = plots,legend = 'right')
fig_heart_3


#"SR-B1"
VlnPlot(heart_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
#fig_1_3 <- VlnPlot(heart_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"),group.by = 'celltype')
DotPlot(heart_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
fig_heart_4 <- DotPlot(heart_1,features = c("ACE2","TMPRSS2","NRP1","AXL","FURIN","CTSL"),group.by = 'celltype')
fig_heart_4
fig_3 <- CombinePlots(plots = list(fig_heart_1,fig_heart_2,fig_heart_3,fig_heart_4))

ggsave("./fig_heart/A.pdf",plot=fig_heart_1)
ggsave("./fig_heart/B.pdf",plot=fig_heart_2)
ggsave("./fig_heart/C.pdf",plot=fig_heart_3)
ggsave("./fig_heart/D.pdf",plot=fig_heart_4)

ggsave("./fig_3/fig_3.pdf",plot=fig_3)
ggsave("./GSE109816_heart/pdf/fig_3.pdf", plot = fig_3, width = 15, height = 18)



####enrich_process####
sce = heart_1
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
                       subset.ident = c("SMC","FB"))
top100_NRP1 <- rownames(markers_NRP1[1:100,])


markers_AXL <- FindMarkers(sce, ident.1 = "high", 
                            group.by = 'highORlow_AXL', 
                            subset.ident = c("SMC","FB"))

top100_AXL <- rownames(markers_AXL[1:100,])

markers_CTSL <- FindMarkers(sce, ident.1 = "high", 
                           group.by = 'highORlow_CTSL', 
                           subset.ident = c("SMC","FB"))

top100_CTSL <- rownames(markers_CTSL[1:100,])


markers_FURIN <- FindMarkers(sce, ident.1 = "high", 
                             group.by = 'highORlow_FURIN', 
                             subset.ident = c("SMC","FB"))

top100_FURIN <- rownames(markers_FURIN[1:100,])

#write.table(top100_NRP1,file = "top100_NRP1.txt",quote = F,row.names = F)
#write.table(top100_AXL,file = "top100_AXL.txt",quote = F,row.names = F)
#write.table(top100_CTSL,file = "top100_CTSL.txt",quote = F,row.names = F)
#write.table(top100_FURIN,file = "top100_FURIN.txt",quote = F,row.names = F)

####GO, NRP1####
go_NRP1 <- read.table("./top100_NRP1.txt",sep=" ")

go_NRP1 <- t(go_NRP1)
keytypes(org.Hs.eg.db)

go_NRP1_id_trance <- bitr(go_NRP1, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)
#write.table(go_NRP1_id_trance,"./go_NRP1_id_trance.txt")

f <- read.table("./genelist_NRP1.txt")
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
#write.table(NRP1_go,"./NRP1_go.txt",row.names = F,col.names = F)


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
write.table(NRP1_kegg,"./NRP1_kegg.txt",row.names = F,col.names = F)
####GO, AXL####
go_AXL <- read.table("./top100_AXL.txt",sep=" ")

go_AXL <- t(go_AXL)
keytypes(org.Hs.eg.db)

go_AXL_id_trance <- bitr(go_AXL, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_AXL_id_trance$ENSEMBL,"./go_AXL_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./genelist_AXL.txt")
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

####GO, FURIN####
go_FURIN <- read.table("./top100_FURIN.txt",sep=" ")

go_FURIN <- t(go_FURIN)
keytypes(org.Hs.eg.db)

go_FURIN_id_trance <- bitr(go_FURIN, fromType = "SYMBOL", toType = "ENSEMBL",
                         OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_FURIN_id_trance$ENSEMBL,"./go_FURIN_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./go_FURIN_id_trance.txt")
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
go_CTSL <- read.table("./top100_CTSL.txt",sep=" ")

go_CTSL <- t(go_CTSL)
keytypes(org.Hs.eg.db)

go_CTSL_id_trance <- bitr(go_CTSL, fromType = "SYMBOL", toType = "ENSEMBL",
                           OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_CTSL_id_trance$ENSEMBL,"./go_CTSL_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./go_CTSL_id_trance.txt")
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

####enrichment result####
load("./GSE109816_heart/enrichment_heart.Rdata")
Reduce(intersect,list(NRP1_go,AXL_go,FURIN_go,CTSL_go))
Reduce(intersect,list(NRP1_kegg,AXL_kegg,FURIN_kegg,CTSL_kegg))

#save(NRP1_go,AXL_go,FURIN_go,CTSL_go,NRP1_kegg,AXL_kegg,FURIN_kegg,CTSL_kegg,
#     file = "./GSE109816_heart/enrichment_heart.Rdata")

load("./GSE109816_heart/enrichment_heart.Rdata")
NRP1_heart_go <- NRP1_go
AXL_heart_go <- AXL_go
FURIN_heart_go <- FURIN_go
CTSL_heart_go <- CTSL_go 

NRP1_heart_kegg <- NRP1_kegg
AXL_heart_kegg <- AXL_kegg
FURIN_heart_kegg <- FURIN_kegg
CTSL_heart_kegg <- CTSL_kegg

save(NRP1_heart_go,AXL_heart_go,FURIN_heart_go,CTSL_heart_go,NRP1_heart_kegg,AXL_heart_kegg,FURIN_heart_kegg,CTSL_heart_kegg,
     file = "./GSE109816_heart/enrichment_heart.Rdata")

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