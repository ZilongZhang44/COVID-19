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
library(cowplot)


setwd("D:/准备文章/covid19/data/")

kidney_1<- Read10X(data.dir = "./GSE131685_kidney/1/")
kidney_2<- Read10X(data.dir = "./GSE131685_kidney/2/")
kidney_3<- Read10X(data.dir = "./GSE131685_kidney/3/")

kidney_1<-CreateSeuratObject(counts = kidney_1, project = "kidney_1",min.cells = 8, min.features = 200)
kidney_2<-CreateSeuratObject(counts = kidney_2, project = "kidney_2",min.cells = 6, min.features = 200)
kidney_3<-CreateSeuratObject(counts = kidney_3, project = "kidney_3",min.cells = 10, min.features = 200)

kid <- merge(x = kidney_1, y = list(kidney_2, kidney_3))

kid[["percent.mt"]] <- PercentageFeatureSet(kid, pattern = "^MT-")
VlnPlot(kid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
kid <- subset(kid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 50)
kid <- NormalizeData(kid, normalization.method = "LogNormalize", scale.factor = 10000)
#kid <- NormalizeData(kid)
kid <- FindVariableFeatures(kid, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(kid), 10)
plot1 <- VariableFeaturePlot(kid)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
kid <- CellCycleScoring(kid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all.genes <- rownames(kid)
kid <- ScaleData(kid, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)

#Eliminate batch effects with harmony and cell classification
kid <- RunPCA(kid, pc.genes = kid@var.genes, npcs = 20, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
kid <- kid %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(kid, 'harmony')
harmony_embeddings[1:5, 1:5]
kid <- kid %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.25) %>% 
  identity()
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7,8,9,10,11)
names(new.cluster.ids) <- levels(kid)
kid <- RenameIdents(kid, new.cluster.ids)

#Calculating differentially expressed genes (DEGs) and Save rds file
kid.markers <- FindAllMarkers(kid, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#write.table(kid.markers,sep="\t",file="./GSE131685_kidney/seurat/0.2_20.xls")
saveRDS(kid,file="./GSE131685_kidney/kid.rds")

#Some visual figure generation
#load("./GSE131685_kidney/kid.rds")
DimPlot(kid, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
DimPlot(kid, reduction = "umap", group.by = "Phase", pt.size = .1)
DimPlot(kid, reduction = "umap", label = TRUE, pt.size = .1)
p_kidney <- DimPlot(kid, reduction = "umap", label = TRUE, pt.size = .1)
p_kidney
fig_kidney_1 <- p_kidney
ggsave("./GSE131685_kidney/kidney_pdf/clusting.pdf",plot = p_kidney)

DoHeatmap(kid, features = c("SLC13A3","SLC34A1","GPX3","DCXR","SLC17A3","SLC22A8","SLC22A7","GNLY","NKG7","CD3D","CD3E","LYZ","CD14","KRT8","KRT18","CD24","VCAM1","UMOD","DEFB1","CLDN8","AQP2","CD79A","CD79B","ATP6V1G3","ATP6V0D2","TMEM213"))
VlnPlot(kid, pt.size =0, idents= c(1,2,3), features = c("GPX3", "DCXR","SLC13A3","SLC34A1","SLC22A8","SLC22A7"))
VlnPlot(kid, idents= c(8,10), features = c("AQP2", "ATP6V1B1","ATP6V0D2","ATP6V1G3"))

####annotation####
load("./GSE131685_kidney/kid.rds")
#Proximal tubule
FeaturePlot(kid,features = "SLC13A3")
FeaturePlot(kid,features = "SLC34A1")
FeaturePlot(kid,features = "GPX3")
FeaturePlot(kid,features = "DCXR")

#Proximal convoluted tubule
FeaturePlot(kid,features = "SLC22A8")

#Proximal straight tubule
FeaturePlot(kid,features = "SLC22A7")

#Glomerular parietal epithelial cells
FeaturePlot(kid,features = "KRT8")
FeaturePlot(kid,features = "KRT18")
FeaturePlot(kid,features = "CD24")
FeaturePlot(kid,features = "VCAM1")

#Monocytes
FeaturePlot(kid,features = "LYZ")
FeaturePlot(kid,features = "CD14")

#NK cells
FeaturePlot(kid,features = "GNLY")
FeaturePlot(kid,features = "NKG7")

#T cells
FeaturePlot(kid,features = "CD3D")
FeaturePlot(kid,features = "CD3E")
FeaturePlot(kid,features = "IL7R")

#B cells
FeaturePlot(kid,features = "CD79A")
FeaturePlot(kid,features = "CD79B")

#Distal tubule
FeaturePlot(kid,features = "UMOD")
FeaturePlot(kid,features = "DEFB1")

#Collecting duct
FeaturePlot(kid,features = "CLDN8")

#Collecting duct principal cells
FeaturePlot(kid,features = "AQP2")

#Collecting duct intercalated cells
FeaturePlot(kid,features = "ATP6V1G3")
FeaturePlot(kid,features = "ATP6V0D2")
FeaturePlot(kid,features = "TMEM213")

kid_1 <- RenameIdents(object = kid,
                      '1'='Proximal convoluted tubule cells',
                      '2'='Proximal tubule cells',
                      '3'='Proximal tubule cells',
                      '4'='Proximal straight tubule cells',
                      '5'='NK-T cells',
                      '6'='Monocytes',
                      '7'='Glomerular parietal epithelial cells',
                      '8'='Distal tubule cells',
                      '9'='Collecting duct cells',
                      '10'='B cells',
                      '11'='Collecting duct intercalated cells')
saveRDS(kid_1,file="./GSE131685_kidney/kid_1.rds")

fig_kidney_2<- DimPlot(object = kid_1, reduction = 'umap',label = T,repel = T)
fig_kidney_2
ggsave(filename = "./GSE131685_kidney/kidney_pdf/kidney_annotation.pdf",plot = fig_kidney_2)


markers <- c( "LYZ","GNLY","CD3D","CD79A","UMOD","AQP2",
              "ATP6V1G3","KRT8","SLC22A7","SLC22A8","DCXR")


if (FALSE){
markers = c('SLC13A3',
            'SLC34A1',
            'GPX3',
            'DCXR',
            'SLC22A8',
            'SLC22A7',
            'KRT8',
            'KRT18',
            'CD24',
            'VCAM1',
            'LYZ',
            'CD14',
            'GNLY',
            'NKG7',
            'CD3D',
            'CD3E',
            'IL7R',
            'CD79A',
            'CD79B',
            'UMOD',
            'DEFB1',
            'CLDN8',
            'AQP2',
            'ATP6V1G3',
            'ATP6V0D2',
            'TMEM213'
            
)
}

pp_temp = FeaturePlot(object = kid_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(family = 'sans',face='italic',size=20),legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
pp
fig_kidney_3 <- pp
fig_kidney_3

pp_dot = DotPlot(kid_1, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp_dot = pp_dot + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + theme(axis.line = element_line(size = 0.6))

pp_dot
fig_1_2 <- pp_dot

supplyment_kid_1 <- DoHeatmap(kid, features = markers)
#save(kid_1,file = "./GSE131685_kidney/annotation.Rdata")

####downstream####
load("./GSE131685_kidney/kid.rds")
load("./GSE131685_kidney/kid_1.rds")

kid_1$celltype = Idents(kid_1)

p_kidney <- DimPlot(kid, reduction = "umap", label = TRUE, pt.size = .1)
p_kidney
fig_kidney_1 <- p_kidney


fig_kidney_2<- DimPlot(object = kid_1, reduction = 'umap',label = T,repel = T)
fig_kidney_2


markers <- c( "LYZ","GNLY","CD3D","CD79A","UMOD","AQP2",
              "ATP6V1G3","KRT8","SLC22A7","SLC22A8","DCXR")

pp_temp = FeaturePlot(object = kid_1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 12),
                        plot.title = element_text(family = 'sans',face='italic',size=16),
                        legend.text = element_text(size = 16),legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
pp
fig_kidney_3 <- pp
fig_kidney_3


fig_kidney_4 <- DotPlot(kid_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"),group.by = 'celltype')
fig_kidney_4


fig_supplymentary_1 <- plot_grid(VlnPlot(kid_1, features=c("ACE2"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
                          VlnPlot(kid_1, features=c("TMPRSS2"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
                          VlnPlot(kid_1, features=c("NRP1"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
                          VlnPlot(kid_1, features=c("FURIN"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
                          VlnPlot(kid_1, features=c("CTSL"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip())
fig_supplymentary_1




#p1 <- VlnPlot(kid_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"))
#p1+coord_flip()
#p1

#plot_grid(VlnPlot(kid_1, features=c("ACE2"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
#          VlnPlot(kid_1, features=c("TMPRSS2"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
#          VlnPlot(kid_1, features=c("NRP1"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
#          VlnPlot(kid_1, features=c("FURIN"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip(),
#          VlnPlot(kid_1, features=c("CTSL"))+theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)+coord_flip())



#RidgePlot(kid_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"))


fig_kidney <- CombinePlots(plots = list(fig_kidney_1,fig_kidney_2,fig_kidney_3,fig_kidney_4))

ggsave("./fig_kidney/A.pdf",plot=fig_kidney_1)
ggsave("./fig_kidney/B.pdf",plot=fig_kidney_2)
ggsave("./fig_kidney/C.pdf",plot=fig_kidney_3)
ggsave("./fig_kidney/D.pdf",plot=fig_kidney_4)
ggsave("./fig_kidney/E.pdf",plot=fig_kidney_5,width = 15, height = 18)
ggsave("./fig_kidney/fig_kidney.pdf",plot=fig_kidney)
####enrich_process####
sce = kid_1
mat <- as.matrix(sce@assays$RNA@counts)



high_NRP1=colnames(subset(x = sce, subset = NRP1 > 0, slot = 'counts'))

high_FURIN=colnames(subset(x = sce, subset = FURIN > 0, slot = 'counts'))
high_CTSL=colnames(subset(x = sce, subset = CTSL > 0, slot = 'counts'))



highORlow_NRP1=ifelse(colnames(sce) %in% high_NRP1,'high','low')

highORlow_FURIN=ifelse(colnames(sce) %in% high_FURIN,'high','low')
highORlow_CTSL=ifelse(colnames(sce) %in% high_CTSL,'high','low')

table(highORlow_NRP1)

table(highORlow_FURIN)
table(highORlow_CTSL)

sce@meta.data$highORlow_NRP1=highORlow_NRP1

sce@meta.data$highORlow_FURIN=highORlow_FURIN
sce@meta.data$highORlow_CTSL=highORlow_CTSL



markers_NRP1 <- FindMarkers(sce, ident.1 = "high", 
                            group.by = 'highORlow_NRP1', 
                            subset.ident = c("Proximal straight tubule cells"))
top100_NRP1 <- rownames(markers_NRP1[1:100,])





markers_FURIN <- FindMarkers(sce, ident.1 = "high", 
                             group.by = 'highORlow_FURIN', 
                             subset.ident = c("Collecting duct cells"))

top100_FURIN <- rownames(markers_FURIN[1:100,])


markers_CTSL <- FindMarkers(sce, ident.1 = "high", 
                            group.by = 'highORlow_CTSL', 
                            subset.ident = c("Collecting duct cells","Monocytes","Proximal straight tubule cells",
                                             "Proximal tubule cells","Proximal convoluted tubule cells"
                            ))

top100_CTSL <- rownames(markers_CTSL[1:100,])


write.table(top100_NRP1,file = "./GSE131685_kidney/top100_NRP1.txt",quote = F,row.names = F)

write.table(top100_CTSL,file = "./GSE131685_kidney/top100_CTSL.txt",quote = F,row.names = F)
write.table(top100_FURIN,file = "./GSE131685_kidney/top100_FURIN.txt",quote = F,row.names = F)




####GO, NRP1####
go_NRP1 <- read.table("./GSE131685_kidney/top100_NRP1.txt",sep=" ")

go_NRP1 <- t(go_NRP1)
keytypes(org.Hs.eg.db)

go_NRP1_id_trance <- bitr(go_NRP1, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)
write.table(go_NRP1_id_trance$ENSEMBL,"./GSE131685_kidney/go_NRP1_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE131685_kidney/go_NRP1_id_trance.txt")
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
write.table(NRP1_go,"./GSE131685_kidney/NRP1_go.txt",row.names = F,col.names = F)


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
write.table(NRP1_kegg,"./GSE131685_kidney/NRP1_kegg.txt",row.names = F,col.names = F)

####GO, FURIN####
go_FURIN <- read.table("./GSE131685_kidney/top100_FURIN.txt",sep=" ")

go_FURIN <- t(go_FURIN)
keytypes(org.Hs.eg.db)

go_FURIN_id_trance <- bitr(go_FURIN, fromType = "SYMBOL", toType = "ENSEMBL",
                           OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_FURIN_id_trance$ENSEMBL,"./GSE131685_kidney/go_FURIN_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE131685_kidney/go_FURIN_id_trance.txt")
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
go_CTSL <- read.table("./GSE131685_kidney/top100_CTSL.txt",sep=" ")

go_CTSL <- t(go_CTSL)
keytypes(org.Hs.eg.db)

go_CTSL_id_trance <- bitr(go_CTSL, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)

write.table(go_CTSL_id_trance$ENSEMBL,"./GSE131685_kidney/go_CTSL_id_trance.txt",
            row.names = F,col.names = F)

f <- read.table("./GSE131685_kidney/go_CTSL_id_trance.txt")
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

NRP1_kidney_go <- NRP1_go

FURIN_kidney_go <- FURIN_go
CTSL_kidney_go <- CTSL_go 

NRP1_kidney_kegg <- NRP1_kegg

FURIN_kidney_kegg <- FURIN_kegg
CTSL_kidney_kegg <- CTSL_kegg


save(NRP1_kidney_go,FURIN_kidney_go,CTSL_kidney_go,NRP1_kidney_kegg,
     FURIN_kidney_kegg,CTSL_kidney_kegg,file = "./GSE131685_kidney/enrichment_kidney.Rdata")


####enrichment result####
#load("./GSE109816_kidney/enrichment_kidney.Rdata")

#Reduce(intersect,list(NRP1_go,AXL_go,FURIN_go,CTSL_go))
#Reduce(intersect,list(NRP1_kegg,AXL_kegg,FURIN_kegg,CTSL_kegg))

#save(NRP1_go,AXL_go,FURIN_go,CTSL_go,NRP1_kegg,AXL_kegg,FURIN_kegg,CTSL_kegg,
#     file = "./GSE109816_kidney/enrichment_kidney.Rdata")

load("./GSE109816_heart/enrichment_heart.Rdata")
load("./GSE115469_liver/enrichment_liver.Rdata")
load("./GSE129845_bladder/enrichment_bladder.Rdata")

Reduce(intersect,list(NRP1_heart_go,NRP1_liver_go,NRP1_bladder_go,NRP1_bladder_go))
Reduce(intersect,list(AXL_heart_go,AXL_liver_go,AXL_bladder_go))
Reduce(intersect,list(FURIN_heart_go,FURIN_bladder_go,FURIN_kidney_go))
Reduce(intersect,list(CTSL_heart_go,CTSL_liver_go,CTSL_bladder_go,CTSL_kidney_go))

Reduce(intersect,list(NRP1_heart_kegg,NRP1_liver_kegg,NRP1_bladder_kegg))

Reduce(intersect,list(NRP1_heart_go,NRP1_liver_go,NRP1_bladder_go,NRP1_bladder_go,CTSL_bladder_go,
                      CTSL_heart_go,CTSL_liver_go,CTSL_kidney_go))
Reduce(intersect,list(NRP1_heart_go,AXL_heart_go,CTSL_heart_go,FURIN_heart_go))

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





#load("./GSE131685_kidney/annotation.Rdata")
kid_1$celltype <- Idents(kid_1)

VlnPlot(kid_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"),group.by = 'celltype')
#fig_1_3 <- VlnPlot(kid_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"),group.by = 'celltype')
DotPlot(kid_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"),group.by = 'celltype')
fig_kidney_4 <- DotPlot(kid_1,features = c("ACE2","TMPRSS2","NRP1","FURIN","CTSL"),group.by = 'celltype')
fig_kidney_4

fig_1 <- CombinePlots(plots = list(fig_1_1,fig_1_2,fig_1_3))
ggsave("./GSE131685_kidney/kidney_pdf/kidney_fig.pdf", plot = fig_1, width = 15, height = 18)





####~####
if (FALSE){
  ##tSNE Plot
  kid <-RunTSNE(kid, reduction = "harmony", dims = 1:20)
  TSNEPlot(kid, do.label = T, label = TRUE, do.return = T, pt.size = 1)
  TSNEPlot(kid, do.return = T, pt.size = 1, group.by = "orig.ident", split.by = 'orig.ident')
  TSNEPlot(kid, do.return = T, pt.size = 1, group.by = "Phase")
  
  #Select a subset of PT cells
  PT <- SubsetData(kid, ident.use = c(0,1,2), subset.raw = T)
  saveRDS(PT,file="/staging/yuzhenyuan/PT.rds")
}