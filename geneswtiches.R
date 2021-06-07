library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(future)
'''
list.of.packages <- c("SingleCellExperiment", "Biobase", "fastglm", "ggplot2", "monocle",
                      "plyr", "RColorBrewer", "ggrepel", "ggridges", "gridExtra", "devtools",
                      "mixtools")

## for package "fastglm", "ggplot2", "plyr", "RColorBrewer", "ggrepel", "ggridges", "gridExtra", "mixtools"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## for package "SingleCellExperiment", "Biobase"
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)


devtools::install_github("SGDDNB/GeneSwitches")
'''
library(GeneSwitches)

setwd("/home/rstudio/covid/monocle/geneswitches/")

load(file='../../macro_p.Rdata')
load("../monocle_p_4.Rdata")


plot_monocle_State(mycds_p_4)


logexpdata <- as.matrix(macro_p@assays$RNA@data)
sce_p1 <- convert_monocle2(monocle2_obj = mycds_p_4, 
                           states = c(4,3), expdata = logexpdata)
sce_p2 <- convert_monocle2(monocle2_obj = mycds_p_4, 
                           states = c(4,2,1), expdata = logexpdata)
####path 1####
sce_p1 <- binarize_exp(sce_p1, ncores = 1)

#save(sce_p1,file = "sce_p1.Rdata")

###########
#load(file = "sce_p1.Rdata")
sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE)


## filter top 15 best fitting switching genes among all the genes
sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, topnum = 15)
## filter top 15 best fitting switching genes among surface proteins and TFs only
sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, topnum = 20,
                                genelists = gs_genelists, genetype = c("Surface proteins", "TFs"))
## combine switching genes and remove duplicated genes from sg_allgenes
sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])

plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3)
p1_1 <- plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3)
ggsave("p1_1.pdf", plot = p1_1)

plot_gene_exp(sce_p1, gene = "FCN1", reduction = "monocleRD", downsample = F)


## filter genes for pathway analysis using r^2 cutoff 0.1
sg_pw <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.1)
## apply hypergeometric test and determine the switching time
switch_pw <- find_switch_pathway(rowData(sce_p1), sig_FDR = 0.05,
                                 pathways = msigdb_h_c2_c5, sg_pw)
## remove redundant pathways
switch_pw_reduce <- reduce_pathways(switch_pw, pathways = msigdb_h_c2_c5, 
                                    redundant_pw_rate = 0.8)

plot_pathway_density(switch_pw_reduce[1:10,], sg_pw, orderbytime = TRUE)
p1_2 <- plot_pathway_density(switch_pw_reduce[1:10,], sg_pw, orderbytime = TRUE)
ggsave("p1_2.pdf", plot = p1_2)

sg_vis <- filter_switchgenes(sce_p1, topnum = 50, pathway_name = c("GO_DEFENSE_RESPONSE"))
plot_timeline_ggplot(sg_vis, timedata=sce_p1$Pseudotime, txtsize=3)


######path_2####
sce_p2 <- binarize_exp(sce_p2,ncores = 4)
#save(sce_p2,file = "sce_p2.Rdata")
#load(file = "sce_p2.Rdata")
sce_p2 <- find_switch_logistic_fastglm(sce_p2, downsample = TRUE, show_warning = FALSE)


## filter top 15 best fitting switching genes among all the genes
sg_allgenes <- filter_switchgenes(sce_p2, allgenes = TRUE, topnum = 15)
## filter top 15 best fitting switching genes among surface proteins and TFs only
sg_gtypes <- filter_switchgenes(sce_p2, allgenes = FALSE, topnum = 20,
                                genelists = gs_genelists, genetype = c("Surface proteins", "TFs"))
## combine switching genes and remove duplicated genes from sg_allgenes
sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])

plot_timeline_ggplot(sg_vis, timedata = sce_p2$Pseudotime, txtsize = 3)
p2_1 <- plot_timeline_ggplot(sg_vis, timedata = sce_p2$Pseudotime, txtsize = 3)
p2_1
ggsave("p2_1.pdf", plot = p2_1)

plot_gene_exp(sce_p2, gene = "FCN1", reduction = "monocleRD", downsample = F)


## filter genes for pathway analysis using r^2 cutoff 0.1
sg_pw <- filter_switchgenes(sce_p2, allgenes = TRUE, r2cutoff = 0.1)
## apply hypergeometric test and determine the switching time
switch_pw <- find_switch_pathway(rowData(sce_p2), sig_FDR = 0.05,
                                 pathways = msigdb_h_c2_c5, sg_pw)
## remove redundant pathways
switch_pw_reduce <- reduce_pathways(switch_pw, pathways = msigdb_h_c2_c5, 
                                    redundant_pw_rate = 0.8)

plot_pathway_density(switch_pw_reduce[1:10,], sg_pw, orderbytime = TRUE)
p2_2 <- plot_pathway_density(switch_pw_reduce[1:10,], sg_pw, orderbytime = TRUE)
ggsave("p2_2.pdf", plot = p2_2)

sg_vis <- filter_switchgenes(sce_p2, topnum = 50, pathway_name = c("HALLMARK_MYC_TARGETS_V1"))
plot_timeline_ggplot(sg_vis, timedata=sce_p1$Pseudotime, txtsize=3)






######path_commom
sce_p2 <- find_switch_logistic_fastglm(sce_p2, downsample = TRUE, show_warnings = FALSE)

sg_p1 <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.05)
sg_p2 <- filter_switchgenes(sce_p2, allgenes = TRUE, r2cutoff = 0.05)

sg_com <- common_genes(sg_p1, sg_p2, r2cutoff = 0.2,
                       path1name = "Definitive CM", path2name = "non-contractile")
common_genes_plot(sg_com, timedata = sce_p1$Pseudotime)
p_c_1 <- common_genes_plot(sg_com, timedata = sce_p1$Pseudotime)
ggsave("p_c_1.pdf", plot = p_c_1)

sg_disgs <- distinct_genes(sg_p1, sg_p2, r2cutoff = 0.4,
                           path1name = "Definitive CM", path2name = "non-contractile",
                           path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime)
plot_timeline_ggplot(sg_disgs, timedata = sce_p1$Pseudotime, color_by = "Paths", 
                     iffulltml = FALSE, txtsize = 3)

sg_disgs_scale <- distinct_genes(sg_p1, sg_p2, r2cutoff = 0.4, 
                                 path1name = "Definitive CM", path2name = "non-contractile",
                                 path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime, 
                                 scale_timeline = T, bin = 100)
# timedata need to be 1 to (number of bins)
plot_timeline_ggplot(sg_disgs_scale, timedata = 1:100, color_by = "Paths", 
                     iffulltml = FALSE, txtsize = 3)
p_c_2 <- plot_timeline_ggplot(sg_disgs_scale, timedata = 1:100, color_by = "Paths", 
                              iffulltml = FALSE, txtsize = 3)
ggsave("p_c_2.pdf", plot = p_c_2)

gn <- "PARAL1"
p <- plot_gene_exp(sce_p1, gene = gn, reduction = "monocleRD", 
                   downsample = FALSE, fitting = TRUE)
#
p <- plot_gene_exp(sce_p2, gene = gn, reduction = "monocleRD", 
                   downsample = FALSE, fitting = TRUE)

p


