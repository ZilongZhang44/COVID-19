library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
load("nCoV.integrated2.Rdata")

h <- nCoV.integrated2[,nCoV.integrated2@meta.data$group=='HC']
table(h$celltype_new)

m <- nCoV.integrated2[,nCoV.integrated2@meta.data$group=='M']
table(m$celltype_new)
s <- nCoV.integrated2[,nCoV.integrated2@meta.data$group=='S']
table(s$celltype_new)

table(h$celltype_new)

####all_celltype####
healthy <- data.frame(cell_type = h$celltype_new)
medium <- data.frame(cell_type = m$celltype_new)
severe <- data.frame(cell_type = s$celltype_new)

healthy <- healthy %>% mutate(condition="Healthy")
medium <- medium %>% mutate(condition="Medium")
severe <- severe %>% mutate(condition="Severe")

DF <- bind_rows(healthy,medium,severe)

DF <- DF %>% group_by(condition, cell_type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = round(Nb/C*100,1))

p <- ggplot(DF, aes(x = condition, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = paste(percent,"%")), position = position_stack(vjust = 0.5))

#ggsave("celltype_all.pdf", p, width = 12, height = 15)

####macrophage####

h_m <- subset(h,celltype_new=="M1" | celltype_new=="M2")
m_m <- subset(m,celltype_new=="M1" | celltype_new=="M2") 
s_m <- subset(s,celltype_new=="M1" | celltype_new=="M2")

healthy <- data.frame(cell_type = h_m$celltype_new)
medium <- data.frame(cell_type = m_m$celltype_new)
severe <- data.frame(cell_type = s_m$celltype_new)

healthy <- healthy %>% mutate(condition="Healthy")
medium <- medium %>% mutate(condition="Medium")
severe <- severe %>% mutate(condition="Severe")

DF <- bind_rows(healthy,medium,severe)

DF <- DF %>% group_by(condition, cell_type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = round(Nb/C*100,1))

p <- ggplot(DF, aes(x = condition, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = paste(percent,"%")), position = position_stack(vjust = 0.5))

ggsave("celltype_macrophage-1.pdf", p, width = 8, height = 6)  ###C

####example####
set.seed(123)
Untreated <- data.frame(Cell_Type = sample(LETTERS[1:4],10, replace = TRUE))
Treated <- data.frame(Cell_Type =sample(LETTERS[1:4],25, replace = TRUE))

library(dplyr)
Untreated <- Untreated %>% mutate(Condition = "Untreated")
Treated <- Treated %>% mutate(Condition = "Treated")
DF <- bind_rows(Untreated, Treated)

DF <- DF %>% group_by(Condition, Cell_Type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

library(ggplot2)
ggplot(DF, aes(x = Condition, y = percent, fill = Cell_Type))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = paste(percent,"%")), position = position_stack(vjust = 0.5))


#####################################  2C  ###############################

markers <- c('FCN1','CCL2','CCL3','TREM2')

pp_temp = FeaturePlot(object =mac1, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 10),axis.text = element_text(size = 10),plot.title = element_text(family = 'sans',face='italic',size=12),legend.text = element_text(size = 10),legend.key.height = unit(1.2,"line"),legend.key.width = unit(0.8,"line") ))
pp = CombinePlots(plots = plots,ncol = 2,legend = 'right')

pp

ggsave("macro_markers.pdf",
       plot = pp,width = 8, height = 12)


