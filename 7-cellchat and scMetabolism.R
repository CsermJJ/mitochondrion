####
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(cowplot)
library(tidyverse)
library(CellChat)
Mac.ann$Timm23group=ifelse(Mac.ann@assays$SCT@data["Timm23",]>0,"Mac_Timm23_pos","Mac_Timm23_neg")
scRNA_Mac_merge=scRNA.ann
Idents(Mac.ann)=Mac.ann$Timm23group
Mac.ann$Timm23group=factor(Mac.ann$Timm23group,levels = c("Mac_Timm23_pos","Mac_Timm23_neg"))
scRNA_Mac_merge$Celltype = as.character(Idents(scRNA_Mac_merge))
Mac.ann$Celltype = as.character(Idents(Mac.ann))
scRNA_Mac_merge$Celltype[match(colnames(Mac.ann),colnames(scRNA_Mac_merge))] =  Mac.ann$Celltype
Idents(scRNA_Mac_merge)=scRNA_Mac_merge$Celltype
table(scRNA_Mac_merge$Celltype)
scRNA_Mac_merge$type <-Idents(scRNA_Mac_merge)
metadata=scRNA_Mac_merge@meta.data
cellchat = createCellChat(object = scRNA_Mac_merge,meta =metadata, group.by = "type")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB = CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use = CellChatDB
cellchat@DB = CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(object=cellchat,raw.use = TRUE,population.size=T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

pdf('Cellchat-Interaction number.pdf', width=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf('Cellchat-Interaction weights.pdf', width=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

LRpro=netVisual_bubble(cellchat)
LRpro=LRpro$data[1:8]
write.table(LRpro,"2.3_all LR_data.txt",sep = "\t",row.names = F,quote = F)

library(CCPlotR)
input=LRpro[,c(1:5)]
colnames(input)=c("source","target","ligand","receptor","score")
pdf("cellchat_num.pdf",width = 8.5,height = 6.5)
cc_heatmap(input)
dev.off()

####
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat)
library(nichenetr)
library(dplyr)
Mac.ann <- Mouse.sc.metabolism(Mac.ann, metabolism.type = 'REACTOME')
df = data.frame(t(Mac.ann@assays[["METABOLISM"]][["score"]]))
names(Mac.ann$orig.ident)
rownames(df) <- gsub(".", "-", rownames(df), fixed = TRUE)
df = df[names(Mac.ann$Timm23group),]
df$Cluster <- Mac.ann$Timm23group
avg_df =aggregate(df[,1:ncol(df)-1],list(df$Cluster),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]

avg_df <- as.data.frame(t(avg_df))
write.table(avg_df,"Rectome metabolism result.txt",sep = "\t",quote = F)
library(pheatmap)
pdf(file = "Rectome metabolism.pdf",width = 15,height = 8)
pheatmap(t(avg_df), show_colnames = T,scale='column', cluster_rows = F,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = T)
dev.off()
