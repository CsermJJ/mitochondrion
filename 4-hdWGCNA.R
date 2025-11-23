###
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())
set.seed(123123)
enableWGCNAThreads(nThreads = 32)
seurat_obj <- SetupForWGCNA(
  Mac.ann,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "DPN"
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Celltype"),  
  reduction = "harmony", 
  k = 25, assay = "SCT",
  max_shared = 10, 
  ident.group = "Celltype"  
)
seurat_obj <- NormalizeMetacells(seurat_obj)
metacell_obj=GetMetacellObject(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Mac_Cd83+"), 
  group.by='Celltype', 
  assay = 'SCT', 
  slot = 'data' 
)
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed')
plot_list <- PlotSoftPowers(seurat_obj)
pdf("1.1_SoftPowers.pdf",width = 9.5,height = 8)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(seurat_obj)
write.table(power_table,"1.1_power table.txt",sep = "\t",quote = F)

seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=14,minModuleSize=100,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'Mac_Cd83+' 
)
pdf("1.2_hdWGCNA Dendrogram.pdf",width = 8,height = 4)
PlotDendrogram(seurat_obj, main='Mac_Cd83+ hdWGCNA Dendrogram')
dev.off()

TOM <- GetTOM(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj,group.by.vars = "orig.ident",assay = "SCT")

hMEs <- GetMEs(seurat_obj)
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Celltype', group_name = 'Mac_Cd83+'
)

p <- PlotKMEs(seurat_obj, ncol=5)
pdf("1.3_kME.pdf",width = 10,height = 5)
p
dev.off()

modules <- GetModules(seurat_obj)
write.table(modules,"1.4_module genes.txt",sep = "\t",quote = F)
hub_df <- GetHubGenes(seurat_obj = seurat_obj)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25, # topN hub genes
  method = "UCell")
plot_list <- ModuleFeaturePlot(
  seurat_obj,point_size = 0.1,
  features = 'hMEs', 
  order = T # order so the points with highest hMEs are on top
)
pdf("1.5_module featureplot.pdf",width = 15,height = 7.5)
wrap_plots(plot_list, ncol = 5)
dev.off()

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p <- DotPlot(seurat_obj, features = mods, group.by = "Celltype")
p <- p + 
  #     coord_flip() + 
  RotatedAxis() + 
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")  
pdf("1.6_module dotplot.pdf",width = 15,height = 7.5)
p
dev.off()

pdf("1.7_module red hubgene.pdf",width = 6,height = 6)
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "red")
dev.off()