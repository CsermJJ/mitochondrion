#####DPN_DRG_GSE176017
dir='~/DPN_DRG_GSE176017' 
samples=list.files( dir ,pattern = 'gz')
samples 
library(data.table)
ctList = lapply(samples,function(pro){ 
  print(pro)
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  colnames(ct) = paste(gsub('.txt.gz','',pro),
                       colnames(ct) ,sep = '_')
  ct=ct[,-1] 
  return(ct)
})
lapply(ctList, dim)
tmp =table(unlist(lapply(ctList, rownames)))
cg = names(tmp)[tmp==length(samples)]
bigct = do.call(cbind,
                lapply(ctList,function(ct){ 
                  ct = ct[cg,] 
                  return(ct)
                }))
scRNA=CreateSeuratObject(counts =bigct, 
                         min.cells = 3, 
                         min.features = 200,
                         project = "DRG" )
scRNA=SCTransform(scRNA)
scRNA=RunPCA(scRNA, npcs = 50, verbose = FALSE)
scRNA=RunUMAP(scRNA, reduction = 'pca', dims = 1:15)%>% 
  FindNeighbors(., reduction = 'pca', dims = 1:15)%>% 
  FindClusters(., resolution = 0.5)
My_DoubletFinder <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = T,num.cores = 64)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),]$pK))
  DoubletRate = 0.008
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi <- round(DoubletRate*ncol(data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  colnames(data@meta.data)[ncol(data@meta.data)] = "DF.classifications"
  return(data)
}
scRNA <- My_DoubletFinder(scRNA) 
DPN_GSE176017 <- subset(scRNA, subset = (DF.classifications == "Singlet"))

################DPN_DRG_GSE248328
setwd("~/DPN_DRG_GSE248328")
folders=list.files('./')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
for (i in seq_along(sceList)) {
  sceList[[i]][["percent.mt"]] <- PercentageFeatureSet(sceList[[i]], pattern = "^mt-")
  sceList[[i]][["percent.rp"]] <- PercentageFeatureSet(sceList[[i]],pattern = "^Rp[sl]")
  sceList[[i]] <- subset(sceList[[i]], 
                         subset = nCount_RNA > 500 & 
                           nFeature_RNA < 5000 & 
                           percent.mt < 10)
}
for (i in seq_along(sceList)) {
  sceList[[i]] <-  SCTransform(sceList[[i]])
}
for (i in seq_along(sceList)) {
  sceList[[i]] <- RunPCA(sceList[[i]], npcs = 50, verbose = FALSE)
}
for (i in seq_along(sceList)) {
  sceList[[i]] <- RunUMAP(sceList[[i]], reduction = 'pca', dims = 1:15)%>% 
    FindNeighbors(., reduction = 'pca', dims = 1:15)%>% 
    FindClusters(., resolution = 0.5) 
}
My_DoubletFinder <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = T,num.cores = 32)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),]$pK))
  DoubletRate = ncol(data)*8*1e-6
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi <- round(DoubletRate*ncol(data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  colnames(data@meta.data)[ncol(data@meta.data)] = "DF.classifications"
  return(data)
}

library(DoubletFinder)
for (i in seq_along(sceList)) {
  sceList[[i]] <- My_DoubletFinder(sceList[[i]]) 
}

duPlot <- list()
for (i in seq_along(sceList)) {
  p = DimPlot(sceList[[i]], group.by = "DF.classifications",reduction = "umap")+ggtitle(unique(sceList[[i]]$orig.ident))
  duPlot[[i]] <- p
  
}
library(cowplot)
pdf(file = "1.1_DoubletFinder.pdf",width = 10,height = 4)
plot_grid(duPlot[[1]],duPlot[[2]],ncol =2)
dev.off()

for (i in seq_along(sceList)) {
  sceList[[i]] <- subset(sceList[[i]], subset = (DF.classifications == "Singlet"))
}

DPN_GSE248328 <- merge(sceList[[1]], 
                       y = c(sceList[[2]]), 
                       add.cell.ids = folders,
                       project = "DRG")

##########################merge
scRNA <- merge(DPN_GSE176017, 
               y = c(DPN_GSE248328))

scRNA= SCTransform(scRNA)
scRNA=RunPCA(scRNA, features = VariableFeatures(scRNA),npcs = 50)
scRNA <- FindNeighbors(scRNA, dims = 1:15)
scRNA <- FindClusters(scRNA, resolution = 0.5)
scRNA <- RunUMAP(scRNA, dims = 1:15)


scRNA.harm <- RunHarmony(scRNA,c("Dataset"))
scRNA.harm <- FindNeighbors(scRNA.harm, reduction = "harmony",dims = 1:15)
scRNA.harm=FindClusters(object=scRNA.harm,resolution=c(seq(0.2,1,0.1)))
scRNA.harm <- FindClusters(scRNA.harm, resolution = 0.6) 
scRNA.harm<- RunUMAP(scRNA.harm, reduction = "harmony",dims = 1:15)


####
clumarkers <- FindAllMarkers(scRNA.harm, 
                             only.pos = T, 
                             min.pct = 0.2, 
                             logfc.threshold = 0.5) 
logFCfilter=0.5
adjPvalFilter=0.05
sig.clumarkers=clumarkers[(abs(as.numeric(as.vector(clumarkers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(clumarkers$p_val_adj))<adjPvalFilter),]
write.table(sig.clumarkers,file="cluMarkers.txt",sep="\t",row.names=F,quote=F)

###
library(scRNAtoolVis)
library(Seurat)
library(jjPlot)
reduc <-
  data.frame(Seurat::Embeddings(scRNA.harm, reduction = "umap"))
meta <-scRNA.harm@meta.data
pc12 <- cbind(reduc, meta)
label_id <- pc12 %>% group_by(seurat_clusters) %>%
  summarise(x_m = median(UMAP_1),y_m = median(UMAP_2))
label_id$facet <- "umap"
cell_num <- pc12 %>% group_by(seurat_clusters) %>%
  summarise(num = n())
cell_num$y <- 1:nrow(cell_num)
pc12$facet <- "umap"
cell_num$facet <- "numbers"
cell_num$facet <- factor(cell_num$facet,levels = c("umap","numbers"))
pc12$facet <- factor(pc12$facet,levels = c("umap","numbers"))
label_id$facet <- factor(label_id$facet,levels = c("umap","numbers"))
col=c("#54990F","#843C39","#31A354", "#8C6D31", "#C7E9C0", "#E6550D","#1B9E77","#66A61E", "#3182BD",  "#BD9E39", "#E7BA52",  "#E41A1C", 
      "#6BAED6", "#9ECAE1","#74C476", "#E7CB94", "#FDAE6B",  "#A1D99B", "#AD494A", "#99600F", 
      "#E7298A", "#C3BC3F", "#D6616B", "#FF7F00",  "#B3823E",  
      "#F1788D", "#C6DBEF", "#E6550D", "#E7969C")
pdf("1.5_umap_cluster.pdf",width = 6.5,height = 5)
ggplot() +
  scale_y_continuous(position = "right") +
  geom_point(data = pc12,
             aes(color = seurat_clusters,x = UMAP_1,y = UMAP_2),size=0.3,
             key_glyph = draw_number_circle) +
  geom_label(data = label_id,
             aes(x = x_m,y = y_m,label = seurat_clusters),
             fill = "grey90",
             label.size = NA,label.r = unit(0.25,"cm"),size = 3) +
  geom_markArrow(data = pc12,rel.pos = 0.06,corner.pos = "left_b",label.shift = c(0.025,-0.025)) +
  guides(color = guide_legend(override.aes = list(label = 0:16))) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_bw(base_size = 12) +
  theme(panel.spacing.y = unit(0,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank())
dev.off()


######################
new.cluster.ids <- c("Satellite glial cells","Satellite glial cells","Satellite glial cells","Satellite glial cells",
                     "Satellite glial cells","Neurons","Neurons", "Satellite glial cells","Endothelial cells","Endothelial cells",
                     "Schwann cells","Neurons","Satellite glial cells","Mural cells","Mural cells","Macrophages","Fibroblasts") 
scRNA.ann=scRNA.harm
names(new.cluster.ids) <- levels(scRNA.ann)
scRNA.ann <- RenameIdents(scRNA.ann, new.cluster.ids)
scRNA.ann$Celltype=scRNA.ann@active.ident

reduc <-
  data.frame(Seurat::Embeddings(scRNA.ann, reduction = "umap"))
meta <-scRNA.ann@meta.data
pc12 <- cbind(reduc, meta)
label_id <- pc12 %>% group_by(Celltype) %>%
  summarise(x_m = median(UMAP_1),y_m = median(UMAP_2))
label_id$facet <- "umap"
cell_num <- pc12 %>% group_by(Celltype) %>%
  summarise(num = n())
cell_num$y <- 1:nrow(cell_num)
pc12$facet <- "umap"
cell_num$facet <- "numbers"
cell_num$facet <- factor(cell_num$facet,levels = c("umap","numbers"))
pc12$facet <- factor(pc12$facet,levels = c("umap","numbers"))
label_id$facet <- factor(label_id$facet,levels = c("umap","numbers"))
col1=c("#5fa664","#fbbab6","#e1c548","#f0e2a3","#abd0a7","#f9766e","#ca6a6b","#e5b5b5","#4e79a6","#bac4d0","#45337f","#a199be","#aedd2f","#d7ee96")
pdf("2.1_umap_cluster.pdf",width = 8.5,height = 5)
ggplot() +
  scale_y_continuous(position = "right") +
  geom_rect(data = cell_num,aes(xmin = 0,xmax = num,
                                ymin = rev(y - 0.25),
                                ymax = rev(y + 0.25),
                                fill = Celltype),
            show.legend = F) +
  geom_point(data = pc12,
             aes(color = Celltype,x = UMAP_1,y = UMAP_2),size=0.3,
             key_glyph = draw_number_circle) +
  geom_label(data = label_id,
             aes(x = x_m,y = y_m,label = Celltype),
             fill = "grey90",
             label.size = NA,label.r = unit(0.25,"cm"),size = 3) +
  geom_markArrow(data = pc12,rel.pos = 0.06,corner.pos = "left_b",label.shift = c(0.025,-0.025)) +
  guides(color = guide_legend(override.aes = list(label = 0:6))) +
  scale_fill_manual(values = col1) +
  scale_color_manual(values = col1) +
  theme_bw(base_size = 12) +
  theme(panel.spacing.y = unit(0,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank())+
  facet_wrap(~facet, scales = "free",ncol = 2,
             strip.position = "top") +
  ggh4x::force_panelsizes(cols = c(5,2))
dev.off()

###
library(ComplexHeatmap)
library(circlize)
scRNA.markers <- FindAllMarkers(scRNA.ann, 
                                only.pos = T, 
                                min.pct = 0.2, 
                                logfc.threshold = 0.5) 
logFCfilter=0.5
adjPvalFilter=0.05
sig.markers=scRNA.markers[(abs(as.numeric(as.vector(scRNA.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(scRNA.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="2.2_cellMarkers.txt",sep="\t",row.names=F,quote=F)
cm_top20 <- sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
cm_top20_table=unstack(cm_top20, gene ~ cluster)
names(cm_top20_table)=gsub("X","cluster",names(cm_top20_table))
write.csv(file="2.2_top20_cellMarkers.csv",cm_top20_table,row.names=F)

cal2_cols=c("#54990F","#843C39","#31A354", "#8C6D31", "#C7E9C0", "#E6550D","#1B9E77","#66A61E", "#3182BD",  "#BD9E39", "#E7BA52",  "#E41A1C", 
            "#6BAED6", "#9ECAE1","#74C476", "#E7CB94", "#FDAE6B")
avg <- AverageExpression(object =scRNA.ann, group.by = 'seurat_clusters',assays = "SCT", slot = 'data',features = cm_top20$gene) # Return average expression values across cells in each cluster for the selected genes from top3 object
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")) 
cell_order <- c("0","1","2","3","4","7","12","5","6","11","8","9","10","13","14","15","16")
cell_lineages <- c("Satellite glial cells", "Satellite glial cells", "Satellite glial cells","Satellite glial cells",
                   "Satellite glial cells","Satellite glial cells","Satellite glial cells","Neurons","Neurons","Neurons",
                   "Endothelial cells","Endothelial cells","Schwann cells","Mural cells","Mural cells","Macrophages",
                   "Fibroblasts")
cell_lineage_colors <- c("Satellite glial cells" = "#5fa664","Neurons" ="#fbbab6", "Endothelial cells" = "#e1c548","Schwann cells" = "#f0e2a3",
                         "Mural cells" ="#abd0a7","Macrophages" = "#f9766e", "Fibroblasts" ="#ca6a6b")
cell_counts <- as.numeric(table(factor(scRNA.ann$seurat_clusters, levels = cell_order)))
names(cell_order) <- cal2_cols
anno_df <- data.frame(cell_order, cell_lineages)
names(cal2_cols) <- levels(droplevels(as.factor(scRNA.ann$seurat_clusters)))
ha = HeatmapAnnotation(annotation_label = c("Seurat_clusters", "Celltype"), df = anno_df,
                       border = TRUE, which = 'col', col = list(cell_order = cal2_cols,cell_lineages = c(cell_lineage_colors)))
avg <- as.data.frame(avg$SCT)
genes_show =c('Tyrp1','Fabp7','Timp1','Gng3','Nefm',"Cldn5",'Flt1','Eng',
              'Mpz','Mbp','Cldn19','Acta2','Des','Rgs5','C1qc','C1qa',
              'Ccl3','Dcn','Lum','Col1a1') 
pdf(file = "cellmarker heatmap.pdf",width = 15,height =6)
Heatmap(t(scale(t(avg[,cell_order]))), name = "mat", rect_gp = gpar(col = "black", lwd = 0), col = col_fun, column_order = cell_order,cluster_rows =F,column_names_side = c("top"),
        clustering_method_rows = "single", show_row_names = FALSE, row_names_gp = grid::gpar(fontsize = 25),column_names_gp = grid::gpar(fontsize = 15), top_annotation =  ha) +
  rowAnnotation(link=anno_mark(at=which(rownames(avg) %in% genes_show), labels = rownames(avg)[rownames(avg) %in% genes_show], which = "rows"))
dev.off()

#####
meta=scRNA.ann@meta.data
meta$Group=factor(meta$Group,levels = c("Control","DPN"))
summary <- table(meta[,c('Celltype','Group')])
roe <- as.data.frame(ROIE(summary))
roe$Celltype=rownames(roe)
roe1=melt(roe,id.variables="Celltype",variable.name=c("Group"),value.name = c("Ro/e"))
roe1$pstar=ifelse(roe1$`Ro/e`<1,
                  ifelse(roe1$`Ro/e`<0.8,
                         ifelse(roe1$`Ro/e`<0.2,ifelse(roe1$`Ro/e`<=0,"???","+/???"),"+"),"++"),"+++")
library(showtext)
showtext_auto()
pdf(file="Ro-e.pdf",width=6,height=10) 
ggplot(roe1,aes(Group,Celltype))+geom_tile(aes(fill=`Ro/e`),colour="white",size=1)+
  scale_fill_distiller(palette = "Reds",direction =1)+
  geom_text(aes(label=pstar),col="black",size=5)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size=11),
        axis.text.x = element_text(angle=45,hjust = 1,size = 11),
        axis.text.y = element_text(size=11))+
  labs(fill=paste0("+++       Ro/e> 1","\n\n",
                   "++  0.8 <Ro/e¡Ü 1","\n\n",
                   "+    0.2 ¡ÜRo/e¡Ü 0.8","\n\n",
                   "+/???    0 <Ro/e< 0.2","\n\n",
                   "???              Ro/e= 0","\n\n",
                   "Ro/e"))
dev.off()
showtext_auto(FALSE)
