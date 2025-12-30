####
set.seed(123123)
Mac=subset(scRNA.ann,ident="Macrophages")
mito_genes=rownames(Mac)[grep("^mt-", rownames(Mac),ignore.case = T)]
ribo_genes=rownames(Mac)[grep("^Rp[sl]", rownames(Mac),ignore.case = T)]
Hb_genes=rownames(Mac)[grep("^Hb[^(p)]", rownames(Mac),ignore.case = T)]
a=rownames(Mac)
a=as.data.frame(a)
a=a$a
q=unique(c(a[!a %in% Hb_genes], a[duplicated(a)]))
q=unique(c(q[!q %in% mito_genes], q[duplicated(q)]))
q=unique(c(q[!q %in% ribo_genes], q[duplicated(q)]))
q=as.data.frame(q)
Mac <- subset(Mac, features = q$q)

Mac=SCTransform(Mac)
Mac=RunPCA(Mac, features = VariableFeatures(Mac),npcs = 50)
Mac <- RunHarmony(Mac,"Dataset")
Mac <- FindNeighbors(Mac, reduction = "harmony",dims = 1:15)
Mac=FindClusters(object=Mac,resolution=c(seq(0.2,1,0.1)))
Mac <- FindClusters(Mac, resolution = 0.5) 
Mac<- RunUMAP(Mac, reduction = "harmony",dims = 1:15)

clumarkers <- FindAllMarkers(Mac, 
                             only.pos = T, 
                             min.pct = 0.2, 
                             logfc.threshold = 0.25) 
logFCfilter=0.25
adjPvalFilter=0.05
sig.clumarkers=clumarkers[(abs(as.numeric(as.vector(clumarkers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(clumarkers$p_val_adj))<adjPvalFilter),]
write.table(sig.clumarkers,file="cluMarkers.txt",sep="\t",row.names=F,quote=F)

###
markergene1<- c('Cd83','Ccl7','Cd74','Cxcl1',
                'S100a8','S100a9','Ccl5','Il2ra')   
pdf('1.4_dotplot_marker_cluster.pdf', width = 12, height = 8)
DotPlot(object = Mac, features =  markergene1,scale = T,group.by = "seurat_clusters") + 
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) 
dev.off()

###
new.cluster.ids <- c("Mac_Cd83+","Mac_Cd74+", "Mac_S100a8+", "Mac_Ccl5+", "Mac_Cd74+")
Mac.ann=Mac
names(new.cluster.ids) <- levels(Mac.ann)
Mac.ann <- RenameIdents(Mac.ann, new.cluster.ids)
Mac.ann$Celltype=Mac.ann@active.ident

###
library(scatterpie)
library(tidydr)
df <- Mac.ann@reductions$umap@cell.embeddings%>%   as.data.frame() %>%  cbind(Celltype = Mac.ann@meta.data$Celltype)
Cellratio<-prop.table(table(Mac.ann$Group,Idents(Mac.ann)),margin=2)
Cellratio <- as.data.frame(Cellratio)
library(tidyr)
freq <-spread(Cellratio, Var1, Freq)
colnames(freq)[1] <- 'Celltype'
freq <- freq[sort(freq$Celltype),]
label <- df %>%group_by(Celltype) %>%  summarise(UMAP_1 = median(UMAP_1),            
                                                 UMAP_2 = median(UMAP_2))%>%  as.data.frame()
rownames(label) <- label$Celltype
cell_number <- as.data.frame(table(Mac.ann$Celltype))
colnames(cell_number)[2]<-'cellnumber'
cell_number$cellnumber<-log2(cell_number$cellnumber)/20
data = cbind(freq,label[,c(2:3)], cell_number[,c(2)])
colnames(data)[6]<- 'cellnumber'
col2=c("#57B1E7","#009C73","#E39D01","#D25E00")
p1=ggplot()+  
  geom_point(data=df, aes(x= UMAP_1 , y = UMAP_2 ,color =Celltype),size = 2,shape=16) +  
  scale_color_manual(values = col2)+   
  theme_dr()+  
  theme(panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),
        legend.position = "none")+  
  geom_scatterpie(data=data,                  
                  aes(x=UMAP_1,y=UMAP_2,                      
                      group=Celltype,                      
                      r=cellnumber),                  
                  cols=names(freq)[2:3])+  
  scale_fill_manual(values = c("#87B5B2","#CC88B0"),name='Group')

col2=c("#57B1E7","#009C73","#E39D01","#D25E00")
clust_freq<-as.data.frame(table(Mac.ann$Celltype))
colnames(clust_freq)=c('Celltype','cell_num')
clust_freq=clust_freq[order(clust_freq$cell_num,decreasing = T),]
clust_freq$Celltype=factor(clust_freq$Celltype,levels =c("Mac_Cd83+","Mac_Cd74+", "Mac_S100a8+", "Mac_Ccl5+"))
cell_level = c("Mac_Cd83+","Mac_Cd74+", "Mac_S100a8+", "Mac_Ccl5+")
clust_freq=clust_freq %>% 
  mutate(across(Celltype, ~ factor(.x)))
p2<-ggplot(clust_freq,aes(x = order(cell_num,Celltype),y = cell_num,fill=Celltype))+
  geom_bar(stat="identity",width = 0.5)+ggtitle("") +
  theme_classic2() + scale_fill_manual(values = col2)+
  theme(axis.ticks.length = unit(0, 'cm'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +coord_flip() +guides(fill = guide_legend(reverse = F)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0)))+ylim(0,max(clust_freq$cell_num)+10)
pdf(file = "umap_merge.pdf",width = 8,height = 5)
ggarrange(p1,p2,nrow = 1,ncol = 2,widths = c(1.4,1),vjust =0,hjust = 0,align = c("h"))
dev.off()

####
library(Nebulosa)
pdf(file = "density_Cd83.pdf",width = 5.5,height = 5)
plot_density(Mac.ann, "Cd83",size = 2)
dev.off()
pdf(file = "density_Cd74.pdf",width = 5.5,height = 5)
plot_density(Mac.ann, "Cd74",size = 2)
dev.off()
pdf(file = "density_S100a8.pdf",width = 5.5,height = 5)
plot_density(Mac.ann, "S100a8",size = 2)
dev.off()
pdf(file = "density_Ccl5.pdf",width = 5.5,height = 5)
plot_density(Mac.ann, "Ccl5",size = 2)
dev.off()

####
library(ggalluvial)
sample_clust1<-as.matrix(table(scRNA.ann@active.ident,scRNA.ann$orig.ident))
sample_clust1=apply(sample_clust1,1,function(x){return(x/sum(x))})
sample_clust1=reshape2::melt(sample_clust1)
colnames(sample_clust1)<-c("Group","Cell","proportion")
write.table(sample_clust1,"proportion.txt",sep="\t",quote=F)
p1=ggplot(sample_clust1, aes(Cell, y = proportion,
                             fill = Group,
                             stratum = Group, alluvium = Group))+
  geom_stratum(width = 0.7, color='white')+
  geom_alluvium(alpha = 0.4,
                width = 0.7,
                color='white',
                linewidth = 1,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0,0))+
  labs(y="Percentage(%)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=17),
        axis.text.x = element_text(size=17, angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size=17),
        strip.text = element_text(color = "black", size = 17),
        strip.background = element_rect(color = "black", fill="grey90"))+
  theme(legend.position = "top")+
  guides(fill=guide_legend(keywidth = 1.2, keyheight = 1.2)) +
  scale_fill_manual(values = zxcolor1)+coord_flip()
pdf(file = "proportion.pdf",width = 9,height = 5)
p1
dev.off()

####miloR
library(miloR)
miloMac <- as.SingleCellExperiment(scRNA.ann)
meta=scRNA.ann@meta.data
milo.obj <- buildGraph(miloMac, k=20, d=30)
milo.obj <- makeNhoods(milo.obj, k=20, d=30, refined=TRUE, prop=0.2)
milo.obj <- calcNhoodDistance(milo.obj, d=30)
milo.obj <- countCells(milo.obj, samples="orig.ident", meta.data=meta)
milo.design <- as.data.frame(xtabs(~ Group + orig.ident, data=meta))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$orig.ident
milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]
milo.res <- testNhoods(milo.obj, design=~Group, design.df=milo.design)
head(milo.res)
milo.res <- annotateNhoods(milo.obj, milo.res, coldata_col = "Celltype")
milo.res$Celltype=factor(milo.res$Celltype,levels = c("Mac_Cd83+", "Mac_Cd74+", "Mac_S100a8+", "Mac_Ccl5+"))
pdf("miloR.pdf",height = 11,width = 8)
milo.res %>% mutate(is_signif = ifelse(SpatialFDR < 1,1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 1,logFC, NA)) %>% arrange(Celltype) %>% mutate(Nhood = factor(Nhood,levels = unique(Nhood))) %>% ggplot(aes(Celltype, logFC, color = logFC_color,size=logFC)) + scale_y_continuous(limits = c(-6, 6), breaks = c(-4, -2, 0, 2,4))+ 
  scale_color_gradient2(low = muted("blue"), high = muted("red")) + guides(color = "none")+xlab("Mac_subtype") + ylab("Log Fold Change") + geom_quasirandom(alpha = 1) + coord_flip() + theme_bw(base_size = 18) + theme(strip.text.y = element_text(angle = 0))

dev.off()

###
library(AUCell)
library(ggplot2)
library(Seurat)
library(clusterProfiler)
cells_rankings <- AUCell_buildRankings(Mac.ann@assays$SCT@data,  nCores=32, plotStats=TRUE) 
c2 <- read.gmt("Mito_metabolism.gmt")
c2=c2[!c2$gene=="",]
geneSets <- lapply(unique(c2$term), function(x){print(x);c2$gene[c2$term == x]})
names(geneSets) <- unique(c2$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores =32, aucMaxRank=nrow(cells_rankings)*0.1)
geneSet <- "Mito_metabolism"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
Mac.ann$AUC <- aucs
df<- data.frame(Mac.ann@meta.data, Mac.ann@reductions$umap@cell.embeddings)
df$Celltype=as.factor(df$Celltype)
class_avg <- df %>%
  group_by(Celltype) %>%
  summarise(
    umap_1 = median(UMAP_1),
    umap_2 = median(UMAP_2)
  )
df <- arrange(df, AUC)
pdf(file = "Mitochondrial metabolism_celltype.pdf",width = 8,height = 4.5)
ggviolin(df,"Celltype","AUC",
         color = 'Celltype',add = 'mean_sd',fill = 'Celltype',
         add.params = list(color = 'black'))+
  scale_color_manual(values = col2)+
  scale_fill_manual(values = col2)+
  theme(axis.text.x.bottom = element_text(angle = 30,vjust = 0.5,hjust = 0.5))+
  NoLegend()+labs(x = '',y='AUCell score',title = "Mitochondrial metabolism")
dev.off()

pdf(file = "Mitochondrial metabolism_Group.pdf",width = 5,height = 5)
ggplot(df,aes(Group,AUC,fill=Group))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.25,size=0.8,outlier.color =NA)+
  geom_signif(comparisons = list(c("Control","DPN")),
              map_signif_level = T,
              test = t.test,
              tip_length = c(0,0,0,0,0,0),vjust = 0.6,
              size=1.5,color="black",textsize = 8)+
  theme_bw()+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5,size = 20),
        panel.border = element_rect(size = 1),axis.title = element_text(size = 15),
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black",size = 15),
        legend.position = "none",
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="AUC score")+
  scale_fill_manual(values = c("#87B5B2","#CC88B0"))
dev.off()




