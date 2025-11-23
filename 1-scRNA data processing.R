############################
seu_data_count <- GetAssayData(scRNA.harm, slot = 'counts')
seumeta=scRNA.harm@meta.data
seumeta=seumeta[,-c(2:9,12:23)]
scRNA<- CreateSeuratObject(counts = seu_data_count,  
                           min.cells = 3, 
                           min.features = 200)
q=intersect(rownames(seumeta),colnames(scRNA))
seumeta=seumeta[q,]
scRNA=AddMetaData(scRNA,seumeta)

scRNA= SCTransform(scRNA)
scRNA=RunPCA(scRNA, features = VariableFeatures(scRNA),npcs = 50)
pdf('1.1_ElbowPlot.pdf', width = 9)
ElbowPlot(scRNA,ndims = 30)
dev.off()
scRNA.harm <- RunHarmony(scRNA,"Dataset")
scRNA.harm <- FindNeighbors(scRNA.harm, reduction = "harmony",dims = 1:10)
scRNA.harm=FindClusters(object=scRNA.harm,resolution=c(seq(0.2,1,0.1)))
pdf(file = "1.2_clustree_harm.pdf",width = 12,height = 10)
clustree(scRNA.harm@meta.data,prefix="SCT_snn_res.")
dev.off()

scRNA.harm <- FindClusters(scRNA.harm, resolution = 0.5) 
scRNA.harm<- RunUMAP(scRNA.harm, reduction = "harmony",dims = 1:10)

markergene1<- c('Avil','Gap43','Nefm','Nefl',        #Neurons
                'Fabp7','Tyrp1','Timp3','Adam17',    #Satellite glial cells
                "Cldn5",'Flt1','Emcn','Prom1',       #Endothelial cells
                'Tagln','Acta2','Rgs5','Des',        #Mural cells  
                'Pdgfra','Lum','Col1a1','Dcn',       #Fibroblasts
                'Mpz','Mag','Mbp','Cntf',            #Schwann cells
                'Lyz2','Cd68','Mrc1','Csf1r',        #Macs
                'Top2a','Mki67','Cdkn3','Cdk1')      #Proliferating cells  
pdf('1.3_dotplot_marker_cluster.pdf', width = 12, height = 8)
DotPlot(object = scRNA.harm, features =  markergene1,scale = T,group.by = "seurat_clusters") + 
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) 
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



