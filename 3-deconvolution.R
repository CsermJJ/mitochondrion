###
library(IOBR)
sm<- generateRef_seurat(sce = Mac.ann, celltype = "Celltype",assay_deg = "SCT",slot_deg = "data",assay_out = "SCT", slot_out = "data")
sm=as.data.frame(sm)
library(homologene)
a=rownames(sm)
a=as.data.frame(a)
b=mouse2human(a$a)
colnames(a)="mouseGene"
q=merge(a,b,by="mouseGene")
q=q[!duplicated(q$humanGene),]
q=q[!duplicated(q$mouseGene),]
rownames(q)=q[,1]
w=intersect(rownames(q),rownames(sm))
sm=sm[w,]
q=q[w,]
all=cbind(q,sm)
all=all[!duplicated(all$humanGene),]
rownames(all)=all[,2]
all=all[,-c(1:4)]

svr<- deconvo_tme(eset = eset, reference  = all,  method = "svr", arrays  = T, absolute.mode = FALSE, perm = 1000)

library(gghalves)
iris=read.table("deconvo.txt",sep = "\t",header = T,check.names = F)
ggplot(iris, aes(x =Group, y = Mac_Cd83)) +
  geom_half_violin(
    aes(fill = Group), 
    side = 'r',
    position = position_nudge(x = .25, y = 0), 
    adjust = 2/3) + 
  geom_boxplot(aes(fill = Group),width = 0.1,
               position = position_nudge(x = .25, y = 0)) +
  geom_point(aes(color = Group),
             position = position_jitter(width = 0.15), 
             size = 2) + 
  geom_signif(comparisons = list(c("Non-progressive", "Progressive")), 
              map_signif_level = T, 
              step_increase = 0.05) +
  stat_compare_means(method = "t.test",paired = F,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "")),label = "p.signif",
                     size = 8, 
                     comparisons=list(c("Non-progressive", "Progressive")))+
  scale_fill_manual(values = c("#A0C1D4","#91A3BB")) +scale_colour_manual(values = c("#A0C1D4","#91A3BB"))+
  xlab('') +
  theme_classic() +
  coord_flip()