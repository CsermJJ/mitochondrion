###CIBERSORT
library(tidyr)
library(ggpubr)
library(ggsci)
library(introdataviz)
df=read.table("data.txt",sep = "\t",header = T,check.names = F)
df$Group=factor(df$Group,levels = c("Low PDRM_score","High PDRM_score"))

pdf("Cibersort.pdf",width = 14,height = 5)
ggplot(df,aes(x=gene,y=expression,fill=Group)) +
  geom_boxplot(position=position_dodge(1))+
  theme_classic()+ scale_fill_brewer(palette = "Accent")+
  stat_boxplot(geom = "errorbar",width=0.3,position=position_dodge(1))+
  labs(x="",y="Cibersort score")+rotate_x_text(60)+
  stat_compare_means(aes(group=Group),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")+
  theme(legend.position = "top")
dev.off()  

####
library(limma)
library(estimate)
filterCommonGenes(input.f="symbolGSE24290.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina") 
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
out=rbind(ID=colnames(scores),scores)
write.table(out,file="Estimate scores.txt",sep="\t",quote=F,col.names=F)

data=read.table("clipboard",sep="\t",header=T,check.names=F)
length=length(levels(factor(data$PDRM_group)))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(PDRM_score, StromalScore)) + 
  xlab("PDRM_score")+ylab("Stromal Score")+
  geom_point(aes(colour=PDRM_group))+
  scale_color_manual(values=bioCol[1:length])+ 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =PDRM_score, y =StromalScore))
pdf(file="StromalScore cor.pdf", width=6, height=4.5)
print(p1)
dev.off()
p1=ggplot(data, aes(PDRM_score, ImmuneScore)) + 
  xlab("PDRM_score")+ylab("Immune Score")+
  geom_point(aes(colour=PDRM_group))+
  scale_color_manual(values=bioCol[1:length])+ 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =PDRM_score, y =ImmuneScore))
pdf(file="ImmuneScore cor.pdf", width=6, height=4.5)
print(p1)
dev.off()
p1=ggplot(data, aes(PDRM_score, ESTIMATEScore)) + 
  xlab("PDRM_score")+ylab("ESTIMATE Score")+
  geom_point(aes(colour=PDRM_group))+
  scale_color_manual(values=bioCol[1:length])+ 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =PDRM_score, y =ESTIMATEScore))
pdf(file="ESTIMATEScore cor.pdf", width=6, height=4.5)
print(p1)
dev.off()

####
library(ggplot2)
library(ggsignif)
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
data$PDRM_group=factor(data$PDRM_group,levels =c("Low PDRM_score","High PDRM_score") )
pdf(file = "ImmuCellAI infiltrationscore.pdf",width = 7.5,height = 6)
ggplot(data=data,aes(x=PDRM_group,y=InfiltrationScore,colour = PDRM_group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=PDRM_group,y=InfiltrationScore,colour=PDRM_group,fill=PDRM_group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=PDRM_group,y=InfiltrationScore,colour = PDRM_group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Low PDRM_score","High PDRM_score"), 
                    values =c("#7FC97F","#BEAED4"))+
  scale_color_manual(limits=c("Low PDRM_score","High PDRM_score"), 
                     values=c("#7FC97F","#BEAED4"))+ 
  geom_signif(mapping=aes(x=PDRM_group,y=InfiltrationScore), 
              comparisons = list(c("Low PDRM_score","High PDRM_score")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0),
              size=1, 
              textsize = 7, 
              test = "wilcox.test")+ 
  theme_bw()+
  labs(x="",y="InfiltrationScore")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()


library(ggplot2) 
library(ggsignif) 
library(gghalves)
df <- read.table("clipboard",header = T,sep = "\t")
df$PDRM_group <- factor(df$PDRM_group,levels = c("Low PDRM_score","High PDRM_score"))
ggplot(df,aes(PDRM_group,Score,fill=PDRM_group))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill=PDRM_group),shape=21,size=2.5,width=0.2)+
  geom_signif(comparisons = list(c("Low PDRM_score","High PDRM_score")),
              map_signif_level = T, 
              test = wilcox.test,
              tip_length = c(0,0,0,0,0,0),vjust = 0.5,
              size=1.5,color="black",textsize = 10)+
  theme_bw()+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5,size = 20),
        panel.border = element_rect(size = 1),axis.title = element_text(size = 15),
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black",size = 15),
        legend.position = "none",
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="ImmuCellAI score")+
  scale_fill_manual(values = c("#7FC97F","#BEAED4"))


library(ggalluvial)
sample_clust1=read.table("clipboard",sep = "\t",header = T,check.names = F)
sample_clust1$Group=factor(sample_clust1$Group,levels = c("Low PDRM_score","High PDRM_score"))
p1=ggplot(sample_clust1, aes(Group, y = proportion,
                             fill = Cell,
                             stratum = Cell, alluvium = Cell))+
  geom_stratum(width = 0.7, color='white')+
  geom_alluvium(alpha = 0.4,
                width = 0.7,
                color='white',
                linewidth = 1,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0,0))+
  labs(y="Percentage(%)",x="")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=14, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 12),
        strip.text = element_text(color = "black", size = 12),
        strip.background = element_rect(color = "black", fill="grey90"))+
  theme(legend.position = "right")+
  guides(fill=guide_legend(keywidth = 1.2, keyheight = 1.2)) +
  scale_fill_manual(values = c("#FDC086","#386CB0"))
pdf(file = "ImmuCellAI response proportion.pdf",width = 7,height = 6.5)
p1
dev.off()


###
library(tidyverse)
library(ggpubr)
library(ggsci)
library(introdataviz)
Type=read.table("Estimate scores.txt",sep="\t",check.names=F,row.names=1,header=T)
rt=read.table("ICI.txt",sep="\t",check.names=F,row.names=1,header=T)
Type=Type[colnames(rt),]
rt=t(rt)
data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(expression=(rt[,i]),gene=i,Group=Type[,5]))
}
write.table(data,file="data.txt",sep="\t",quote=F)

df=read.table("data.txt",sep = "\t",header = T,check.names = F)
df$Group=factor(df$Group,levels = c("Low PDRM_score","High PDRM_score"))
pdf(file = "ICI.pdf",width = 15,height = 5)
ggplot(df,aes(x = gene,y = expression,fill = Group)) +
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1.5) +
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1),
        legend.position = 'top') +
  scale_fill_manual(values = c("#7FC97F","#BEAED4"))+
  stat_compare_means(aes(group=Group),method = "wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "")),label = "p.signif",
                     size = 5)
dev.off()