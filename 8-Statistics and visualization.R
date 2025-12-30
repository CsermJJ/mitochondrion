##################Figure1
######
library(ggplot2)
library(ggsignif)
cols=c("#FEBF6E","#FF7F00","#FB9B98","#E31A20","#C7E2E4","#8FC2C7",
       "#A6CEE4","#1F78B4","#B2E08A","#30A12C",
       "#FEBF6E","#FF7F00","#FB9B98","#E31A20","#C7E2E4","#8FC2C7")
df <- read.table("clipboard",header = T,sep = "\t")
df$Group <- factor(df$Group,levels = c("Non progressive","Progressive"))
df$Glu=as.numeric(df$Glu)
pdf(file = "Figure1B-Blood glucose.pdf",width = 7.5,height = 6)
ggplot(data=df,aes(x=Day,y=Glu,colour = Group))+ 
  geom_violin(alpha = 0.8,
              scale = 'width',trim = TRUE)+
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA)+
  geom_jitter(alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Non progressive","Progressive"), 
                    values =c("#9392BE","#51B1B7"))+
  scale_color_manual(limits=c("Non progressive","Progressive"), 
                     values=c("#9392BE","#51B1B7"))+ 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  stat_compare_means(method = "t.test",paired = F,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "")),label = "p.signif",
                     size = 5, 
                     comparisons=list(c("Day 0","Day 7")))+
  theme_bw()+
  labs(x="",y="Blood glucose")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))+
  facet_wrap2(vars(Group), ncol = 2,
              strip  = strip_nested(background_x = elem_list_rect(fill = cols,color = NA)))

dev.off()


####MNCV
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "Figure1C-MNCV.pdf",width = 5.5,height = 6)
ggplot(data=data,aes(x=Group,y=MNCV,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=MNCV,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=MNCV,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Non progressive","Progressive"), 
                    values =c("#9392BE","#51B1B7"))+
  scale_color_manual(limits=c("Non progressive","Progressive"), 
                     values=c("#9392BE","#51B1B7"))+ 
  geom_signif(mapping=aes(x=Group,y=MNCV), 
              comparisons = list(c("Non progressive","Progressive")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(60), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="MNCV (m/s)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15),
        legend.position = "top")
dev.off()

####SNCV
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "Figure1C-SNCV.pdf",width = 5.5,height = 6)
ggplot(data=data,aes(x=Group,y=SNCV,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=SNCV,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=SNCV,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Non progressive","Progressive"), 
                    values =c("#9392BE","#51B1B7"))+
  scale_color_manual(limits=c("Non progressive","Progressive"), 
                     values=c("#9392BE","#51B1B7"))+ 
  geom_signif(mapping=aes(x=Group,y=SNCV), 
              comparisons = list(c("Non progressive","Progressive")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(40), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="SNCV (m/s)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15),
        legend.position = "top")
dev.off()

#######Hot plate test 
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "Figure1C-Hot plate test .pdf",width =5.5,height = 6)
ggplot(data=data,aes(x=Group,y=latency,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=latency,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=latency,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Non progressive","Progressive"), 
                    values =c("#9392BE","#51B1B7"))+
  scale_color_manual(limits=c("Non progressive","Progressive"), 
                     values=c("#9392BE","#51B1B7"))+ 
  geom_signif(mapping=aes(x=Group,y=latency), 
              comparisons = list(c("Non progressive","Progressive")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(30), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="Withdrawal latency (s)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15),
        legend.position = "top")
dev.off()


#######von Frey test
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "Figure1C-von Frey test.pdf",width = 5.5,height = 6)
ggplot(data=data,aes(x=Group,y=Threshold,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=Threshold,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Threshold,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Non progressive","Progressive"), 
                    values =c("#9392BE","#51B1B7"))+
  scale_color_manual(limits=c("Non progressive","Progressive"), 
                     values=c("#9392BE","#51B1B7"))+ 
  geom_signif(mapping=aes(x=Group,y=Threshold), 
              comparisons = list(c("Non progressive","Progressive")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(6), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="Threshold (g)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15),
        legend.position = "top")
dev.off()


#####################Figure 2
####glucose
library(ggplot2)
library(ggsignif)
cols=c("#A6CEE4","#1F78B4","#B2E08A","#30A12C",
       "#FEBF6E","#FF7F00","#FB9B98","#E31A20","#C7E2E4","#8FC2C7",
       "#A6CEE4","#1F78B4","#B2E08A","#30A12C",
       "#FEBF6E","#FF7F00","#FB9B98","#E31A20","#C7E2E4","#8FC2C7")
df <- read.table("clipboard",header = T,sep = "\t")
df$Group <- factor(df$Group,levels = c("Control","DPN"))
df$Glu=as.numeric(df$Glu)
pdf(file = "Blood glucose.pdf",width = 7.5,height = 6)
ggplot(data=df,aes(x=Group,y=Glu,colour = Group))+ 
  geom_violin(alpha = 0.8,
              scale = 'width',trim = TRUE)+
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA)+
  geom_jitter(alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Control","DPN"), 
                    values =c( "#f0b87f","#e59698"))+
  scale_color_manual(limits=c("Control","DPN"), 
                     values=c("#f0b87f","#e59698"))+ 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  stat_compare_means(method = "t.test",paired = F,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "")),label = "p.signif",
                     size = 5, 
                     comparisons=list(c("Control","DPN")))+
  theme_bw()+
  labs(x="",y="Blood glucose")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))+
  facet_wrap2(vaMNCV(Day), ncol = 2,
              strip  = strip_nested(background_x = elem_list_rect(fill = cols,color = NA)))

dev.off()

####MNCV
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "MNCV.pdf",width = 6.5,height = 6)
ggplot(data=data,aes(x=Group,y=MNCV,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=MNCV,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=MNCV,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Control","DPN"), 
                    values =c( "#f0b87f","#e59698","#abd0cd"))+
  scale_color_manual(limits=c("Control","DPN"), 
                     values=c("#f0b87f","#e59698","#abd0cd"))+ 
  geom_signif(mapping=aes(x=Group,y=MNCV), 
              comparisons = list(c("Control","DPN")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(60), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="MNCV (m/s)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()

####SNCV
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "SNCV.pdf",width = 6.5,height = 6)
ggplot(data=data,aes(x=Group,y=SNCV,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=SNCV,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=SNCV,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Control","DPN"), 
                    values =c( "#f0b87f","#e59698","#abd0cd"))+
  scale_color_manual(limits=c("Control","DPN"), 
                     values=c("#f0b87f","#e59698","#abd0cd"))+ 
  geom_signif(mapping=aes(x=Group,y=SNCV), 
              comparisons = list(c("Control","DPN")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(40), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="SNCV (m/s)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()

#######Hot plate test 
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "Hot plate test .pdf",width = 6.5,height = 6)
ggplot(data=data,aes(x=Group,y=latency,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=latency,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=latency,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Control","DPN"), 
                    values =c( "#f0b87f","#e59698","#abd0cd"))+
  scale_color_manual(limits=c("Control","DPN"), 
                     values=c("#f0b87f","#e59698","#abd0cd"))+ 
  geom_signif(mapping=aes(x=Group,y=latency), 
              comparisons = list(c("Control","DPN")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(30), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="Withdrawal latency (s)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()


#######von Frey test
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "von Frey test.pdf",width = 6.5,height = 6)
ggplot(data=data,aes(x=Group,y=Threshold,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=Threshold,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Threshold,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Control","DPN"), 
                    values =c( "#f0b87f","#e59698","#abd0cd"))+
  scale_color_manual(limits=c("Control","DPN"), 
                     values=c("#f0b87f","#e59698","#abd0cd"))+ 
  geom_signif(mapping=aes(x=Group,y=Threshold), 
              comparisons = list(c("Control","DPN")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(6), 
              size=1, 
              textsize = 7, 
              test = "t.test")+ 
  theme_bw()+
  labs(x="",y="Threshold (g)")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()


##########Figure8
###PCR
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
data$Group=factor(data$Group,levels = c("Control","Low MNCV","High MNCV"))
pdf(file = "DRG PCR.pdf",width = 7.5,height = 6)
ggplot(data=data,aes(x=Group,y=Exp,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=Exp,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Exp,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Control","Low MNCV","High MNCV"), 
                    values =c( "#f0b87f","#9392BE","#51B1B7"))+
  scale_color_manual(limits=c("Control","Low MNCV","High MNCV"), 
                     values=c("#f0b87f","#9392BE","#51B1B7"))+ 
  geom_signif(mapping=aes(x=Group,y=Exp), 
              comparisons = list(c("Control","Low MNCV"),c("Low MNCV","High MNCV"),c("Control","High MNCV")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(1.5,1.4,1.3), 
              size=1, 
              textsize = 7, 
              test = "t.test")+
  theme_bw()+
  labs(x="",y="Relative expression")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()

###
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
data$Group=factor(data$Group,levels = c("Control","Low MNCV","High MNCV"))
pdf(file = "SN PCR.pdf",width = 7.5,height = 6)
ggplot(data=data,aes(x=Group,y=Exp,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=Exp,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Exp,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Control","Low MNCV","High MNCV"), 
                    values =c( "#f0b87f","#9392BE","#51B1B7"))+
  scale_color_manual(limits=c("Control","Low MNCV","High MNCV"), 
                     values=c("#f0b87f","#9392BE","#51B1B7"))+ 
  geom_signif(mapping=aes(x=Group,y=Exp), 
              comparisons = list(c("Control","Low MNCV"),c("Low MNCV","High MNCV"),c("Control","High MNCV")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(1.5,1.4,1.3), 
              size=1, 
              textsize = 7, 
              test = "t.test")+
  theme_bw()+
  labs(x="",y="Relative expression")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()

########Tissue Timm23 IF
###
library(ggplot2)
library(ggpubr)
dt=read.table("clipboard",sep = "\t",header = T,check.names = F)
dt$Group=factor(dt$Group,levels = c("Control","Low MNCV","High MNCV"))
dt$Timm23=as.numeric(dt$Timm23)
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
mytheme2<-theme_classic()+
  theme(text=element_text(size = 15),
        legend.position = "top",
        axis.line = element_line(size = 0.6),
        axis.ticks = element_line(size = 0.6),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
mycol=c("#f0b87f","#9392BE","#51B1B7")
p1<-ggbarplot(dt, x = "Group",
              y= "Timm23",
              color= "Group",
              fill= "Group",
              alpha=0.9,
              width=0.5,
              position= position_dodge(0.65),
              add= "mean_se",add.params = list(color = "black"))+
  geom_jitter(width=0.1,shape = 1,size = 2.5)+
  geom_signif(comparisons = list(c("Control","Low MNCV"),c("Low MNCV","High MNCV")),
              map_signif_level = T, 
              test = t.test,vjust = 0.7,
              tip_length = c(0,0,0,0,0,0),y_position = c(85,90),
              size=1,color="black",textsize = 8)+
  xlab("")+ylab("MFI of Timm23")+
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,100),
                     breaks=seq(0,100,10))+
  scale_colour_manual(values=mycol)+
  scale_fill_manual(values=mycol)+
  mytheme2
pdf(file = "DRG IF.pdf",width = 7.5,height = 7.5)
p1
dev.off()

###
dt=read.table("clipboard",sep = "\t",header = T,check.names = F)
dt$Group=factor(dt$Group,levels = c("Control","Low MNCV","High MNCV"))
dt$Timm23=as.numeric(dt$Timm23)
p1<-ggbarplot(dt, x = "Group",
              y= "Timm23",
              color= "Group",
              fill= "Group",
              alpha=0.9,
              width=0.5,
              position= position_dodge(0.65),
              add= "mean_se",add.params = list(color = "black"))+
  geom_jitter(width=0.1,shape = 1,size = 2.5)+
  geom_signif(comparisons = list(c("Control","Low MNCV"),c("Low MNCV","High MNCV")),
              map_signif_level = T, 
              test = t.test,vjust = 0.7,
              tip_length = c(0,0,0,0,0,0),y_position = c(130,135),
              size=1,color="black",textsize = 8)+
  xlab("")+ylab("MFI of Timm23")+
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,150),
                     breaks=seq(0,150,15))+
  scale_colour_manual(values=mycol)+
  scale_fill_manual(values=mycol)+
  mytheme2
pdf(file = "SN IF.pdf",width = 7.5,height = 7.5)
p1
dev.off()


###############Figure9
####PCR
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
data$Group=factor(data$Group,levels = c("Mock","Vector","oe-Timm23","oe-Timm23-"))
pdf(file = "oe-Timm23 PCR.pdf",width = 8.5,height = 6)
ggplot(data=data,aes(x=Group,y=Exp,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=Exp,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Exp,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Mock","Vector","oe-Timm23","oe-Timm23-"), 
                    values =c( "#7ABFE2","#FFDE9C","#B1ADD6","#C9CBE6"))+
  scale_color_manual(limits=c("Mock","Vector","oe-Timm23","oe-Timm23-"), 
                     values=c( "#7ABFE2","#FFDE9C","#B1ADD6","#C9CBE6"))+ 
  geom_signif(mapping=aes(x=Group,y=Exp), 
              comparisons = list(c("Mock","Vector"),c("Mock","oe-Timm23"),c("Mock","oe-Timm23-")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(2.3,2.4,2.5), 
              size=1, 
              textsize = 7, 
              test = "t.test")+
  theme_bw()+
  labs(x="",y="Relative expression")+ 
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.text.x= element_text(size = 15),axis.text.y=element_text(size=15),strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
dev.off()




####ATP
dt=read.table("clipboard",sep = "\t",header = T,check.names = F)
dt$Group=factor(dt$Group,levels = c("Mock","Vector","oe-Timm23","oe-Timm23-"))
dt$ATP=as.numeric(dt$ATP)
mycol=c("#7ABFE2","#FFDE9C","#B1ADD6","#C9CBE6")
p1<-ggbarplot(dt, x = "Group",
              y= "ATP",
              color= "Group",
              fill= "Group",
              alpha=0.9,
              width=0.5,
              position= position_dodge(0.65),
              add= "mean_se",add.params = list(color = "black"))+
  geom_jitter(width=0.1,shape = 1,size = 2.5)+
  geom_signif(comparisons = list(c("Mock","Vector"),c("Mock","oe-Timm23"),c("Mock","oe-Timm23-")),
              map_signif_level = T, 
              test = t.test,vjust = 0.7,
              tip_length = c(0,0,0,0,0,0),y_position = c(2,2.1,2.2),
              size=1,color="black",textsize = 8)+
  xlab("")+ylab("ATP (Fold change of Control)")+
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,2.5),
                     breaks=seq(0,2.5,0.5))+
  scale_colour_manual(values=mycol)+
  scale_fill_manual(values=mycol)+
  mytheme2
pdf("ATP-PA.pdf",width = 11,height = 9)
p1
dev.off()


######mitotracker
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "Number of mitochondrion.pdf",width = 10,height = 10)
ggplot(data=data,aes(x=Group,y=Num,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=Num,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Num,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Mock","Vector","oe-Timm23","oe-Timm23-"), 
                    values =c( "#7ABFE2","#FFDE9C","#B1ADD6","#C9CBE6"))+
  scale_color_manual(limits=c("Mock","Vector","oe-Timm23","oe-Timm23-"), 
                     values=c( "#7ABFE2","#FFDE9C","#B1ADD6","#C9CBE6"))+ 
  geom_signif(mapping=aes(x=Group,y=Num), 
              comparisons = list(c("Mock","Vector"),c("Mock","oe-Timm23"),c("Mock","oe-Timm23-")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(23,24,25), 
              size=1, 
              textsize = 7, 
              test = "t.test")+
  theme_bw()+
  labs(x="",y="Number of mitochondrion/per cell")+ 
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x= element_text(size = 22),axis.text.y=element_text(size=22),strip.text = element_text(size = 22),
        legend.title = element_text(size = 22),legend.text = element_text(size = 22),legend.position = "top")
dev.off()

######mitoSOX
data=read.table("clipboard",sep = "\t",header = T,check.names = F)
pdf(file = "MFI of mitoSOX.pdf",width = 10,height = 10)
ggplot(data=data,aes(x=Group,y=MFI,colour = Group))+ 
  geom_violin(
    alpha = 0.8,
    scale = 'width',
    
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Group,y=MFI,colour=Group,fill=Group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=MFI,colour = Group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Mock","Vector","oe-Timm23","oe-Timm23-"), 
                    values =c( "#7ABFE2","#FFDE9C","#B1ADD6","#C9CBE6"))+
  scale_color_manual(limits=c("Mock","Vector","oe-Timm23","oe-Timm23-"), 
                     values=c( "#7ABFE2","#FFDE9C","#B1ADD6","#C9CBE6"))+ 
  geom_signif(mapping=aes(x=Group,y=MFI), 
              comparisons = list(c("Mock","Vector"),c("Mock","oe-Timm23"),c("Mock","oe-Timm23-")),
              map_signif_level=T,
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(108,114,120), 
              size=1, 
              textsize = 7, 
              test = "t.test")+
  theme_bw()+
  labs(x="",y="MFI of mitoSOX")+ 
  theme(axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),axis.text.x= element_text(size = 22),axis.text.y=element_text(size=22),strip.text = element_text(size = 22),
        legend.title = element_text(size = 22),legend.text = element_text(size = 22),legend.position = "top")
dev.off()



