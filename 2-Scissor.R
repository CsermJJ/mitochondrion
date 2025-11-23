######################
library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
genehuman2mouse = getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",values = bulk_dataset$gene,mart = human,attributesL=c("mgi_symbol"),martL = mouse,uniqueRows = T)
genehuman2mouse=genehuman2mouse[!duplicated(genehuman2mouse$HGNC.symbol),]
rownames(genehuman2mouse)=genehuman2mouse[,1]
rownames(bulk_dataset)=bulk_dataset[,1]
q=intersect(rownames(bulk_dataset),rownames(genehuman2mouse))
bulk_dataset=bulk_dataset[q,]
genehuman2mouse=genehuman2mouse[q,]
all=cbind(genehuman2mouse,bulk_dataset)
all=all[!duplicated(all$MGI.symbol),]
rownames(all)=all[,2]
all=all[,-c(1:3)]
library(Scissor)
library(preprocessCore)
library(scAB)
bulk_dataset=all
bulk_dataset=as.data.frame(bulk_dataset)
bulk_dataset[is.na(bulk_dataset)]<-0
phenotype=read.table("group.txt",sep="\t",header = T,check.names = F)
phenotype1 <- as.numeric(phenotype[,2])
names(phenotype1)<-phenotype[,1]
sc_dataset=scRNA.ann
tag<-c("0","1")
infos1<-Scissor(bulk_dataset,sc_dataset,phenotype1,alpha=0.05,tag=tag,
                family="binomial",Save_file='./Scissor_progressive.RData')

Scissor_select<-rep(0,ncol(sc_dataset))
names(Scissor_select)<-colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos]<-"Scissor+"
Scissor_select[infos1$Scissor_neg]<-"Scissor-"
sc_dataset<-AddMetaData(sc_dataset,metadata=Scissor_select,col.name="scissor")
pdf("3.1_scissor.pdf",width = 6.5,height = 5.5)
DimPlot(sc_dataset,reduction='umap',
        group.by='scissor',order = T,
        cols=c('#D9F5F4','#1C6AB1','#A04294'),
        pt.size=0.001)
dev.off()


library(gplots)
pdf("1.2_scissor_balloon.pdf",width = 8,height = 7)
balloonplot(table(sc_dataset$scissor,sc_dataset$Celltype))
dev.off()