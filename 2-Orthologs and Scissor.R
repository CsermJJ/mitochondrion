####################one-to-one orthologs

library(biomaRt)
human <- useMart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
rat <- useMart(biomart = "ensembl",dataset = "rnorvegicus_gene_ensembl")
human_genes <- bulk_dataset[, 1]
orthologs <- getLDS(
  attributes  = c("hgnc_symbol"),
  filters     = "hgnc_symbol",
  values      = human_genes,
  mart        = human,
  attributesL = c("rgd_symbol"),
  martL       = rat
)
colnames(orthologs) <- c("human_gene", "rat_gene")
orthologs <- orthologs %>%
  filter(human_gene != "", rat_gene != "") %>%
  distinct(human_gene, .keep_all = TRUE)

bulk_df <- as.data.frame(bulk_dataset)
colnames(bulk_df)[1] <- "human_gene"
bulk_rat <- bulk_df %>%
  inner_join(orthologs, by = "human_gene")
bulk_rat <- bulk_rat %>%
  select(rat_gene, everything(), -human_gene)
bulk_rat_mat <- as.data.frame(bulk_rat)
rownames(bulk_rat_mat) <- bulk_rat_mat$rat_gene
bulk_rat_mat$rat_gene <- NULL
bulk_dataset=bulk_rat_mat

###############Scissor
library(Scissor)
library(preprocessCore)
library(scAB)
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
pdf("Scissor.pdf",width = 6.5,height = 5.5)
DimPlot(sc_dataset,reduction='umap',
        group.by='scissor',order = T,
        cols=c('#D9F5F4','#1C6AB1','#A04294'),
        pt.size=0.001)
dev.off()


library(gplots)
pdf("Scissor_balloon.pdf",width = 8,height = 7)
balloonplot(table(sc_dataset$scissor,sc_dataset$Celltype))
dev.off()


