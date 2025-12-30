####
library(DESeq2)
library(BiocParallel)
library(reshape2)
library(tximport)
library(readr)
rt=read.table("DRG exp.txt",header = T,sep = "\t",row.names =1,check.names=F )
rt=as.matrix(rt)
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=round(as.matrix(rt))
data=rt                                      
data=data[rowMeans(rt)>1,]                        
group=read.table("DRG group.txt",header = T,sep = "\t")
condition=group$condition
colData=data.frame(condition=as.factor(condition))
row.names(colData)=colnames(data)
dds <- DESeqDataSetFromMatrix(countData=data, colData=colData, design=~condition)
dds_norm <- DESeq(dds)
normalized_count=counts(dds_norm,normalized=T)
normalized_counts_mad=apply(normalized_count,1,mad)
normalized_counts=normalized_count[order(normalized_counts_mad,decreasing = T),]
write.table(normalized_counts,file = "DRG_normalized.txt",quote = F,sep = "\t",row.names = T,col.names = T)
rld=rlog(dds_norm,blind = F)
rlogMat=assay(rld)
rlogMat=rlogMat[order(normalized_counts_mad,decreasing=T),]
write.table(rlogMat,file = "DRG_nromalized_vst.txt",quote=F,sep = "\t",row.names = T,col.names = T)
res=results(dds_norm,contrast = c("condition","NP_DPN","P_DPN"))
res=res[order(res$pvalue),]
summary(res)
write.table(res,file="DRG diff-All.txt",sep = "\t")



###
library('GSEABase')
library(GSVA)
library(msigdbr)
library("org.Mm.eg.db")
library(clusterProfiler)
library(GseaVis)
exp=read.table("DRG diff-All.txt",sep = "\t",check.names = F,header = T)
exp1 <- arrange(exp,desc(logFC))
gene_list <- exp1$logFC
names(gene_list) <- exp1$SYMBOL
m_df = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
geneset = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
gg <- GSEA(gene_list, TERM2GENE=geneset,verbose=F,
           pvalueCutoff=0.1, pAdjustMethod = "BH")
sortgg<- gg[order(gg$NES, decreasing = T),]
sortgg<- sortgg[sortgg$p.adjust <0.05,]
write.table(sortgg,"GSEA_gobp_DRG.txt",sep = "\t",quote = F,row.names = F)

setid <- c("GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY","GOBP_ATP_METABOLIC_PROCESS",
           "GOBP_OXIDATIVE_PHOSPHORYLATION",
           "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT","GOBP_MITOCHONDRIAL_GENE_EXPRESSION")
pdf(file="GSEA_gobp_DRG.pdf",width=11,height=8)
gseaNb(object = gg,
       geneSetID = setid,
       curveCol = ggsci::pal_npg()(7),
       subPlot = 2,
       addPval = T,
       pvalX = 1,
       pvalY = 1)
dev.off()


#####
rt=read.table("SN exp.txt",header = T,sep = "\t",row.names =1,check.names=F )
rt=as.matrix(rt)
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=round(as.matrix(rt))
data=rt                                      
data=data[rowMeans(rt)>1,]                        
group=read.table("SN group.txt",header = T,sep = "\t")
condition=group$condition
colData=data.frame(condition=as.factor(condition))
row.names(colData)=colnames(data)
dds <- DESeqDataSetFromMatrix(countData=data, colData=colData, design=~condition)
dds_norm <- DESeq(dds)
normalized_count=counts(dds_norm,normalized=T)
normalized_counts_mad=apply(normalized_count,1,mad)
normalized_counts=normalized_count[order(normalized_counts_mad,decreasing = T),]
write.table(normalized_counts,file = "SN_normalized.txt",quote = F,sep = "\t",row.names = T,col.names = T)
rld=rlog(dds_norm,blind = F)
rlogMat=assay(rld)
rlogMat=rlogMat[order(normalized_counts_mad,decreasing=T),]
write.table(rlogMat,file = "SN_nromalized_vst.txt",quote=F,sep = "\t",row.names = T,col.names = T)
res=results(dds_norm,contrast = c("condition","NP_DPN","P_DPN"))
res=res[order(res$pvalue),]
summary(res)
write.table(res,file="SN diff-All.txt",sep = "\t")

###
exp=read.table("SN diff-All.txt",sep = "\t",check.names = F,header = T)
exp1 <- arrange(exp,desc(logFC))
gene_list <- exp1$logFC
names(gene_list) <- exp1$SYMBOL
gg <- GSEA(gene_list, TERM2GENE=geneset,verbose=F,
           pvalueCutoff=0.5, pAdjustMethod = "BH")
sortgg<- gg[order(gg$NES, decreasing = T),]
sortgg<- sortgg[sortgg$p.adjust <0.05,]
write.table(sortgg,"GSEA_gobp_SN.txt",sep = "\t",quote = F,row.names = F)

setid <- c("GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY","GOBP_MITOCHONDRIAL_TRANSLATION",
           "GOBP_MITOCHONDRIAL_GENE_EXPRESSION",
           "GOBP_OXIDATIVE_PHOSPHORYLATION","GOBP_CELLULAR_RESPIRATION")
pdf(file="GSEA_gobp_SN.pdf",width=11,height=8)
gseaNb(object = gg,
       geneSetID = setid,
       curveCol = ggsci::pal_npg()(7),
       subPlot = 2,
       addPval = T,
       pvalX = 1,
       pvalY = 1)
dev.off()
