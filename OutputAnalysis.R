Cancer_Type<-c(rep("PRAD",2),rep("BRCA",2),rep("UCEC",2),rep("COAD",2),rep("READ",2))
Reg_Type<-c(rep(c("Up","Down"),5))
Value<-c(34,27,71,87,142,103,186,153,165,114)
Dataset1<-data.frame(Cancer_Type,Reg_Type,Value)

Table_ceOutput_ALL<-data.frame(table(ceOutput_ALL$lncRNA_Gene_miRNA))
Sig_5times<-Table_ceOutput_ALL[Table_ceOutput_ALL$Freq>4,]
Sig_4times<-Table_ceOutput_ALL[Table_ceOutput_ALL$Freq>3&&Table_ceOutput_ALL$Freq<5,]
Sig_3times<-Table_ceOutput_ALL[Table_ceOutput_ALL$Freq>2&&Table_ceOutput_ALL$Freq<4,]
Sig_2times<-Table_ceOutput_ALL[Table_ceOutput_ALL$Freq>1&&Table_ceOutput_ALL$Freq<3,]
###Convert Ensembl ID into Gene Symbol
library("biomaRt")
mart<-useDataset("hsapiens_gene_ensembl",useMart("ensembl"))
genes_PRAD<-ceOutput_PRAD$Genes
Gene_Symbol_PRAD<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_PRAD,mart=mart)
genes_BRCA<-ceOutput_BRCA$Genes
Gene_Symbol_BRCA<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_BRCA,mart=mart)
genes_UCEC<-ceOutput_UCEC$Genes
Gene_Symbol_UCEC<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_UCEC,mart=mart)
genes_COAD<-ceOutput_COAD$Genes
Gene_Symbol_COAD<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_COAD,mart=mart)
genes_READ<-ceOutput_READ$Genes
Gene_Symbol_READ<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_READ,mart=mart)

###Convert ncRNA ID into Symbol
genes_PRAD<-ceOutput_PRAD$lncRNAs
LncRNA_Gene_Symbol_PRAD<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_PRAD,mart=mart)
genes_BRCA<-ceOutput_BRCA$lncRNAs
LncRNA_Gene_Symbol_BRCA<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_BRCA,mart=mart)
genes_UCEC<-ceOutput_UCEC$lncRNAs
LncRNA_Gene_Symbol_UCEC<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_UCEC,mart=mart)
genes_COAD<-ceOutput_COAD$lncRNAs
LncRNA_Gene_Symbol_COAD<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_COAD,mart=mart)
genes_READ<-ceOutput_READ$lncRNAs
LncRNA_Gene_Symbol_READ<-getBM(filters = "ensembl_gene_id",attributes=c("ensembl_gene_id","entrezgene_accession","description"),values=genes_READ,mart=mart)







install.packages("SuperExactTest")
gene_list<-list(PRAD=Gene_Symbol_PRAD$entrezgene_accession,BRCA=Gene_Symbol_BRCA$entrezgene_accession,UCEC=Gene_Symbol_UCEC$entrezgene_accession,COAD=Gene_Symbol_COAD$entrezgene_accession,READ=Gene_Symbol_READ$entrezgene_accession)
(length.gene.sets=sapply(gene_list, length))
#PRAD BRCA UCEC COAD READ 
#42  196  212   71   44 
#Total 565
total=565
num.expected.overlap<-total*do.call(prod,as.list(length.gene.sets/total))
p<-sapply(0:42,function(i) dpsets(i,length.gene.sets,n=total))
