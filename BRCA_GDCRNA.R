project_BRCA<-'TCGA-BRCA'
rnadir_BRCA<-paste(project_BRCA,'RNAseq',sep="/")
mirdir_BRCA<-paste(project_BRCA,"miRNAs",sep="/")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-BRCA',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir_BRCA)
gdcRNADownload(project.id = "TCGA-BRCA",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir_BRCA)
#Download clinical data
clinicaldir_BRCA<-paste(project_BRCA,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-BRCA",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir_BRCA)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_BRCA<-gdcParseMetadata(project.id = "TCGA-BRCA",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_BRCA<-gdcFilterDuplicate(metaMatrix.RNA_BRCA)
metaMatrix.RNA_BRCA<-gdcFilterSampleType(metaMatrix.RNA_BRCA)

####Parse miRNA metadata
metaMatrix.MIR_BRCA<-gdcParseMetadata(project.id = "TCGA-BRCA",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_BRCA<-gdcFilterDuplicate(metaMatrix.MIR_BRCA)
metaMatrix.MIR_BRCA<-gdcFilterSampleType(metaMatrix.MIR_BRCA)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_BRCA<-gdcRNAMerge(metadata = metaMatrix.RNA_BRCA,path = rnadir_BRCA,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_BRCA<-gdcRNAMerge(metadata = metaMatrix.MIR_BRCA,path = mirdir_BRCA,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_BRCA<-gdcClinicalMerge(path = clinicaldir_BRCA,key.info = TRUE)

ClinicalDa_BRCA[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_BRCA<-gdcVoomNormalization(counts = rnaCounts_BRCA,filter = FALSE)
####Normalization of miRNAs data
mirExpr_BRCA<-gdcVoomNormalization(counts = mirCounts_BRCA,filter = FALSE)
####Differential gene expression analysis
DEGAll_BRCA<-gdcDEAnalysis(counts = rnaCounts_BRCA,group = metaMatrix.RNA_BRCA$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
DEGMIR_BRCA<-gdcDEAnalysis(counts = mirCounts_BRCA,group = metaMatrix.MIR_BRCA$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma")
####As I need log2FC for all set of genes I added a modification for gdcDEAnlaysis
DEGMIR_BRCA2<-gdcDEAnalysis(counts = mirCounts_BRCA,group = metaMatrix.MIR_BRCA$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
setDT(DEGMIR_BRCA2, keep.rownames = TRUE)[]
demiR_BRCA2<-gdcDEReport(DEGMIR_BRCA2,fc=0.05)
data("DEGAll_BRCA")
####All DEGs
deALL_BRCA<-gdcDEReport(deg = DEGAll_BRCA,gene.type = "all")
####DE long-noncoding
deLNC_BRCA<-gdcDEReport(deg = DEGAll_BRCA,gene.type = "long_non_coding")
####DE protein coding genes
dePC_BRCA<-gdcDEReport(deg=DEGAll_BRCA,gene.type = "protein_coding")
##DE pseudogenes
dePseudo_BRCA<-gdcDEReport(deg=DEGAll_BRCA,gene.type = "pseudogene")
demiR_BRCA<-gdcDEReport(deg = DEGMIR_BRCA)
################################################################################################
DEGALL_BRCA_MIR<-gdcDEAnalysis(counts = mirCounts_BRCA,group = metaMatrix.MIR_BRCA$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
gdcVolcanoPlot(DEGALL_BRCA_MIR)
gdcHeatmap(deg.id = rownames(DEGALL_BRCA_MIR),metadata = metaMatrix.MIR_BRCA,rna.expr = mirExpr_BRCA)

####DEG visualisation
gdcVolcanoPlot(DEGMIR_BRCA)
demiR_BRCA<-gdcDEReport(deg=DEGMIR_BRCA,gene.type = "all")
gdcBarPlot(deg=demiR_BRCA,angle = 45,data.type = "miRNAs")
gdcBarPlot(deg = deALL_BRCA,angle = 45,data.type = "RNAseq")
gdcHeatmap(deg.id = rownames(demiR_BRCA),metadata = metaMatrix.MIR_BRCA,rna.expr = mirExpr_BRCA)
degName=rownames(deALL_BRCA)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_BRCA,rna.expr = rnaExpr_BRCA)
gdcBarPlot(deg=demiR_BRCA,angle=45,data.type="miRNAs")
t1<-table(which(demiR_BRCA$logFC<0))
dim(t1)
t2<-table(which(demiR_BRCA$logFC>0))
dim(t2)
##ceRNAs network analysis

ceOutput_BRCA<- gdcCEAnalysis(lnc= rownames(deLNC_BRCA),pc= rownames(dePC_BRCA),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_BRCA,mir.expr= mirExpr_BRCA)
ceOutput_BRCA2<- gdcCEAnalysis(lnc= rownames(deLNC_BRCA),pc= rownames(dePC_BRCA),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_BRCA,mir.expr= mirExpr_BRCA)
####Network Visualization in Cytoscape

ceOutput2_BRCA<-ceOutput_BRCA[ceOutput_BRCA$hyperPValue<0.01 & ceOutput_BRCA$corPValue<0.01 & ceOutput_BRCA$regSim !=0,]
ceOutput2_BRCA2<-ceOutput_BRCA2[ceOutput_BRCA2$hyperPValue<0.01 & ceOutput_BRCA2$corPValue<0.01 & ceOutput_BRCA2$regSim !=0,]
edges_BRCA<-gdcExportNetwork(ceNetwork = ceOutput2_BRCA,net = "edges")
nodes_BRCA<-gdcExportNetwork(ceNetwork = ceOutput2_BRCA,net = "nodes")

write.table(edges_BRCA,file = "edges_BRCA.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_BRCA,file="nodes_BRCA.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(ceOutput2_BRCA,file="ceOutput_BRCA.txt",sep="\t",quote = FALSE,row.names=FALSE)
head(edges)

shinyCorPlot(gene1    = rownames(dePC), 
             gene2    = rownames(demiR), 
             rna.expr = mirExpr, 
             metadata = metaMatrix.MIR)


shinyCorPlot(gene1    = 'hsa-miR-7a-5p', 
             gene2    = 'hsa-miR-7a-3p', 
             rna.expr = mirExpr, 
             metadata = metaMatrix.MIR)

p_ENSG00000137726_BRCA<-gdcKMPlot(gene="ENSG00000137726",rna.expr = rnaExpr_BRCA,metadata = metaMatrix.RNA_BRCA,sep= 'median')
p_ENSG00000137726_BRCA<-p_ENSG00000137726_BRCA$plot+theme(legend.position = c(0.2,0.2))
#p_ENSG00000065534<-gdcKMPlot(gene = "ENSG00000065534",rna.expr = rnaExpr_BRCA,metadata = metaMatrix.RNA_BRCA,sep = 'median')
p_ENSG00000135272<-gdcKMPlot(gene="ENSG00000135272",rna.expr = rnaExpr_BRCA,metadata = metaMatrix.RNA_BRCA,sep= 'median')
p_ENSG00000135272<-p_ENSG00000135272$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000137726$plot|p_ENSG00000135272
#p_ENSG00000095637<-gdcKMPlot(gene = "ENSG00000095637",rna.expr = rnaExpr_BRCA,metadata = metaMatrix.RNA_BRCA,sep = 'median')

#i was able to change legend postion
library(ggplot2)
#par(mfrow=c(2,2))
layout(matrix(c(1,2), 1,1, byrow = TRUE))
p_ENSG00000137726$plot+theme(legend.position = c(0.2,0.2))
#p_ENSG00000065534$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000135272$plot+theme(legend.position = c(0.2,0.2))
#p_ENSG00000095637$plot+theme(legend.position = c(0.2,0.2))

####Functional Enrichment Analysis
###Enrichment analysis
library(GDCRNATools)
enrichOutput_BRCA<- gdcEnrichAnalysis(gene = rownames(deALL_BRCA), simplify = TRUE)
enrichoutput_BRCA<-enrichOutput_BRCA %>% separate(Terms,c("Terms2","Terms"),sep="~")# I removed Term code as it disturbs for plot space
BRCA_GO<-gdcEnrichPlot(enrichoutput_BRCA, type = 'bar', category = 'GO',num.terms = 10)
BRCA_GO<-BRCA_GO+theme_update(text=element_text(size = 10))+theme_bw()+ylab("-log10(FDR)-BRCA")+theme(legend.position ="none")
BRCA_KEGG<-gdcEnrichPlot(enrichoutput_BRCA, type = 'bubble', category = 'KEGG',num.terms = 10)+ylab("Fold enrichment-BRCA")+theme(legend.position = "none")+theme_update(text=element_text(size = 10))+theme_bw()
BRCA_DO<-gdcEnrichPlot(enrichoutput_BRCA, type = 'bar', category = 'DO',num.terms = 10)
BRCA_DO<-BRCA_DO+ylab("-log10(FDR)-BRCA")+geom_bar(stat = "identity",fill="#1C9C5A")+theme_update(text=element_text(size = 15))+theme_bw()









