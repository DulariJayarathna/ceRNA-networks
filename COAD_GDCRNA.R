project<-'TCGA-COAD'
rnadir<-paste(project,"RNAseq",sep="/")
mirdir<-paste(project,"miRNAs",sep="/")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-COAD',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir)
gdcRNADownload(project.id = "TCGA-COAD",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir)
#Download clinical data
clinicaldir<-paste(project,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-COAD",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_COAD<-gdcParseMetadata(project.id = "TCGA-COAD",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_COAD<-gdcFilterDuplicate(metaMatrix.RNA_COAD)
metaMatrix.RNA_COAD<-gdcFilterSampleType(metaMatrix.RNA_COAD)

####Parse miRNA metadata
metaMatrix.MIR_COAD<-gdcParseMetadata(project.id = "TCGA-COAD",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_COAD<-gdcFilterDuplicate(metaMatrix.MIR_COAD)
metaMatrix.MIR_COAD<-gdcFilterSampleType(metaMatrix.MIR_COAD)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_COAD<-gdcRNAMerge(metadata = metaMatrix.RNA_COAD,path = rnadir,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_COAD<-gdcRNAMerge(metadata = metaMatrix.MIR_COAD,path = mirdir,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_COAD<-gdcClinicalMerge(path = clinicaldir,key.info = TRUE)
ClinicalDa_COAD[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_COAD<-gdcVoomNormalization(counts = rnaCounts_COAD,filter = FALSE)
####Normalization of miRNAs data
mirExpr_COAD<-gdcVoomNormalization(counts = mirCounts_COAD,filter = FALSE)
####Differential gene expression analysis
DEGAll_COAD<-gdcDEAnalysis(counts = rnaCounts_COAD,group = metaMatrix.RNA_COAD$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")

DEGMIR_COAD<-gdcDEAnalysis(counts = mirCounts_COAD,group = metaMatrix.MIR_COAD$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma")
####As I need log2FC for all set of genes I added a modification for gdcDEAnlaysis
DEGMIR_COAD2<-gdcDEAnalysis(counts = mirCounts_COAD,group = metaMatrix.MIR_COAD$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
setDT(DEGMIR_COAD2, keep.rownames = TRUE)[]
demiR_COAD2<-gdcDEReport(DEGMIR_COAD2,fc=0.05)
data("DEGAll")
####All DEGs
deALL_COAD<-gdcDEReport(deg = DEGAll_COAD,gene.type = "all")
####DE long-noncoding
deLNC_COAD<-gdcDEReport(deg = DEGAll_COAD,gene.type = "long_non_coding")
####DE protein coding genes
dePC_COAD<-gdcDEReport(deg=DEGAll_COAD,gene.type = "protein_coding")
##DE pseudogenes
dePseudo_COAD<-gdcDEReport(deg=DEGAll_COAD,gene.type = "pseudogene")
demiR_COAD<-gdcDEReport(DEGMIR_COAD)
################################################################################################
DEGALL_COAD_MIR<-gdcDEAnalysis(counts = mirCounts_COAD,group = metaMatrix.MIR_COAD$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
gdcVolcanoPlot(DEGALL_COAD_MIR)
gdcHeatmap(deg.id = rownames(DEGALL_COAD_MIR),metadata = metaMatrix.MIR_COAD,rna.expr = mirExpr_COAD)
demiR_COAD<-gdcDEReport(deg=DEGALL_COAD_MIR,gene.type = "miRNAs")
####DEG visualisation
gdcVolcanoPlot(DEGAll_COAD)
gdcBarPlot(deg = deALL_COAD,angle = 45,data.type = "RNAseq")
degName=rownames(deALL_COAD)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_COAD,rna.expr = rnaExpr_COAD)
gdcBarPlot(deg=demiR_COAD,angle=45,data.type="miRNAs")
t1<-table(which(demiR_COAD$logFC<0))
dim(t1)
t2<-table(which(demiR_COAD$logFC>0))
dim(t2)
##ceRNAs network analysis

ceOutput_COAD<- gdcCEAnalysis(lnc= rownames(deLNC_COAD),pc= rownames(dePC_COAD),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_COAD,mir.expr    = mirExpr_COAD)
ceOutput_COAD2<- gdcCEAnalysis(lnc= rownames(deLNC_COAD),pc= rownames(dePC_COAD),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_COAD,mir.expr    = mirExpr_COAD)

####Network Visualization in Cytoscape

ceOutput2_COAD<-ceOutput_COAD[ceOutput_COAD$hyperPValue<0.01 & ceOutput_COAD$corPValue<0.01 & ceOutput_COAD$regSim !=0,]
ceOutput2_COAD2<-ceOutput_COAD2[ceOutput_COAD2$hyperPValue<0.01 & ceOutput_COAD2$corPValue<0.01 & ceOutput_COAD2$regSim !=0,]

edges_COAD<-gdcExportNetwork(ceNetwork = ceOutput2_COAD,net = "edges")
nodes_COAD<-gdcExportNetwork(ceNetwork = ceOutput2_COAD,net = "nodes")

write.table(edges_COAD,file = "edges_COAD.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_COAD,file="nodes_COAD.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(ceOutput2_COAD,file="ceOutput_COAD.txt",sep="\t",quote = FALSE,row.names=FALSE)
head(edges_COAD)

shinyCorPlot(gene1    = rownames(deLNC_COAD), 
             gene2    = rownames(dePC_COAD), 
             rna.expr = rnaExpr_COAD, 
             metadata = metaMatrix.RNA_COAD)
####Survival Analysis
#####Cox-PH Model
survOutput_COAD_COXPH<- gdcSurvivalAnalysis(gene     = rownames(deALL_COAD), 
                                            method   = 'coxph', 
                                            rna.expr = rnaExpr_COAD, 
                                            metadata = metaMatrix.RNA_COAD)

#####Kaplan Meir Plot
survOutput_COAD_KM<- gdcSurvivalAnalysis(gene     = rownames(deALL_COAD), 
                                         method   = 'KM', 
                                         rna.expr = rnaExpr_COAD, 
                                         metadata = metaMatrix.RNA_COAD, 
                                         sep      = 'median')

gdcKMPlot(gene     = 'ENSG00000137331',
          rna.expr = rnaExpr_COAD,
          metadata = metaMatrix.RNA_COAD,
          sep      = 'median')

p_ENSG00000101955_COAD<-gdcKMPlot(gene="ENSG00000101955",rna.expr = rnaExpr_COAD,metadata = metaMatrix.RNA_COAD,sep= 'median')
p_ENSG00000131471<-gdcKMPlot(gene = "ENSG00000131471",rna.expr = rnaExpr_COAD,metadata = metaMatrix.RNA_COAD,sep = 'median')
p_ENSG00000173641<-gdcKMPlot(gene="ENSG00000173641",rna.expr = rnaExpr_COAD,metadata = metaMatrix.RNA_COAD,sep= 'median')
p_ENSG00000136160_COAD<-gdcKMPlot(gene = "ENSG00000136160",rna.expr = rnaExpr_COAD,metadata = metaMatrix.RNA_COAD,sep = 'median')
#p_ENSG00000173641<-gdcKMPlot(gene = "ENSG00000173641",rna.expr = rnaExpr_COAD,metadata = metaMatrix.RNA_COAD,sep = 'median')
#i was able to change legend postion
library(ggplot2)
#par(mfrow=c(2,2))
layout(matrix(c(1,2), 1,1, byrow = TRUE))
#p_ENSG00000101955$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000131471<-p_ENSG00000131471$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000173641<-p_ENSG00000173641$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000136160_COAD<-p_ENSG00000136160_COAD$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000136160+p_ENSG00000144642+p_ENSG00000137726+p_ENSG00000135272+p_ENSG00000131471+p_ENSG00000173641+p_ENSG00000136160+p_ENSG00000145936
#p_ENSG00000173641$plot+theme(legend.position = c(0.2,0.2))

library(GDCRNATools)
enrichOutput_COAD<- gdcEnrichAnalysis(gene = rownames(deALL_COAD), simplify = TRUE)
enrichoutput_COAD<-enrichOutput_COAD %>% separate(Terms,c("Terms2","Terms"),sep="~")# I removed Term code as it disturbs for plot space
COAD_GO<-gdcEnrichPlot(enrichoutput_COAD, type = 'bar', category = 'GO',num.terms = 10)
COAD_GO<-COAD_GO+theme_update(text=element_text(size = 15))+theme_bw()+theme(legend.position = "none")+ylab("-log10(FDR)-COAD")
COAD_KEGG<-gdcEnrichPlot(enrichoutput_COAD, type = 'bubble', category = 'KEGG',num.terms = 10)+ylab("Fold enrichment-COAD")+theme(legend.position = "none")+theme_update(text=element_text(size = 10))+theme_bw()
COAD_DO<-gdcEnrichPlot(enrichoutput_COAD, type = 'bar', category = 'DO',num.terms = 10)
COAD_DO<-COAD_DO+ylab("-log10(FDR)-COAD")+geom_bar(stat = "identity",fill="#1C9C5A")+theme_update(text=element_text(size = 15))+theme_bw()



