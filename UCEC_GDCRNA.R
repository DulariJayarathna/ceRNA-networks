project<-'TCGA-UCEC'
rnadir<-paste(project,"RNAseq",sep="/")
mirdir<-paste(project,"miRNAs",sep="/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("GDCRNATools",lib="C:/Users/n10136142/AppData/Local/Temp/RtmpMhnjtL/downloaded_packages")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-UCEC',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir)
gdcRNADownload(project.id = "TCGA-UCEC",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir)
#Download clinical data
clinicaldir<-paste(project,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-UCEC",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_UCEC<-gdcParseMetadata(project.id = "TCGA-UCEC",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_UCEC<-gdcFilterDuplicate(metaMatrix.RNA_UCEC)
metaMatrix.RNA_UCEC<-gdcFilterSampleType(metaMatrix.RNA_UCEC)

####Parse miRNA metadata
metaMatrix.MIR_UCEC<-gdcParseMetadata(project.id = "TCGA-UCEC",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_UCEC<-gdcFilterDuplicate(metaMatrix.MIR_UCEC)
metaMatrix.MIR_UCEC<-gdcFilterSampleType(metaMatrix.MIR_UCEC)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_UCEC<-gdcRNAMerge(metadata = metaMatrix.RNA_UCEC,path = rnadir,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_UCEC<-gdcRNAMerge(metadata = metaMatrix.MIR_UCEC,path = mirdir,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_UCEC<-gdcClinicalMerge(path = clinicaldir,key.info = TRUE)
ClinicalDa_UCEC[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_UCEC<-gdcVoomNormalization(counts = rnaCounts_UCEC,filter = FALSE)
####Normalization of miRNAs data
mirExpr_UCEC<-gdcVoomNormalization(counts = mirCounts_UCEC,filter = FALSE)
####Differential gene expression analysis
DEGAll_UCEC<-gdcDEAnalysis(counts = rnaCounts_UCEC,group = metaMatrix.RNA_UCEC$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
DEGMIR_UCEC<-gdcDEAnalysis(counts = mirCounts_UCEC,group = metaMatrix.MIR_UCEC$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma")
####As I need log2FC for all set of genes I added a modification for gdcDEAnlaysis
DEGMIR_UCEC2<-gdcDEAnalysis(counts = mirCounts_UCEC,group = metaMatrix.MIR_UCEC$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
setDT(DEGMIR_UCEC2, keep.rownames = TRUE)[]
demiR_UCEC2<-gdcDEReport(DEGMIR_UCEC2,fc=0.05)

data("DEGAll")
####All DEGs
deALL_UCEC<-gdcDEReport(deg = DEGAll_UCEC,gene.type = "all")
####DE long-noncoding
deLNC_UCEC<-gdcDEReport(deg = DEGAll_UCEC,gene.type = "long_non_coding")
####DE protein coding genes
dePC_UCEC<-gdcDEReport(deg=DEGAll_UCEC,gene.type = "protein_coding")
##DE pseudogenes
dePseudo_UCEC<-gdcDEReport(deg=DEGAll_UCEC,gene.type = "pseudogene")
demiR_UCEC<-gdcDEReport(DEGMIR_UCEC)
################################################################################################
DEGALL_UCEC_MIR<-gdcDEAnalysis(counts = mirCounts_UCEC,group = metaMatrix.MIR_UCEC$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
gdcVolcanoPlot(DEGALL_UCEC_MIR)
gdcHeatmap(deg.id = rownames(DEGALL_UCEC_MIR),metadata = metaMatrix.MIR_UCEC,rna.expr = mirExpr_UCEC)
demiR_UCEC<-gdcDEReport(deg=DEGALL_UCEC_MIR,gene.type = "miRNAs")
####DEG visualisation
gdcVolcanoPlot(DEGAll_UCEC)
gdcBarPlot(deg = deALL_UCEC,angle = 45,data.type = "RNAseq")
degName=rownames(deALL_UCEC)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_UCEC,rna.expr = rnaExpr_UCEC)
gdcBarPlot(deg=demiR_UCEC,angle=45,data.type="miRNAs")
t1<-table(which(demiR_UCEC$logFC<0))
dim(t1)
t2<-table(which(demiR_UCEC$logFC>0))
dim(t2)
##ceRNAs network analysis

ceOutput_UCEC<- gdcCEAnalysis(lnc= rownames(deLNC_UCEC),pc= rownames(dePC_UCEC),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_UCEC,mir.expr    = mirExpr_UCEC)
ceOutput_UCEC2<- gdcCEAnalysis(lnc= rownames(deLNC_UCEC),pc= rownames(dePC_UCEC),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_UCEC,mir.expr    = mirExpr_UCEC)

####Network Visualization in Cytoscape

ceOutput2_UCEC<-ceOutput_UCEC[ceOutput_UCEC$hyperPValue<0.01 & ceOutput_UCEC$corPValue<0.01 & ceOutput_UCEC$regSim !=0,]
ceOutput2_UCEC2<-ceOutput_UCEC2[ceOutput_UCEC2$hyperPValue<0.01 & ceOutput_UCEC2$corPValue<0.01 & ceOutput_UCEC2$regSim !=0,]
edges_UCEC<-gdcExportNetwork(ceNetwork = ceOutput2_UCEC,net = "edges")
nodes_UCEC<-gdcExportNetwork(ceNetwork = ceOutput2_UCEC,net = "nodes")

write.table(edges_UCEC,file = "edges_UCEC.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_UCEC,file="nodes_UCEC.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(ceOutput2_UCEC,file="ceOutput_UCEC.txt",sep="\t",quote = FALSE,row.names=FALSE)
head(edges_UCEC)

shinyCorPlot(gene1    = rownames(deLNC_UCEC), 
             gene2    = rownames(dePC_UCEC), 
             rna.expr = rnaExpr_UCEC, 
             metadata = metaMatrix.RNA_UCEC)
####Survival Analysis
#####Cox-PH Model
survOutput_UCEC_COXPH<- gdcSurvivalAnalysis(gene     = rownames(deALL_UCEC), 
                                          method   = 'coxph', 
                                          rna.expr = rnaExpr_UCEC, 
                                          metadata = metaMatrix.RNA_UCEC)

#####Kaplan Meir Plot
survOutput_UCEC_KM<- gdcSurvivalAnalysis(gene     = rownames(deALL_UCEC), 
                                       method   = 'KM', 
                                       rna.expr = rnaExpr_UCEC, 
                                       metadata = metaMatrix.RNA_UCEC, 
                                       sep      = 'median')

survOutput_UCEC_KM_miRNA<- gdcSurvivalAnalysis(gene     = rownames(demiR_UCEC), 
                                         method   = 'KM', 
                                         rna.expr = mirExpr_UCEC, 
                                         metadata = metaMatrix.MIR_UCEC, 
                                         sep      = 'median')
gdcKMPlot(gene     = 'ENSG00000137331',
          rna.expr = rnaExpr_UCEC,
          metadata = metaMatrix.RNA_UCEC,
          sep      = 'median')

p_ENSG00000101955_UCEC<-gdcKMPlot(gene="ENSG00000101955",rna.expr = rnaExpr_UCEC,metadata = metaMatrix.RNA_UCEC,sep= 'median')
p_ENSG00000137726_UCEC<-gdcKMPlot(gene = "ENSG00000137726",rna.expr = rnaExpr_UCEC,metadata = metaMatrix.RNA_UCEC,sep = 'median')
p_ENSG00000141404<-gdcKMPlot(gene="ENSG00000141404",rna.expr = rnaExpr_UCEC,metadata = metaMatrix.RNA_UCEC,sep= 'median')
p_ENSG00000134198<-gdcKMPlot(gene = "ENSG00000134198",rna.expr = rnaExpr_UCEC,metadata = metaMatrix.RNA_UCEC,sep = 'median')

#i was able to change legend postion
library(ggplot2)
#par(mfrow=c(2,2))
layout(matrix(c(1,2), 1,1, byrow = TRUE))
p_ENSG00000101955_UCEC<-p_ENSG00000101955_UCEC$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000137726_UCEC<-p_ENSG00000137726_UCEC$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000141404<-p_ENSG00000141404$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000134198<-p_ENSG00000134198$plot+theme(legend.position = c(0.2,0.2))
p_ENSG00000136160+p_ENSG00000144642+p_ENSG00000137726+p_ENSG00000135272+p_ENSG00000131471+p_ENSG00000173641+p_ENSG00000136160_1+p_ENSG00000145936+p_ENSG00000101955+p_ENSG00000137726+p_ENSG00000141404+p_ENSG00000134198
####Functional Enrichment Analysis
###Enrichment analysis
library(GDCRNATools)
enrichOutput_UCEC<- gdcEnrichAnalysis(gene = rownames(deALL_UCEC), simplify = TRUE)
enrichoutput_UCEC<-enrichOutput_UCEC %>% separate(Terms,c("Terms2","Terms"),sep="~")# I removed Term code as it disturbs for plot space
UCEC_GO<-gdcEnrichPlot(enrichoutput_UCEC, type = 'bar', category = 'GO',num.terms = 10)
UCEC_GO<-UCEC_GO+theme_update(text=element_text(size = 10))+theme_bw()+theme(legend.position =c(0.85,0.4))+ylab("-log10(FDR)-UCEC")
UCEC_KEGG<-gdcEnrichPlot(enrichoutput_UCEC, type = 'bubble', category = 'KEGG',num.terms = 10)+ylab("Fold enrichment-UCEC")+theme(legend.position = "none")+theme_update(text=element_text(size = 10))+theme_bw()
UCEC_DO<-gdcEnrichPlot(enrichoutput_UCEC, type = 'bar', category = 'DO',num.terms = 10)
UCEC_DO<-UCEC_DO+ylab("-log10(FDR)-UCEC")+geom_bar(stat = "identity",fill="#1C9C5A")+theme_update(text=element_text(size = 15))+theme_bw()




######Analysis for AMSI Abstract
dim(hyperout_corr_UCEC)
#784  14
dim(MC_HC_ceOutput_UCEC)
#2160   14
intersect(ceOutput2_UCEC_new$Genes,hyperout_corr_UCEC$Gene)
#"ENSG00000111913" "ENSG00000177311" "ENSG00000198853" "ENSG00000159256" "ENSG00000131018" "ENSG00000079308" "ENSG00000082397"
hyperout_corr_UCEC$miRNAs
write.table(hyperout_corr_UCEC,file = "hyperout_corr_UCEC_miR_added.txt",quote = FALSE,row.names = FALSE,sep = "\t")
hyperout_corr_UCEC_miR_added <- read.delim("C:/Users/n10136142/OneDrive - Queensland University of Technology/Desktop/AMSI_Bioinformatics/hyperout_corr_UCEC_miR_added.txt")
miRs_Pseudo_UCEC<-unlist(strsplit(as.character(hyperout_corr_UCEC_miR_added$miRNAs),","))
miRs_lncRNA_UCEC<-unlist(strsplit(as.character(MC_HC_ceOutput_UCEC$miRNAs),","))
intersect(miRs_lncRNA_UCEC,miRs_Pseudo_UCEC)
#"hsa-miR-490-3p"  "hsa-miR-217"     "hsa-miR-455-5p"  "hsa-miR-137"     "hsa-miR-129-5p"  "hsa-miR-142-3p"  "hsa-miR-184"     "hsa-miR-375"    
#"hsa-miR-125a-5p" "hsa-miR-139-5p"  "hsa-miR-126-3p" 
survOutput_UCEC_COXPH_PC<- gdcSurvivalAnalysis(gene     = rownames(dePC_UCEC), 
                                            method   = 'coxph', 
                                            rna.expr = rnaExpr_UCEC, 
                                            metadata = metaMatrix.RNA_UCEC)
survOutput_UCEC_COXPH_MIR<- gdcSurvivalAnalysis(gene     = rownames(demiR_UCEC), 
                                               method   = 'coxph', 
                                               rna.expr = mirExpr_UCEC, 
                                               metadata = metaMatrix.MIR_UCEC)
survOutput_UCEC_COXPH_LNC<- gdcSurvivalAnalysis(gene     = rownames(deLNC_UCEC), 
                                               method   = 'coxph', 
                                               rna.expr = rnaExpr_UCEC, 
                                               metadata = metaMatrix.RNA_UCEC)
survOutput_UCEC_COXPH_Pseudo<- gdcSurvivalAnalysis(gene     = rownames(dePseudo_UCEC), 
                                                method   = 'coxph', 
                                                rna.expr = rnaExpr_UCEC, 
                                                metadata = metaMatrix.RNA_UCEC)

Surv_Sig_PC<-subset(survOutput_UCEC_COXPH_PC,survOutput_UCEC_COXPH_PC$pValue<0.00005)
dim(Surv_Sig_PC)
Surv_Sig_LNC<-subset(survOutput_UCEC_COXPH_LNC,survOutput_UCEC_COXPH_LNC$pValue<0.0005)
dim(Surv_Sig_LNC)
Surv_Sig_miR<-subset(survOutput_UCEC_COXPH_MIR,survOutput_UCEC_COXPH_MIR$pValue<0.005)
dim(Surv_Sig_miR)
Surv_Sig_Pseudo<-subset(survOutput_UCEC_COXPH_Pseudo,survOutput_UCEC_COXPH_Pseudo$pValue<0.005)
dim(Surv_Sig_Pseudo)

ENSG00000048540_UCEC<-gdcKMPlot(gene="ENSG00000048540",rna.expr = rnaExpr_UCEC,metadata = metaMatrix.RNA_UCEC,sep= 'median')
ENSG00000238105_UCEC<-gdcKMPlot(gene="ENSG00000238105",rna.expr = rnaExpr_UCEC,metadata = metaMatrix.RNA_UCEC,sep= 'median')
miR_UCEC<-gdcKMPlot(gene="hsa-miR-190b",rna.expr = mirExpr_UCEC,metadata = metaMatrix.MIR_UCEC,sep= 'median')
