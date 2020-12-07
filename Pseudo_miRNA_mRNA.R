lnc<-rownames(dePseudo_PRAD)
pc<-rownames(dePC_PRAD)
mir<-rownames(demiR_PRAD)
###############################

#######################################
library(dplyr)
lncDa <- unlist(rnaExpr_PRAD[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_PRAD[pc,]) %>% t()
mirDa<-unlist(mirExpr_PRAD[mir,]) %>% t()
lncmir_intsect<-intersect(rownames(lncDa),rownames(mirDa))
lncDa<-subset(lncDa,rownames(lncDa) %in% lncmir_intsect) 
mirDa<-subset(mirDa,rownames(mirDa) %in% lncmir_intsect) 
correlation_lm<-list()
i<-0
for( pseudo in lnc) {
  for(miRs in mir) {
    i=i+1
    
    corlm<-cor.test(lncDa[,pseudo], mirDa[,miRs],alternative = "less")
    reglm<-corlm$estimate
    plm<-corlm$p.value
    correlation_lm[[i]]<-c(pseudo,miRs,reglm,plm)
    
  }}

correlation_lm<-do.call(rbind,correlation_lm)
colnames(correlation_lm)<-c('Pseudogene','miRNA','reglm','CorPval')
correlation_lm<-correlation_lm %>% as.numeric()
correlation_lm<-correlation_lm[correlation_lm[,4]<0.05,]
correlation_lm2<-as.data.frame(correlation_lm) %>% group_by(Pseudogene) %>% summarise(Val=paste(miRNA,collapse=",")) %>% as.data.frame()

###########################################
pcmir_intsect<-intersect(rownames(pcDa),rownames(mirDa))
pcDa<-subset(pcDa,rownames(pcDa) %in% pcmir_intsect) 
mirDa<-subset(mirDa,rownames(mirDa) %in% pcmir_intsect) 
correlation_pm<-list()
i<-0
for( gene in pc) {
  for(miRs in mir) {
    i=i+1
    
    corpm<-cor.test(pcDa[,gene], mirDa[,miRs],alternative = "less")
    regpm<-corpm$estimate
    ppm<-corpm$p.value
    correlation_pm[[i]]<-c(gene,miRs,regpm,ppm)
    
  }}

correlation_pm<-do.call(rbind,correlation_pm)
colnames(correlation_pm)<-c('Gene','miRNA','regpm','CorPval')
correlation_pm<-correlation_pm[correlation_pm[,4]<0.05,]
correlation_pm2<-as.data.frame(correlation_pm) %>% group_by(Gene) %>% summarise(Val=paste(miRNA,collapse=",")) %>% as.data.frame()
##############################################
correlation_pl<-list()
i<-0
for( gene in pc) {
  for(pseudo in lnc) {
    i=i+1
    
    corpc<-cor.test(pcDa[,gene], lncDa[,pseudo],alternative = "greater")
    regpc<-corpc$estimate
    ppc<-corpc$p.value
    correlation_pl[[i]]<-c(gene,pseudo,regpc,ppc)
    
  }}

correlation_pl<-do.call(rbind,correlation_pl)
colnames(correlation_pl)<-c('Gene','Pseudo','regpc','CorPval')
correlation_pl<-as.data.frame(as.matrix(correlation_pl))
correlation_pl<-correlation_pl[correlation_pl[,4]<0.05,]
################Correlation PRAD####################
lnc<-rownames(dePseudo_PRAD)
pc<-rownames(dePC_PRAD)
mir<-rownames(demiR_PRAD)
lncDa <- unlist(rnaExpr_PRAD[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_PRAD[pc,]) %>% t()
mirDa<-unlist(mirExpr_PRAD[mir,]) %>% t()
df_ENSG00000179277<-data.frame(label=paste(rownames(dePC_PRAD)),estimate="",p.value="")
estimates=numeric(1513)
pvalues=numeric(1513)
for (i in 1:1513)
{
    
  test<-cor.test(pcDa[,i],lncDa[,8])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000179277$estimate<-estimates
df_ENSG00000179277$p.value<-pvalues
df_ENSG00000179277

df_ENSG00000214796<-data.frame(label=paste(rownames(dePC_PRAD)),estimate="",p.value="")
estimates=numeric(1513)
pvalues=numeric(1513)
for (i in 1:1513)
{
  
  test<-cor.test(pcDa[,i],lncDa[,12])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000214796$estimate<-estimates
df_ENSG00000214796$p.value<-pvalues
df_ENSG00000214796

df_PseudoPC_PRAD<-data.frame(df_ENSG00000179277,df_ENSG00000214796)
###########################################################Correlation BRCA####################
lnc<-rownames(dePseudo_BRCA)
pc<-rownames(dePC_BRCA)
n<-length(pc)
lncDa <- unlist(rnaExpr_BRCA[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_BRCA[pc,]) %>% t()

df_ENSG00000179277<-data.frame(label=paste(rownames(dePC_BRCA)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,8])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000179277$estimate<-estimates
df_ENSG00000179277$p.value<-pvalues
df_ENSG00000179277

df_ENSG00000214796<-data.frame(label=paste(rownames(dePC_BRCA)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,12])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000214796$estimate<-estimates
df_ENSG00000214796$p.value<-pvalues
df_ENSG00000214796

df_PseudoPC_BRCA<-data.frame(df_ENSG00000179277,df_ENSG00000214796)
##########################Correlation COAD###########################################
lnc<-rownames(dePseudo_COAD)
pc<-rownames(dePC_COAD)
n<-length(pc)
lncDa <- unlist(rnaExpr_COAD[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_COAD[pc,]) %>% t()

df_ENSG00000179277<-data.frame(label=paste(rownames(dePC_COAD)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,8])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000179277$estimate<-estimates
df_ENSG00000179277$p.value<-pvalues
df_ENSG00000179277

df_ENSG00000214796<-data.frame(label=paste(rownames(dePC_COAD)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,12])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000214796$estimate<-estimates
df_ENSG00000214796$p.value<-pvalues
df_ENSG00000214796

df_PseudoPC_COAD<-data.frame(df_ENSG00000179277,df_ENSG00000214796)
##################################correlation READ##################
lnc<-rownames(dePseudo_READ)
pc<-rownames(dePC_READ)
n<-length(pc)
lncDa <- unlist(rnaExpr_READ[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_READ[pc,]) %>% t()

df_ENSG00000179277<-data.frame(label=paste(rownames(dePC_READ)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,8])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000179277$estimate<-estimates
df_ENSG00000179277$p.value<-pvalues
df_ENSG00000179277

df_ENSG00000214796<-data.frame(label=paste(rownames(dePC_READ)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,12])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000214796$estimate<-estimates
df_ENSG00000214796$p.value<-pvalues
df_ENSG00000214796

df_PseudoPC_READ<-data.frame(df_ENSG00000179277,df_ENSG00000214796)
#############################Correlation UCEC#########################
lnc<-rownames(dePseudo_UCEC)
pc<-rownames(dePC_UCEC)
n<-length(pc)
lncDa <- unlist(rnaExpr_UCEC[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_UCEC[pc,]) %>% t()

df_ENSG00000179277<-data.frame(label=paste(rownames(dePC_UCEC)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,8])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000179277$estimate<-estimates
df_ENSG00000179277$p.value<-pvalues
df_ENSG00000179277

df_ENSG00000214796<-data.frame(label=paste(rownames(dePC_UCEC)),estimate="",p.value="")
estimates=numeric(n)
pvalues=numeric(n)
for (i in 1:n)
{
  
  test<-cor.test(pcDa[,i],lncDa[,12])
  estimates[i]=test$estimate
  pvalues[i]=test$p.value
}

df_ENSG00000214796$estimate<-estimates
df_ENSG00000214796$p.value<-pvalues
df_ENSG00000214796

df_PseudoPC_UCEC<-data.frame(df_ENSG00000179277,df_ENSG00000214796)

#################################Hypergeometric Test#############################################################################
#########################Hyper Test for PRAD#########################
dePC_PRAD_ENSGID<-data.frame(rownames(dePC_PRAD))
dePseudo_PRAD_ENSGID<-data.frame(rownames(dePseudo_PRAD))
miRcode_Latest_PCPRAD<-dplyr::left_join(dePC_PRAD_ENSGID,miRcode_Latest,by=c("rownames.dePC_PRAD."="ID"))
miRcode_Latest_PseudoPRAD<-dplyr::left_join(dePseudo_PRAD_ENSGID,miRcode_Latest,by=c("rownames.dePseudo_PRAD."="ID"))
library(dplyr)
miRcode_UniqID_PCPRAD<-miRcode_Latest_PCPRAD %>% group_by(rownames.dePC_PRAD.) %>% summarise(Val=paste(microRNAs,collapse = ","))
miRcode_UniqID_PseudoPRAD<-miRcode_Latest_PseudoPRAD %>% group_by(rownames.dePseudo_PRAD.) %>% dplyr::summarise(Val=paste(microRNAs,collapse=","))
head(miRcode_UniqID_PCPRAD)
MIRCODE_PCPRAD<-strsplit(as.character(miRcode_UniqID_PCPRAD$Val),',')
names(MIRCODE_PCPRAD)<-miRcode_UniqID_PCPRAD$rownames.dePC_PRAD.
View(MIRCODE_PCPRAD)
MIRCODE_PseudoPRAD<-strsplit(as.character(miRcode_UniqID_PseudoPRAD$Val),',')
names(MIRCODE_PseudoPRAD)<-miRcode_UniqID_PseudoPRAD$rownames.dePseudo_PRAD.
View(MIRCODE_PseudoPRAD)
####################################
lnc<-rownames(dePseudo_PRAD)
pc<-rownames(dePC_PRAD)
ceLNC<-lnc[lnc %in% names(MIRCODE_PseudoPRAD)]
cePC<-pc[pc %in% names(MIRCODE_PCPRAD)]
popTotal<-length(union(unique(unlist(MIRCODE_PseudoPRAD)),unique(unlist(MIRCODE_PCPRAD))))
hyperOutput_PRAD <- list()
i <- 0
for (lncID in ceLNC) {
  listTotal <- length(MIRCODE_PseudoPRAD[[lncID]])
  for (gene in cePC) {
    i = i + 1
    ovlp <- intersect(MIRCODE_PseudoPRAD[[lncID]], MIRCODE_PCPRAD[[gene]])
    
    popHits <- length(MIRCODE_PCPRAD[[gene]])
    Counts <- length(ovlp)
    
    ovlpMIRs <- paste(ovlp, collapse = ',')
    foldEnrichment <- Counts/listTotal*popTotal/popHits
    pValue <- phyper(Counts-1, popHits, popTotal-popHits, 
                     listTotal, lower.tail=FALSE, log.p=FALSE)
    
    ceMIR <- Reduce(intersect, list(ovlp, demiR_PRAD))
    deMIRs <- paste(ceMIR, collapse = ',')
    deMIRCounts <- length(ceMIR)
    
    hyperOutput_PRAD[[i]] <- c(lncID, gene, Counts, listTotal,
                          popHits,popTotal,foldEnrichment,pValue,ovlpMIRs,
                          deMIRCounts, deMIRs)
    
  }
}

#hyperOutput_PRAD <- Reduce(rbind, hyperOutput_PRAD)  ## slower
hyperOutput_PRAD <- do.call(rbind, hyperOutput_PRAD)
#hyperOutput_PRAD <- rbind_list(hyperOutput_PRAD) ## not test

colnames(hyperOutput_PRAD) <- c('lncRNAs','Genes','Counts','listTotal',
                           'popHits','popTotal','foldEnrichment','hyperPValue','miRNAs',
                           'deMIRCounts','deMIRs')
hyperOutput_PRAD <- as.data.frame(as.matrix(hyperOutput_PRAD), 
                             stringsAsFactors=FALSE)
hyperOutput_PRAD <- hyperOutput_PRAD[as.numeric(hyperOutput_PRAD$Counts)>0,]

#hyperOutput_PRAD$FDR <- p.adjust(as.numeric(as.character(hyperOutput_PRAD$pValue)),
#method = 'fdr')
#hyperOutput_PRAD <- hyperOutput_PRAD[hyperOutput_PRAD$Counts>0,]
#hyperOutput_PRAD$lncRNAs <- ensembl2symbolFun(hyperOutput_PRAD$lncRNAs)
#hyperOutput_PRAD$gene <- ensembl2symbolFun(hyperOutput_PRAD$gene)

if (is.null(demiR_PRAD)) {
  hyperOutput_PRAD <- hyperOutput_PRAD[,! colnames(hyperOutput_PRAD) %in% 
                               c('deMIRCounts','deMIRs')]
}

return (hyperOutput_PRAD)
}
################################
library(tidyr)
hyperOutput_PRAD$pseudogene_Gene<-paste(hyperOutput_PRAD$lncRNAs,hyperOutput_PRAD$Genes,sep = "_")
names(hyperOutput_PRAD)[10]<-"Pseudo_Gene"
hyperOutput_PRAD<-hyperOutput_PRAD[hyperOutput_PRAD$hyperPValue<0.01,1:10] 
write.table(hyperOutput_PRAD,file = "hyperOutput_PRAD.txt",sep = "\t",quote = FALSE,row.names = FALSE)
##############################HyperTest for BRCA#####################################################
dePC_BRCA_ENSGID<-data.frame(rownames(dePC_BRCA))
dePseudo_BRCA_ENSGID<-data.frame(rownames(dePseudo_BRCA))
miRcode_Latest_PCBRCA<-dplyr::left_join(dePC_BRCA_ENSGID,miRcode_Latest,by=c("rownames.dePC_BRCA."="ID"))
miRcode_Latest_PseudoBRCA<-dplyr::left_join(dePseudo_BRCA_ENSGID,miRcode_Latest,by=c("rownames.dePseudo_BRCA."="ID"))
library(dplyr)
miRcode_UniqID_PCBRCA<-miRcode_Latest_PCBRCA %>% group_by(rownames.dePC_BRCA.) %>% summarise(Val=paste(microRNAs,collapse = ","))
miRcode_UniqID_PseudoBRCA<-miRcode_Latest_PseudoBRCA %>% group_by(rownames.dePseudo_BRCA.) %>% dplyr::summarise(Val=paste(microRNAs,collapse=","))
head(miRcode_UniqID_PCBRCA)
MIRCODE_PCBRCA<-strsplit(as.character(miRcode_UniqID_PCBRCA$Val),',')
names(MIRCODE_PCBRCA)<-miRcode_UniqID_PCBRCA$rownames.dePC_BRCA.
View(MIRCODE_PCBRCA)
MIRCODE_PseudoBRCA<-strsplit(as.character(miRcode_UniqID_PseudoBRCA$Val),',')
names(MIRCODE_PseudoBRCA)<-miRcode_UniqID_PseudoBRCA$rownames.dePseudo_BRCA.
View(MIRCODE_PseudoBRCA)
########
lnc<-rownames(dePseudo_BRCA)
pc<-rownames(dePC_BRCA)
ceLNC<-lnc[lnc %in% names(MIRCODE_PseudoBRCA)]
cePC<-pc[pc %in% names(MIRCODE_PCBRCA)]
popTotal<-length(union(unique(unlist(MIRCODE_PseudoBRCA)),unique(unlist(MIRCODE_PCBRCA))))

hyperOutput_BRCA <- list()
i <- 0
for (lncID in ceLNC) {
  listTotal <- length(MIRCODE_PseudoBRCA[[lncID]])
  for (gene in cePC) {
    i = i + 1
    ovlp <- intersect(MIRCODE_PseudoBRCA[[lncID]], MIRCODE_PCBRCA[[gene]])
    
    popHits <- length(MIRCODE_PCBRCA[[gene]])
    Counts <- length(ovlp)
    
    ovlpMIRs <- paste(ovlp, collapse = ',')
    foldEnrichment <- Counts/listTotal*popTotal/popHits
    pValue <- phyper(Counts-1, popHits, popTotal-popHits, 
                     listTotal, lower.tail=FALSE, log.p=FALSE)
    
    ceMIR <- Reduce(intersect, list(ovlp, demiR_BRCA))
    deMIRs <- paste(ceMIR, collapse = ',')
    deMIRCounts <- length(ceMIR)
    
    hyperOutput_BRCA[[i]] <- c(lncID, gene, Counts, listTotal,
                               popHits,popTotal,foldEnrichment,pValue,ovlpMIRs,
                               deMIRCounts, deMIRs)
    
  }
}

#hyperOutput_BRCA <- Reduce(rbind, hyperOutput_BRCA)  ## slower
hyperOutput_BRCA <- do.call(rbind, hyperOutput_BRCA)
#hyperOutput_BRCA <- rbind_list(hyperOutput_BRCA) ## not test

colnames(hyperOutput_BRCA) <- c('lncRNAs','Genes','Counts','listTotal',
                                'popHits','popTotal','foldEnrichment','hyperPValue','miRNAs',
                                'deMIRCounts','deMIRs')
hyperOutput_BRCA <- as.data.frame(as.matrix(hyperOutput_BRCA), 
                                  stringsAsFactors=FALSE)
hyperOutput_BRCA <- hyperOutput_BRCA[as.numeric(hyperOutput_BRCA$Counts)>0,]

#hyperOutput_BRCA$FDR <- p.adjust(as.numeric(as.character(hyperOutput_BRCA$pValue)),
#method = 'fdr')
#hyperOutput_BRCA <- hyperOutput_BRCA[hyperOutput_BRCA$Counts>0,]
#hyperOutput_BRCA$lncRNAs <- ensembl2symbolFun(hyperOutput_BRCA$lncRNAs)
#hyperOutput_BRCA$gene <- ensembl2symbolFun(hyperOutput_BRCA$gene)

if (is.null(demiR_BRCA)) {
  hyperOutput_BRCA <- hyperOutput_BRCA[,! colnames(hyperOutput_BRCA) %in% 
                                         c('deMIRCounts','deMIRs')]
}

return (hyperOutput_BRCA)
}
hyperOutput_BRCA$Pseudo_Gene<-paste(hyperOutput_BRCA$lncRNAs,hyperOutput_BRCA$Genes,sep = "_")
head(hyperOutput_BRCA)
hyperOutput_BRCA<-hyperOutput_BRCA[hyperOutput_BRCA$hyperPValue<0.01,1:10]
write.table(hyperOutput_BRCA,file = "hyperOutput_BRCA.txt",sep = "\t",quote = FALSE,row.names = FALSE)

####################################Hypertest for COAD#######################################
dePC_COAD_ENSGID<-data.frame(rownames(dePC_COAD))
dePseudo_COAD_ENSGID<-data.frame(rownames(dePseudo_COAD))
miRcode_Latest_PCCOAD<-dplyr::left_join(dePC_COAD_ENSGID,miRcode_Latest,by=c("rownames.dePC_COAD."="ID"))
miRcode_Latest_PseudoCOAD<-dplyr::left_join(dePseudo_COAD_ENSGID,miRcode_Latest,by=c("rownames.dePseudo_COAD."="ID"))
library(dplyr)
miRcode_UniqID_PCCOAD<-miRcode_Latest_PCCOAD %>% group_by(rownames.dePC_COAD.) %>% summarise(Val=paste(microRNAs,collapse = ","))
miRcode_UniqID_PseudoCOAD<-miRcode_Latest_PseudoCOAD %>% group_by(rownames.dePseudo_COAD.) %>% dplyr::summarise(Val=paste(microRNAs,collapse=","))
head(miRcode_UniqID_PCCOAD)
MIRCODE_PCCOAD<-strsplit(as.character(miRcode_UniqID_PCCOAD$Val),',')
names(MIRCODE_PCCOAD)<-miRcode_UniqID_PCCOAD$rownames.dePC_COAD.
View(MIRCODE_PCCOAD)
MIRCODE_PseudoCOAD<-strsplit(as.character(miRcode_UniqID_PseudoCOAD$Val),',')
names(MIRCODE_PseudoCOAD)<-miRcode_UniqID_PseudoCOAD$rownames.dePseudo_COAD.
View(MIRCODE_PseudoCOAD)
######################################
lnc<-rownames(dePseudo_COAD)
pc<-rownames(dePC_COAD)
ceLNC<-lnc[lnc %in% names(MIRCODE_PseudoCOAD)]
cePC<-pc[pc %in% names(MIRCODE_PCCOAD)]
popTotal<-length(union(unique(unlist(MIRCODE_PseudoCOAD)),unique(unlist(MIRCODE_PCCOAD))))
hyperOutput_COAD <- list()
i <- 0
for (lncID in ceLNC) {
  listTotal <- length(MIRCODE_PseudoCOAD[[lncID]])
  for (gene in cePC) {
    i = i + 1
    ovlp <- intersect(MIRCODE_PseudoCOAD[[lncID]], MIRCODE_PCCOAD[[gene]])
    
    popHits <- length(MIRCODE_PCCOAD[[gene]])
    Counts <- length(ovlp)
    
    ovlpMIRs <- paste(ovlp, collapse = ',')
    foldEnrichment <- Counts/listTotal*popTotal/popHits
    pValue <- phyper(Counts-1, popHits, popTotal-popHits, 
                     listTotal, lower.tail=FALSE, log.p=FALSE)
    
    ceMIR <- Reduce(intersect, list(ovlp, demiR_COAD))
    deMIRs <- paste(ceMIR, collapse = ',')
    deMIRCounts <- length(ceMIR)
    
    hyperOutput_COAD[[i]] <- c(lncID, gene, Counts, listTotal,
                               popHits,popTotal,foldEnrichment,pValue,ovlpMIRs,
                               deMIRCounts, deMIRs)
    
  }
}

#hyperOutput_COAD <- Reduce(rbind, hyperOutput_COAD)  ## slower
hyperOutput_COAD <- do.call(rbind, hyperOutput_COAD)
#hyperOutput_COAD <- rbind_list(hyperOutput_COAD) ## not test

colnames(hyperOutput_COAD) <- c('lncRNAs','Genes','Counts','listTotal',
                                'popHits','popTotal','foldEnrichment','hyperPValue','miRNAs',
                                'deMIRCounts','deMIRs')
hyperOutput_COAD <- as.data.frame(as.matrix(hyperOutput_COAD), 
                                  stringsAsFactors=FALSE)
hyperOutput_COAD <- hyperOutput_COAD[as.numeric(hyperOutput_COAD$Counts)>0,]

#hyperOutput_COAD$FDR <- p.adjust(as.numeric(as.character(hyperOutput_COAD$pValue)),
#method = 'fdr')
#hyperOutput_COAD <- hyperOutput_COAD[hyperOutput_COAD$Counts>0,]
#hyperOutput_COAD$lncRNAs <- ensembl2symbolFun(hyperOutput_COAD$lncRNAs)
#hyperOutput_COAD$gene <- ensembl2symbolFun(hyperOutput_COAD$gene)

if (is.null(demiR_COAD)) {
  hyperOutput_COAD <- hyperOutput_COAD[,! colnames(hyperOutput_COAD) %in% 
                                         c('deMIRCounts','deMIRs')]
}

return (hyperOutput_COAD)
hyperOutput_COAD$Pseudo_Gene<-paste(hyperOutput_COAD$lncRNAs,hyperOutput_COAD$Genes,sep = "_")
head(hyperOutput_COAD)
hyperOutput_COAD<-hyperOutput_COAD[hyperOutput_COAD$hyperPValue<0.01,1:10]
write.table(hyperOutput_COAD,file = "hyperOutput_COAD.txt",sep = "\t",quote = FALSE,row.names = FALSE)
########################HyperTest for READ###########################
dePC_READ_ENSGID<-data.frame(rownames(dePC_READ))
dePseudo_READ_ENSGID<-data.frame(rownames(dePseudo_READ))
miRcode_Latest_PCREAD<-dplyr::left_join(dePC_READ_ENSGID,miRcode_Latest,by=c("rownames.dePC_READ."="ID"))
miRcode_Latest_PseudoREAD<-dplyr::left_join(dePseudo_READ_ENSGID,miRcode_Latest,by=c("rownames.dePseudo_READ."="ID"))
library(dplyr)
miRcode_UniqID_PCREAD<-miRcode_Latest_PCREAD %>% group_by(rownames.dePC_READ.) %>% summarise(Val=paste(microRNAs,collapse = ","))
miRcode_UniqID_PseudoREAD<-miRcode_Latest_PseudoREAD %>% group_by(rownames.dePseudo_READ.) %>% dplyr::summarise(Val=paste(microRNAs,collapse=","))
head(miRcode_UniqID_PCREAD)
MIRCODE_PCREAD<-strsplit(as.character(miRcode_UniqID_PCREAD$Val),',')
names(MIRCODE_PCREAD)<-miRcode_UniqID_PCREAD$rownames.dePC_READ.
View(MIRCODE_PCREAD)
MIRCODE_PseudoREAD<-strsplit(as.character(miRcode_UniqID_PseudoREAD$Val),',')
names(MIRCODE_PseudoREAD)<-miRcode_UniqID_PseudoREAD$rownames.dePseudo_READ.
View(MIRCODE_PseudoREAD)
#########
lnc<-rownames(dePseudo_READ)
pc<-rownames(dePC_READ)
ceLNC<-lnc[lnc %in% names(MIRCODE_PseudoREAD)]
cePC<-pc[pc %in% names(MIRCODE_PCREAD)]
popTotal<-length(union(unique(unlist(MIRCODE_PseudoREAD)),unique(unlist(MIRCODE_PCREAD))))
hyperOutput_READ <- list()
i <- 0
for (lncID in ceLNC) {
  listTotal <- length(MIRCODE_PseudoREAD[[lncID]])
  for (gene in cePC) {
    i = i + 1
    ovlp <- intersect(MIRCODE_PseudoREAD[[lncID]], MIRCODE_PCREAD[[gene]])
    
    popHits <- length(MIRCODE_PCREAD[[gene]])
    Counts <- length(ovlp)
    
    ovlpMIRs <- paste(ovlp, collapse = ',')
    foldEnrichment <- Counts/listTotal*popTotal/popHits
    pValue <- phyper(Counts-1, popHits, popTotal-popHits, 
                     listTotal, lower.tail=FALSE, log.p=FALSE)
    
    ceMIR <- Reduce(intersect, list(ovlp, demiR_READ))
    deMIRs <- paste(ceMIR, collapse = ',')
    deMIRCounts <- length(ceMIR)
    
    hyperOutput_READ[[i]] <- c(lncID, gene, Counts, listTotal,
                               popHits,popTotal,foldEnrichment,pValue,ovlpMIRs,
                               deMIRCounts, deMIRs)
    
  }
}

#hyperOutput_READ <- Reduce(rbind, hyperOutput_READ)  ## slower
hyperOutput_READ <- do.call(rbind, hyperOutput_READ)
#hyperOutput_READ <- rbind_list(hyperOutput_READ) ## not test

colnames(hyperOutput_READ) <- c('lncRNAs','Genes','Counts','listTotal',
                                'popHits','popTotal','foldEnrichment','hyperPValue','miRNAs',
                                'deMIRCounts','deMIRs')
hyperOutput_READ <- as.data.frame(as.matrix(hyperOutput_READ), 
                                  stringsAsFactors=FALSE)
hyperOutput_READ <- hyperOutput_READ[as.numeric(hyperOutput_READ$Counts)>0,]

#hyperOutput_READ$FDR <- p.adjust(as.numeric(as.character(hyperOutput_READ$pValue)),
#method = 'fdr')
#hyperOutput_READ <- hyperOutput_READ[hyperOutput_READ$Counts>0,]
#hyperOutput_READ$lncRNAs <- ensembl2symbolFun(hyperOutput_READ$lncRNAs)
#hyperOutput_READ$gene <- ensembl2symbolFun(hyperOutput_READ$gene)

if (is.null(demiR_READ)) {
  hyperOutput_READ <- hyperOutput_READ[,! colnames(hyperOutput_READ) %in% 
                                         c('deMIRCounts','deMIRs')]
}

return (hyperOutput_READ)
}
hyperOutput_READ$Pseudo_Gene<-paste(hyperOutput_READ$lncRNAs,hyperOutput_READ$Genes,sep = "_")
head(hyperOutput_READ)
  hyperOutput_READ<-hyperOutput_READ[hyperOutput_READ$hyperPValue<0.01,1:10]
  write.table(hyperOutput_READ,file = "hyperOutput_READ.txt",sep = "\t",quote = FALSE,row.names = FALSE)  

######################################################HyperTest UCEC##################################
  dePC_UCEC_ENSGID<-data.frame(rownames(dePC_UCEC))
  dePseudo_UCEC_ENSGID<-data.frame(rownames(dePseudo_UCEC))
  miRcode_Latest_PCUCEC<-dplyr::left_join(dePC_UCEC_ENSGID,miRcode_Latest,by=c("rownames.dePC_UCEC."="ID"))
  miRcode_Latest_PseudoUCEC<-dplyr::left_join(dePseudo_UCEC_ENSGID,miRcode_Latest,by=c("rownames.dePseudo_UCEC."="ID"))
  library(dplyr)
  miRcode_UniqID_PCUCEC<-miRcode_Latest_PCUCEC %>% group_by(rownames.dePC_UCEC.) %>% summarise(Val=paste(microRNAs,collapse = ","))
  miRcode_UniqID_PseudoUCEC<-miRcode_Latest_PseudoUCEC %>% group_by(rownames.dePseudo_UCEC.) %>% dplyr::summarise(Val=paste(microRNAs,collapse=","))
  head(miRcode_UniqID_PCUCEC)
  MIRCODE_PCUCEC<-strsplit(as.character(miRcode_UniqID_PCUCEC$Val),',')
  names(MIRCODE_PCUCEC)<-miRcode_UniqID_PCUCEC$rownames.dePC_UCEC.
  View(MIRCODE_PCUCEC)
  MIRCODE_PseudoUCEC<-strsplit(as.character(miRcode_UniqID_PseudoUCEC$Val),',')
  names(MIRCODE_PseudoUCEC)<-miRcode_UniqID_PseudoUCEC$rownames.dePseudo_UCEC.
  View(MIRCODE_PseudoUCEC)
  ########
  lnc<-rownames(dePseudo_UCEC)
  pc<-rownames(dePC_UCEC)
  ceLNC<-lnc[lnc %in% names(MIRCODE_PseudoUCEC)]
  cePC<-pc[pc %in% names(MIRCODE_PCUCEC)]
  popTotal<-length(union(unique(unlist(MIRCODE_PseudoUCEC)),unique(unlist(MIRCODE_PCUCEC))))
  hyperOutput_UCEC <- list()
  i <- 0
  for (lncID in ceLNC) {
    listTotal <- length(MIRCODE_PseudoUCEC[[lncID]])
    for (gene in cePC) {
      i = i + 1
      ovlp <- intersect(MIRCODE_PseudoUCEC[[lncID]], MIRCODE_PCUCEC[[gene]])
      
      popHits <- length(MIRCODE_PCUCEC[[gene]])
      Counts <- length(ovlp)
      
      ovlpMIRs <- paste(ovlp, collapse = ',')
      foldEnrichment <- Counts/listTotal*popTotal/popHits
      pValue <- phyper(Counts-1, popHits, popTotal-popHits, 
                       listTotal, lower.tail=FALSE, log.p=FALSE)
      
      ceMIR <- Reduce(intersect, list(ovlp, demiR_UCEC))
      deMIRs <- paste(ceMIR, collapse = ',')
      deMIRCounts <- length(ceMIR)
      
      hyperOutput_UCEC[[i]] <- c(lncID, gene, Counts, listTotal,
                                 popHits,popTotal,foldEnrichment,pValue,ovlpMIRs,
                                 deMIRCounts, deMIRs)
      
    }
  }
  
  #hyperOutput_UCEC <- Reduce(rbind, hyperOutput_UCEC)  ## slower
  hyperOutput_UCEC <- do.call(rbind, hyperOutput_UCEC)
  #hyperOutput_UCEC <- rbind_list(hyperOutput_UCEC) ## not test
  
  colnames(hyperOutput_UCEC) <- c('lncRNAs','Genes','Counts','listTotal',
                                  'popHits','popTotal','foldEnrichment','hyperPValue','miRNAs',
                                  'deMIRCounts','deMIRs')
  hyperOutput_UCEC <- as.data.frame(as.matrix(hyperOutput_UCEC), 
                                    stringsAsFactors=FALSE)
  hyperOutput_UCEC <- hyperOutput_UCEC[as.numeric(hyperOutput_UCEC$Counts)>0,]
  
  #hyperOutput_UCEC$FDR <- p.adjust(as.numeric(as.character(hyperOutput_UCEC$pValue)),
  #method = 'fdr')
  #hyperOutput_UCEC <- hyperOutput_UCEC[hyperOutput_UCEC$Counts>0,]
  #hyperOutput_UCEC$lncRNAs <- ensembl2symbolFun(hyperOutput_UCEC$lncRNAs)
  #hyperOutput_UCEC$gene <- ensembl2symbolFun(hyperOutput_UCEC$gene)
  
  if (is.null(demiR_UCEC)) {
    hyperOutput_UCEC <- hyperOutput_UCEC[,! colnames(hyperOutput_UCEC) %in% 
                                           c('deMIRCounts','deMIRs')]
  }
  
  return (hyperOutput_UCEC)
  }
hyperOutput_UCEC$Pseudo_Gene<-paste(hyperOutput_UCEC$lncRNAs,hyperOutput_UCEC$Genes,sep = "_")
head(hyperOutput_UCEC)
  hyperOutput_UCEC<-hyperOutput_UCEC[hyperOutput_UCEC$hyperPValue<0.01,1:10]
  write.table(hyperOutput_UCEC,file = "hyperOutput_UCEC.txt",sep = "\t",quote = FALSE,row.names = FALSE)  
  
  ############################DONE!!!!!!#############################
  
hyperOutput_BRCA$Cancer<-rep("BRCA",length(hyperOutput_BRCA$lncRNAs))
hyperOutput_COAD$Cancer<-rep("COAD",length(hyperOutput_COAD$lncRNAs))
hyperOutput_READ$Cancer<-rep("READ",length(hyperOutput_READ$lncRNAs))
hyperOutput_PRAD$Cancer<-rep("PRAD",length(hyperOutput_PRAD$lncRNAs))
hyperOutput_UCEC$Cancer<-rep("UCEC",length(hyperOutput_UCEC$lncRNAs))

hyperOutput_ALL<-rbind(hyperOutput_BRCA,hyperOutput_COAD,hyperOutput_PRAD,hyperOutput_READ,hyperOutput_UCEC)
##############################PRAD-pseudogene-gene correlation########################################
lnc<-rownames(dePseudo_PRAD)
pc<-rownames(dePC_PRAD)
mir<-rownames(demiR_PRAD)
library(dplyr)
lncDa <- unlist(rnaExpr_PRAD[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_PRAD[pc,]) %>% t()
mirDa<-unlist(mirExpr_PRAD[mir,]) %>% t()
lncpc_intsect<-intersect(rownames(lncDa),rownames(pcDa))
lncDa<-subset(lncDa,rownames(lncDa) %in% lncpc_intsect) 
pcDa<-subset(pcDa,rownames(pcDa) %in% lncpc_intsect) 

correlation_pl<-list()
i<-0
for( gene in pc) {
  for(pseudo in lnc) {
    i=i+1
    
    corpc<-cor.test(pcDa[,gene], lncDa[,pseudo],alternative = "greater")
    regpc<-corpc$estimate
    ppc<-corpc$p.value
    correlation_pl[[i]]<-c(gene,pseudo,regpc,ppc)
    
  }}

correlation_pl_PRAD<-do.call(rbind,correlation_pl)
colnames(correlation_pl_PRAD)<-c('Gene','Pseudo','regpc','CorPval')
correlation_pl_PRAD<-as.data.frame(correlation_pl_PRAD)
correlation_pl_PRAD$Pseudo_Gene<-paste(correlation_pl_PRAD$Pseudo,correlation_pl_PRAD$Gene,sep = "_")
write.table(correlation_pl_PRAD,file = "correlation_pl_PRAD.txt",quote = FALSE,sep = "\t")
########################BRCA##########################################
lnc<-rownames(dePseudo_BRCA)
pc<-rownames(dePC_BRCA)
mir<-rownames(demiR_BRCA)
library(dplyr)
lncDa <- unlist(rnaExpr_BRCA[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_BRCA[pc,]) %>% t()
mirDa<-unlist(mirExpr_BRCA[mir,]) %>% t()
lncpc_intsect<-intersect(rownames(lncDa),rownames(pcDa))
lncDa<-subset(lncDa,rownames(lncDa) %in% lncpc_intsect) 
pcDa<-subset(pcDa,rownames(pcDa) %in% lncpc_intsect) 

correlation_pl<-list()
i<-0
for( gene in pc) {
  for(pseudo in lnc) {
    i=i+1
    
    corpc<-cor.test(pcDa[,gene], lncDa[,pseudo],alternative = "greater")
    regpc<-corpc$estimate
    ppc<-corpc$p.value
    correlation_pl[[i]]<-c(gene,pseudo,regpc,ppc)
    
  }}

correlation_pl_BRCA<-do.call(rbind,correlation_pl)
colnames(correlation_pl_BRCA)<-c('Gene','Pseudo','regpc','CorPval')
correlation_pl_BRCA<-as.data.frame(correlation_pl_BRCA)
correlation_pl_BRCA$Pseudo_Gene<-paste(correlation_pl_BRCA$Pseudo,correlation_pl_BRCA$Gene,sep = "_")
write.table(correlation_pl_BRCA,file = "correlation_pl_BRCA.txt",quote = FALSE)
###########################COAD####################################
lnc<-rownames(dePseudo_COAD)
pc<-rownames(dePC_COAD)
mir<-rownames(demiR_COAD)
library(dplyr)
lncDa <- unlist(rnaExpr_COAD[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_COAD[pc,]) %>% t()
mirDa<-unlist(mirExpr_COAD[mir,]) %>% t()
lncpc_intsect<-intersect(rownames(lncDa),rownames(pcDa))
lncDa<-subset(lncDa,rownames(lncDa) %in% lncpc_intsect) 
pcDa<-subset(pcDa,rownames(pcDa) %in% lncpc_intsect) 

correlation_pl<-list()
i<-0
for( gene in pc) {
  for(pseudo in lnc) {
    i=i+1
    
    corpc<-cor.test(pcDa[,gene], lncDa[,pseudo],alternative = "greater")
    regpc<-corpc$estimate
    ppc<-corpc$p.value
    correlation_pl[[i]]<-c(gene,pseudo,regpc,ppc)
    
  }}

correlation_pl_COAD<-do.call(rbind,correlation_pl)
colnames(correlation_pl_COAD)<-c('Gene','Pseudo','regpc','CorPval')
correlation_pl_COAD<-as.data.frame(correlation_pl_COAD)
correlation_pl_COAD$Pseudo_Gene<-paste(correlation_pl_COAD$Pseudo,correlation_pl_COAD$Gene,sep = "_")
write.table(correlation_pl_COAD,file = "correlation_pl_COAD.txt",quote = FALSE)
################################READ#####################################
lnc<-rownames(dePseudo_READ)
pc<-rownames(dePC_READ)
mir<-rownames(demiR_READ)
library(dplyr)
lncDa <- unlist(rnaExpr_READ[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_READ[pc,]) %>% t()
mirDa<-unlist(mirExpr_READ[mir,]) %>% t()
lncpc_intsect<-intersect(rownames(lncDa),rownames(pcDa))
lncDa<-subset(lncDa,rownames(lncDa) %in% lncpc_intsect) 
pcDa<-subset(pcDa,rownames(pcDa) %in% lncpc_intsect) 

correlation_pl<-list()
i<-0
for( gene in pc) {
  for(pseudo in lnc) {
    i=i+1
    
    corpc<-cor.test(pcDa[,gene], lncDa[,pseudo],alternative = "greater")
    regpc<-corpc$estimate
    ppc<-corpc$p.value
    correlation_pl[[i]]<-c(gene,pseudo,regpc,ppc)
    
  }}

correlation_pl_READ<-do.call(rbind,correlation_pl)
colnames(correlation_pl_READ)<-c('Gene','Pseudo','regpc','CorPval')
correlation_pl_READ<-as.data.frame(correlation_pl_READ)
correlation_pl_READ$Pseudo_Gene<-paste(correlation_pl_READ$Pseudo,correlation_pl_READ$Gene,sep = "_")
write.table(correlation_pl_READ,file = "correlation_pl_READ.txt",quote = FALSE)
##########################UCEC################################
lnc<-rownames(dePseudo_UCEC)
pc<-rownames(dePC_UCEC)
mir<-rownames(demiR_UCEC)
library(dplyr)
lncDa <- unlist(rnaExpr_UCEC[lnc,]) %>% t()
pcDa <- unlist(rnaExpr_UCEC[pc,]) %>% t()
mirDa<-unlist(mirExpr_UCEC[mir,]) %>% t()
lncpc_intsect<-intersect(rownames(lncDa),rownames(pcDa))
lncDa<-subset(lncDa,rownames(lncDa) %in% lncpc_intsect) 
pcDa<-subset(pcDa,rownames(pcDa) %in% lncpc_intsect) 

correlation_pl<-list()
i<-0
for( gene in pc) {
  for(pseudo in lnc) {
    i=i+1
    
    corpc<-cor.test(pcDa[,gene], lncDa[,pseudo],alternative = "greater")
    regpc<-corpc$estimate
    ppc<-corpc$p.value
    correlation_pl[[i]]<-c(gene,pseudo,regpc,ppc)
    
  }}

correlation_pl_UCEC<-do.call(rbind,correlation_pl)
colnames(correlation_pl_UCEC)<-c('Gene','Pseudo','regpc','CorPval')
correlation_pl_UCEC<-as.data.frame(correlation_pl_UCEC)
correlation_pl_UCEC$Pseudo_Gene<-paste(correlation_pl_UCEC$Pseudo,correlation_pl_UCEC$Gene,sep = "_")
write.table(correlation_pl_UCEC,file = "correlation_pl_UCEC.txt",quote = FALSE)
######################Merging Hypergeometeric test and correlation test results#################
hyper_cor_PRAD<-left_join(hyperOutput_PRAD,correlation_pl_PRAD,by="Pseudo_Gene")
hyper_cor_BRCA<-left_join(hyperOutput_BRCA,correlation_pl_BRCA,by="Pseudo_Gene")
hyper_cor_COAD<-left_join(hyperOutput_COAD,correlation_pl_COAD,by="Pseudo_Gene")
hyper_cor_READ<-left_join(hyperOutput_READ,correlation_pl_READ,by="Pseudo_Gene")
hyper_cor_UCEC<-left_join(hyperOutput_UCEC,correlation_pl_UCEC,by="Pseudo_Gene")
hyper_cor_ALL<-rbind(hyper_cor_BRCA,hyper_cor_PRAD,hyper_cor_COAD,hyper_cor_READ,hyper_cor_UCEC)
t1<-as.data.frame(table(hyper_cor_ALL$Pseudo_Gene))
t2<-t1[t1$Freq>4,]
hyper_cor_ALL_common<-subset(hyper_cor_ALL,hyper_cor_ALL$Pseudo_Gene %in% t2$Var1)
write.table(hyper_cor_ALL_common,file = "hyper_cor_ALL_common.txt",sep = "\t",quote = FALSE)
#########################################

###I merged all correlation results and hypergeometric results
hyperout_corr_UCEC<-dplyr::left_join(hyperOutput_UCEC,correlation_pl_UCEC,by="Pseudo_Gene")
hyperout_corr_PRAD<-dplyr::left_join(hyperOutput_PRAD,correlation_pl_PRAD,by=c("lnc_Gene"="Pseudo_Gene"))
hyperout_corr_COAD<-dplyr::left_join(hyperOutput_COAD,correlation_pl_COAD,by="Pseudo_Gene")
hyperout_corr_READ<-dplyr::left_join(hyperOutput_READ,correlation_pl_READ,by="Pseudo_Gene")
hyperout_corr_BRCA<-dplyr::left_join(hyperOutput_BRCA,correlation_pl_BRCA,by="Pseudo_Gene")


###sent them to my desktop :)
hyperout_corr_UCEC<-subset(hyperout_corr_UCEC,hyperout_corr_UCEC$regpc>0.4)
hyperout_corr_PRAD<-subset(hyperout_corr_PRAD,hyperout_corr_PRAD$regpc>0.4)
hyperout_corr_COAD<-subset(hyperout_corr_COAD,hyperout_corr_COAD$regpc>0.4)
hyperout_corr_READ<-subset(hyperout_corr_READ,hyperout_corr_READ$regpc>0.4)
hyperout_corr_BRCA<-subset(hyperout_corr_BRCA,hyperout_corr_BRCA$regpc>0.4)
write.table(hyperout_corr_UCEC,file = "hyperout_corr_UCEC.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(hyperout_corr_PRAD,file = "hyperout_corr_PRAD.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(hyperout_corr_COAD,file = "hyperout_corr_COAD.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(hyperout_corr_READ,file = "hyperout_corr_READ.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(hyperout_corr_BRCA,file = "hyperout_corr_BRCA.txt",sep = "\t",quote = FALSE,row.names = FALSE)

##Counting unique miRNA,lncRNA, pseudoegens, and mRNAs in ceRNAs

miRNAs_lncRNA_PRAD<-strsplit(MC_HC_ceOutput_PRAD$miRNAs,",")
length(unique(unlist(miRNAs_lncRNA_PRAD)))
length(unique(MC_HC_ceOutput_PRAD$lncRNAs))
length(unique(MC_HC_ceOutput_PRAD$Genes))

miRNAs_lncRNA_BRCA<-strsplit(MC_HC_ceOutput_BRCA$miRNAs,",")
length(unique(unlist(miRNAs_lncRNA_BRCA)))
length(unique(MC_HC_ceOutput_BRCA$lncRNAs))
length(unique(MC_HC_ceOutput_BRCA$Genes))

miRNAs_lncRNA_COAD<-strsplit(MC_HC_ceOutput_COAD$miRNAs,",")
length(unique(unlist(miRNAs_lncRNA_COAD)))
length(unique(MC_HC_ceOutput_COAD$lncRNAs))
length(unique(MC_HC_ceOutput_COAD$Genes))

miRNAs_lncRNA_READ<-strsplit(MC_HC_ceOutput_READ$miRNAs,",")
length(unique(unlist(miRNAs_lncRNA_READ)))
length(unique(MC_HC_ceOutput_READ$lncRNAs))
length(unique(MC_HC_ceOutput_READ$Genes))

miRNAs_lncRNA_UCEC<-strsplit(MC_HC_ceOutput_UCEC$miRNAs,",")
length(unique(unlist(miRNAs_lncRNA_UCEC)))
length(unique(MC_HC_ceOutput_UCEC$lncRNAs))
length(unique(MC_HC_ceOutput_UCEC$Genes))

miRNAs_Pseudogene_UCEC<-strsplit(hyperout_corr_UCEC$miRNAs,",")
length(unique(unlist(miRNAs_Pseudogene_UCEC)))
length(unique(hyperout_corr_UCEC$lncRNAs))
length(unique(hyperout_corr_UCEC$Genes))

