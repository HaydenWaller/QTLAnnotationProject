library(topGO)
library(Rgraphviz)
library(tidyverse)


# LK

GOuniverse<-read_csv("C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lk_all_GOs.recips.csv")
geneID2GO<-readMappings(file = "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lk_all_GOs.recips.csv", sep = ",", IDsep = ";")
geneUniverse <-names(geneID2GO)
QTLers<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lk_QTL_GOs.recips.csv", sep = ",", header = TRUE)
geneList<-factor(as.integer(geneUniverse %in% as.vector(QTLers$transcript)))
names(geneList)<-geneUniverse
myGOdata_BP<-new("topGOdata", description= "lk", ontology = "BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata_MF<-new("topGOdata", description= "lk", ontology = "MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
myGOdata_CC<-new("topGOdata", description= "lk", ontology = "CC", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
resultFisher_BP <- runTest(myGOdata_BP, algorithm="weight01", statistic="fisher")
resultFisher_MF <- runTest(myGOdata_MF, algorithm="weight01", statistic="fisher")
resultFisher_CC <- runTest(myGOdata_CC, algorithm="weight01", statistic="fisher")

allRes_BP <- GenTable(myGOdata_BP, Fisher = resultFisher_BP, topNodes=20)
allRes_BP <- allRes_BP[which(allRes_BP$Significant > 0),]
allRes_BP$fdr_pvalue<-p.adjust(allRes_BP$Fisher,"fdr",nrow(allRes_BP))
allRes_MF <- GenTable(myGOdata_MF, Fisher = resultFisher_MF, topNodes=20)
allRes_MF <- allRes_MF[which(allRes_MF$Significant > 0),]
allRes_MF$fdr_pvalue<-p.adjust(allRes_MF$Fisher,"fdr",nrow(allRes_MF))
allRes_CC <- GenTable(myGOdata_CC, Fisher = resultFisher_CC, topNodes=20)
allRes_CC <- allRes_CC[which(allRes_CC$Significant > 0),]
allRes_CC$fdr_pvalue<-p.adjust(allRes_CC$Fisher,"fdr",nrow(allRes_CC))

write_tsv(allRes_MF, "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lk_MF.tsv")
write_tsv(allRes_BP, "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lk_BP.tsv")
write_tsv(allRes_CC, "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lk_CC.tsv")


showSigOfNodes(myGOdata_BP, score(resultFisher_BP), firstSigNodes = 5, useInfo = 'all')

# LP

GOuniverse<-read_csv("C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lp_all_GOs.recips.csv")
geneID2GO<-readMappings(file = "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lp_all_GOs.recips.csv", sep = ",", IDsep = ";")
geneUniverse <-names(geneID2GO)
QTLers<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lp_QTL_GOs.recips.csv", sep = ",", header = TRUE)
geneList<-factor(as.integer(geneUniverse %in% as.vector(QTLers$transcript)))
names(geneList)<-geneUniverse
myGOdata_BP<-new("topGOdata", description= "lp", ontology = "BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata_MF<-new("topGOdata", description= "lp", ontology = "MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
myGOdata_CC<-new("topGOdata", description= "lp", ontology = "CC", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
resultFisher_BP <- runTest(myGOdata_BP, algorithm="weight01", statistic="fisher")
resultFisher_MF <- runTest(myGOdata_MF, algorithm="weight01", statistic="fisher")
resultFisher_CC <- runTest(myGOdata_CC, algorithm="weight01", statistic="fisher")

allRes_BP <- GenTable(myGOdata_BP, Fisher = resultFisher_BP, topNodes=20)
allRes_BP <- allRes_BP[which(allRes_BP$Significant > 0),]
allRes_BP$fdr_pvalue<-p.adjust(allRes_BP$Fisher,"fdr",nrow(allRes_BP))
allRes_MF <- GenTable(myGOdata_MF, Fisher = resultFisher_MF, topNodes=20)
allRes_MF <- allRes_MF[which(allRes_MF$Significant > 0),]
allRes_MF$fdr_pvalue<-p.adjust(allRes_MF$Fisher,"fdr",nrow(allRes_MF))
allRes_CC <- GenTable(myGOdata_CC, Fisher = resultFisher_CC, topNodes=20)
allRes_CC <- allRes_CC[which(allRes_CC$Significant > 0),]
allRes_CC$fdr_pvalue<-p.adjust(allRes_CC$Fisher,"fdr",nrow(allRes_CC))

write_tsv(allRes_MF, "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lp_MF.tsv")
write_tsv(allRes_BP, "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lp_BP.tsv")
write_tsv(allRes_CC, "C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/annotations/recips/lp_CC.tsv")


pvalFis <- score(resultFisher_BP)
hist(pvalFis, 50, xlab = "p-values")

showSigOfNodes(myGOdata_BP, score(resultFisher_BP), firstSigNodes = 5, useInfo = 'all')
