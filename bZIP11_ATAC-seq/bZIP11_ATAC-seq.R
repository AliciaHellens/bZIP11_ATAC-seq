#if DiffBind not installed 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind")

library(tidyverse)
library(DiffBind)
library(rtracklayer)
library(edgeR)

library(readr)
coldata_bZIP45 <- read_csv("coldata_bZIP45.csv")
View(coldata_bZIP45)

#full data

#sampleTable.df = read.csv("C:/Users/uqahell2/Dropbox/R/ATAC_bZIP_45/coldata_bZIP45.csv")

sampleTable.df = coldata_bZIP45

DiffBindObj <- dba(sampleSheet=sampleTable.df, peakCaller="bed")
plot(DiffBindObj)
DiffBindObj1 <- dba.count(DiffBindObj)
# takes a long time
DiffBindObj2 <- dba.contrast(DiffBindObj1, categories = DBA_CONDITION, minMembers=3)
DiffBindObjAll<- dba.analyze(DiffBindObj2, method=DBA_ALL_METHODS)
DiffBindReportAll <- dba.report(DiffBindObjAll, contrast = 1, file = 'bZIP1145all', initString = "diffBoundPeaks", ext = "csv")

#plots
dba.plotPCA(DiffBindObjAll,attributes = c(DBA_TISSUE, DBA_CONDITION), label = c(DBA_REPLICATE, DBA_TREATMENT) )

pvals <- dba.plotBox(DiffBindObjAll, contrast = 1,method = DBA_EDGER)
dba.plotVolcano(DiffBindObjAll, contrast = 1,)
plot(DiffBindObjAll, contrast = 1,method = DBA_EDGER)
dba.plotPCA(DiffBindObjAll, contrast=1,label=DBA_TISSUE,method = DBA_DESEQ2)
corvals <- dba.plotHeatmap(DiffBindObjAll, contrast=1, correlations=FALSE,method = DBA_DESEQ2)

#bZIP-DEX v WT 

bZIPDEXvWT = read.csv("C:/Users/uqahell2/Dropbox/R/ATAC_bZIP_45/coldata_bZIPDEXvWT.csv")
DiffBindObja <- dba(sampleSheet=bZIPDEXvWT, peakCaller="bed")
plot(DiffBindObja)
DiffBindObjb <- dba.count(DiffBindObja)
DiffBindObjc <- dba.contrast(DiffBindObjb, categories = DBA_TISSUE, minMembers=3) 
DiffBindObjBDWT<- dba.analyze(DiffBindObjc, method=DBA_ALL_METHODS)
DiffBindReportBDWT <- dba.report(DiffBindObjBDWT, contrast = 1, file = 'WTvBZD', initString = "diffBoundPeaks", ext = "csv")

coldiffbindpeaks <- DiffBindReportBDWT
export.bed(coldiffbindpeaks,"bZIP_DEX_v_WT_DAR.bed")

dba.plotPCA(DiffBindObjBDWT,attributes = c(DBA_TISSUE, DBA_CONDITION), label = DBA_REPLICATE)
dba.plotVolcano(DiffBindObjBDWT, contrast = 1,)

#bZIP-DEX v WT - CHX 

CHX_WT_BZD = read.csv("C:/Users/uqahell2/Dropbox/R/ATAC_bZIP_45/coldataCHX_WTvBZD.csv")
DiffBindObjd <- dba(sampleSheet=CHX_WT_BZD, peakCaller="bed")
plot(DiffBindObjd)
DiffBindObje <- dba.count(DiffBindObjd)
DiffBindObjf <- dba.contrast(DiffBindObje, categories = DBA_TISSUE, minMembers=3) 
DiffBindObjCWB<- dba.analyze(DiffBindObjf, method=DBA_ALL_METHODS)
DiffBindReportCWB <- dba.report(DiffBindObjCWB, contrast = 1, file = 'CHX_WTvBZD', initString = "diffBoundPeaks", ext = "csv")

coldiffbindpeaks <- DiffBindReportCWB
export.bed(coldiffbindpeaks,"CHX_bZIP_DEX_v_WT_DAR.bed")

dba.plotPCA(DiffBindObjCWB,attributes = c(DBA_TISSUE, DBA_CONDITION), label = DBA_REPLICATE)
dba.plotVolcano(DiffBindObjCWB, contrast = 1,)

#bZIP-DEX v MOCK 

Mock_BZD = read.csv("C:/Users/uqahell2/Dropbox/R/ATAC_bZIP_45/coldataMOCKvBZD.csv")
DiffBindObjg <- dba(sampleSheet=MOCK_BZD, peakCaller="bed")
plot(DiffBindObjg)
DiffBindObjh <- dba.count(DiffBindObjg)
DiffBindObji <- dba.contrast(DiffBindObjg, categories = DBA_CONDITION, minMembers=3) 
DiffBindObjMOD<- dba.analyze(DiffBindObji, method=DBA_ALL_METHODS)
DiffBindReportMOD <- dba.report(DiffBindObjMOD, contrast = 1, file = 'MOCKvBZD', initString = "diffBoundPeaks", ext = "csv")

coldiffbindpeaks <- DiffBindReportMOD
export.bed(coldiffbindpeaks,"bZIP_DEX_v_MOCK_DAR.bed")

dba.plotPCA(DiffBindObjMOD,attributes = c(DBA_TISSUE, DBA_CONDITION), label = DBA_REPLICATE)
dba.plotVolcano(DiffBindObjMOD, contrast = 1,)


#bZIP-DEX v MOCK - CHX 

CHX_Mock_BZD = read.csv("C:/Users/uqahell2/Dropbox/R/ATAC_bZIP_45/coldataCHX_MOCKvBZD.csv")
DiffBindObjj <- dba(sampleSheet=CHX_MOCK_BZD, peakCaller="bed")
plot(DiffBindObjj)
DiffBindObjk <- dba.count(DiffBindObjj)
DiffBindObjl <- dba.contrast(DiffBindObjk, categories = DBA_CONDITION, minMembers=3) 
DiffBindObjCMOD<- dba.analyze(DiffBindObjl, method=DBA_ALL_METHODS)
DiffBindReportCMOD <- dba.report(DiffBindObjCMOD, contrast = 1, file = 'CHX_MOCKvBZD', initString = "diffBoundPeaks", ext = "csv")

coldiffbindpeaks <- DiffBindReportCMOD
export.bed(coldiffbindpeaks,"CHX_bZIP_DEX_v_MOCK_DAR.bed")

dba.plotPCA(DiffBindObjCMOD,attributes = c(DBA_TISSUE, DBA_CONDITION), label = DBA_REPLICATE)
dba.plotVolcano(DiffBindObjCMOD, contrast = 1,)

