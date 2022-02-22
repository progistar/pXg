
install.packages("readxl")
library("readxl")
library(ggplot2)
library(pheatmap)

redC <- "#E41A1C" ## FOR Novel things
blueC <- "#377EB8" ## FOR Protein coding things
greenC <- "#4DAF4A" ## FOR ELSE
purpleC <- "#984EA3" ##
orangeC <- "#FF7F00"
yellowC <- "#FFFF33"
brownC <- "#A65628"
pinkC <- "#F781BF"
grayC <- "#999999"

setwd("C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results/2.withCalibration/")

## Function List
getBindingRatio <- function (data, fdr = NULL) {
  cPSMs <- data[data$IsCanonical == T, ]
  ncPSMs <- data[data$IsCanonical == F, ]
  
  if(!is.null(fdr)) {
    ncPSMs <- ncPSMs[ncPSMs$`Denovo score` >= fdr, ]
  }
  
  cMAPSMs <- cPSMs[cPSMs$BestScore < 2, ]
  ncMAPSMs <- ncPSMs[ncPSMs$BestScore < 2, ]
  
  
  cPSM_Ratio <- nrow(cMAPSMs) / nrow(cPSMs) * 100
  ncPSM_Ratio <- nrow(ncMAPSMs) / nrow(ncPSMs) * 100
  
  cPSMs <- unique(cPSMs$InferredPeptide)
  ncPSMs <- unique(ncPSMs$InferredPeptide)
  cMAPSMs <- unique(cMAPSMs$InferredPeptide)
  ncMAPSMs <- unique(ncMAPSMs$InferredPeptide)
  
  cMAP_Ratio <- length(cMAPSMs) / length(cPSMs) * 100
  ncMAP_Ratio <- length(ncMAPSMs) / length(ncPSMs) * 100
  
  res <- list()
  
  res$cPSM_Ratio <- cPSM_Ratio
  res$ncPSM_Ratio <- ncPSM_Ratio
  res$cMAP_Ratio <- cMAP_Ratio
  res$ncMAP_Ratio <- ncMAP_Ratio
  
  res
}

## Read Laumont Results
subject1Data <- read_excel(path = "Laumont_Results.xlsx", sheet = "subject1")
subject2Data <- read_excel(path = "Laumont_Results.xlsx", sheet = "subject2")
subject3Data <- read_excel(path = "Laumont_Results.xlsx", sheet = "subject3")
subject4Data <- read_excel(path = "Laumont_Results.xlsx", sheet = "subject4")

subjectNames <- c("subject1", "subject2", "subject3", "subject4")
subjectFDRScore <- c(80, 82, 77, 84)
#subjectFDRScore <- c(0, 0, 0, 0)
allData <- list(subject1Data, subject2Data, subject3Data, subject4Data)

## IDed PSMs
idx <- 1
for(data in allData){
  cPSMs <- data[data$IsCanonical == T, ]
  ncPSMs <- data[data$IsCanonical == F, ]
  ncPSMs <- ncPSMs[ncPSMs$`Denovo score` >= subjectFDRScore[idx], ]
  
  print("PSMs")
  print(nrow(cPSMs))
  print(nrow(ncPSMs))
  print("Peptides")
  print(length(unique(cPSMs$InferredPeptide)))
  print(length(unique(ncPSMs$InferredPeptide)))
  print("========")
  
  idx <- idx+1
}
idx <- 1
## print ratio
for(data in allData){
  #ratio <- getBindingRatio(data)
  ratio <- getBindingRatio(data, subjectFDRScore[idx])
  cstr <- paste("cPSM_Ratio:",ratio$cPSM_Ratio,sep =" ")
  ncstr <- paste("ncPSM_Ratio:",ratio$ncPSM_Ratio,sep =" ")
  str <- paste(cstr, ncstr, sep=" ")
  print(str)
  
  cstr <- paste("cMAP_Ratio:",ratio$cMAP_Ratio,sep =" ")
  ncstr <- paste("ncMAP_Ratio:",ratio$ncMAP_Ratio,sep =" ")
  str <- paste(cstr, ncstr, sep=" ")
  print(str)
  
  idx <- idx+1
}

## print cMAPs
counts <- data.frame(row.names = subjectNames, stringsAsFactors = as.double())
idx <- 1
for(data in allData){
  cPSMs <- data[data$IsCanonical == T, ]
  #cPSMs <- cPSMs[cPSMs$BestScore < 2, ]
  
  MAPGenes <- unique(cPSMs$GeneNames)
  
  for(MAPGene in MAPGenes) {
    count <- nrow(cPSMs[cPSMs$GeneNames == MAPGene, ])
    if(count > 0) {
      counts[subjectNames[idx], MAPGene] <- log2(count+1)
    }
  }
  idx <- idx+1
}
counts[is.na(counts)] <- 0
pheatmap(t(as.matrix(counts)), color = colorRampPalette(c("white", blueC))(50))

## print ncMAPs
counts <- data.frame(row.names = subjectNames, stringsAsFactors = as.double())
idx <- 1
for(data in allData){
  ncPSMs <- data[data$IsCanonical == F, ]
  #ncPSMs <- ncPSMs[ncPSMs$BestScore < 2, ]
  ncPSMs <- ncPSMs[ncPSMs$`Denovo score` >= subjectFDRScore[idx], ]
  ## unique
  #ncPSMs <- ncPSMs[ncPSMs$GenomicLociCount == 1, ]
  #ncPSMs <- ncPSMs[ncPSMs$EventCount == 1, ]
  #ncPSMs <- ncPSMs[ncPSMs$GeneIDCount == 1, ]
  geneCnt <- 0
  
  # for gene level
  #ncMAPGenes <- unique(ncPSMs$GeneNames)
  #for(ncMAPGene in ncMAPGenes) {
  #  count <- nrow(ncPSMs[ncPSMs$GeneNames == ncMAPGene, ])
  #  if(count > 0) {
  #    counts[subjectNames[idx], ncMAPGene] <- log2(count+1)
  #  }
  #}
  
  # for peptide-gene pair level
  #ncMAPGenes <- unique(paste(ncPSMs$InferredPeptide,ncPSMs$Events, sep="@"))
  #for(ncMAPGene in ncMAPGenes) {
  #  count <- nrow(ncPSMs[paste(ncPSMs$InferredPeptide,ncPSMs$Events, sep="@") == ncMAPGene, ])
  #  if(count > 0) {
  #    counts[subjectNames[idx], ncMAPGene] <- log2(count+1)
  #    geneCnt <- geneCnt+1
  #  }
  #}
  
  ncMAPGenes <- unique(ncPSMs$InferredPeptide)
  for(ncMAPGene in ncMAPGenes) {
    count <- nrow(ncPSMs[ncPSMs$InferredPeptide == ncMAPGene, ])
    if(count > 0) {
      counts[subjectNames[idx], ncMAPGene] <- log2(count+1)
      geneCnt <- geneCnt+1
    }
  }
  
  idx <- idx+1
  
  print(geneCnt)
}
counts[is.na(counts)] <- 0
pheatmap(t(as.matrix(counts)), 
         color = colorRampPalette(c("white", blueC))(50),
         clustering_method = "average", kmeans_k = 6)
?pheatmap
         