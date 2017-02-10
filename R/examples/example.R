#loading necessary libraries
library(data.table)
library(ggplot2)

#reading in annotation


#reading in region data
#reading in bisulfite data
#will generalize/replace these later
load("C:/cygwin64/Comp_Epigenetics/dump/BSDTListBcellEpithelial.RData")
load("C:/cygwin64/Comp_Epigenetics/dump/encodeRegionSampling.RData")
BSDTList=BSDTList[c(1:3,11:13)]
save(BSDTList,file="C:/cygwin64/Comp_Epigenetics/MIRA_Projects/BSDTListBcellepithelial3_3.RData")
someFeatures=someFeatures[1:6]

#converting to proper format
BSDTList=addMethCol(BSDTList = BSDTList)


#doing MIRA analysis
bigBin=lapply(X = BSDTList,FUN = returnMIRABins,GRList=someFeatures,binNum=11,sampleNameInBSDT=TRUE)

bigBinDT=rbindlist(bigBin)#need to make sure that sample names have a column in here

#adding sampleType from annotation object


#visualizing MIRA signatures
plotMIRARegions(bigBinDT)

#visualizing MIRA scores



