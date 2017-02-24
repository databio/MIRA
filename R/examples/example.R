
#initial parameter
binNum=11

#loading sample annotation
load(system.file("data","annoDF.RData",package="MIRA"))
annoDF=as.data.table(annoDF)

#reading in region data and bisulfite data
someFeatures=get(load(system.file("data","Gm12878Nrf1_Subset.RData",package="MIRA")))
BSDT=get(load(system.file("data","GM06990_1_ExampleSet.RData",package = "MIRA")))

#converting to lists 
BSDTList=list(BSDT)
someFeatures=GRangesList(someFeatures)
names(someFeatures) <- "GM12878Nrf1"

#converting to proper format
BSDTList=addMethCol(BSDTList = BSDTList)


#doing MIRA analysis
bigBin=lapply(X = BSDTList,FUN = returnMIRABins,GRList=someFeatures,binNum=11,sampleNameInBSDT=TRUE)

bigBinDT=rbindlist(bigBin)#need to make sure that sample names have a column in here

#adding sampleType from annotation object
setkey(bigBinDT,sampleName)
setkey(annoDF,sampleName)
bigBinDT=merge(bigBinDT,annoDF, all.x=TRUE)

#scoring samples
sampleScores=bigBinDT[,.(score = scoreDip(methyl,binNum)),by=.(featureID,sampleName)]
setkey(sampleScores,sampleName)
setkey(annoDF,sampleName)
sampleScores=merge(sampleScores,annoDF,all.x=TRUE)

#visualizing MIRA signature
plotMIRARegions(binnedRegDT = bigBinDT)


