
#loading sample annotation
#annoDF object, classes: data.table, data.frame
load(system.file("data","annoDF.RData",package="MIRA")) 

#reading in region data and bisulfite data
someFeatures=get(load(system.file("data","Gm12878Nrf1_Subset.RData",package="MIRA")))
BSDT=get(load(system.file("data","GM06990_1_ExampleSet.RData",package = "MIRA")))

#converting to lists 
BSDTList=list(BSDT)

#doing MIRA analysis
bigBin=lapply(X = BSDTList,FUN = returnMIRABins,GRList=someFeatures,binNum=11,sampleNameInBSDT=TRUE)

bigBinDT=bigBin[[1]]

#scoring samples
sampleScores=bigBinDT[,.(score = scoreDip(methyl,binCount=11)),by=.(featureID,sampleName)]



