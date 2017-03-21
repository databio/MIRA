#Unit tests
library(MIRA)
library(GenomicRanges)
library(data.table)

# context("Testthat context...")
# 
# #test that input to BSBinAggregate and wrappers is a data.table not a list with a data.table in it
# #test that BSBinAggregate can deal with 0 methylation scores and low methylation in general 
# #test that there is no division by zero during scoring
# # and that it is not possible to get a negative proportion of methylation
# #test that nonunique regionSet gives result for the unique version of that region set
# #check what happens if annotation is messed up


######Testing binning and aggregation functions#######
#making test data
#making fake BSDT data to test

devtools::load_all(paste0(getOption("CODEBASE"),"MIRA"))#UPDATE, remove later
numBins=13 #don't change this (code below depends on it)
numCpG=numBins
regionNum=2 #changing will break code
chr=rep("chr1",numCpG*regionNum)#numCpG times number of separate regions in start 
#CpGs are 20 apart within region and first CpGs of each region are 1000 apart
start=c((1:numCpG)*20,((1:numCpG)*20+1000))
readCount=rep(10000,numCpG*regionNum)
#linear increase
methyl=c(rep(1:numCpG/numCpG,regionNum))
hitCount=readCount*methyl
sampleName=rep("TestData",numCpG*regionNum)
testBSDT=data.table(chr,start,hitCount,readCount,methyl,sampleName)


#making fake region data to test
chr=rep("chr1",regionNum)
start=c(10,10+1000)
end=c(numCpG*20+10,numCpG*20+1010)
strand=c("+","-")
testGR=GRanges(seqnames = chr,ranges =IRanges(start, end),strand=strand)
testGRDT=grToDt(testGR,includeStrand = TRUE)
#cleaning up variables so they are not used by data.table j expressions later
rm(list=c("chr","start","end","strand","methyl","hitCount","readCount"))
rm("sampleName")

test_that("binRegion is working and flipping - strand"){
    #testing that binRegion is reversing - strands during binning if given strand info
    
    #Using two regions, first has "+", second "-" strand
    #testing binregion as a "j command"
    testBins=testGRDT[,binRegion(start,end,numBins,chr,strand)]
    expect_equal(testBins$id,rep(c(1,2),each=13))#testing id col
    expect_equal(testBins$binID,rep(1:numBins,2))#testing binID col
    expect_equal(testBins[17,start],1190) #spot check
    #testing range of binned region
    s2=testGRDT[2,start]
    e2=testGRDT[2,end]
    #lowest value is original start
    expect_equal(min(testBins[1:numBins+numBins,start]),s2)
    expect_true(max(testBins[1:numBins+numBins,start])<e2)#all starts are in range
    #highest value is original end
    expect_equal(max(testBins[1:numBins+numBins,end]),e2)
    expect_true(min(testBins[1:numBins+numBins,end])>s2)#all ends are in range
    #check that size of output is as expected (each Region should have been split up into binNum rows)
    expect_equal(nrow(testBins),(numBins * regionNum) )
    #last start val should be less than first if orientation was flipped
    #due to "-" strand on second region as it should have been
    expect_true(testBins[numBins*2,start]<testBins[numBins+1,start])
    
    #testing binRegion with two "+" strands
    testGRDT[,strand:=c("+","+")]
    testBins=testGRDT[,binRegion(start,end,numBins,chr,strand)]
    expect_equal(testBins[17,start],1070) #spot check
    #last start val should be greater than first if orientation was not flipped
    expect_true(testBins[numBins*2,start]>testBins[numBins+1,start])
    
    
    #testing binRegion with no strand (should be equivalent to "*")
    testBins=testGRDT[,binRegion(start,end,numBins,chr)]
    expect_equal(testBins[17,start],1070) #spot check
    #last start val should be greater than first if orientation was not flipped
    expect_true(testBins[numBins*2,start]>testBins[numBins+1,start])
    #testing range of binned region
    s2=testGRDT[2,start]
    e2=testGRDT[2,end]
    #lowest value is original start
    expect_equal(min(testBins[1:numBins+numBins,start]),s2)
    expect_true(max(testBins[1:numBins+numBins,start])<e2)#all starts are in range
    #highest value is original end
    expect_equal(max(testBins[1:numBins+numBins,end]),e2)
    expect_true(min(testBins[1:numBins+numBins,end])>s2)#all ends are in range
    
}

test_that("BSAggregate"){
    #setting up prereqs for BSAggregate
    binnedDT=testGRDT[, binRegion(start, end, numBins, chr,strand)]
    binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
    #
    binnedBSDT = BSAggregate(BSDT=testBSDT, regionsGRL=GRangesList(binnedGR), jCommand=buildJ(c("methyl", "readCount"), c("mean", "sum")), byRegionGroup=TRUE, splitFactor=NULL)
    #all should be the same since having opposite strands resulted in
    #flipping the second region
    expect_equal(round(binnedBSDT[,methyl],2),rep(.54,nrow(binnedBSDT)))
    
    #now with no flipping of the second region (changing - to +)
    testGRDT[,strand := c("+","+")]
    #adding to readCount so that not all readCount values will be the same
    #should not affect methylation aggregation
    testBSDT[numBins*2,readCount:= readCount*2]
    binnedDT=testGRDT[, binRegion(start, end, numBins, chr,strand)]
    binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
    binnedBSDT = BSAggregate(BSDT=testBSDT, regionsGRL=GRangesList(binnedGR), jCommand=buildJ(c("methyl", "readCount"), c("mean", "sum")), byRegionGroup=TRUE, splitFactor=NULL)
    #all should be the same since having opposite strands resulted in
    #flipping the second region
    expect_equal(round(binnedBSDT[,methyl],2),round(1:numBins/numBins,2))
    #testing that expected readcount in each bin is obtained
    expect_equal(binnedBSDT[,readCount],c(rep(20000,numBins-1),30000))
    
    #now with no strand info given
    #tests that no strand will be symmetrically averaged 
    #changing first value so symmetry comes from averaging, not from input shape
    testBSDT[1,methyl := methyl*2]
    binnedDT=testGRDT[, binRegion(start, end, numBins, chr)] #no strand
    binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
    binnedBSDT = BSAggregate(BSDT=testBSDT, regionsGRL=GRangesList(binnedGR), jCommand=buildJ(c("methyl", "readCount"), c("mean", "sum")), byRegionGroup=TRUE, splitFactor=NULL)
    #hard coded so may need to be updated later
    expect_equal(round(binnedBSDT[,methyl],2),c(.56,rep(.54,numBins-2),.56))
    
    #test more than two CpG sites in one bin
    testBSDT[24,start:=(testBSDT[25,start]-2)]#putting 24th CpG very close to 25th CpG
    binnedDT=testGRDT[, binRegion(start, end, numBins, chr,strand)] 
    binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
    binnedBSDT = BSAggregate(BSDT=testBSDT, regionsGRL=GRangesList(binnedGR), jCommand=buildJ(c("methyl", "readCount"), c("mean", "sum")), byRegionGroup=TRUE, splitFactor=NULL)
    expect_equal(round(binnedBSDT[12,methyl],2),0.9)
    
    #test that error is given if BSDT input does not have a "methyl"...
    #...column with proportion of reads methylated
    testBSDT[,methyl:=NULL]
    expect_error(BSAggregate(BSDT=testBSDT, regionsGRL=GRangesList(binnedGR), jCommand=buildJ(c("methyl", "readCount"), c("mean", "sum")), byRegionGroup=TRUE, splitFactor=NULL))

}

###############testing scoring and wrapper functions#################

test_that("returnMIRABins and MIRAScore"){
    #testing returnMIRABins
    binnedBSDT=returnMIRABins(BSDT = testBSDT,GRList=testGR,binNum=numBins,minReads = 0)
    expect_equal(round(binnedBSDT[,methyl],2),rep(0.54,numBins))
    #testing that expected readcount in each bin is obtained
    expect_equal(binnedBSDT[,readCount],c(rep(20000,numBins)))
    
    
    #testing MIRAScore
    scoreDT=MIRAScore(BSDT = testBSDT,GRList=testGR,binNum=numBins,minReads = 0,scoringMethod = "logRatio")
    #since the MIRA signature was flat, score=0
    expect_equal(scoreDT$score,0)
    #region set name should have been given automatically
    expect_equal(scoreDT$featureID,"RegionSet1")
}

test_that("scoreDip"){
    #test with odd bin number
    x=-10:10
    y=x^2
    binNumber=length(y) #21
    #default shoulderShift is based on length of input (number of bins)
    testScore=round(scoreDip(values = y,binCount = binNumber),2)
    expScore=round(log((mean(y[c(1,binNumber)]))/(mean(y[10:12]))),2)
    expect_equal(testScore,expScore)
    
    testScore=scoreDip(values = y,binCount = length(y),shoulderShift = 9.5)
    #check by hand, expect_equal(round(testScore,2),3.62)
    
    #test with even bin number
    x=-4:5
    y=x^2
    binNumber=length(y) #10
    testScore=round(scoreDip(values = y,binCount = binNumber),2)
    expScore=round(log(mean(y[c(1,binNumber)])/((0.5*y[4]+y[5]+y[6]+.5*y[7])/3)),2)
    expect_equal(testScore,expScore)
    #test with non default shoulderShift
    testScore=round(scoreDip(values = y,binCount = binNumber,shoulderShift=3),2)
    expScore=round(log(mean(y[c(2,9)])/((0.5*y[4]+y[5]+y[6]+.5*y[7])/3)),2)
    expect_equal(testScore,expScore)
    
    
}

##########testing smaller functions##########

test_that("addMethCol"){
    #Note: this function assigns methyl column by reference (:=)
    #still returns object from function call though
    
    #testing if it returns expected column names, only adds "methyl"
    testBSDT=testBSDT[,methyl:=NULL]#removing methyl column
    oldColNames=colnames(copy(testBSDT))#copy avoids updating by reference
    addMethCol(testBSDT)#by reference, still returns object though
    newColNames=colnames(testBSDT)
    colDiff=base::setdiff(newColNames,oldColNames)
    expect_equal(colDiff,"methyl")
    #testing that methyl values are as expected
    expect_equal(round(testBSDT$methyl,2),rep(round(1:numCpG/numCpG,2),2))
}

#############cleaning up the environment after finishing tests##############




