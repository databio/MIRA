#Unit tests
library(MIRA)
library(GenomicRanges)
library(data.table)
library(bsseq)

context("Testthat context...")
# 
# #test that BSBinAggregate can deal with 0 methylation scores and 
# #low methylation in general 
# # and that it is not possible to get a negative proportion of methylation
# #test that nonunique regionSet gives result for the unique version of that region set
# #check what happens if annotation is messed up


######Testing binning and aggregation functions#######
#making test data
#making fake BSDT data to test

numBins <- 13 #don't change this (code below depends on it)
numCpG <- numBins
regionNum <- 2 #changing will break code
chr <- rep("chr1", numCpG * regionNum)#numCpG times number of separate regions in start 
#CpGs are 20 apart within region and first CpGs of each region are 1000 apart
start <- c((1:numCpG) * 20, ((1:numCpG) * 20 + 1000))
coverage <- rep(10000, numCpG * regionNum)
#linear increase
methylProp <- c(rep(1:numCpG / numCpG, regionNum))
methylCount <- coverage * methylProp
sampleName <- rep("TestData", numCpG * regionNum)
testBSDT <- data.table(chr, start, methylCount, coverage, methylProp, sampleName)
origtestBSDT <- copy(testBSDT) #so it can be the same for each test section

#making fake region data to test
chr <- rep("chr1", regionNum)
start <- c(10, 10 + 1000)
end <- c(numCpG * 20 + 10, numCpG * 20 + 1010)
strand <- c("+", "-")
testGR <- GRanges(seqnames = chr, ranges = IRanges(start, end), strand = strand)
testGRDT <- MIRA:::grToDt(testGR, includeStrand = TRUE)
origtestGRDT <- copy(testGRDT)
origtestGR <- copy(testGR)
#cleaning up variables so they are not used by data.table j expressions later
rm(list = c("chr", "start", "end", "strand", "methylProp", "methylCount", "coverage"))
rm("sampleName")

#(re)setting test objects
testBSDT <- copy(origtestBSDT)
testGR <- copy(origtestGR)
testGRDT <- copy(origtestGRDT)
test_that("binRegion is working and flipping - strand", {
    #making sure input is as I expect
    expect_equal(testGRDT[, strand], c("+", "-"))
    
    #testing that binRegion is reversing - strands during binning if given strand info
    
    #Using two regions, first has "+", second "-" strand
    #testing binregion as a "j command"
    testBins <- testGRDT[, binRegion(start, end, numBins, chr, strand)]
    expect_equal(testBins$id, rep(c(1, 2), each = 13))#testing id col
    expect_equal(testBins$binID, rep(1:numBins, 2))#testing binID col
    expect_equal(testBins[17, start], 1190) #spot check
    #testing range of binned region
    s2 <- testGRDT[2, start]
    e2 <- testGRDT[2, end]
    #lowest value is original start
    expect_equal(min(testBins[1:numBins + numBins, start]), s2)
    expect_true(max(testBins[1:numBins + numBins, start]) < e2)#all starts are in range
    #highest value is original end
    expect_equal(max(testBins[1:numBins + numBins, end]), e2)
    expect_true(min(testBins[1:numBins + numBins, end]) > s2)#all ends are in range
    #check that size of output is as expected 
    #(each Region should have been split up into binNum rows)
    expect_equal(nrow(testBins), (numBins * regionNum) )
    #last start val should be less than first if orientation was flipped
    #due to "-" strand on second region as it should have been
    expect_true(testBins[numBins * 2, start] < testBins[numBins + 1, start])
    
    #testing binRegion with two "+" strands
    testGRDT[, strand := c("+", "+")]
    testBins <- testGRDT[, binRegion(start, end, numBins, chr, strand)]
    expect_equal(testBins[17, start], 1070) #spot check
    #last start val should be greater than first if orientation was not flipped
    expect_true(testBins[numBins * 2, start] > testBins[numBins + 1, start])
    
    
    #testing binRegion with no strand (should be equivalent to "*")
    testBins <- testGRDT[, binRegion(start, end, numBins, chr)]
    expect_equal(testBins[17, start], 1070) #spot check
    #last start val should be greater than first if orientation was not flipped
    expect_true(testBins[numBins * 2, start] > testBins[numBins + 1, start])
    #testing range of binned region
    s2 <- testGRDT[2, start]
    e2 <- testGRDT[2, end]
    #lowest value is original start
    expect_equal(min(testBins[1:numBins + numBins, start]), s2)
    #all starts are in range
    expect_true(max(testBins[1:numBins + numBins, start]) < e2)
    #highest value is original end
    expect_equal(max(testBins[1:numBins + numBins, end]), e2)
    expect_true(min(testBins[1:numBins + numBins, end]) > s2)#all ends are in range
    
})

#resetting test objects
testBSDT <- copy(origtestBSDT)
testGR <- copy(origtestGR)
testGRDT <- copy(origtestGRDT)
test_that("BSAggregate", {
    # NOTE: BSAggregate still can output a coverage column, (not sumCoverage)
    #making sure input is as I expect
    expect_equal(testGRDT[, strand], c("+", "-"))
    
    #setting up prereqs for BSAggregate
    binnedDT <- testGRDT[, binRegion(start, end, numBins, chr, strand)]
    binnedGR <- sapply(split(binnedDT, binnedDT$binID), dtToGr)
    #
    binnedBSDT <- 
        BSAggregate(BSDT = testBSDT, regionsGRL = GRangesList(binnedGR), 
                    jCommand = buildJ(c("methylProp", "coverage"), 
                                      c("mean", "sum"),
                                      newColNames=c("methylProp", "coverage")), 
                    byRegionGroup = TRUE, splitFactor = NULL)
    #all should be the same since having opposite strands resulted in
    #flipping the second region
    expect_equal(round(binnedBSDT[, methylProp], 2), rep(.54, nrow(binnedBSDT)))
    
    #now with no flipping of the second region (changing - to + )
    testGRDT[, strand := c("+", "+")]
    #adding to coverage so that not all coverage values will be the same
    #should not affect methylation aggregation
    testBSDT[numBins * 2, coverage := coverage * 2]
    binnedDT <- testGRDT[, binRegion(start, end, numBins, chr, strand)]
    binnedGR <- sapply(split(binnedDT, binnedDT$binID), dtToGr)
    binnedBSDT <-
        BSAggregate(BSDT = testBSDT, regionsGRL = GRangesList(binnedGR), 
                    jCommand = buildJ(c("methylProp", "coverage"), 
                                      c("mean", "sum"),
                                      newColNames=c("methylProp", "coverage")), 
                    byRegionGroup = TRUE, splitFactor = NULL)
    #all should be the same since having opposite strands resulted in
    #flipping the second region
    expect_equal(round(binnedBSDT[, methylProp], 2), round(1:numBins / numBins, 2))
    #testing that expected coverage in each bin is obtained
    expect_equal(binnedBSDT[, coverage], c(rep(20000, numBins - 1), 30000))
    
    #now with no strand info given
    #tests that no strand will be symmetrically averaged 
    #changing first value so symmetry comes from averaging, not from input shape
    testBSDT[1, methylProp := methylProp * 2]
    binnedDT <- testGRDT[, binRegion(start, end, numBins, chr)] #no strand
    binnedGR <- sapply(split(binnedDT, binnedDT$binID), dtToGr)
    binnedBSDT <- 
        BSAggregate(BSDT = testBSDT, regionsGRL = GRangesList(binnedGR), 
                    jCommand = buildJ(c("methylProp", "coverage"), 
                                      c("mean", "sum"),
                                      newColNames=c("methylProp", "coverage")), 
                    byRegionGroup = TRUE, splitFactor = NULL)
    #hard coded so may need to be updated later
    expect_equal(round(binnedBSDT[, methylProp], 2), 
                 c(.56, rep(.54, numBins - 2), .56))
    
    #test more than two CpG sites in one bin
    #putting 24th CpG very close to 25th CpG
    testBSDT[24, start := (testBSDT[25, start] - 2)]
    binnedDT <- testGRDT[, binRegion(start, end, numBins, chr, strand)] 
    binnedGR <- sapply(split(binnedDT, binnedDT$binID), dtToGr)
    binnedBSDT <-
        BSAggregate(BSDT = testBSDT, regionsGRL = GRangesList(binnedGR), 
                    jCommand = buildJ(c("methylProp", "coverage"), 
                                      c("mean", "sum"),
                                      newColNames=c("methylProp", "coverage")), 
                    byRegionGroup = TRUE, splitFactor = NULL)
    expect_equal(round(binnedBSDT[12, methylProp], 2), 0.9)
    
    #test that error is given if BSDT input does not have a "methylProp"...
    #...column with proportion of reads methylated
    testBSDT[, methylProp := NULL]
    expect_error(BSAggregate(BSDT = testBSDT, 
                             regionsGRL = GRangesList(binnedGR), 
                             jCommand = buildJ(c("methylProp", "coverage"), 
                                               c("mean", "sum"),
                                               newColNames=c("methylProp", "coverage")), 
                             byRegionGroup = TRUE, splitFactor = NULL))

})

###############testing scoring and wrapper functions#################
#resetting test objects
testBSDT <- copy(origtestBSDT)
testGR <- copy(origtestGR)
testGRDT <- copy(origtestGRDT)
test_that("aggregateMethyl and MIRAScore", {
    # NOTE: aggregateMethyl does not output a coverage column, now sumCoverage
    #making sure input is as I expect
    expect_equal(testGRDT[, strand], c("+", "-"))
    
    #testing aggregateMethyl, warning about names is expected
    # ignore warning about needing a named list/GrangesList
    binnedBSDT <- aggregateMethyl(BSDT = testBSDT, GRList = testGR, 
                              binNum = numBins, minBaseCovPerBin = 0)
    expect_equal(round(binnedBSDT[, methylProp], 2), rep(0.54, numBins))
    #testing that expected sumCoverage in each bin is obtained
    expect_equal(binnedBSDT[, sumCoverage], c(rep(20000, numBins)))
    
    # making sure the output is same for input of data.table or bsseq
    # making bsseq version of testBSDT (there should only be one sampleName)
    testBSDTbsseq <- BSseq(M=as.matrix(testBSDT$methylCount),
                          Cov=as.matrix(testBSDT$coverage),
                          chr=testBSDT$chr,
                          pos=testBSDT$start,
                          sampleNames=unique(testBSDT$sampleName))
    binnedBSDTbsseq <- aggregateMethyl(BSDT = testBSDTbsseq, GRList = testGR,
                                    binNum = numBins, minBaseCovPerBin = 0)
    names(binnedBSDTbsseq) <- NULL # taking off name for the following comparison
    expect_equal(binnedBSDTbsseq, list(binnedBSDT))
    # test whether aggregateMethyl deals with cases where BSDT doesn't overlap 
    # with region set/s or minBaseCovPerBin causes region set to be screened out
    # screened out
    screenedBSDT = suppressWarnings(aggregateMethyl(BSDT = testBSDT, GRList = testGR,
                                    binNum = numBins, minBaseCovPerBin = 1000000))
    expect_equal(nrow(screenedBSDT), 0)
    # for no overlap
    noOLGR <- GRanges(seqnames = c("chrN", "chrN"), 
                      ranges = IRanges(c(1, 100), c(10, 110)), 
                      strand = c("*", "*"))
    noOLBinnedBSDT= suppressWarnings(aggregateMethyl(BSDT = testBSDT, GRList = noOLGR,
                                     binNum = numBins, minBaseCovPerBin = 0))
    expect_equal(nrow(noOLBinnedBSDT), 0)

    
    
    # testing MIRAScore, warning about names is expected
    # ignore warning about needing a named list/GrangesList
    scoreDT <- MIRAScore(BSDT = testBSDT, GRList = testGR, binNum = numBins, minBaseCovPerBin = 0, 
                      scoringMethod = "logRatio")
    #since the MIRA signature was flat, score = 0
    expect_equal(scoreDT$score, 0)
    #region set name should have been given automatically
    expect_equal(scoreDT$featureID, "RegionSet1")
})

#resetting test objects
testBSDT <- copy(origtestBSDT)
testGR <- copy(origtestGR)
testGRDT <- copy(origtestGRDT)
test_that("calcMIRAScore, findShoulder, and isProfileConcaveUp", {
    
    #test with odd bin number
    x <- -10:10
    y <- x^2 + 1
    binNumber <- length(y) #21
    #default shoulderShift is based on length of input (number of bins)
    # ignore warning about "essentially perfect fit"
    testScore <- round(calcMIRAScore(binnedDT = y), 2)
    expScore <- round(log((mean(y[c(1, binNumber)])) / (y[11])), 2)
    expect_equal(testScore, expScore)
    
    #testing that averaging will happen if middle is not lowest value
    y[11] <- 5
    testScore <- round(calcMIRAScore(binnedDT = y), 2)
    expScore <- round(log((mean(y[c(1, binNumber)])) / ((y[10] + y[11] + y[12]) / 3)), 2)
    expect_equal(testScore, expScore)
    
    
    
    #testScore <- calcMIRAScore(binnedDT = y, shoulderShift = 9.5)
    #check by hand, expect_equal(round(testScore, 2), 3.62)
    
    #test with even bin number
    x <- -4:5
    y <- x^2 + 1
    binNumber <- length(y) #10
    # ignore warning about "essentially perfect fit"
    testScore <- round(calcMIRAScore(binnedDT = y), 2)
    expScore <- round(log(mean(y[c(1, binNumber)])
                        / y[5]), 2)
    expect_equal(testScore, expScore)
    #test with non default shoulderShift
    # ignore warning about "essentially perfect fit"
    testScore <- round(calcMIRAScore(binnedDT = y, shoulderShift = 3), 2)
    expScore <- round(log(mean(y[c(2, 9)]) / y[5]), 2)
    expect_equal(testScore, expScore)
    
    #testing results with automatic shoulder detection, symmetrical
    #chooses pos. 12 because (14+16)/2 is not > 15
    #this test will fail if standard is changed to >=
    jagged <- c(16, 14, 15, 11, 12, 7, 4, 4, 7, 12, 11, 15, 14, 16)
    testScore <- round(calcMIRAScore(binnedDT = jagged, shoulderShift = "auto"), 2)
    expScore <- round(log(mean(jagged[c(3, 12)]) / jagged[7]), 2)
    expect_equal(testScore, expScore)
    
    #testing that averaging will happen if middle point is not lowest
    jagged[c(7, 8)] <- 8 
    testScore <- round(calcMIRAScore(binnedDT = jagged, shoulderShift = "auto"), 2)
    expMidpoint <- (.5 * jagged[6] + jagged[7] + jagged[8] + .5 * jagged[9]) / 3
    expScore <- round(log(mean(jagged[c(3, 12)]) / expMidpoint), 2)
    expect_equal(testScore, expScore)
    
    #testing the findShoulder function with unsymmetrical input
    unsymmetric <- jagged
    unsymmetric[1] <- 17 #now the first spot should be counted as the shoulder
    centerSpot <- (1 + length(unsymmetric)) / 2
    rShoulderShift <- findShoulder(unsymmetric, length(unsymmetric), 
                                  centerSpot, whichSide = "right")
    expect_equal(rShoulderShift, 4.5)
    #middle is 7.5 outer edge is 1, shoulderShift should be 6.5
    lShoulderShift <- findShoulder(unsymmetric, length(unsymmetric), 
                                  centerSpot, whichSide = "left")
    expect_equal(lShoulderShift, 6.5)
    
    #testing calcMIRAScore function with unsymmetrical input
    unSymScore <- round(calcMIRAScore(binnedDT = unsymmetric, 
                                shoulderShift = "auto", usedStrand = TRUE), 2)
    expMidpoint <- (.5 * unsymmetric[6] + unsymmetric[7] + 
                       unsymmetric[8] + .5 * unsymmetric[9]) / 3
    #if symmetrical, you would use values 1 and 14 or 3 and 12
    expScore <- round(log(mean(unsymmetric[c(1, 12)]) / expMidpoint), 2)
    expect_equal(expScore, unSymScore)
    
    # testing calcMIRAScore with input of a data.table and no data.table syntax
    # when calling calcMIRAScore
    jagged <- c(16, 14, 15, 11, 12, 7, 4, 4, 7, 12, 11, 15, 14, 16)
    jagLen <- length(jagged)
    binnedDT <- data.table(bin=seq_along(jagged), methylProp=jagged, 
                          sampleName=rep("Sample1", jagLen), 
                          featureID=rep("RegionSet1", jagLen))
    newSyntaxScoreDT <- calcMIRAScore(binnedDT)
    newSyntaxScore <- round(newSyntaxScoreDT$score, 2)
    expScore <- round(log(mean(jagged[c(3, 12)]) / jagged[7]), 2)
    oldSyntaxScore <- round(binnedDT[, .(score=calcMIRAScore(methylProp)), 
                                    by=.(sampleName, featureID)]$score, 2)
    expect_equal(expScore, newSyntaxScore)
    expect_equal(oldSyntaxScore, newSyntaxScore)
    
    # testing isProfileConcaveUp
    jagged <- c(16, 14, 15, 11, 12, 7, 4, 4, 7, 12, 11, 15, 14, 16)
    expect_true(isProfileConcaveUp(jagged, 14))
    invJagged <- 1 - jagged
    expect_true(!isProfileConcaveUp(invJagged, 14))
    miniDip <- c(1, 3, 5, 8, 7, 5, 4, 5, 7, 8, 5, 3, 1)
    expect_true(!isProfileConcaveUp(miniDip, 13))
    miniDip2 <- c(1, miniDip, 1)
    expect_true(isProfileConcaveUp(miniDip2, 15))
    
    # testing that calcMIRAScore will give expected score
    testScore <- round(calcMIRAScore(miniDip2), 3)
    expScore <- round(log(mean(miniDip2[c(5,11)])/miniDip2[8]), 3)
    expect_equal(testScore, expScore)
    invMiniDip2 <- 10 - miniDip2
    testScore <- round(calcMIRAScore(invMiniDip2), 3)
    expScore <- round(log(mean(invMiniDip2[c(5,11)])/invMiniDip2[8]), 3)
    expect_equal(testScore, expScore)
})

##########testing miscellaneous functions##########
#resetting test objects
testBSDT <- copy(origtestBSDT)
testGR <- copy(origtestGR)
testGRDT <- copy(origtestGRDT)
test_that("addMethPropCol", {
    #Note: this function assigns methylProp column by reference ( := )
    #still returns object from function call though
    
    #making sure input is as I expect
    expect_equal(testGRDT[, strand], c("+", "-"))
    
    #testing if it returns expected column names, only adds "methylProp"
    testBSDT <- testBSDT[, methylProp := NULL]#removing methylProp column
    oldColNames <- colnames(copy(testBSDT))#copy avoids updating by reference
    addMethPropCol(testBSDT)#by reference, still returns object though
    newColNames <- colnames(testBSDT)
    colDiff <- base::setdiff(newColNames, oldColNames)
    expect_equal(colDiff, "methylProp")
    #testing that methylProp values are as expected
    expect_equal(round(testBSDT$methylProp, 2), rep(round(1:numCpG / numCpG, 2), 2))
})

#resetting test objects
testBSDT <- copy(origtestBSDT)
testGR <- copy(origtestGR)
testGRDT <- copy(origtestGRDT)
test_that("BSreadBiSeq",{
    
    tRRBS <- BSreadBiSeq(system.file("extdata", "shortRRBS.bed", package = "MIRA"))
    tRRBSCols <- c("chr", "start", "methylCount", "coverage")
    expect_equal(tRRBSCols, colnames(tRRBS))
    # coverage should be greater than or equal to methylCount, true for all
    compareCounts <- (tRRBS$coverage >= tRRBS$methylCount)
    expect_true(all(compareCounts))
})

#resetting test objects
testBSDT <- copy(origtestBSDT)
testGR <- copy(origtestGR)
testGRDT <- copy(origtestGRDT)
test_that("SummarizedExperimentToDataTable and bsseqToDataTable",{
    data("exampleBSseqObj")
    MIRAFormatBSDTList <- SummarizedExperimentToDataTable(coordinates = bsseq::granges(exampleBSseqObj), 
                                      methylCountDF = bsseq::getCoverage(BSseq = exampleBSseqObj, type = "M"), 
                                      coverageDF = bsseq::getCoverage(BSseq = exampleBSseqObj, type = "Cov"), 
                                      methylPropDF = bsseq::getMeth(BSseq = exampleBSseqObj, type = "raw"), 
                                      sample_names = bsseq::sampleNames(exampleBSseqObj))
    numLoci <- nrow(exampleBSseqObj)
    # checking that screening out noncovered C's is working
    # also accomplishes basic check that general length of output is as expected
    numNotCovered1 <- sum(getCoverage(exampleBSseqObj[, 1]) == 0)
    numNotCovered2 <- sum(getCoverage(exampleBSseqObj[, 2]) == 0)
    expect_equal(nrow(MIRAFormatBSDTList[[1]]), numLoci - numNotCovered1)
    expect_equal(nrow(MIRAFormatBSDTList[[2]]), numLoci - numNotCovered2)
    withinRange <- all((MIRAFormatBSDTList[[1]][, methylProp] <= 1) &
                          (MIRAFormatBSDTList[[1]][, methylProp] >= 0))
    # checking that methylProp values are within expected range
    expect_true(withinRange)
    # making sure coverage is greater than/equal to methylCount
    isCovGreater <- all(MIRAFormatBSDTList[[1]][, coverage] >= 
                           MIRAFormatBSDTList[[1]][, methylCount])
    expect_true(isCovGreater)
    
    
    # quick check of bsseqToDataTable
    MIRAFormatBSDTList2 <- bsseqToDataTable(exampleBSseqObj)
    expect_equal(MIRAFormatBSDTList, MIRAFormatBSDTList2)
})

#############cleaning up the environment after finishing tests##############




