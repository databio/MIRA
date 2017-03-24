# PACKAGE DOCUMENTATION
#' Methylation-based Inference of Regulatory Activity (MIRA)
#' MIRA is a score that measures the degree of dip in methylation
#' level surrounding a regulatory site of interest, such as a 
#' transcription factor binding sites.
#' This script provides functions for aggregating methylation 
#' data across region sets, in bins.
#'
#' @docType package
#' @name MIRA
#' @author Nathan Sheffield
#'
#' @references \url{http://github.com/databio}
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand seqnames 
#'             granges
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_jitter geom_line
#' @import BiocGenerics S4Vectors IRanges
#' @importFrom data.table ":=" setDT data.table setkey fread setnames 
#'             setcolorder melt setkeyv rbindlist setattr
#' @importFrom Biobase sampleNames
#' @importFrom bsseq getBSseq hasBeenSmoothed
NULL


# Because of some issues with CRAN submission, 
#(see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation, 
# in order to pass some R CMD check NOTES.
if(getRversion() >= "2.15.1"){
    utils::globalVariables(c(
    ".", "V1", "chr", "featureID", "hitCount", "meth", "methyl", "readCount", 
    "regionGroupID", "regionID", "sampleName", "sampleType"))
}


#' Function to aggregate methylation data into bins 
#' over all regions in each region set;
#' 
#'
#' @param BSDT A single data table that has DNA methylation data on individual 
#' sites including a "chr" column with chromosome, a "start" column with the 
#' coordinate number for the cytosine, a "methyl" column with proportion of 
#' methylation (0 to 1), a "hitCount" column with number of methylated reads for
#' each site, and a "readCount" column with total number of reads for each site.
#' A "sampleName" column is preferred.
#' @param GRList A GRangesList object containing region sets, each set 
#' corresponding to a regulatory element. Each regionSet in the list should be 
#' named. A named list of data.tables also works. 
#' @param binNum How many bins each region should be split into for aggregation 
#' of the DNA methylation data.
#' @param minReads Filter out bins with fewer than X reads before returning.
#' @param sampleNameInBSDT boolean for whether the BSDT has a sampleName column
#' 
#' @return a data.table with binNum rows for each region set containing
#' aggregated methylation data.
#' Each region was split into bins; methylation was put in these bins; 
#' Output contains sum of the all corresponding bins for the regions of each 
#' region set, ie for all regions in each region set: first bins summed, second 
#' bins summed, etc. Columns of the output should be "regionGroupID", "methyl", 
#' "readCount", "featureID", and possibly "sampleName".
#' 
#' @export
#' @examples
#' data("GM06990_1_ExampleSet", package = "MIRA") #exampleBSDT
#' data("Gm12878Nrf1_Subset", package = "MIRA") #exampleRegionSet
#' exBinDT = returnMIRABins(exampleBSDT, exampleRegionSet)
returnMIRABins = function(BSDT, GRList, binNum = 11, minReads = 500, 
                          sampleNameInBSDT = TRUE){
  
    #########returnMIRABins:Preprocessing and formatting###############
    #converting to list format if GRList is a data.table or GRanges object
    if (class(GRList) %in% "GRanges"){
      GRList = GRangesList(GRList)
      message("Converting to GRangesList...")
    }
    if (class(GRList) %in% "data.table"){
        GRList = list(GRList)
        message("Converting to list...")
    }

    #checking that input is in list format
    if (!class(GRList) %in% c("list", "GRangesList")){
        stop("GRList should be a named list/GRangesList.")
    }

    #checking if region sets have names
    if (is.null(names(GRList))){
        warning("GRList should be a named list/GRangesList. 
                Sequential names given according to order in object.")
        names(GRList)<-paste0(rep("RegionSet", length(GRList)), 1:length(GRList))
    }

    #checking that all objects in GRList are the same type 
    #and converting to data.tables
    if (all(sapply(X = GRList, FUN = class) %in% "GRanges")){
        #GRanges to data.tables
        GRDTList = lapply(X = GRList, FUN = grToDt, includeStrand = TRUE)
        #below statement will be true if all objects in the list are of 
        #class data.table
        #necessary since data.tables also include data.frame as a class
    }else if (all(sapply(
        X = lapply(X = GRList, FUN = function(x) class(x) %in% "data.table"), 
        FUN = any))){
        
        GRDTList = GRList #this case is okay
    }else{
        stop("GRList should be a GRangesList or a list of data.tables")
    }
    
    if (sampleNameInBSDT){
        if(!"sampleName" %in% colnames(BSDT)){
            stop("BSDT should have sampleName col if sampleNameInBSDT = TRUE")
        }
    }

    #adding a methyl column if it is not already in the BSDT
    if (!"methyl" %in% names(BSDT)){
        BSDTList = addMethCol(list(BSDT))
        BSDT = BSDTList[[1]] 
    }

    ########returnMIRABins:Binning and processing output####################


    methylByBin = lapply(X = GRDTList, 
                       FUN = function(x) BSBinAggregate(BSDT = BSDT, rangeDT = x, 
                                                        binCount = binNum, 
                                                        splitFactor = NULL, 
                                                        minReads = minReads))
    names(methylByBin) = names(GRList)#preserving names
    #adding a feature ID column to each data.table that 
    #should identify what region set was used
    for (i in 1:length(methylByBin)){
        methylByBin[[i]][, featureID := rep(names(methylByBin)[i], 
                                           nrow(methylByBin[[i]]))][]
    }
    #screening out region sets that had incomplete binning
    binNumScreen = sapply(X = methylByBin, FUN = nrow)
    #taking out incomplete region sets
    methylByBin = methylByBin[!(binNumScreen<binNum)]

    bigMethylByBin = rbindlist(methylByBin)
    if (sampleNameInBSDT){
        #creating new sampleName column
        bigMethylByBin[, sampleName := rep(BSDT[1, sampleName])][] 
    }


    return(bigMethylByBin)
}

#' Function to take DNA methylation and region data and return a MIRA score;
#' a wrapper for returnMIRABins and scoreDip
#' 
#'
#' @param BSDT A single data table that has DNA methylation data on individual 
#' sites including a "chr" column with chromosome, a "start" column with the 
#' coordinate number for the cytosine, a "methyl" column with proportion of 
#' methylation (0 to 1), a "hitCount" column with number of methylated reads 
#' for each site, and a "readCount" column with total number of reads for each 
#' site. A "sampleName" column is preferred.
#' @param GRList A GRangesList object containing region sets, each set 
#' corresponding to a regulatory element (or having regions with the 
#' same biological annotation).
#' Each regionSet in the list should be named. 
#' @param binNum How many bins each region should be split into for aggregation 
#' of the DNA methylation data.
#' @param scoringMethod Method to calculate MIRA score after binning. 
#' "logRatio" is currently the only option. See scoreDip function.
#' @param sampleNameInBSDT boolean for whether the BSDT has a sampleName column
#' @param minReads Filter out bins with fewer than X reads before returning.
#' 
#' @return A MIRA score for each region set in GRList. See ?scoreDip. 
#' @examples 
#' data("GM06990_1_ExampleSet", package = "MIRA") #exampleBSDT
#' data("Gm12878Nrf1_Subset", package = "MIRA") #exampleRegionSet
#' MIRAScore(BSDT = exampleBSDT, GRList = exampleRegionSet)
#' 
#' @export
MIRAScore = function(BSDT, GRList, binNum = 11, scoringMethod = "logRatio", 
                     sampleNameInBSDT = TRUE, minReads = 500){

    #making sure methyl column is part of input BSDT
    if (!"methyl" %in% names(BSDT)){
        stop("BSDT must have a methyl column with proportion of methylation. 
             addMethCol() will add this.")
    }

    MIRAresults = list()


    bigBin = returnMIRABins(BSDT = BSDT, GRList = GRList, binNum = binNum, 
                          sampleNameInBSDT = sampleNameInBSDT, 
                          minReads = minReads)
  
    #using binned methylation data to calculate MIRA score
    scoreDT = bigBin[, .(score = scoreDip(methyl, binNum, method = scoringMethod)), 
                   by = .(featureID, sampleName)]

    return(scoreDT)
}

#' My dip scoring function - for MIRA scores;
#' That's Methylation-based Inference of Regulatory Activity
#' 
#' @param values A vector with proportion of methylation values for each bin. 
#'  Between 0 and 1.
#' @param binCount Number of bins, also length of "values" vector.
#' @param shoulderShift The number of bins away from the center to use as the
#' shoulders. I have used 5 or 3.
#' @param method The scoring method. "logRatio" is the log of the ratio of
#' the average of outside values (shoulders) divided by 
#' the average of the middle values. 
#' A higher score with "logRatio" corresponds to a deeper dip. "logRatio" is the
#' only scoring method currently but more methods may be added in the future.
#' 
#' @return A MIRA score. The MIRA score quantifies the "dip" of 
#' the MIRA signature which is an aggregation of methylation 
#' over all regions in a region set. 
#' 
#' @export
#' @examples
#' data("exampleBins")
#' binCount = 11 #bin number for exampleBins 
#' exampleBins[, .(score = scoreDip(methyl, binCount)), by = .(featureID, sampleName)]
scoreDip = function(values, binCount, 
                    shoulderShift = floor((binCount-1)/2), method = "logRatio") {
    if (!method %in% "logRatio"){ #add new methods eventually
        stop("Invalid scoring method. Check spelling/capitalization.")
    }
    if (method == "logRatio"){
        centerSpot = (binCount + 1)/2 # X.5 for even binCount
        #floor and ceiling are only relevant when binCount is even 
        #(which means centerSpot is X.5)
        leftSide = floor(centerSpot - shoulderShift)  
        rightSide = ceiling(centerSpot + shoulderShift)
        if ((binCount %% 2) == 0){ #if binCount is even, centerSpot is X.5
            #includes 4 bins but outer two bins are weighted by half
            #approximates having 3 middle bins
            midpoint = (.5*values[centerSpot-1.5] + values[centerSpot-.5]
                      + values[centerSpot + .5] + 0.5*values[centerSpot + 1.5])/3
        }else{#if binCount is odd, centerSpot is X.0
            #three middle bins
            midpoint = (values[centerSpot] + values[centerSpot + 1] 
                       + values[centerSpot-1] ) /3
        }
        shoulders = ((values[leftSide] + values[rightSide])/2)
        if (midpoint<.000001){
            warning("Division by zero. Consider adding a small constant 
                    to all bins for this region set.")
        }
        if (shoulders<.000001){
            warning("Taking log of zero. Consider adding a small constant 
                    to all bins for this region set.")
        }
        # log ratio...
        score = log( shoulders / midpoint )
    }

    # #alternate way of scoring by the area in the dip
    # if (method == "area"){
    #     maxMethyl = max(values)
    #     score = maxMethyl*binCount-sum(values)
    # }

    # #another alternate method
    # if (method == "parabola"){
    #     #fit2 <- lm(y~poly(x, 2, raw = TRUE))
    #     #lines(xx, predict(fit2, data.frame(x = xx)), col = "green")
    # }
    return(score)
}


#' Aggregating signals in bins across a set of regions
#'
#' given a start, end, and number of bins, to divide, 
#' this function will split the regions into bins.
#' Bins will be only approximately the same size, due to rounding.
#' (they should not be more than 1 different).
#'
#' Use case: take a set of regions, like CG islands, and bin them; now you can
#' aggregate signal scores across the bins, giving you an aggregate signal
#' in bins across many regions of the same type.
#'
#' In theory, this just runs on 3 values, but you can run it inside a 
#' data.table j expression to divide a bunch of regions in the same way.
#' @param start Coordinate for beginning of range/range.
#' @param end Coordinate for end of range/region.
#' @param bins How many bins to divide this range/region.
#' @param idDF A string/vector of strings that has chromosome (e.g. "chr1") 
#' for given start and end values
#'
#' @return
#' A data.table, expanded to nrow = number of bins, with these id columns:
#'      id: region ID
#'      binID: repeating ID (this is the value to aggregate across)
#'      ubinID: unique bin IDs
#' @export
#' @examples
#' start = c(100, 1000, 3000)
#' end = c(500, 1400, 3400)
#' chr = c("chr1", "chr1", "chr2")
#' strand = c("*", "*", "*")
#' #strand not included in GRanges object 
#' #since MIRA assumes "*" already unless given something else
#' regionsToBin = GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
#' regionsDT = grToDt(regionsToBin, includeStrand = FALSE)
#' numberOfBins = 15
#' #data.table "j command" using column names and numberOfBins variable
#' binnedRegionDT = regionsDT[, binRegion(start, end, numberOfBins, chr)]
binRegion = function(start, end, bins, idDF = NULL, strand = "*") {
    #if (!is.null(idDF) & ( ! "data.frame"  %in% class(idDF))) {
    #   stop("idDF should be a data.frame")
    #}
    
    #conditionally altered later
    finalColNames = c("chr", "start", "end", "id", "binID", "ubinID")
    if (!"*" %in% strand){
        if ("+" %in% strand){
            plusIndex = which(strand == "+")
            #once for plus strand coordinates
            plusStart = start[plusIndex]
            plusEnd = end[plusIndex]
            
            binSize = (plusEnd-plusStart)/(bins)
            breaks = round(rep(plusStart, each = (bins + 1)) 
                          + (0:(bins)) * rep(binSize, each = (bins + 1)))
            
            endpoints = (bins + 1) * (1:(length(plusStart)))
            startpoints = 1 + (bins + 1)  * (0:(length(plusStart)-1))
            
            plusDT = data.table(start = breaks[-endpoints], 
                                end = breaks[-startpoints], 
                                id = rep(plusIndex, each = bins), 
                                binID = 1:bins, 
                                strand = "+", 
                                key = "id")
            dt = plusDT #placeholder but may be returned
        }
        if ("-" %in% strand){
            minusIndex = which(strand == "-")
        
            minusStart = start[minusIndex]
            minusEnd = end[minusIndex]
        
            binSize = (minusEnd-minusStart)/(bins)
            breaks = round(rep(minusStart, each = (bins + 1)) 
                          + (0:(bins)) * rep(binSize, each = (bins + 1)))
        
            endpoints = (bins + 1) * (1:(length(minusStart)))
            startpoints = 1 + (bins + 1)  * (0:(length(minusStart)-1))
        
            minusDT = data.table(start = breaks[-endpoints], 
                                 end = breaks[-startpoints], 
                                 id = rep(minusIndex, each = bins), 
                                 binID = bins:1, 
                                 strand = "-", 
                                 key = "id")
            dt = minusDT #placeholder but may be returned
        }
        
        #if there are both + and - strands
        #combining and sorting plus and minus data.tables 
        if (("+" %in% strand) && ("-" %in% strand)){
            dt = rbindlist(list(plusDT, minusDT))
            setorder(x = dt, id, binID)#setorder(dt, id) might also work?
        }
        
        finalColNames = c(finalColNames, "strand")#included only if + /- are present
        dt[, ubinID := 1:nrow(dt)] 
        
    }else { #some strand information is "*", don't flip bin directions
    
        binSize = (end-start)/(bins)
        breaks = round(rep(start, each = (bins + 1)) 
                      + (0:(bins)) * rep(binSize, each = (bins + 1)))

        endpoints = (bins + 1) * (1:(length(start)))
        startpoints = 1 + (bins + 1)  * (0:(length(start)-1))
        #do all regions in the same order
        dt = data.table(start = breaks[-endpoints], 
                        end = breaks[-startpoints], 
                        id = rep((1:length(start)), each = bins), 
                        binID = 1:bins, 
                        ubinID = 1:length(breaks[-startpoints]), 
                        key = "id")
    
    }
    
    if (!(is.null(idDF))){
        chr = rep(idDF, each = bins)
        dt = dt[, chr := chr]
        setcolorder(dt, finalColNames)#putting chr first, does not copy
    }
    
    return(dt[])
}

#' A wrapper of BSAggregate that first bins regions and then aggregates
#' each bin across a set of regions, individually.
#'
#' Produced originally for binning Ewing RRBS data across various region sets
#' 
#' @param BSDT A single data table that has DNA methylation data 
#' on individual sites including a "chr" column with chromosome, 
#' a "start" column with the coordinate number for the cytosine, 
#' a "methyl" column with proportion of methylation (0 to 1), 
#' a "hitCount" column with number of methylated reads for each site, and 
#' a "readCount" column with total number of reads for each site.
#' @param rangeDT A data table with the sets of regions to be binned, 
#' with columns named "start", "end".
#' @param binCount Number of bins across the region.
#' @param byRegionGroup Default TRUE will aggregate methylation over 
#' corresponding bins for each region (all bin1's aggregated, all bin2's, etc).
#' byRegionGroup = FALSE is deprecated.
#' @param minReads Filter out bins with fewer than X reads before returning.
#' @param splitFactor With default NULL, aggregation will be done 
#' separately/individually for each sample.
#' 
#' @return With splitFactor = NULL, it will return a data.table 
#' with binCount rows, 
#' containing aggregated methylation data over regions in region set "rangeDT".
#' Each region was split into bins; methylation was put in these bins; 
#' Output contains sum of the all corresponding bins 
#' for the regions of each region set ie for all regions in each region set: 
#' first bins summed, second bins summed, etc.
#' Columns of the output should be "regionGroupID", "methyl", and "readCount"
#' @examples
#' data("GM06990_1_ExampleSet") #exampleBSDT
#' data("Gm12878Nrf1_Subset") #exampleRegionSet
#' exampleRegionSet = grToDt(exampleRegionSet, includeStrand = TRUE)
#' aggregateBins = BSBinAggregate(BSDT = exampleBSDT, rangeDT = exampleRegionSet, 
#'                              binCount = 11, splitFactor = NULL)
#' 
#' @export
BSBinAggregate = function(BSDT, rangeDT, binCount, minReads = 500, 
                          byRegionGroup = TRUE, splitFactor = NULL) {
    if (! "data.table" %in% class(rangeDT)) {
        stop("rangeDT must be a data.table")
    }
    seqnamesColName = "seqnames"  # default column name
        if ("chr" %in% colnames(rangeDT)) {
            message("seqnames column name set to: chr")
            seqnamesColName = "chr"
        } else {
            # Got neither.
            stop("rangeDT must have a seqnames column")
        }
    
    if(! "strand" %in% colnames(rangeDT)){
        rangeDT[, strand := "*"]
        warning("Since strand not given, set to '*' ")
    }

    ##if(!silent){
    #message("Binning...")
    ##}
    binnedDT = rangeDT[, binRegion(start, end, 
                                   binCount, get(seqnamesColName), strand)]
    binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
    #message("Aggregating...")
    binnedBSDT = BSAggregate(BSDT = BSDT, 
                             regionsGRL = GRangesList(binnedGR), 
                             jCommand = buildJ(c("methyl", "readCount"), 
                                             c("mean", "sum")), 
                             byRegionGroup = byRegionGroup, 
                             splitFactor = splitFactor)
    # If we aren't aggregating by bin, then don't restrict to min reads!
    if (byRegionGroup) {
        binnedBSDT = binnedBSDT[readCount >= minReads, ]
        if (nrow(binnedBSDT)<binCount){
            stop("Less than minReads. Unable to give MIRA score.")
        }
    }
    
    
    return(binnedBSDT)
}

#' BSaggregate -- Aggregate a BSDT across regions or region groups, 
#' for multiple samples at a time.
#' This function is as BScombineByRegion, but can handle not only multiple
#' samples in BSDT, but also simultaneously multiple region sets by passing
#' a regionsGRL (GRangesList object).However currently code for symmetrical 
#' averaging will cause it to only work with one region set (which may be
#' split up into multiple GRanges in a GRangesList).
#' you can use jCommand to do other functions.

#' Given a bisulfite data table as input, with an identifier column for
#' different samples; plus a GRanges objects with regions to aggregate.
#'
#' @param BSDT The bisulfite data.table (output from one of the parsing
#' functions for methylation calls) that you wish to aggregate. It can
#' be a combined table, with individual samples identified by column passed
#' to splitFactor. To be safe, "chr", "start", "hitCount", "readCount", and 
#' "methyl" columns should be in BSDT.
#' @param regionsGRL Regions across which you want to aggregate.
#' @param excludeGR A GenomicRanges object with regions you want to 
#' exclude from the aggregation function. These regions will be eliminated
#' from the input table and not counted.
#' @param regionsGRL.length Vector with number of regions in each bin.
#' From bin1 to binN. With default NULL value, it will be auto assigned.
#' @param splitFactor Used to make "by string" to be plugged into a data.table
#' "by = " statemnt. With default NULL value, by string will be "list(regionID)"
#' @param keepCols Deprecated, NULL value should be used for MIRA aggregation.
#' @param sumCols Deprecated, NULL value should be used for MIRA aggregation.
#' @param jCommand You can pass a custom command in the j slot to data.table
#' specifying which columns to aggregate, and which functions to use. You
#' can use buildJ() to build a jCommand argument easily.
#' @param byRegionGroup You can aggregate by regionID or by regionGroupID; 
#' this reflects the regionsGRL that you pass; by default, BSAggregate will
#' aggregate each region individually -- scores will then be contiguous, and
#' the output is 1 row per region.
#' Turn on this flag to aggregate across all region groups, making the result
#' uncontiguous, and resulting in 1 row per *region group*.
#' @param keep.na Not used in general MIRA context.
#' 
#' @return In context of MIRA, with byRegionGroup = TRUE and jCommand = 
#' list( methyl = mean(methyl), readCount = sum(readCount) )", this function
#' will return a data.table with binCount rows (parameter for BSBinAggregate)
#' containing aggregated methylation from BSDT over binned regions from a region
#' set.
#'
BSAggregate = function(BSDT, regionsGRL, excludeGR = NULL, 
                       regionsGRL.length = NULL, splitFactor = NULL, 
                       keepCols = NULL, sumCols = NULL, jCommand = NULL, 
                       byRegionGroup = FALSE, keep.na = FALSE) {

    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if( "GRanges" %in% class(regionsGRL)) {
        regionsGRL = GRangesList(regionsGRL);
    } else if (! "GRangesList" %in% class(regionsGRL)) {
        stop("regionsGRL is not a GRanges or GRangesList object");
    }

    #make sure methyl column is present
    if(!"methyl" %in% colnames(BSDT)){
        stop("BSDT must have a methyl column.")
    }
    
    if(! is.null(excludeGR)) {
        BSDT = BSFilter(BSDT, minReads = 0, excludeGR)
    }

    bsgr = BSdtToGRanges(list(BSDT));
    
    colModes = sapply(BSDT, mode);
    if (is.null(sumCols)) {
        sumCols = setdiff(colnames(BSDT), c("chr", "start", "end", 
                                           "strand", splitFactor, keepCols))
        # Restrict to numeric columns.      
        sumCols = intersect(sumCols, 
                            names(colModes[which(colModes == "numeric")]))

    }
    # It's required to do a findoverlaps on each region individually, 
    # Not on a GRL, because of the way overlaps with GRLs work. So, 
    # we must convert the GRL to a GR, but we must keep track of which
    # regions came from which group.
    regionsGR = unlist(regionsGRL)
    
    if(is.null(regionsGRL.length)) {
        if (length(regionsGRL) > 100) {
            message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF = FALSE)
        }
        regionsGRL.length = sapply(regionsGRL, length)
        #message("Done counting regionsGRL lengths.");
    }

    # Build a table to keep track of which regions belong to which group
    region2group = data.table(
        regionID = 1:length(regionsGR), 
        chr = as.vector(seqnames(regionsGR)), 
        start = as.vector(start(regionsGR)), 
        end = as.vector(end(regionsGR)), 
        withinGroupID = as.vector(unlist(sapply(regionsGRL.length, seq))), 
        regionGroupID = rep(1:length(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)


    #message("Finding overlaps...");
    fo = findOverlaps(bsgr[[1]], regionsGR)

    setkey(BSDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))

    #message("Setting regionIDs...");
    BSDT = BSDT[queryHits(fo), ] #restrict the table to CpGs in any region.

    if (NROW(BSDT) < 1) {
        warning("No BSDT sites in the given region list. 
                Please expand your regionsGRL")
        return(NULL)
    }

    BSDT[, regionID := subjectHits(fo)] #record which region they overlapped.
    #if (!keep.na) {
    # BSDT = BSDT[queryHits(fo), ]
    #}

    if (is.null(jCommand)) {
        cols = c(sumCols, keepCols)
        funcs = c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
        jCommand = buildJ(cols, funcs)
    }
    #message("jCommand: ", jCommand)

    # Build the by string
    if (is.null(splitFactor)) {
        byString = paste0("list(regionID)");
    } else {
        byString = paste0("list(", paste("regionID", paste0(splitFactor, ""), 
                                         collapse = ", ", sep = ", "), ")")
    }

    # Now actually do the aggregate:
    #message("Combining...");
    bsCombined = BSDT[, eval(parse(text = jCommand)), by = eval(parse(text = byString))]
    setkey(bsCombined, regionID)
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups, 
    # which takes awhile. I think this preserves memory and is faster.

    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        # must set allow = TRUE here in case there are multiple IDs (splitCol)
        #adds regionGroupID column from region2group to bsCombined
        bsCombined[region2group, regionGroupID := regionGroupID, allow = TRUE]
        if (! is.null(splitFactor) ) { 
            byStringGroup = paste0("list(", 
                                   paste("regionGroupID", 
                                         paste0(splitFactor, collapse = ", "), 
                                         sep = ", "), 
                                   ")")
        } else {
            byStringGroup = "list(regionGroupID)"
        }
        #actual aggregation operation
        bsCombined = bsCombined[, eval(parse(text = jCommand)), 
                              by = eval(parse(text = byStringGroup))]
        
        #if any strand information was not given, averaging the signatures 
        #about the center to account for unknown strand orientation, 
        #also averaging readCount about center
        #ie if any "*" are present then average
        if ("*" %in% unique(as.character(strand(regionsGR)))){
            bsCombined[, methyl := (methyl + rev(methyl))/2]
            bsCombined[, readCount := (readCount + rev(readCount))/2]
        }
        
        return(bsCombined[]);
    } else {
    warning("Using byRegionGroup = FALSE may result in missing functionalities 
            such as symmetrical averaging")
        e = region2group[bsCombined, ]
        setkey(e, regionID);
        return(e);
    }
    # WARNING: There are now 2^2 ways to aggregate, sum vs mean
    # at each level: across regions, then across region sets. THis
    # doesn't give you a choice at this point. 
}


#' A function to plot the binned methylation of samples over a region.
#' 
#' @param binnedRegDT A datatable with specific column names containing:
#' bin numbers(binnedRegionDT column), 
#' aggregated methylation values (methyl column), 
#' name of the region set (featureID column), 
#' case/control column (sampleType column), 
#' sample name (sampleName column).
#' @param featID Region set names in a single string or vector of strings.
#' @param plotType Line or jitter (ggplot2). 
#' 
#' @return A plot of class "gg"/ "ggplot" that shows MIRA signatures
#' @examples
#' data("exampleBins", package = "MIRA")
#' MIRAplot = plotMIRARegions(binnedRegDT = exampleBins)
#' 
#' @export
plotMIRARegions <- function(binnedRegDT, 
                            featID = unique(binnedRegDT[, featureID]), 
                            plotType = "line"){
    
    setkey(binnedRegDT, featureID)
    binPlot = ggplot(data = binnedRegDT[featID], 
                   mapping = aes(x = regionGroupID, y = methyl))
    
    if (!("sampleType" %in% names(binnedRegDT))){
        sampleType = "All samples" 
        #if no sampleType column then all lines/points will be black
        warning("sampleType column must exist if it is desired to split up sample types by color")
    }
    if (plotType == "line"){
        binPlot = binPlot + geom_line(aes(col = sampleType, group = sampleName)) + 
            facet_wrap(~featureID)
    }else if (plotType == "jitter"){
        binPlot = binPlot + geom_jitter(aes(col = sampleType)) + facet_wrap(~featureID)
    }else {
        stop('The only supported values for plotType are "line" and "jitter"')
    }
    return(binPlot)
}

#' A function to plot MIRA scores and compare case/control.
#' 
#' @param scoreDT A datatable with the following columns: 
#' score, featureID (names of regions), sampleType.
#' @param featID Region set name/names in a single string or vector of strings.
#' @return a plot of class "gg" / "ggplot" that shows MIRA scores 
#' with geom_boxplot and geom_jitter.
#' @example
#' ##UPDATE #plotMIRAScores(scoreDT = EwingExample)
#' @export
plotMIRAScores <- function(scoreDT, featID = unique(scoreDT[, featureID])){
    setkey(scoreDT, featureID)
    scorePlot = ggplot(data = scoreDT[featID], mapping = aes(x = sampleType, y = score)) + 
              geom_boxplot() + geom_jitter() + facet_wrap(~featureID)  
    return(scorePlot)
}


#' Adding methyl column that has proportion of reads that were methylated for 
#' each site.
#' Note: Assigns methyl column by reference with " := "
#' 
#' @param BSDTList A bisulfite datatable or list of datatables with a column for
#' number of methylated reads (hitCount) and a column for number of total reads 
#' (readCount) for each cytosine that was measured.
#' @return The BSDTList but with extra methyl column on each data.table in list.
#' @export
#' @examples 
#' data("GM06990_1_ExampleSet", package = "MIRA") #exampleBSDT
#' exampleBSDT[, methyl := NULL] #removing methyl column
#' addMethCol(list(exampleBSDT))
addMethCol <- function(BSDTList){

    #converting to a data.table list if it was a single data.table
    if ("data.table" %in% class(BSDTList)){
        BSDTList = list(BSDTList)
    }

    #stopping the function if the input was not data.table originally
    if (!"data.table" %in% class(BSDTList[[1]])){
        stop('Input must be a single data.table object 
             or list of data.table objects')
    }

    #using anonymous function to apply operation that adds methyl column to each 
    #element of list
    #extra [] on the end is necessary for proper display/printing of the object
    BSDTList = lapply(X = BSDTList, 
                    FUN = function(x) x[, methyl := round(hitCount/readCount, 3)][])

    return(BSDTList)
}

#UPDATE:should normalizeMIRA really be exported? might need to be rehauled first
#' Function to normalize case/experimental samples to the controls
#' 
#' It finds the median of the controls for each bin for each region set
#' then divides by the median (oldBinVal/medianBinVal) for each bin.
#' @param binnedDT A datatable containing bins for each region set for each 
#' sample;bins contain aggregated methylation across regions for that sample; 
#' it should have a column with sample annotation 
#' so cases and controls can be split up: sampleType column (case/control)
#' as well as annotation of the region sets (featureID column). 
#' @return binnedDT The input DT but with normalized "methyl" values
#' 
#'
normalizeMIRA = function(binnedDT){
    if ("list" %in% class(binnedDT)){
        binnedDT = rbindlist(binnedDT)
    }
    if (!"data.table" %in% class(binnedDT)){
        stop("binnedDT must be a data.table (or list of data.tables)")
    }

    features = unique(binnedDT[, featureID]) #get a set of all features

    setkey(binnedDT, featureID, sampleType)

    #getting regionGroupIDs for the loop
    regGroupID = unique(binnedDT[, regionGroupID])
    #normalize for each feature
    for (i in 1:length(features)){ #get median methylation for each feature
        medMeth = binnedDT[.(features[i], "control"), median(methyl), 
                         by = regionGroupID]
    for (j in regGroupID){
        #divide each methylation values for each region group by 
        #the appropriate value
        #do operation on the rows for this feature and 
        #only one regionGroupID at a time
        binnedDT[featureID == features[i] & regionGroupID == j, 
                 methyl := methyl/medMeth[regionGroupID == j, V1] ]
    }

    } 
    return(binnedDT)
}

#' helper function
#' given a vector of columns, and the equally-sized vector of functions
#' to apply to those columns, constructs a j-expression for use in
#' a data.table 
#' (functions applied to columns in corresponding spot in "cols" string).
#' One function may be given to be applied to multiple columns.
#' use it in a DT[, eval(parse(text = buildJ(cols, funcs)))]
#' @param cols A string/vector of strings containing columns 
#' on which to use functions.
#' @param funcs Functions to use on columns.
#' @return A jcommand string. After performing function on column, column 
#' is reassigned the same name.
buildJ = function(cols, funcs) {
    r = paste("list(", paste(paste0(cols, " = ", funcs, "(", cols, ")"), collapse = ", "), ")")
    return(r);
}

# Internal part of a utility to convert data.tables into GRanges objects
# genes = dtToGR(gModels, "chr", "txStart", "txEnd", "strand", "geneId").
# 
# @param DT A data.table with at least "chr" and "start" columns
# @return gr A genomic ranges object derived from DT
dtToGrInternal = function(DT, chr, start, 
                          end = NA, strand = NA, name = NA, metaCols = NA) {
    
    if (is.na(end)) {
        if ("end" %in% colnames(DT)) {
            end = "end"
        } else {
            end = start;
        }
    }
    if (is.na(strand)){
        if ("strand" %in% colnames(DT)){ #checking if strand info is in DT
            strand = "strand"
        }
    }
    if (is.na(strand)) {
        gr = GRanges(seqnames = DT[[`chr`]], 
                   ranges = IRanges(start = DT[[`start`]], end = DT[[`end`]]), 
                   strand = "*")
    } else {
        # GRanges can only handle '*' for no strand, so replace any non-accepted
        # characters with '*'
        DT[, strand := as.character(strand)]
        DT[strand == "1", strand := "+"]
        DT[strand == "-1", strand := "-"]
        DT[[`strand`]] = gsub("[^+-]", "*", DT[[`strand`]])
        gr = GRanges(seqnames = DT[[`chr`]], 
                   ranges = IRanges(start = DT[[`start`]], end = DT[[`end`]]), 
                   strand = DT[[`strand`]])
    }
    if (! is.na(name) ) {
        names(gr) = DT[[`name`]];
    } else {
        names(gr) = 1:length(gr);
    }
    if(! is.na(metaCols)) {
        for(x in metaCols) {
            elementMetadata(gr)[[`x`]] = DT[[`x`]]
        }
    }
    gr;
}

# Convert a data.table to GRanges object.
# 
# @param DT a data.table with at least "chr" and "start" columns
# 
# @return gr A genomic ranges object derived from DT
dtToGr = function(DT, chr = "chr", start = "start", 
                  end = NA, strand = NA, name = NA, splitFactor = NA, metaCols = NA) {
    
    if(is.na(splitFactor)) {
        return(dtToGrInternal(DT, chr, start, end, strand, name, metaCols));
    }
    if ( length(splitFactor) == 1 ) { 
        if( splitFactor %in% colnames(DT) ) {
            splitFactor = DT[, get(splitFactor)];
        }
    }
    lapply(split(1:nrow(DT), splitFactor), 
        function(x) { 
            dtToGrInternal(DT[x, ], chr, start, end, strand, name, metaCols)
        }
    )
}

dtToGR = dtToGr;

#' Converts a list of data.tables into GRanges.
#' @param dtList A list of data.tables, 
#' Each should have "chr", "start", "hitCount", and "readCount" columns.
#' Error results if missing "chr", "start" but if hitCount and readCount are
#' missing, it will still work, just not have that info in the output.

#' @return a list of GRanges objects, strand has been set to "*", 
#' "start" and "end" have both been set to "start" of the DT.
#' hitCount and readCount info is preserved in GRanges object.
BSdtToGRanges = function(dtList) {
    
    gList = list();
    for (i in 1:length(dtList)) {
        #dt = dtList[[i]];
        setkey(dtList[[i]], chr, start)
        #convert the data into granges object
        gList[[i]] = GRanges(seqnames = dtList[[i]]$chr, 
                             ranges = IRanges(start = dtList[[i]]$start, 
                                            end = dtList[[i]]$start), 
                             strand = rep("*", nrow(dtList[[i]])), 
                             hitCount = dtList[[i]]$hitCount, 
                             readCount = dtList[[i]]$readCount)
        #I used to use end = start + 1, but this targets CG instead of just a C, and
        #it's causing edge-effects problems when I assign Cs to tiled windows 
        #using (within). Aug 2014 I'm changing to start/end at the same 
        #coordinate.
    }
    return(gList);
}

# This can run into memory problems if there are too many files...
# because of the way parallel lacks long vector support. The solution is
# to just use a single core; or to pass mc.preschedule = FALSE; This
# makes it so that each file is processed as a separate job. Much better.
#' Read in files from biseq meth caller
#' @param files a list of filenames (use parseInputArg if necessary)
#' @param contrastList  a list of named character vectors, 
#' each with length equal to the number of items in files. 
#' These will translate into column names in the final table.
#' @param sampleNames   a vector of length length(files), name for each file. 
#' You can also just use contrastList to implement the same thing 
#' so this is really unnecessary...
#' @param cores number of processors.
#' @param returnAsList Whether to return the output as a list 
#' or as one big data.table.
#' @return Data from each input file joined together into one big data.table.
#' If returnAsList = TRUE, then input from each file will be in its own data.table
#' in a list.
#'
#' @export
BSreadBiSeq = function(files, contrastList = NULL, 
                       sampleNames = extractSampleName(files), 
                       cores = 4, returnAsList = FALSE) {
    
    cores = min(length(files), cores); #not more cores than files!
    setLapplyAlias(cores);
    if (!is.null(contrastList)) {
        if( any(sapply(contrastList, length) ! = length(files))) {
            stop("contrastList must be a list, 
                 with each value having the same number of elements as files.");
        }
    }
    message("Reading ", length(files), " files..");
    freadList = lapplyAlias(files, fread, mc.preschedule = FALSE);
    colNames = names(contrastList)
    message("File reading finished (", length(files), " files). Parsing Biseq format...", appendLF = FALSE);
    # TODO: This parsing takes awhile, and could be done in parallel.
    freadListParsed = lapplyAlias(freadList, parseBiseq, mc.preschedule = FALSE)

    message("Parsing complete, building final tables and cleaning up...")
    numberOfFiles = length(freadListParsed);
    for (i in 1:numberOfFiles) {
        if (numberOfFiles > 1) {
            message(i, ": ", sampleNames[i], "; ", appendLF = FALSE)
        }
        if (numberOfFiles > 1 && i == numberOfFiles){
            message("", appendLF = TRUE)
        }
        DT = freadListParsed[[i]]; #convenience alias.
        if(!is.null(contrastList)) {
            DT[, get("colNames") := as.list(sapply(contrastList, "[[", i))]
        }
        if (!is.null(sampleNames)) {
            DT[, sampleName := sampleNames[i]]
        }
        freadListParsed[[i]] = DT
    }

    #filteredList = do.call(rbind, freadListParsed)
    #gc(); #rbind call is memory inefficient; this helps.
    # rbindlist supposedly does the same thing as do.call(rbind, list) but 
    # faster
    #default (returnAsList = FALSE) is to return as 
    #one combined data.table/data.frame
    if (!returnAsList){ 
        filteredList = rbindlist(freadListParsed)
    }else{
        filteredList = freadListParsed
    }

    return(filteredList);
}

#' Takes a data.table from BSreadBiSeq and parses the strange x/y format
#' of methylation calls, splitting them into individual columns
#' @param DT data.table to parse
#' "chr", "start", "end", "meth", "rate", "strand" columns expected
#' in that order.
#' @return data.table with separate methylated and unmethylated columns.
#' Specific col names are set
#' 
#' 
parseBiseq = function(DT) {
    message(".", appendLF = FALSE);
    setnames(DT, paste0("V", 1:6), 
             c("chr", "start", "end", "meth", "rate", "strand"))
    DT[, meth := gsub("'", "", meth)]
    #split the '12/12' format of meth calls
    ll = unlist(strsplit(DT$meth, "/", fixed = TRUE))
    idx = seq(1, length(ll), by = 2)
    DT[, ` := `(hitCount = as.integer(ll[idx]), 
              readCount = as.integer(ll[idx + 1]))]
    DT[, start := as.integer(start + 1)] #re-index
    DT[, c("rate", "end", "meth" ) := NULL] #remove unnecessary columns
    DT[, strand := NULL]
    DT = DT[, list(hitCount = sum(hitCount), readCount = sum(readCount)), 
          by = list(chr, start)] #smash measurements
    setcolorder(DT, c("chr", "start", "hitCount", "readCount"));
    DT = DT[ !grep("_", chr), ]; #clean Chrs
    return(DT)
}


#' convert a GenomicRanges into a data.table
#' 
#' @param GR A GRanges object
#' @param includeStrand Boolean, whether to include strand from GR in output DT
#' @return A data.table object with columns:
#' "chr", "start", and "end" (possibly strand)
grToDt = function(GR, includeStrand = FALSE) {
    DF = as.data.frame(elementMetadata(GR))
    if( ncol(DF) > 0) {
        if(includeStrand){
            DT = data.table(chr = as.vector(seqnames(GR)), 
                            start = start(GR), 
                            end = end(GR), 
                            strand = as.vector(strand(GR), mode = "character"), 
                            DF)    
        } else{
            DT = data.table(chr = as.vector(seqnames(GR)), 
                            start = start(GR), 
                            end = end(GR), 
                            DF)    
        }
    } else {
        if(includeStrand){
            DT = data.table(chr = as.vector(seqnames(GR)), 
                            start = start(GR), 
                            end = end(GR), 
                            strand = as.vector(strand(GR), mode = "character")) 
        } else{
            DT = data.table(chr = as.vector(seqnames(GR)), 
                            start = start(GR), 
                            end = end(GR))    
        }
    }
    return(DT)
}

#' This function is a drop-in replacement for the base list() function, 
#' which automatically names your list according to the names of the 
#' variables used to construct it.
#' It seemlessly handles lists with some names and others absent, 
#' not overwriting specified names while naming any unnamed parameters.
#'
#' @param ...   arguments passed to list()
#' @return A named list object.
#' @export
#' @examples
#' x = 5
#' y = 10
#' nlist(x, y) # returns list(x = 5, y = 10)
#' list(x, y) # returns unnamed list(5, 10)
nlist = function(...) {
    fcall = match.call(expand.dots = FALSE)
    l = list(...);
    if(!is.null(names(list(...)))) { 
        names(l)[names(l) == ""] = fcall[[2]][names(l) == ""]
    } else {  
        names(l) = fcall[[2]];
    }
    return(l)
}

#' To make parallel processing a possibility but not required, 
#' I use an lapply alias which can point at either the base lapply
#' (for no multicore), or it can point to mclapply, 
#' and set the options for the number of cores (what mclapply uses).
#' With no argument given, returns intead the number of cpus currently selected.
#'
#' @param cores Number of cpus
#' @return None
setLapplyAlias = function(cores = 0) {
    if (cores < 1) {
        return(getOption("mc.cores"))
    }
    if(cores > 1) { #use multicore?
        if (requireNamespace("parallel", quietly = TRUE)) {
            options(mc.cores = cores)
        } else {
            warning("You don't have package parallel installed. Setting cores to 1.")
            options(mc.cores = 1) #reset cores option.
        }
    } else {
        options(mc.cores = 1) #reset cores option.
    }
}

#' Function to run lapply or mclapply, depending on the option set in
#' getOption("mc.cores"), which can be set with setLapplyAlias().
#'
#' @param ... Arguments passed lapply() or mclapply()
#' @param mc.preschedule Argument passed to mclapply
#' @return Result from lapply or parallel::mclapply
lapplyAlias = function(..., mc.preschedule = TRUE) {
    if (is.null(getOption("mc.cores"))) { setLapplyAlias(1) }
    if(getOption("mc.cores") > 1) {
        return(parallel::mclapply(..., mc.preschedule = mc.preschedule))
    } else {
        return(lapply(...))
    }
}

#' Extract sample names from file names as the first part of the file name 
#' (before any suffix).
#' @param fileNames A string/vector containing names of files.
#' fileNames can contain the full path. The actual file name should be the 
#' sample name followed by a character (suffixSep) that is not in 
#' the sample name then whatever else part of the file name 
#' (like file extension).
#' @param suffixSep The character in between the sample name 
#' and the rest of the file name.
#' @param pathSep The character in between directories in the file path.
#' @return The sample name portion of a file name as a string/vector of strings.
extractSampleName = function(fileNames, suffixSep = "\\.", pathSep = "/") {
    sapply(strsplit(fileNames, pathSep), 
           function(x) strsplit(rev(x)[1], suffixSep)[[1]][1])
}

#' Given a BSDT (bisulfite data.table), remove any entries that overlap
#' regions given in the excludeGR argument and/or filter out sites
#' that have lower than a minimum number of reads.
#' @param BSDT Bisulfite data.table to filter
#' @param minReads Require at least this level of coverage at a cpg.
#' @param excludeGR GRanges object with regions to filter.
#' 
#' @return The BSDT with appropriate regions removed.
BSFilter = function(BSDT, minReads = 10, excludeGR = NULL) {
    # First, filter for minimum reads.
    if (minReads > 0) {
        BSDT = BSDT[readCount >= minReads, ]
    }
    if (NROW(BSDT) == 0) { return(data.table(NULL)) }
    # Now, filter entries overlapping a region in excludeGR.
    if ( !is.null(excludeGR) ) {
        gr = dtToGr(BSDT)
        fo = findOverlaps(gr, excludeGR)
        qh = unique(queryHits(fo))
        length(qh)
        nrow(BSDT)
        BSDT = BSDT[-qh, ]
    }
    return(BSDT)
}

#check whether object is smoothed
#check for names in phenoData
#fix names in data.table
#add methyl column?, use addMethCol or bsseq built in getMeth()?
# @param bsseqObj An object of class bsseq
# @return MIRAFormatBSDTList A list of data.tables in MIRA format
# One data table for each sample column of the bsseq object.
bsseqToMIRA <-function(bsseqObj){
    # if (hasBeenSmoothed(bsseqObj)){
    #   warning("Raw (not smoothed) methylation and coverage values are being used.")
    # }
    MIRAFormatBSDTList = list() #to store output
    #obtaining coordinates as GRanges obj. and changing to data.table
    coordinates = grToDt(granges(bsseqObj))
    for (i in 1:ncol(bsseqObj)){ #each column is a different sample
        hitCount = getBSseq(BSseq = bsseqObj[, i], type = "M")
        readCount = getBSseq(BSseq = bsseqObj[, i], type = "Cov")
        #index for taking out rows with 0 coverage
        notCovered = which(readCount == 0)
        warning("Taking out rows with no coverage. Genomic coordinates may not have identical row numbers in different samples now.")
        MIRAFormatBSDTList[[i]] = data.table(chr = coordinates[, chr], 
                                           start = coordinates[, start], 
                                           hitCount = hitCount, 
                                           readCount = readCount
                                           )[!notCovered]#filtering
        setnames(MIRAFormatBSDTList[[i]], 
                 c("chr", "start", "hitCount", "readCount"))
    }
    #names for list (by reference)
    setattr(MIRAFormatBSDTList, "names", Biobase::sampleNames(bsseqObj))
    
    return(MIRAFormatBSDTList)
}
