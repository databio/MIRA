# PACKAGE DOCUMENTATION
#' Methylation-based Inference of Regulatory Activity (MIRA)
#' 
#' MIRA is a score that infers regulatory activity of genomic elements
#' based on DNA methylation data. It assess the degree of dip in methylation
#' level surrounding a regulatory site of interest, such as 
#' transcription factor binding sites.
#' This package provides functions for aggregating methylation 
#' data across region sets, in bins.
#'
#' @docType package
#' @name MIRA
#' @author Nathan Sheffield
#' @author John Lawson
#'
#' @references \url{http://github.com/databio}
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand
#'             seqnames granges
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_jitter geom_line
#'             theme_classic xlab ylab geom_hline ylim scale_color_discrete
#'             scale_x_discrete scale_fill_brewer scale_color_manual
#'             scale_color_brewer
#' @import BiocGenerics S4Vectors IRanges
#' @importFrom data.table ":=" setDT data.table setkey fread setnames 
#'             setcolorder rbindlist setattr setorder copy
#' @importFrom Biobase sampleNames
#' @importFrom stats lm coefficients poly
#' @importFrom bsseq getCoverage getMeth
#' @importFrom methods is
NULL


# Because of some issues with CRAN submission, 
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation, 
# in order to pass some R CMD check NOTES.
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
    ".", "bin", "binID", "chr", "element_blank", "featureID", 
    "geom_violin", "methylCount", "id", "meth", 
    "methylProp", "coverage", "regionGroupID", "regionID", 
    "sampleName", "sampleType", "theme", "ubinID", "V1"))
}

##########################################################################

#' Aggregate methylation data to get
#' a summary methylation profile for each region set
#' 
#' The main function for aggregating methylation data in MIRA analysis. 
#' Aggregates methylation across all regions in a given region set to
#' give a summary methylation profile for each region set. 
#' 
#' 
#' Each region is split into bins. For a given set of regions, 
#' methylation is first aggregated (averaged) within each 
#' bin in each region. Then methylation from corresponding bins from each
#' region are aggregated (averaged) across all regions 
#' (all first bins together, all second bins together, etc.), 
#' giving a summary methylation profile.
#' This process is done for each region set.
#'
#' @param BSDT A single data.table that has DNA methylation data on individual 
#' sites. Alternatively a BSseq object is allowed which will be converted 
#' internally to data.tables. The data.table input should have columns:
#' "chr" for chromosome, "start" for 
#' cytosine coordinate, "methylProp" for proportion of 
#' methylation (0 to 1), optionally "methylCount" 
#' for number of methylated reads, and
#' optionally "coverage" for total number of reads.
#' In addition, a "sampleName" column is strongly preferred (and required later
#' for scoring multiple samples at the same time using the 
#' "scoreDip" function in a MIRA workflow).

#' @param GRList A GRangesList object containing region sets, each set 
#' corresponding to a type of regulatory element. 
#' Each regionSet in the list should
#' be named. A named list of data.tables also works. 
#' @param binNum How many bins each region should be split into for aggregation 
#' of the DNA methylation data.
#' @param minBaseCovPerBin Filter out bins with fewer than minBaseCovPerBin reads. Only
#' used if there is a "coverage" column
#' 
#' @return a data.table with binNum rows for each region set containing
#' aggregated methylation data. If the input was a BSseq object
#' with multiple samples, a list of data.tables will be returned with
#' one data.table for each sample.
#' Each region was split into bins; methylation was put in these bins; 
#' Output contains sum of the all corresponding bins for the regions of each 
#' region set, ie for all regions in each region set: first bins summed, second 
#' bins summed, etc. Columns of the output should be "bin", "methylProp", 
#' "coverage" (if coverage was an input column), "featureID", 
#' and possibly "sampleName".
#' For information on symmetry of bins and output when a region set has
#' strand info, see ?BSBinAggregate.
#' 
#' @export
#' @examples
#' data("exampleBSDT", package = "MIRA")
#' data("exampleRegionSet", package = "MIRA")
#' exBinDT <- aggregateMethyl(exampleBSDT, exampleRegionSet)
aggregateMethyl <- function(BSDT, GRList, binNum = 11, minBaseCovPerBin = 500){
    
    if (is(BSDT, "BSseq")) {
        # if input is not a data.table but rather a BSseq object
        # bsseq objects can include multiple samples so make a list
        
        ## not adding sample name to BSDT because of memory it would require
        ## so sample names are no longer required in BSseq object
        # # checking for sample names
        # if (length(sampleNames(BSDT)) != ncol(BSDT)) {
        #     stop(cleanws("BSseq object must have sample name for each sample.
        #                 Check output of bsseq::sampleNames(BSDT)"))
        # }
        BSDTList <- bsseqToDataTable(BSDT)
        
        ## not adding sample name to BSDT because of memory it would require
        # # for each data.table, make a sampleName column with the appropriate
        # # name (by reference)
        # mapply(FUN = function(x, y) x[, sampleName := y], 
        #        BSDTList,
        #        names(BSDTList))
        # do the aggregation step on each data.table
        bigMethylByBin <- lapply(X = BSDTList, 
                                FUN = function(x) aggregateMethylInt(BSDT = x, 
                                                                     GRList = GRList, 
                                                                     binNum = binNum, 
                                                                     minBaseCovPerBin = minBaseCovPerBin))
    } else {
        bigMethylByBin <- aggregateMethylInt(BSDT = BSDT, 
                                            GRList = GRList, 
                                            binNum = binNum, 
                                            minBaseCovPerBin = minBaseCovPerBin)
    }
    
    

    return(bigMethylByBin)
}

aggregateMethylInt <- function(BSDT, GRList, binNum = 11, minBaseCovPerBin = 500) {
    ######### aggregateMethyl:Preprocessing and formatting###############
    # BSDT should not be a list but can be converted
    if (is(BSDT, "list")) {
        if (length(BSDT) == 1) {
            BSDT <- BSDT[[1]]
        } else {
            stop(cleanws("Only one BSDT may be given to function. 
                         BSDT should not be a list."))
        }
        }
    if (!is(BSDT, "data.table")) {
        stop("BSDT must be a data.table")
    }
    
    # converting to list format if GRList is a data.table or GRanges object
    if (is(GRList, "GRanges")) {
        GRList <- GRangesList(GRList)
        message("Converting to GRangesList...")
    }
    if (is(GRList, "data.table")) {
        GRList <- list(GRList)
        message("Converting to list...")
    }
    
    # checking that input is in list format
    # should be either list or "GRangesList"
    if (!(is(GRList, "list") || is(GRList, "GRangesList"))) {
        stop("GRList should be a named list/GRangesList.")
    }
    
    # checking if region sets have names
    if (is.null(names(GRList))) {
        warning(cleanws("GRList should be a named list/GRangesList. 
                        The region sets were assigned sequential names 
                        based on their order in the list."))
        names(GRList) <- paste0(rep("RegionSet", length(GRList)), 
                                seq_along(GRList))
    }
    
    # checking that all objects in GRList are the same type 
    # and converting to data.tables
    # first check if all objects are GRanges
    if (all(sapply(X = GRList, FUN = function(x) is(x, "GRanges")))) {
        # GRanges to data.tables
        GRDTList <- lapply(X = GRList, FUN = grToDt, includeStrand = TRUE)
        
        # next check if all objects are data.tables 
    }else if (all(sapply(X = GRList, FUN = function(x) is(x, "data.table")))) {
        
        GRDTList <- GRList # this case is okay
    }else{
        stop("GRList should be a GRangesList or a list of data.tables")
    }
    
    
    # adding a methylProp column if it is not already in the BSDT
    if (!("methylProp" %in% names(BSDT))) {
        BSDTList <- addMethPropCol(list(BSDT))
        BSDT <- BSDTList[[1]] 
    }
    
    # checking for a coverage column in BSDT
    hasCoverage <- "coverage" %in% colnames(BSDT)
    
    
    ######## aggregateMethyl:Binning and processing output####################
    
    
    methylByBin <- lapply(X = GRDTList, 
                         FUN = function(x) BSBinAggregate(BSDT = BSDT, 
                                                          rangeDT = x, 
                                                          binNum = binNum, 
                                                          splitFactor = NULL, 
                                                          minBaseCovPerBin = minBaseCovPerBin,
                                                          hasCoverage = hasCoverage))
    names(methylByBin) <- names(GRList)# preserving names
    # adding a feature ID column to each data.table that 
    # should identify what region set was used
    for (i in seq_along(methylByBin)) {
        methylByBin[[i]][, featureID := rep(names(methylByBin)[i], 
                                            nrow(methylByBin[[i]]))][]
    }
    # screening out region sets that had incomplete binning
    binNumScreen <- sapply(X = methylByBin, FUN = nrow)
    # taking out incomplete region sets
    methylByBin <- methylByBin[!(binNumScreen < binNum)]
    
    bigMethylByBin <- rbindlist(methylByBin)
    ## letting user take care of sampleName stuff outside aggregateMethyl()
    # sampleNameInBSDT <- "sampleName" %in% colnames(BSDT)
    # if (sampleNameInBSDT && (ncol(bigMethylByBin) != 0)) {
    #     # creating new sampleName column
    #     bigMethylByBin[, sampleName := rep(BSDT[1, sampleName])][] 
    # }
    
    return(bigMethylByBin)
}


#' Get MIRA score for each sample/region set combination
#' 
#' Takes methylation data and sets of regions then aggregates
#' methylation for each region set and scores the resulting profile.
#' A wrapper for aggregateMethyl and scoreDip but it does not return
#' the summary methylation profiles, just the scores. This function is 
#' given only for convenience in working with small numbers of samples/region
#' sets. For large analyses, "aggregateMethyl" and "scoreDip" are recommended.
#' See vignettes for recommended use of MIRA.
#'
#' @param BSDT A single data table that has DNA methylation data on individual 
#' sites. #' Alternatively a BSseq object may be input which will be converted 
#' internally to data.table/s (sample names must be in BSseq object). 
#' For the data.table input, it should 
#' include a "chr" column with chromosome, a "start" column with the 
#' coordinate number for the cytosine, a "methylProp" column with proportion of 
#' methylation (0 to 1), optionally a "methylCount" column with 
#' number of methylated reads 
#' for each site, optionally a "coverage" column with 
#' total number of reads for each 
#' site, and a "sampleName" column with a sample identifier/name (required). 

#' @param GRList A GRangesList object containing region sets, each set 
#' corresponding to a regulatory element (or having regions with the 
#' same biological annotation).
#' Each regionSet in the list should be named. 
#' @param binNum How many bins each region should be split into for aggregation 
#' of the DNA methylation data.
#' @param scoringMethod Method to calculate MIRA score after binning. 
#' "logRatio" is currently the only option. See scoreDip function.
#' @param minBaseCovPerBin Filter out bins with fewer than minBaseCovPerBin reads. Only used
#' if there is a 'coverage' column in BSDT
#' 
#' @return A data.table with a MIRA score for each region set in GRList. 
#' See ?scoreDip. 
#' If input for "BSDT" is a BSseq object, output will be a list of 
#' data.tables if there were multiple samples in the BSseq object. 
#' @examples 
#' data("exampleBSDT", package = "MIRA") 
#' data("exampleRegionSet", package = "MIRA") 
#' MIRAScore(BSDT = exampleBSDT, GRList = exampleRegionSet)
#' 
#' @export
MIRAScore <- function(BSDT, GRList, binNum = 11, scoringMethod = "logRatio", 
                     minBaseCovPerBin = 500){

    ## not requiring sampleName column to save memory
    # # checking for sampleName column
    # if (is(BSDT, "data.table")) {
    #     if (!("sampleName" %in% colnames(BSDT))) {
    #         stop("sampleName column must be present in BSDT")
    #     }
    # } 
    
    # if BSseq object was given, convert to a list of data.tables
    # then run aggregation on each with lapply
    if (is(BSDT, "BSseq")) {
        ## not requiring sampleName to save on memory
        # # checking for sample names
        # if (length(sampleNames(BSDT)) != ncol(BSDT)) {
        #     stop(cleanws("BSseq object must have sample name for each sample.
        #                 Check output of bsseq::sampleNames(BSDT)"))
        # }
        BSDTList <- bsseqToDataTable(BSDT)
        ## not adding sampleName column in order to save on memory
        # # for each data.table, make a sampleName column with the appropriate
        # # name (by reference)
        # mapply(FUN = function(x, y) x[, sampleName := y], 
        #        BSDTList,
        #        names(BSDTList))
        sampleNames = names(BSDTList) # could be NULL
        bigBinList <- lapply(X = BSDTList, 
               FUN = function(x) aggregateMethylInt(BSDT = x, 
                                                    GRList = GRList, 
                                                    binNum = binNum, 
                                                    minBaseCovPerBin = minBaseCovPerBin))
        
        # using binned methylation data to calculate MIRA score
        scoreDT <- lapply(X = bigBinList, FUN = function(x) x[, .(score = scoreDip(methylProp, 
                                                                        method = scoringMethod)), 
                                                   by = .(featureID)])
        if (!is.null(sampleNames)) {
            setnames(scoreDT, sampleNames)
        }
    } else {
        
        bigBin <- aggregateMethyl(BSDT = BSDT, GRList = GRList, binNum = binNum, 
                             minBaseCovPerBin = minBaseCovPerBin)
        # using binned methylation data to calculate MIRA score
        scoreDT <- bigBin[, .(score = scoreDip(methylProp, 
                                              method = scoringMethod)), 
                         by = .(featureID)]
        
    }
  
    
    return(scoreDT)
}

#' Score methylation profile based on its shape
#' 
#' The dip scoring function for MIRA scores. This will take a data.table
#' that has the methylation level in each bin in the MIRA profile and
#' return a single score summarizing how large the 'dip' 
#' in methylation is at the center of that methylation profile.
#' A column for sample ID/name and a column for region set ID/name
#' should be included in the data.table because a separate score will be given 
#' for each sample/region set combination.
#'  
#' See `method` parameter for details on scoring calculations.
#' 
#' @param binnedDT A data.table with columns for:
#' bin ("bin"), methylation level ("methylProp"), 
#' region set ID/name (default expected column name is "featureID" 
#' but this is configurable via a parameter), sample name (default expected column name is 
#' "sampleName" but this is configurable via a parameter). 
#' The bin column is not used for calculations since it is assumed by
#' the function that the rows will be in the order of the bins (so the
#' function will work without a bin column although the bin column assists
#' in human readability of the input data.table) 
#' @param shoulderShift Used to determine the number of bins away from the 
#' center to use as the shoulders. Default value "auto" optimizes the 
#' shoulderShift variable for each sample/region set combination to try find the 
#' outside edges of the dip. shoulderShift may be manually set as an integer
#' that will be used for all sample/region set combinations. "auto" does not 
#' currently work with region sets that include strand info.
#' @param method The scoring method. "logRatio" is the log of the ratio of outer
#' edges to the middle. This ratio is the average of outside values 
#' of the dip (shoulders) divided by the center value if 
#' it is lower than the two surrounding values (lower for concave up profiles or
#' higher for concave down profiles) or if it is not lower (higher for
#' concave down profiles), an 
#' average of the three middle values. For an even binNum, the middle four
#' values would be averaged with the 1st and 4th being weighted by half (as
#' if there were 3 values). 
#' A higher score with "logRatio" corresponds to a deeper dip. "logRatio" is the
#' only scoring method currently but more methods may be added in the future.
#' @param usedStrand If strand information is included as part of an
#' input region set when aggregating methylation, 
#' the MIRA signature will probably not be 
#' symmetrical. In this case, the automatic 
#' shoulderShift sensing (done when shoulderShift="auto") needs to
#' be done for both sides of the dip instead of just one side so set
#' usedStrand=TRUE if strand was included for a region set.
#' usedStrand=TRUE only has an effect on the function when shoulderShift="auto".
#' @param regionSetIDColName A character object. The name of the column
#' that has region set names/identifiers.
#' @param sampleIDColName A character object. The name of 
#' the column that has sample names/identifiers.
#' 
#' @return A data.table with a column for region set ID 
#' (default name is featureID), sample ID (default name is sampleName), 
#' and MIRA score (with name "score"). There will
#' be one row and MIRA score for each sample/region set combination.
#' The MIRA score quantifies the "dip" of 
#' the MIRA signature which is an aggregation of methylation 
#' over all regions in a region set. 
#' 
#' @export
#' @examples
#' data("exampleBins")
#' scoreDip(exampleBins)
#' 
scoreDip <- function(binnedDT, 
                    shoulderShift = "auto", 
                    method = "logRatio",
                    usedStrand = FALSE,
                    regionSetIDColName = "featureID",
                    sampleIDColName = "sampleName"){
    
    if (is(binnedDT, "data.table")) {
    
        # the expected/necessary columns
        expectedCols <- c("methylProp", regionSetIDColName,
                         sampleIDColName)
        expectedColsPresent <- expectedCols %in% colnames(binnedDT)
        
        if (!all(expectedColsPresent)) {
            stop(paste("Missing the following expected columns:",
                       paste(expectedCols[!expectedColsPresent], collapse=" "), 
                       sep=" "))
        }
            
        scoreDT <- binnedDT[, .(score = scoreDipInt(methylProp, 
                                                   shoulderShift=shoulderShift,
                                                   method=method,
                                                   usedStrand=usedStrand)), 
                           by = .(get(regionSetIDColName), get(sampleIDColName))] 
        #fixing names (not assigned normally because of using "get")
        setnames(scoreDT, c(regionSetIDColName, sampleIDColName, "score"))
    } else if (is(binnedDT, "numeric")) {
        # if binnedDT is actually a vector as was original behaviour of function
        # preserving ability to use scoreDip on a single vector or in a 
        # data.table j expression
        # scoring by sample and region set will be accomplished outside
        # this function in the "by" part of the data.table expression
        # output is called scoreDT but is actually a single number object
        scoreDT <- scoreDipInt(values=binnedDT,
                    shoulderShift=shoulderShift,
                    method=method,
                    usedStrand=usedStrand)
    } else {
        # the preferred input is data.table
        stop("Input should be a data.table object")
    }
        
    return(scoreDT)
}


# This function was the original scoreDip but it was a wrapper was
# created to add some features and simplify the user interface.
# This will take a vector
# describing the methylation pattern and return a single score summarizing
# that vector for how large the 'dip' in methylation is at the center of
# the vector.
# @param values A vector with proportion of methylation values for each bin. 
# Between 0 and 1.
# @return A MIRA score. The MIRA score quantifies the "dip" of 
# the MIRA signature which is an aggregation of methylation 
# over all regions in a region set. 
# @examples
# data("exampleBins")
# exampleBins[, .(score = scoreDip(methylProp)), 
# 
scoreDipInt <- function(values, 
                    shoulderShift = "auto", 
                    method = "logRatio",
                    usedStrand = FALSE){
    
    if (!(method %in% "logRatio")) { # add new methods eventually
        stop("Invalid scoring method. Check spelling/capitalization.")
    }
    
    # determining number of bins
    binNum <- length(values)
    
    # determining whether signature is concave up or down in general
    # because values for finding shoulder need to be altered if concave
    # down ('logratio scoring') 
    # also it matters for getting middle value for 'logratio' scoring
    concaveUp <- isProfileConcaveUp(values, binNum)
    
    if (method == "logRatio") {
        centerSpot <- (binNum + 1) / 2 # X.5 for even binNum
        
        
        if ((binNum %% 2) == 0) { # if binNum is even, centerSpot is X.5
            # if one of middle 2 vals is lowest, use it, otherwise average
            # order of midVals vector matters because of which.min
            midVals=values[c(centerSpot - .5, centerSpot + .5, centerSpot - 1.5,
                             centerSpot + 1.5)]
            if (which.min(midVals) == 1) { # first middle bin
                midpoint <- midVals[1]
            } else if (which.min(midVals) == 2) { # second middle bin
                midpoint <- midVals[2]
            } else { # otherwise take average
                # includes 4 bins but outer two bins are weighted by half
                # approximates having 3 middle bins
                midpoint <- (.5 * values[centerSpot - 1.5] 
                            + values[centerSpot - .5]
                            + values[centerSpot + .5] 
                            + 0.5 * values[centerSpot + 1.5]) / 3    
            }
        }else {# if binNum is odd, centerSpot is X.0
            
            # if the middle is lowest for concave up/highest for concave down
            # use it, otherwise average will be taken
            if (concaveUp) {
                # order matters in which.min (for ties) so centerSpot is first
                useMiddlePoint <- (which.min(values[c(centerSpot, 
                                                      centerSpot - 1, 
                                                      centerSpot + 1)]) == 1)
            } else {
                # if concave down, see if middle point is the max of 
                # the three middle points
                useMiddlePoint <- (which.max(values[c(centerSpot, 
                                                      centerSpot - 1, 
                                                      centerSpot + 1)]) == 1)
            }
            
            if (useMiddlePoint) {
                midpoint <- values[centerSpot]
            } else { # if centerSpot does not have lowest/highest value of 
                # the three, average the three middle bins
                midpoint <- (values[centerSpot] + values[centerSpot + 1] 
                            + values[centerSpot - 1] ) / 3
            }
        }
        # automatically figuring out shoulderShift based on each signature
        if (shoulderShift == "auto") {
            
            
            
            if (concaveUp) {
                values2 <- values 
            } else {
                # converts to concave up so same `findShoulder` algorithm 
                # can be used for both cases
                # values2 is not used for scoring, just for finding shoulders
                values2 <- 1 - values 
            }
            
            
            
            # signature will probably not be symmetrical if strand was used
            if (usedStrand) { # probably not common but still an option
                shoulderShiftL <- findShoulder(values2, binNum, centerSpot,
                                             whichSide = "left")
                shoulderShiftR <- findShoulder(values2, binNum, centerSpot, 
                                             whichSide = "right")
            } else { # most common use case, strand was not used
                # either side would work since signature is symmetrical
                shoulderShift <- findShoulder(values2, binNum, centerSpot, 
                                             whichSide = "right")
            }
        }
        
        
        
        # floor and ceiling are only relevant when binNum is even 
        #(which means centerSpot is X.5)
        if (usedStrand && (shoulderShift == "auto")) { # probably uncommon case
            # floor and ceiling should not matter when "auto" is used but
            # I kept them just in case
            leftSide <- floor(centerSpot - shoulderShiftL)  
            rightSide <- ceiling(centerSpot + shoulderShiftR)
        } else { # most common case, strand was not used,
            #"auto" may or may not have been used
            leftSide <- floor(centerSpot - shoulderShift)  
            rightSide <- ceiling(centerSpot + shoulderShift)
        }
        shoulders <- ((values[leftSide] + values[rightSide]) / 2)
        if (midpoint < .000001) {
            warning("Division by zero. Consider adding a small constant 
                    to all bins for this region set.")
        }
        if (shoulders < .000001) {
            warning("Taking log of zero. Consider adding a small constant 
                    to all bins for this region set.")
        }
        # log ratio...
        score <- log(shoulders / midpoint)
    }

    # # alternate way of scoring by the area in the dip
    # if (method == "area") {
    #     maxMethyl <- max(values)
    #     score <- maxMethyl * binNum - sum(values)
    # }

    # # another alternate method
    # if (method == "parabola") {
    #     # fit2 <- lm(y~poly(x, 2, raw = TRUE))
    #     # lines(xx, predict(fit2, data.frame(x = xx)), col = "green")
    # }
    return(score)
}

# helper function for determining concavity of MIRA signature
# fits a parabola to the values then looks at the x^2 coefficient to 
# determine concavity (+ is concave up, - is concave down)
# if number of bins is high enough, a wide parabola is fit using all the data
# and then a narrow parabola is fit using the middle 50% or so values. The
# fit with the best adjusted R^2 value is used. This could help when the
# dip is mostly in the middle half but is somewhat arbitrary and 
# there might be a better way to do this.
# @return TRUE if concave up, FALSE if concave down
isProfileConcaveUp <- function(values, binNum) {
    wideFit <- lm(values ~ poly(seq_along(values), 2))
    wideX2Coef <- coefficients(wideFit)[3]
    wideR2 <- summary(wideFit)$adj.r.squared
    # if x2 coef is positive, it is concave up
    concaveUp <- (wideX2Coef >= 0) # up is preferred so inclusive of 0
    
    # only do a narrow fit if enough bins are used
    # want at least 7 bins to be a part of the narrow fit (smallerHalf)
    # floor(15 / 2) = 7
    # because 3 bins can be used for center in scoring (with odd binNum) 
    # and we want at least 2 bins on each side of those center bins
    if (binNum >= 15) {
        # using middle 50% or so for this fit
        smallerHalf <- floor(binNum / 2)
        endFirstQuarter <- ceiling((binNum - smallerHalf) / 2)
        narrowInd <- c(endFirstQuarter + seq(smallerHalf)) # middle ~half
        
        narrowFit <- lm(values[narrowInd] ~ 
                           poly(seq_along(values[narrowInd]), 2))
        narrowX2Coef <- coefficients(narrowFit)[3] # the x^2 coefficient
        narrowR2 <- summary(narrowFit)$adj.r.squared
        
        # use narrow fit if it has a better R2
        if (narrowR2 > wideR2) {
            concaveUp <- (narrowX2Coef >= 0)
        }
        
    }
    
    return(concaveUp)
}

# helper function for automatic shoulder sensing
# for unsymmetrical signatures this function should be run twice, once
# with whichSide="right" and once with whichSide="left"
# only finds the right shoulder so to find the left shoulder, flip the 
# input values "values", then take length(values)-shoulderShift
# used in scoreDip
# @param whichSide the side of the signature to find the shoulder for
# @param for other params see ?scoreDip
# @return shoulderShift The distance from the center of the signature to the shoulder.
#         It may be X.0 (ie an integer) if centerSpot is an integer or 
#         X.5 if centerSpot was X.5

findShoulder <- function(values, binNum, centerSpot, whichSide="right"){
    if (whichSide == "left") {
        values <- rev(values)
    }
    
    # first value that's not part of midpoint calculations
    shoulderStart <- ceiling(centerSpot) + 2 
    shoulderSpot <- shoulderStart
    for (i in shoulderStart:(binNum - 2)) {
        if (values[i + 1] > values[i]) {
            shoulderSpot <- i + 1
        } else if (((values[i + 2] + values[i + 1]) / 2) > values[i]) {
            shoulderSpot <- i + 1
        } else {
            break()
        }
    }
    # testing the last/most outside point if appropriate
    if (shoulderSpot == (binNum - 1)) {
        if (values[shoulderSpot + 1] > values[shoulderSpot]) {
            shoulderSpot <- shoulderSpot + 1
        }
    }
    # should work for X.0 and X.5 values of centerSpot
    # assuming calculations for leftSide and rightSide vars stay the same
    shoulderShift <- shoulderSpot - centerSpot
    
    return(shoulderShift)
}




#' Add column for proportion of methylation
#' 
#' Adding methylProp column that has proportion of reads that were 
#' methylated for each site based on number of methylated reads
#' divided by total number of reads.
#' Note: Assigns methylProp column by reference with ":="
#' 
#' @param BSDTList A bisulfite datatable or list of datatables 
#' with a column for number of methylated reads (methylCount) and 
#' a column for number of total reads 
#' (coverage) for each cytosine that was measured.
#' @return The BSDTList but with extra `methylProp` column on each 
#' data.table in list.
#' @export
#' @examples 
#' data("exampleBSDT", package = "MIRA")
#' exampleBSDT[, methylProp := NULL] # removing methylProp column
#' addMethPropCol(list(exampleBSDT))
addMethPropCol <- function(BSDTList){

    # converting to a data.table list if it was a single data.table
    if (is(BSDTList, "data.table")) {
        BSDTList <- list(BSDTList)
    }

    # stopping the function if the input was not data.table originally
    if (!is(BSDTList[[1]], "data.table")) {
        stop('Input must be a single data.table object 
             or list of data.table objects')
    }

    # using anonymous function to apply operation 
    # that adds methylProp column to each element of list
    # extra [] on the end is necessary for proper display/printing of the object
    BSDTList <- lapply(X = BSDTList, 
                    FUN = function(x) x[, methylProp := 
                                            round(methylCount / coverage, 3)][])

    return(BSDTList)
}




#' Read in files from biseq meth caller
#' 
#' Parses the x/y format of methylation calls, 
#' splitting them into individual columns: "methylCount" column for
#' number of methylated reads for site and "coverage" column for total
#' number of reads covering that site. Input files should have the following
#' columns: "chr", "start", "end", "meth", "rate", "strand".
#' 
#' 
#' This can run into memory problems if there are too many files...
#' because of the way parallel lacks long vector support. The solution is
#' to just use a single core; or to pass mc.preschedule = FALSE; This
#' makes it so that each file is processed as a separate job. Much better.
#' 
#' @param files a list of filenames (use parseInputArg if necessary)
#' @param contrastList Generally not needed for MIRA. 
#' A list of named character vectors, 
#' each with length equal to the number of items in files. 
#' These will translate into column names in the final table.
#' @param sampleNames   a vector of length length(files), name for each file. 
#' @param cores number of processors.
#' @param returnAsList Whether to return the output as a list 
#' or as one big data.table.
#' @return Data from each input file joined together into one big data.table.
#' If returnAsList = TRUE, then input from each file will be 
#' in its own data.table in a list.
#' @examples
#' shortBSDTFile <- system.file("extdata", "shortRRBS.bed", package = "MIRA") 
#' shortBSDT <- BSreadBiSeq(shortBSDTFile)
#'
#' @export
BSreadBiSeq <- function(files, contrastList = NULL, 
                       sampleNames = tools::file_path_sans_ext(basename(files)), 
                       cores = 4, returnAsList = FALSE) {
    
    cores <- min(length(files), cores); # not more cores than files!
    setLapplyAlias(cores);
    if (!is.null(contrastList)) {
        if (any(sapply(contrastList, length) != length(files))) {
            stop("contrastList must be a list, 
                 with each value having the same number of elements as files.");
        }
    }
    message("Reading ", length(files), " files..");
    freadList <- lapplyAlias(files, fread, mc.preschedule = FALSE);
    colNames <- names(contrastList)
    message("File reading finished (",
        length(files),
        " files). Parsing Biseq format...",
        appendLF = FALSE);
    # TODO: This parsing takes awhile, and could be done in parallel.
    freadListParsed <- lapplyAlias(freadList, parseBiseq, mc.preschedule = FALSE)

    message("Parsing complete, building final tables and cleaning up...")
    numberOfFiles <- length(freadListParsed);
    for (i in 1:numberOfFiles) {
        if (numberOfFiles > 1) {
            message(i, ": ", sampleNames[i], "; ", appendLF = FALSE)
        }
        if (numberOfFiles > 1 && i == numberOfFiles) {
            message("", appendLF = TRUE)
        }
        DT <- freadListParsed[[i]]; # convenience alias.
        if (!is.null(contrastList)) {
            DT[, get("colNames") := as.list(sapply(contrastList, "[[", i))]
        }
        if (!is.null(sampleNames)) {
            DT[, sampleName := sampleNames[i]]
        }
        freadListParsed[[i]] <- DT
    }

    # filteredList <- do.call(rbind, freadListParsed)
    # gc(); # rbind call is memory inefficient; this helps.
    # rbindlist supposedly does the same thing as do.call(rbind, list) but 
    # faster
    # default (returnAsList = FALSE) is to return as 
    # one combined data.table/data.frame
    if (!returnAsList) { 
        filteredList <- rbindlist(freadListParsed)
    }else{
        filteredList <- freadListParsed
    }

    return(filteredList);
}

# Takes a data.table from BSreadBiSeq and parses the strange x/y format
# of methylation calls, splitting them into individual columns
# @param DT data.table to parse
# "chr", "start", "end", "meth", "rate", "strand" columns expected
# in that order.
# @return data.table with separate methylated and unmethylated columns.
# Specific col names are set
parseBiseq <- function(DT) {
    message(".", appendLF = FALSE);
    setnames(DT, paste0("V", 1:6), 
             c("chr", "start", "end", "meth", "rate", "strand"))
    DT[, meth := gsub("'", "", meth)]
    # split the '12/12' format of meth calls
    ll <- unlist(strsplit(DT$meth, "/", fixed = TRUE))
    idx <- seq(1, length(ll), by = 2)
    DT[, `:=` (methylCount = as.integer(ll[idx]), 
              coverage = as.integer(ll[idx + 1]))]
    DT[, start := as.integer(start + 1)] # re-index
    DT[, c("rate", "end", "meth" ) := NULL] # remove unnecessary columns
    DT[, strand := NULL]
    DT <- DT[, list(methylCount = sum(methylCount), coverage = sum(coverage)), 
          by = list(chr, start)] # smash measurements
    setcolorder(DT, c("chr", "start", "methylCount", "coverage"));
    DT <- DT[ !grep("_", chr), ]; # clean Chrs
    return(DT)
}



# Given a BSDT (bisulfite data.table), remove any entries that overlap
# regions given in the excludeGR argument and/or filter out sites
# that have lower than a minimum number of reads.
# @param BSDT Bisulfite data.table to filter
# @param minReads Require at least this level of coverage at a cpg.
# @param excludeGR GRanges object with regions to filter.
# 
# @return The BSDT with appropriate regions removed.
BSFilter <- function(BSDT, minReads = 10, excludeGR = NULL) {
    # First, filter for minimum reads.
    if (minReads > 0) {
        BSDT <- BSDT[coverage >= minReads, ]
    }
    if (NROW(BSDT) == 0) { return(data.table(NULL)) }
    # Now, filter entries overlapping a region in excludeGR.
    if (!is.null(excludeGR)) {
        gr <- dtToGr(BSDT)
        fo <- findOverlaps(gr, excludeGR)
        qh <- unique(queryHits(fo))
        length(qh)
        nrow(BSDT)
        BSDT <- BSDT[-qh, ]
    }
    return(BSDT)
}
