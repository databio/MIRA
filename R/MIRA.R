# PACKAGE DOCUMENTATION
#' Methylation-based Inference of Regulatory Activity (MIRA)
#' MIRA is a score that measures the degree of dip in methylation
#' level surrounding a regulatory site of interest, such as a 
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
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand seqnames 
#'             granges
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_jitter geom_line
#' @import BiocGenerics S4Vectors IRanges
#' @importFrom data.table ":=" setDT data.table setkey fread setnames 
#'             setcolorder rbindlist setattr setorder
#' @importFrom Biobase sampleNames
NULL


# Because of some issues with CRAN submission, 
#(see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation, 
# in order to pass some R CMD check NOTES.
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
    ".", "binID", "chr", "featureID", "hitCount", "id", "meth", 
    "methyl", "readCount", "regionGroupID", "regionID", 
    "sampleName", "sampleType", "ubinID", "V1"))
}

##########################################################################

#' Function to aggregate methylation data into bins 
#' over all regions in each region set.
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
#' For information on symmetry of bins and output when a region set has
#' strand info, see ?BSBinAggregate.
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
    if (class(GRList) %in% "GRanges") {
      GRList = GRangesList(GRList)
      message("Converting to GRangesList...")
    }
    if (class(GRList) %in% "data.table") {
        GRList = list(GRList)
        message("Converting to list...")
    }

    #checking that input is in list format
    if (!(class(GRList) %in% c("list", "GRangesList"))) {
        stop("GRList should be a named list/GRangesList.")
    }

    #checking if region sets have names
    if (is.null(names(GRList))) {
        warning("GRList should be a named list/GRangesList. 
                Sequential names given according to order in object.")
        names(GRList) <- paste0(rep("RegionSet", length(GRList)), 
                              1:length(GRList))
    }

    #checking that all objects in GRList are the same type 
    #and converting to data.tables
    if (all(sapply(X = GRList, FUN = class) %in% "GRanges")) {
        #GRanges to data.tables
        GRDTList = lapply(X = GRList, FUN = grToDt, includeStrand = TRUE)
        #below statement will be true if all objects in the list are of 
        #class data.table
        #necessary since data.tables also include data.frame as a class
    }else if (all(sapply(
        X = lapply(X = GRList, FUN = function(x) class(x) %in% "data.table"), 
        FUN = any))) {
        
        GRDTList = GRList #this case is okay
    }else{
        stop("GRList should be a GRangesList or a list of data.tables")
    }
    
    if (sampleNameInBSDT) {
        if (!("sampleName" %in% colnames(BSDT))) {
            stop("BSDT should have sampleName col if sampleNameInBSDT = TRUE")
        }
    }

    #adding a methyl column if it is not already in the BSDT
    if (!("methyl" %in% names(BSDT))) {
        BSDTList = addMethCol(list(BSDT))
        BSDT = BSDTList[[1]] 
    }

    ########returnMIRABins:Binning and processing output####################


    methylByBin = lapply(X = GRDTList, 
                       FUN = function(x) BSBinAggregate(BSDT = BSDT, 
                                                        rangeDT = x, 
                                                        binCount = binNum, 
                                                        splitFactor = NULL, 
                                                        minReads = minReads))
    names(methylByBin) = names(GRList)#preserving names
    #adding a feature ID column to each data.table that 
    #should identify what region set was used
    for (i in 1:length(methylByBin)) {
        methylByBin[[i]][, featureID := rep(names(methylByBin)[i], 
                                           nrow(methylByBin[[i]]))][]
    }
    #screening out region sets that had incomplete binning
    binNumScreen = sapply(X = methylByBin, FUN = nrow)
    #taking out incomplete region sets
    methylByBin = methylByBin[!(binNumScreen < binNum)]

    bigMethylByBin = rbindlist(methylByBin)
    if (sampleNameInBSDT) {
        #creating new sampleName column
        bigMethylByBin[, sampleName := rep(BSDT[1, sampleName])][] 
    }


    return(bigMethylByBin)
}

#' Function to take DNA methylation and region data and return a MIRA score;
#' a wrapper for returnMIRABins and scoreDip.
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
    if (!("methyl" %in% names(BSDT))) {
        stop("BSDT must have a methyl column with proportion of methylation. 
             addMethCol() will add this.")
    }

    MIRAresults = list()


    bigBin = returnMIRABins(BSDT = BSDT, GRList = GRList, binNum = binNum, 
                          sampleNameInBSDT = sampleNameInBSDT, 
                          minReads = minReads)
  
    #using binned methylation data to calculate MIRA score
    scoreDT = bigBin[, .(score = scoreDip(methyl, binNum, 
                                          method = scoringMethod)), 
                   by = .(featureID, sampleName)]

    return(scoreDT)
}

#' The dip scoring function for MIRA scores.
#' 
#' @param values A vector with proportion of methylation values for each bin. 
#'  Between 0 and 1.
#' @param binCount Number of bins, also length of "values" vector.
#' @param shoulderShift Used to determine the number of bins away from the 
#' center to use as the shoulders. Default value "auto" optimizes the 
#' shoulderShift variable for each sample/region set combination to try find the 
#' outside edges of the dip. shoulderShift may be manually set as an integer
#' that will be used for all sample/region set combinations. "auto" does not 
#' currently work with region sets that include strand info.
#' @param method The scoring method. "logRatio" is the log of the ratio of outer
#' edges to the middle. This ratio is the average of outside values 
#' of the dip (shoulders) divided by the center value if 
#' it is lower than the two surrounding values or if it is not lower, an 
#' average of the three middle values. For an even binCount, the middle four
#' values would be averaged with the 1st and 4th being weighted by half (as
#' if there were 3 values). 
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
#' exampleBins[, .(score = scoreDip(methyl, binCount)), 
#'               by = .(featureID, sampleName)]
scoreDip = function(values, binCount, 
                    shoulderShift = "auto", 
                    method = "logRatio"){
    
    if (!(method %in% "logRatio")) { #add new methods eventually
        stop("Invalid scoring method. Check spelling/capitalization.")
    }
    if (method == "logRatio") {
        centerSpot = (binCount + 1) / 2 # X.5 for even binCount
        
        if ((binCount %% 2) == 0) { #if binCount is even, centerSpot is X.5
            #if one of middle 2 vals is lowest, use it, otherwise average
            #order of midVals vector matters because of which.min
            midVals=values[c(centerSpot - .5, centerSpot + .5, centerSpot - 1.5,
                             centerSpot + 1.5)]
            if (which.min(midVals) == 1) { #first middle bin
                midpoint = midVals[1]
            } else if (which.min(midVals) == 2) { #second middle bin
                midpoint = midVals[2]
            } else { #otherwise take average
                #includes 4 bins but outer two bins are weighted by half
                #approximates having 3 middle bins
                midpoint = (.5 * values[centerSpot - 1.5] 
                            + values[centerSpot - .5]
                            + values[centerSpot + .5] 
                            + 0.5 * values[centerSpot + 1.5]) / 3    
            }
        }else {#if binCount is odd, centerSpot is X.0
            # if the middle is lowest use it, otherwise take average
            # order matters in which.min (for ties) so centerSpot is first
            if (which.min(values[c(centerSpot, centerSpot - 1, centerSpot + 1)])
                == 1) {
                midpoint = values[centerSpot]
            } else { #if centerSpot does not have lowest value of the three
                #average the three middle bins
                midpoint = (values[centerSpot] + values[centerSpot + 1] 
                            + values[centerSpot - 1] ) / 3
            }
        }
        #automatically figuring out shoulderShift based on each signature
        #only works for symmetric MIRA signatures
        if (shoulderShift == "auto") {
            
            #first value that's not part of midpoint calculations
            shoulderStart = ceiling(centerSpot) + 2 
            shoulderSpot = shoulderStart
            for (i in shoulderStart:(binCount - 2)) {
                if (values[i + 1] > values[i]) {
                    shoulderSpot = i + 1
                } else if (((values[i + 2] + values[i + 1]) / 2) > values[i]) {
                    shoulderSpot = i + 1
                } else {
                    break()
                }
            }
            #testing the last/most outside point if appropriate
            if (shoulderSpot == (binCount - 1)) {
                if (values[shoulderSpot + 1] > values[shoulderSpot]) {
                    shoulderSpot = shoulderSpot + 1
                }
            }
            #should work for X.0 and X.5 values of centerSpot
            #assuming calculations for leftSide and rightSide vars stay the same
            shoulderShift = shoulderSpot - centerSpot
        }
        
        
        
        #floor and ceiling are only relevant when binCount is even 
        #(which means centerSpot is X.5)
        leftSide = floor(centerSpot - shoulderShift)  
        rightSide = ceiling(centerSpot + shoulderShift)
        shoulders = ((values[leftSide] + values[rightSide]) / 2)
        if (midpoint < .000001) {
            warning("Division by zero. Consider adding a small constant 
                    to all bins for this region set.")
        }
        if (shoulders < .000001) {
            warning("Taking log of zero. Consider adding a small constant 
                    to all bins for this region set.")
        }
        # log ratio...
        score = log(shoulders / midpoint)
    }

    # #alternate way of scoring by the area in the dip
    # if (method == "area") {
    #     maxMethyl = max(values)
    #     score = maxMethyl * binCount - sum(values)
    # }

    # #another alternate method
    # if (method == "parabola") {
    #     #fit2 <- lm(y~poly(x, 2, raw = TRUE))
    #     #lines(xx, predict(fit2, data.frame(x = xx)), col = "green")
    # }
    return(score)
}



#' Adding methyl column that has proportion of reads that were methylated for 
#' each site.
#' Note: Assigns methyl column by reference with ":="
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
    if ("data.table" %in% class(BSDTList)) {
        BSDTList = list(BSDTList)
    }

    #stopping the function if the input was not data.table originally
    if (!("data.table" %in% class(BSDTList[[1]]))) {
        stop('Input must be a single data.table object 
             or list of data.table objects')
    }

    #using anonymous function to apply operation that adds methyl column to each 
    #element of list
    #extra [] on the end is necessary for proper display/printing of the object
    BSDTList = lapply(X = BSDTList, 
                    FUN = function(x) x[, methyl := 
                                            round(hitCount / readCount, 3)][])

    return(BSDTList)
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
#' If returnAsList = TRUE, then input from each file will be 
#' in its own data.table in a list.
#'
#' @export
BSreadBiSeq = function(files, contrastList = NULL, 
                       sampleNames = tools::file_path_sans_ext(basename(files)), 
                       cores = 4, returnAsList = FALSE) {
    
    cores = min(length(files), cores); #not more cores than files!
    setLapplyAlias(cores);
    if (!is.null(contrastList)) {
        if (any(sapply(contrastList, length) != length(files))) {
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
        if (numberOfFiles > 1 && i == numberOfFiles) {
            message("", appendLF = TRUE)
        }
        DT = freadListParsed[[i]]; #convenience alias.
        if (!is.null(contrastList)) {
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
    if (!returnAsList) { 
        filteredList = rbindlist(freadListParsed)
    }else{
        filteredList = freadListParsed
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
# 
# 
parseBiseq = function(DT) {
    message(".", appendLF = FALSE);
    setnames(DT, paste0("V", 1:6), 
             c("chr", "start", "end", "meth", "rate", "strand"))
    DT[, meth := gsub("'", "", meth)]
    #split the '12/12' format of meth calls
    ll = unlist(strsplit(DT$meth, "/", fixed = TRUE))
    idx = seq(1, length(ll), by = 2)
    DT[, `:=` (hitCount = as.integer(ll[idx]), 
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



# Given a BSDT (bisulfite data.table), remove any entries that overlap
# regions given in the excludeGR argument and/or filter out sites
# that have lower than a minimum number of reads.
# @param BSDT Bisulfite data.table to filter
# @param minReads Require at least this level of coverage at a cpg.
# @param excludeGR GRanges object with regions to filter.
# 
# @return The BSDT with appropriate regions removed.
BSFilter = function(BSDT, minReads = 10, excludeGR = NULL) {
    # First, filter for minimum reads.
    if (minReads > 0) {
        BSDT = BSDT[readCount >= minReads, ]
    }
    if (NROW(BSDT) == 0) { return(data.table(NULL)) }
    # Now, filter entries overlapping a region in excludeGR.
    if (!is.null(excludeGR)) {
        gr = dtToGr(BSDT)
        fo = findOverlaps(gr, excludeGR)
        qh = unique(queryHits(fo))
        length(qh)
        nrow(BSDT)
        BSDT = BSDT[-qh, ]
    }
    return(BSDT)
}
