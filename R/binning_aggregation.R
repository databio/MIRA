#####################################################################
#Binning and aggregation functions
#BSAggregate is primary function for this and workhorse of MIRA
#####################################################################



#' Bins regions and aggregates methylation data across the regions by bin
#' 
#' First bins regions and averages the proportion of methylation for
#' all methylation sites within each bin (ie  the methylation of 
#' all sites within region 1, bin 1 are averaged, then all sites 
#' within region 1, bin 2 are averaged, etc.)
#' Then aggregates methylation across all regions by bin by averaging
#' the proportion of methylation in each corresponding 
#' bin (ie all bin1's together, all bin2's together, etc.). 
#'
#' 
#' @param BSDT A single data table that has DNA methylation data 
#' on individual sites including a "chr" column with chromosome, 
#' a "start" column with the coordinate number for the cytosine, 
#' a "methylProp" column with proportion of methylation (0 to 1), 
#' optionally a "methylCount" column with number of 
#' methylated reads for each site, and 
#' optionally a "coverage" column with total number of reads for each site
#' (hasCoverage param).
#' @param rangeDT A data table with the sets of regions to be binned, 
#' with columns named "start", "end". Strand may also be given and will
#' affect the output. See "Value" section.
#' @param binNum Number of bins across the region.
#' @param byRegionGroup Default TRUE will aggregate methylation over 
#' corresponding bins for each region (all bin1's aggregated, all bin2's, etc).
#' byRegionGroup = FALSE is deprecated.
#' @param minBaseCovPerBin Filter out bins where the sum of coverage values is
#' less than X before returning.
#' @param splitFactor With default NULL, aggregation will be done 
#' separately/individually for each sample.
#' @param hasCoverage Default TRUE. Whether there is a coverage column
#' 
#' @return With splitFactor = NULL, it will return a data.table 
#' with binNum rows, 
#' containing aggregated methylation data over regions in region set "rangeDT".
#' Each region was split into bins; methylation was put in these bins; 
#' Output contains sum of the all corresponding bins 
#' for the regions of each region set ie for all regions in each region set: 
#' first bins summed, second bins summed, etc.
#' Columns of the output should be "bin", "methylProp", and, if
#' coverage was included as input col, "coverage"
#' ###########################################################################
#' Info about how strand of rangeDT affects output:
#' The MIRA signature will be symmetrical if no strand information is given for 
#' the regions (produced by averaging the signature with the reverse of the 
#' signature), because the orientation of the regions is arbitrary with respect 
#' to biological features (like a promoter for instance) that could be 
#' oriented directionally (e.g. 5' to 3'). If strand information is given, 
#' regions on the minus strand will be flipped before being aggregated 
#' with plus strand regions so the MIRA signature will be in 
#' 5' to 3' orientation.
#' ###########################################################################
#' @examples
#' data("exampleBSDT") # exampleBSDT
#' data("exampleRegionSet") # exampleRegionSet
#' exampleBSDT <- addMethPropCol(exampleBSDT)
#' aggregateBins <- BSBinAggregate(BSDT = exampleBSDT, 
#'                              rangeDT = exampleRegionSet, 
#'                              binNum = 11, splitFactor = NULL)
#' 
#' @export
BSBinAggregate <- function(BSDT, rangeDT, binNum, minBaseCovPerBin = 500, 
                          byRegionGroup = TRUE, splitFactor = NULL,
                          hasCoverage = TRUE) {
    
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
    
    
    # if given GRanges object, change to DT
    if (is(rangeDT, "GRanges")) {
        rangeDT <- grToDt(GR = rangeDT, includeStrand = TRUE)
    }
    
    if (!is(rangeDT, "data.table")) {
        stop("rangeDT must be a data.table")
    }
    seqnamesColName <- "seqnames"  # default column name
    if ("chr" %in% colnames(rangeDT)) {
        # message("seqnames column name set to: chr")
        seqnamesColName <- "chr"
    } else {
        # Got neither.
        stop("rangeDT must have a seqnames or chr column")
    }
    
    if (! ("strand" %in% colnames(rangeDT))) {
        rangeDT[, strand := "*"]
        warning("Since strand not given, set to '*' ")
    }
    
    ## if (!silent) {
    # message("Binning...")
    ##}
    binnedDT <- rangeDT[, binRegion(start, end, 
                                   binNum, get(seqnamesColName), strand)]
    binnedGR <- sapply(split(binnedDT, binnedDT$binID), dtToGr)
    # message("Aggregating...")
    if (hasCoverage) {
        # aggregate methylation (mean) and sum coverage values
        aggrJCommand <- buildJ(c("methylProp", "coverage"), 
                              c("mean", "sum"))
    } else {
        # if no coverage only aggregate methylation level
        aggrJCommand <- buildJ("methylProp", "mean")
    }
    binnedBSDT <- BSAggregate(BSDT = BSDT, 
                             regionsGRL = GRangesList(binnedGR), 
                             jCommand = aggrJCommand, 
                             byRegionGroup = byRegionGroup, 
                             splitFactor = splitFactor,
                             hasCoverage = hasCoverage)
    # If we aren't aggregating by bin, then don't restrict to min reads!
    if (byRegionGroup) {
        if (hasCoverage) {
            # only keep rows (bins) with coverage >= minBaseCovPerBin
            binnedBSDT <- binnedBSDT[coverage >= minBaseCovPerBin, ]
        }
        if (nrow(binnedBSDT) < binNum) {
            # telling user what sample failed if sample name is in BSDT
            if ("sampleName" %in% names(BSDT)) {
                thisSample <- BSDT[1, sampleName]
            } else {
                thisSample <- "this sample"
            }
            
            warning(paste0("Less than minBaseCovPerBin. Unable to return bins for ", 
                            thisSample, " for this region set."))
        }
    }
    
    
    return(binnedBSDT)
}

# BSaggregate: Aggregate a BSDT across regions or region groups, 
# for multiple samples at a time.
# This function is as BScombineByRegion, but can handle not only multiple
# samples in BSDT, but also simultaneously multiple region sets by passing
# a regionsGRL (GRangesList object).However currently code for symmetrical 
# averaging will cause it to only work with one region set (which may be
# split up into multiple GRanges in a GRangesList).
# you can use jCommand to do other functions.

# Given a bisulfite data table as input, with an identifier column for
# different samples; plus a GRanges objects with regions to aggregate.
#
# @param BSDT The bisulfite data.table (output from one of the parsing
# functions for methylation calls) that you wish to aggregate. It can
# be a combined table, with individual samples identified by column passed
# to splitFactor. "chr", "start", "methylCount" columns should be in BSDT.
#  "coverage", and "methylProp" columns are optional (hasCoverage
# must be set to FALSE if there is not coverage column)
# @param regionsGRL Regions across which you want to aggregate.
# @param excludeGR A GenomicRanges object with regions you want to 
# exclude from the aggregation function. These regions will be eliminated
# from the input table and not counted.
# @param regionsGRL.length Vector with number of regions in each bin.
# From bin1 to binN. With default NULL value, it will be auto assigned.
# @param splitFactor Used to make "by string" to be plugged into a data.table
# "by=" statemnt. With default NULL value, by string will be "list(regionID)"
# @param keepCols Deprecated, NULL value should be used for MIRA aggregation.
# @param sumCols Deprecated, NULL value should be used for MIRA aggregation.
# @param jCommand You can pass a custom command in the j slot to data.table
# specifying which columns to aggregate, and which functions to use. You
# can use buildJ() to build a jCommand argument easily. jCommand may be
# used twice so changing column names in the jCommand could introduce errors
# @param byRegionGroup You can aggregate by regionID or by regionGroupID; 
# this reflects the regionsGRL that you pass; by default, BSAggregate will
# aggregate each region individually -- scores will then be contiguous, and
# the output is 1 row per region.
# Turn on this flag to aggregate across all region groups, making the result
# uncontiguous, and resulting in 1 row per *region group*. 
# (byRegionGroup=TRUE is used for normal MIRA behaviour that 
# aggregates across regions by bins 
# -- all bin1's together, all bin2's together etc.--, in which case
# "region group" in the description above refers to bin number,
# and you would have 1 row per bin)
# @param keep.na Not used in general MIRA context.
# @param hasCoverage Default TRUE. Assuming there is a coverage column
# unless told otherwise
# 
# @return In context of MIRA, with byRegionGroup = TRUE and jCommand = 
# list( methylProp = mean(methylProp), coverage = sum(coverage) )", 
# this function will return a data.table with 
# binNum rows (parameter for BSBinAggregate)
# containing aggregated methylation from BSDT over binned regions from a region
# set.
#
BSAggregate <- function(BSDT, regionsGRL, excludeGR = NULL, 
                       regionsGRL.length = NULL, splitFactor = NULL, 
                       keepCols = NULL, sumCols = NULL, jCommand = NULL, 
                       byRegionGroup = FALSE, keep.na = FALSE, hasCoverage = TRUE) {
    
    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if (is(regionsGRL, "GRanges")) {
        regionsGRL <- GRangesList(regionsGRL);
    } else if (!is(regionsGRL, "GRangesList")) {
        stop("regionsGRL is not a GRanges or GRangesList object");
    }
    
    # make sure methylProp column is present
    if (!("methylProp" %in% colnames(BSDT))) {
        stop("BSDT must have a methylProp column.")
    }
    
    if (! is.null(excludeGR)) {
        BSDT <- BSFilter(BSDT, minReads = 0, excludeGR)
    }
    
    bsgr <- BSdtToGRanges(list(BSDT));
    
    colModes <- sapply(BSDT, mode);
    if (is.null(sumCols)) {
        sumCols <- setdiff(colnames(BSDT), c("chr", "start", "end", 
                                            "strand", splitFactor, keepCols))
        # Restrict to numeric columns.      
        sumCols <- intersect(sumCols, 
                            names(colModes[which(colModes == "numeric")]))
        
    }
    # It's required to do a findoverlaps on each region individually, 
    # Not on a GRL, because of the way overlaps with GRLs work. So, 
    # we must convert the GRL to a GR, but we must keep track of which
    # regions came from which group.
    regionsGR <- unlist(regionsGRL)
    
    if (is.null(regionsGRL.length)) {
        if (length(regionsGRL) > 100) {
            message(cleanws("BSAggregate: Calculating sizes. You can speed this
                             up by supplying a regionsGRL.length vector..."),
                             appendLF = FALSE)
        }
        regionsGRL.length <- sapply(regionsGRL, length)
        # message("Done counting regionsGRL lengths.");
    }
    
    # Build a table to keep track of which regions belong to which group
    region2group <- data.table(
        regionID = seq_along(regionsGR), 
        chr = as.vector(seqnames(regionsGR)), 
        start = as.vector(start(regionsGR)), 
        end = as.vector(end(regionsGR)), 
        withinGroupID = as.vector(unlist(sapply(regionsGRL.length, seq))), 
        regionGroupID = rep(seq_along(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)
    
    
    # message("Finding overlaps...");
    fo <- findOverlaps(bsgr[[1]], regionsGR)
    
    setkey(BSDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$coverage == BSDT$coverage))
    
    # message("Setting regionIDs...");
    BSDT <- BSDT[queryHits(fo), ] # restrict the table to CpGs in any region.
    
    if (NROW(BSDT) < 1) {
        warning("No BSDT sites in the given region list. 
                Please expand your regionsGRL")
        return(NULL)
    }
    
    BSDT[, regionID := subjectHits(fo)] # record which region they overlapped.
    # if (!keep.na) {
    # BSDT <- BSDT[queryHits(fo), ]
    #}
    
    if (is.null(jCommand)) {
        cols <- c(sumCols, keepCols)
        funcs <- c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
        jCommand <- buildJ(cols, funcs)
    }
    # message("jCommand: ", jCommand)
    
    # Build the by string
    if (is.null(splitFactor)) {
        byString <- paste0("list(regionID)");
    } else {
        byString <- paste0("list(", paste("regionID", paste0(splitFactor, ""), 
                                         collapse = ",", sep = ","), ")")
    }
    
    # Now actually do the aggregate:
    # message("Combining...");
    # for MIRA: average methylProp and sum coverage within each instance of 
    # each bin (ie average methylProp's in bin1 of first region, average 
    # methylProp's in bin1 of second region etc. for all bins and all regions
    # separately)
    bsCombined <- BSDT[, eval(parse(text = jCommand)), 
                      by = eval(parse(text = byString))]
    setkey(bsCombined, regionID)
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups, 
    # which takes awhile. I think this preserves memory and is faster.
    
    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        # must set allow = TRUE here in case there are multiple IDs (splitCol)
        # adds regionGroupID column from region2group to bsCombined
        bsCombined[region2group, regionGroupID := regionGroupID, allow = TRUE]
        if (! is.null(splitFactor)) { 
            byStringGroup <- paste0("list(", 
                                   paste("regionGroupID", 
                                         paste0(splitFactor, collapse = ","), 
                                         sep = ","), 
                                   ")")
        } else {
            byStringGroup <- "list(regionGroupID)"
        }
       
        
        # actual aggregation operation
        # for normal MIRA use: averaging methylProp's and summing
        # coverage for all bins with the same number (ie all
        # bin1's together, all bin2's together, etc.)
        # NOTE: 2nd use of the jCommand so if the first jCommand use changed
        # column names that are required by the jCommand 
        # then this 2nd jCommand use will cause an error
        bsCombined <- bsCombined[, eval(parse(text = jCommand)), 
                                by = eval(parse(text = byStringGroup))]
        
        # if any strand information was not given, averaging the signatures 
        # about the center to account for unknown strand orientation, 
        # also averaging coverage about center
        # ie if any "*" are present then average
        if ("*" %in% unique(as.character(strand(regionsGR)))) {
            bsCombined[, methylProp := (methylProp + rev(methylProp)) / 2]
            if (hasCoverage) {
                bsCombined[, coverage := (coverage + rev(coverage)) / 2]
            }
        }
        
        # changing "regionGroupID" name to "bin" which is less confusing
        # for normal MIRA use cases
        setnames(bsCombined, old = "regionGroupID", new = "bin")
        return(bsCombined[]);
    } else {
        warning(cleanws("Using byRegionGroup = FALSE may 
             result in missing functionalities such as symmetrical averaging"))
        e <- region2group[bsCombined, ]
        setkey(e, regionID);
        return(e);
    }
    # WARNING: There are now 2^2 ways to aggregate, sum vs mean
    # at each level: across regions, then across region sets. THis
    # doesn't give you a choice at this point. 
}



#' Divide region into similarly sized bins.
#'
#' Given a start, end, and number of bins, to divide, 
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
#' @param strand "strand" column of the data.table (or single
#' strand value if binRegion is only used on one region). Default is "*".
#'
#' @return
#' A data.table, expanded to nrow = number of bins, with these id columns:
#'      id: region ID
#'      binID: repeating ID (this is the value to aggregate across)
#'      ubinID: unique bin IDs
#' @export
#' @examples
#' library(data.table)
#' start <- c(100, 1000, 3000)
#' end <- c(500, 1400, 3400)
#' chr <- c("chr1", "chr1", "chr2")
#' strand <- c("*", "*", "*")
#' # strand not included in object 
#' # since MIRA assumes "*" already unless given something else
#' regionsToBinDT <- data.table(chr, start, end)
#' numberOfBins <- 15
#' # data.table "j command" using column names and numberOfBins variable
#' binnedRegionDT <- regionsToBinDT[, binRegion(start, end, numberOfBins, chr)]
binRegion <- function(start, end, bins, idDF = NULL, strand = "*") {
    # if (!is.null(idDF) & (!is(idDF, "data.frame"))) {
    #   stop("idDF should be a data.frame")
    #}
    
    # conditionally altered later
    finalColNames <- c("chr", "start", "end", "id", "binID", "ubinID")
    if (!("*" %in% strand)) {
        if ("+" %in% strand) {
            plusIndex <- which(strand == "+")
            # once for plus strand coordinates
            plusStart <- start[plusIndex]
            plusEnd <- end[plusIndex]
            
            binSize <- (plusEnd - plusStart) / (bins)
            breaks <- round(rep(plusStart, each = (bins + 1)) 
                           + (0:(bins)) * rep(binSize, each = (bins + 1)))
            
            endpoints <- (bins + 1) * (1:(length(plusStart)))
            startpoints <- 1 + (bins + 1)  * (0:(length(plusStart) - 1))
            
            plusDT <- data.table(start = breaks[-endpoints], 
                                end = breaks[-startpoints], 
                                id = rep(plusIndex, each = bins), 
                                binID = 1:bins, 
                                strand = "+", 
                                key = "id")
            dt <- plusDT # placeholder but may be returned
        }
        if ("-" %in% strand) {
            minusIndex <- which(strand == "-")
            
            minusStart <- start[minusIndex]
            minusEnd <- end[minusIndex]
            
            binSize <- (minusEnd - minusStart) / (bins)
            breaks <- round(rep(minusStart, each = (bins + 1)) 
                           + (0:(bins)) * rep(binSize, each = (bins + 1)))
            
            endpoints <- (bins + 1) * (1:(length(minusStart)))
            startpoints <- 1 + (bins + 1)  * (0:(length(minusStart) - 1))
            
            minusDT <- data.table(start = breaks[-endpoints], 
                                 end = breaks[-startpoints], 
                                 id = rep(minusIndex, each = bins), 
                                 binID = bins:1, 
                                 strand = "-", 
                                 key = "id")
            dt <- minusDT # placeholder but may be returned
        }
        
        # if there are both + and - strands
        # combining and sorting plus and minus data.tables 
        if (("+" %in% strand) && ("-" %in% strand)) {
            dt <- rbindlist(list(plusDT, minusDT))
            setorder(x = dt, id, binID)  # setorder(dt, id) might also work?
        }
        
        # included only if + /- are present
        finalColNames <- c(finalColNames, "strand")
        dt[, ubinID := 1:nrow(dt)] 
        
    }else { # some strand information is "*", don't flip bin directions
        
        binSize <- (end - start) / (bins)
        breaks <- round(rep(start, each = (bins + 1)) 
                       + (0:(bins)) * rep(binSize, each = (bins + 1)))
        
        endpoints <- (bins + 1) * (1:(length(start)))
        startpoints <- 1 + (bins + 1)  * (0:(length(start) - 1))
        # do all regions in the same order
        dt <- data.table(start = breaks[-endpoints], 
                        end = breaks[-startpoints], 
                        id = rep((seq_along(start)), each = bins), 
                        binID = 1:bins, 
                        ubinID = seq_along(breaks[-startpoints]), 
                        key = "id")
        
    }
    
    if (!is.null(idDF)) {
        chr <- rep(idDF, each = bins)
        dt <- dt[, chr := chr]
        setcolorder(dt, finalColNames)  # putting chr first, does not copy
    }
    
    return(dt[])
}