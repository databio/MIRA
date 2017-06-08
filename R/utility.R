#########################################################
# Utility functions
# Conversion functions in 2nd half of file
########################################################


# This function is a drop-in replacement for the base list() function, 
# which automatically names your list according to the names of the 
# variables used to construct it.
# It seamlessly handles lists with some names and others absent, 
# not overwriting specified names while naming any unnamed parameters.
#
# @param ...   arguments passed to list()
# @return A named list object.
# @examples
# x = 5
# y = 10
# nlist(x, y) # returns list(x = 5, y = 10)
# list(x, y) # returns unnamed list(5, 10)
nlist = function(...) {
    fcall = match.call(expand.dots = FALSE)
    l = list(...);
    if (!is.null(names(list(...)))) { 
        names(l)[names(l) == ""] = fcall[[2]][names(l) == ""]
    } else {  
        names(l) = fcall[[2]];
    }
    return(l)
}



# Function to run lapply or mclapply, depending on the option set in
# getOption("mc.cores"), which can be set with setLapplyAlias().
#
# @param ... Arguments passed lapply() or mclapply()
# @param mc.preschedule Argument passed to mclapply
# @return Result from lapply or parallel::mclapply
lapplyAlias = function(..., mc.preschedule = TRUE) {
    if (is.null(getOption("mc.cores"))) { setLapplyAlias(1) }
    if (getOption("mc.cores") > 1) {
        return(parallel::mclapply(..., mc.preschedule = mc.preschedule))
    } else {
        return(lapply(...))
    }
}



# To make parallel processing a possibility but not required, 
# I use an lapply alias which can point at either the base lapply
# (for no multicore), or it can point to mclapply, 
# and set the options for the number of cores (what mclapply uses).
# With no argument given, returns intead the number of cpus currently selected.
#
# @param cores Number of cpus
# @return None
setLapplyAlias = function(cores = 0) {
    if (cores < 1) {
        return(getOption("mc.cores"))
    }
    if (cores > 1) { # use multicore?
        if (requireNamespace("parallel", quietly = TRUE)) {
            options(mc.cores = cores)
        } else {
            warning(cleanws("You don't have package parallel installed. 
                    Setting cores to 1."))
            options(mc.cores = 1) # reset cores option.
        }
    } else {
        options(mc.cores = 1) # reset cores option.
    }
}



# helper function
# given a vector of columns, and the equally-sized vector of functions
# to apply to those columns, constructs a j-expression for use in
# a data.table 
# (functions applied to columns in corresponding spot in "cols" string).
# One function may be given to be applied to multiple columns.
# use it in a DT[, eval(parse(text = buildJ(cols, funcs)))]
# @param cols A string/vector of strings containing columns 
# on which to use functions.
# @param funcs Functions to use on columns.
# @return A jcommand string. After performing function on column, column 
# is reassigned the same name.
buildJ = function(cols, funcs) {
    r = paste("list(", paste(paste0(cols, "=", funcs, "(", cols, ")"), collapse = ","), ")")
    return(r);
}


####################  Conversion Functions  ########################
# Convert between different object types

# Convert a data.table to GRanges object.
# 
# @param DT a data.table with at least "chr" and "start" columns
# 
# @return gr A genomic ranges object derived from DT
dtToGr = function(DT, chr = "chr", start = "start", 
                  end = NA, strand = NA, name = NA, 
                  splitFactor = NA, metaCols = NA) {
    
    if (is.na(splitFactor)) {
        return(dtToGrInternal(DT, chr, start, end, strand, name, metaCols));
    }
    if ( length(splitFactor) == 1 ) { 
        if ( splitFactor %in% colnames(DT) ) {
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
    if (is.na(strand)) {
        if ("strand" %in% colnames(DT)) { # checking if strand info is in DT
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
    if (! is.na(name)) {
        names(gr) = DT[[`name`]];
    } else {
        names(gr) = seq_along(gr);
    }
    if (! is.na(metaCols)) {
        for(x in metaCols) {
            elementMetadata(gr)[[`x`]] = DT[[`x`]]
        }
    }
    gr;
}



# Converts a list of data.tables into GRanges.
# @param dtList A list of data.tables, 
# Each should have "chr", "start", "methylCount", and "coverage" columns.
# Error results if missing "chr", "start" but if methylCount and coverage are
# missing, it will still work, just not have that info in the output.
# @return a list of GRanges objects, strand has been set to "*", 
# "start" and "end" have both been set to "start" of the DT.
# methylCount and coverage info is preserved in GRanges object.
BSdtToGRanges = function(dtList) {
    
    gList = list();
    for (i in seq_along(dtList)) {
        # dt = dtList[[i]];
        setkey(dtList[[i]], chr, start)
        # convert the data into granges object
        gList[[i]] = GRanges(seqnames = dtList[[i]]$chr, 
                             ranges = IRanges(start = dtList[[i]]$start, 
                                              end = dtList[[i]]$start), 
                             strand = rep("*", nrow(dtList[[i]])), 
                             methylCount = dtList[[i]]$methylCount, 
                             coverage = dtList[[i]]$coverage)
        # I used to use end = start + 1, but this targets CG instead of just 
        # a C, and it's causing edge-effects problems when I assign Cs to 
        # tiled windows using (within). Aug 2014 I'm changing to start/end at 
        # the same coordinate.
    }
    return(gList);
}



# Convert a GenomicRanges into a data.table
# 
# @param GR A GRanges object
# @param includeStrand Boolean, whether to include strand from GR in output DT
# @return A data.table object with columns:
# "chr", "start", and "end" (possibly strand)
grToDt = function(GR, includeStrand = FALSE) {
    DF = as.data.frame(elementMetadata(GR))
    if ( ncol(DF) > 0) {
        if (includeStrand) {
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
        if (includeStrand) {
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



# cleanws takes multi-line, code formatted strings and just formats them
# as simple strings
# @param string string to clean
# @return A string with all consecutive whitespace characters, including
# tabs and newlines, merged into a single space.
cleanws = function(string) {
    return(gsub('\\s+'," ", string))
}

#' Make MIRA-compatible data.tables using information 
#' from SummarizedExperiment-based classes
#' 
#' Of the methylCountDF, coverageDF, and methylPropDF arguments, either 
#' methylPropDF must be given or both methylCountDF and coverageDF must
#' be given. After that, whichever is not given will be calculated 
#' from the others. If coverageDF and methylCountDF are both not given,
#' coverage will be assumed to be 1 at each site with a methylProp value.
#' Acceptable formats for the three "DF" parameters include:
#' data.frame, data.table, matrix, and DelayedMatrix classes.
#' 
#' @param coordinates Coordinates for the methylation loci as a GRanges object 
#' (in same order as methylCountDF/coverageDF/methylPropDF, whichever is given).
#' The start coordinate should be the coordinate of the cytosine.
#' @param methylCountDF An object of matrix/data.frame/similar format
#' that contains the number of reads with methylated C's for the loci in
#' 'coordinates' argument. It should have one column per sample with
#' rows that correspond to 'coordinates'.
#' @param coverageDF An object of matrix/data.frame/similar format
#' that contains the total number of reads for the loci in
#' 'coordinates' argument. It should have one column per sample with
#' rows that correspond to 'coordinates'.
#' @param methylPropDF An object of matrix/data.frame/similar format
#' that contains the total number of reads for the loci in
#' 'coordinates' argument. It should have one column per sample with
#' rows that correspond to 'coordinates'.
#' @param sample_names A character vector with sample names in the same
#' order as the columns of methylCountDF/coverageDF/methylPropDF
#' samples in the columns
#' @return MIRAFormatBSDTList A list of data.tables containing
#' the methylation data. One data.table per sample with the column
#' names: 'chr', 'start' (methylation loci), 'methylCount' (number of
#' methylated reads), 'coverage' (total number of reads), and 
#' 'methylProp' (proportion of methylated reads). The order of the
#' list is the order of samples in the columns of 
#' methylCountDF/coverageDF/methylPropDF. If sample names 
#' are explicitly given as input ('sample_names' argument)
#' or can be derived from other arguments to the function then
#' a named list will be returned. 
#' @examples  
#' data("exampleBSseqObj")
#' MIRAFormatBSDTList = SumExpToMIRA(coordinates = granges(exampleBSseqObj), 
#'     methylCountDF = getCoverage(BSseq = exampleBSseqObj, type = "M"), 
#'     coverageDF = getCoverage(BSseq = exampleBSseqObj, type = "Cov"),
#'     methylPropDF = getMeth(BSseq = exampleBSseqObj, type = "raw"),
#'     sample_names = bsseq::sampleNames(exampleBSseqObj))
SumExpToMIRA = function(coordinates, methylCountDF=NULL, 
                        coverageDF=NULL, methylPropDF=NULL, 
                        sample_names=NULL) {
    
    haveMCount = !is.null(methylCountDF)
    haveCov = !is.null(coverageDF)
    haveMProp = !is.null(methylPropDF)
    
    # making sure sufficient data has been given
    if (!haveMProp) {
        # if no methylProp given, you need both methylCount and coverage
        if (!haveMCount || !haveCov) {
            stop(cleanws("If methylProp is not given, both 
                         methylCount and coverage must be given."))
        }
        }
    
    # checking that right data types have been given 
    if (!("GRanges" %in% class(coordinates))) {
        stop("'coordinates' argument should be a GRanges object")
    }
    
    # check for accepted formats?: delayed matrix?, matrix, data.frame, data.table
    # not sure if there might be other valid inputs so excluding this for now
    # if (haveMCount) {
    #     if (!("" %in% class(methylCountDF))) {
    #         stop()
    #     }
    # }
    # if (haveCov) {
    #     if (!("" %in% class(coverageDF))) {
    #         stop()
    #     }
    # }
    # if (haveMProp) {
    #     if (!("" %in% class(methylPropDF))) {
    #         stop()
    #     }
    # }
    
    
    MIRAFormatBSDTList = list() # to store output
    
    # changing coordinates from GRanges object to data.table
    coordinates = grToDt(coordinates)
    rowNum = nrow(coordinates)
    
    # setting sampleNum and doing checks on the inputs
    if (haveMProp) {
        sampleNum = ncol(methylPropDF)
        # making sure other inputs also have sampleNum columns
        if (haveMCount) {
            if (sampleNum != ncol(methylCountDF)) {
                stop(cleanws("Input with the count of 
                             methylated reads at each loci and input 
                             with the proportion of methylation should have
                             the same number of columns."))
            }
            }
        if (haveCov) {
            if (sampleNum != ncol(coverageDF)) {
                stop(cleanws("Input with the count of total reads at each loci 
                             and input with the proportion of methylation
                             should have the same number of columns."))
            }
            }
            } else {
                sampleNum = ncol(coverageDF)
                # making sure other input has the same number of columns 
                if (sampleNum != ncol(methylCountDF)) {
                    stop(cleanws("Input with the count of total reads  
                                 and input with the count of methylated reads
                                 at each loci
                                 should have the same number of columns."))
                }
                }
    
    # making sure number of rows are the same in all given inputs
    if (haveMProp) {
        if (rowNum != nrow(methylPropDF)) {
            stop(cleanws("The number of rows for the genomic coordinates and 
                         the number of rows for the proportion of methylation
                         for each loci should be the same."))
        }
        }
    if (haveMCount) {
        if (rowNum != nrow(methylCountDF)) {
            stop(cleanws("The number of rows for the genomic coordinates and 
                         the number of rows for the count of methylated reads
                         for each loci should be the same."))
        }
        }
    if (haveCov) {
        if (rowNum != nrow(coverageDF)) {
            stop(cleanws("The number of rows for the genomic coordinates and 
                         the number of rows for the count of total reads
                         for each loci should be the same."))
        }
        }
    
    # adding sample names if not given as argument but present in other input
    if (is.null(sample_names)) {
        if (haveMProp) {
            sample_names = colnames(methylPropDF)
        } else {
            sample_names = colnames(coverageDF)
        }
    }
    
    
    
    if (!haveCov && !haveMCount) {
        coverage = rep(1, rowNum)
        warning(cleanws("coverage was not included in input 
                        so it was set to 1.
                        It is recommended to set the minreads 
                        argument of aggregateMethyl to 0."))
    }
    
    # looping through columns/samples to make 1 data.table per sample
    for (i in 1:sampleNum) {
        # each column is a different sample
        # order of the below if statements is important
        if (haveMCount) {
            methylCount = methylCountDF[, i]
        }
        if (haveMProp) {
            methylProp = methylPropDF[, i]
        }
        if (haveCov) {
            coverage = coverageDF[, i]
        } else if (haveMCount) {
            # assumed that if !haveCov, you haveMProp
            coverage = round(methylCount / methylProp, 3)
        }
        if (haveCov && !haveMCount) {
            # theoretically this should be a whole number
            methylCount = round(methylProp * coverage, 0)
        }
        
        if (!haveCov && !haveMCount) {
            # this will result in fractional methylCount values (from 0 to 1)
            # but is necessary to be consistent with giving each loci a 
            # coverage of 1 which was done earlier in code 
            # for !haveCov && !haveMCount
            # methylCount = methylProp * coverage = methylProp * 1
            methylCount = methylProp
        }
        
        if (!haveMProp) {
            methylProp = round(methylCount / coverage, 3)
        }
        
        # index for taking out rows with 0 coverage
        notCovered = which(coverage == 0)
        # checking for NAs in methylProp
        # they only should be present if methylPropDF was given
        if (haveMProp) {
            notCoveredNA = which(is.na(methylProp))
            notCovered = sort(c(notCovered, notCoveredNA))
        }
        MIRAFormatBSDTList[[i]] = data.table(chr = coordinates[, chr], 
                                             start = coordinates[, start], 
                                             methylCount = as.vector(methylCount),
                                             coverage = as.vector(coverage),
                                             methylProp = as.vector(methylProp)
        )[!notCovered,] # filtering
        setnames(MIRAFormatBSDTList[[i]], 
                 c("chr", "start", "methylCount", "coverage", "methylProp"))
    }
    
    # add sample names to list (by reference)
    if (!is.null(sample_names)) {
        setattr(MIRAFormatBSDTList, "names", sample_names)
    }
    
    return(MIRAFormatBSDTList)
}

#' Convert a BSseq object to data.table format for MIRA.
#' 
#' Converts a BSseq object to a list of data.tables with one data.table
#' per sample. A wrapper of the SumExpToMIRA function.
#' 
#' @param BSseqObj An object of class BSseq, can have smoothed or raw
#' methylation data.
#' @return MIRAFormatBSDTList A list of data.tables containing
#' the methylation data. One data.table per sample with the column
#' names: 'chr', 'start' (methylation loci), 'methylCount' (number of
#' methylated reads), 'coverage' (total number of reads), and 
#' 'methylProp' (proportion of methylated reads). The order of the
#' list is the order of samples in the columns of the BSseq object.
#' If sample names are in the BSseq object, then a named list
#' will be returned.
#' @examples 
#' data("exampleBSseqObj")
#' MIRAFormatBSDTList = bsseqToMIRA(exampleBSseqObj)
#' 
bsseqToMIRA <- function(BSseqObj) {
    
    if (bsseq::hasBeenSmoothed(BSseqObj)) {
        rawSmooth = "smooth"
    } else {
        rawSmooth = "raw"
    }
    
    # use accessor functions from bsseq to get needed data
    MIRAFormatBSDTList = SumExpToMIRA(coordinates = granges(BSseqObj), 
                 methylCountDF = getCoverage(BSseq = BSseqObj, type = "M"), 
                 coverageDF = getCoverage(BSseq = BSseqObj, type = "Cov"),
                 methylPropDF = getMeth(BSseq = BSseqObj, type = rawSmooth),
                 sample_names = bsseq::sampleNames(BSseqObj))
    
    # # if only one sample, do not return list, just return data.table
    # if (length(MIRAFormatBSDTList) == 1) {
    #     MIRAFormatBSDTList = MIRAFormatBSDTList[[1]]
    # }
    
    return(MIRAFormatBSDTList)
}