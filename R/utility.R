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
            warning("You don't have package parallel installed. Setting cores to 1.")
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



# check whether object is smoothed
# check for names in phenoData
# fix names in data.table
# add methylProp column?, use addMethCol or bsseq built in getMeth()?
# @param bsseqObj An object of class bsseq
# @return MIRAFormatBSDTList A list of data.tables in MIRA format
# One data table for each sample column of the bsseq object.
bsseqToMIRA <- function(bsseqObj){
    # if (bsseq::hasBeenSmoothed(bsseqObj)) {
    #   warning("Raw (unsmoothed) methylation and coverage values are being used.")
    # }
    MIRAFormatBSDTList = list() # to store output
    # obtaining coordinates as GRanges obj. and changing to data.table
    coordinates = grToDt(granges(bsseqObj))
    for (i in 1:ncol(bsseqObj)) { # each column is a different sample
        methylCount = bsseq::getBSseq(BSseq = bsseqObj[, i], type = "M")
        coverage = bsseq::getBSseq(BSseq = bsseqObj[, i], type = "Cov")
        # index for taking out rows with 0 coverage
        notCovered = which(coverage == 0)
        warning(cleanws("Taking out rows with no coverage. Genomic coordinates may 
                         not have identical row numbers in different samples now."))
        MIRAFormatBSDTList[[i]] = data.table(chr = coordinates[, chr], 
                                             start = coordinates[, start], 
                                             methylCount = methylCount, 
                                             coverage = coverage
        )[!notCovered]  # filtering
        setnames(MIRAFormatBSDTList[[i]], 
                 c("chr", "start", "methylCount", "coverage"))
    }
    # names for list (by reference)
    setattr(MIRAFormatBSDTList, "names", Biobase::sampleNames(bsseqObj))
    
    return(MIRAFormatBSDTList)
}


#' cleanws takes multi-line, code formatted strings and just formats them
#' as simple strings
#' @param string string to clean
#' @return A string with all consecutive whitespace characters, including
#' tabs and newlines, merged into a single space.
cleanws = function(string) {
    return(gsub('\\s+'," ", string))
}