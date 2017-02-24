# PACKAGE DOCUMENTATION
#' Methylation-based Inference of Regulatory Activity (MIRA)
#' MIRA is a score that measures the degree of dip in methylation
#' level surrounding a regulatory site of interest, such as a 
#' transcription factor binding sites.
#' This script provides functions for aggregating methylation 
#' data across region sets, in bins
#'
#' @docType package
#' @name MIRA
#' @author Nathan Sheffield
#'
#' @references \url{http://github.com/sheffien}
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand seqnames
#' @importFrom ggplot2 ggplot aes
#' @import BiocGenerics S4Vectors IRanges
#' @importFrom data.table ":=" setDT data.table setkey fread setnames as.data.table setcolorder melt setkeyv rbindlist
NULL

#' Function to aggregate methylation data into bins over all regions in each region set;
#' 
#' 
#'
#' @param BSDT A single data table that has DNA methylation data on individual sites including a "chr" column with chromosome, 
#'  a "start" column with the coordinate number for the cytosine, a "methyl" column with proportion of methylation (0 to 1),
#'  a "hitCount" column with number of methylated reads for each site, and a "readCount" column with total number of reads for each site. A "sampleName" column is preferred.
#' @param GRList A GRangesList object containing region sets, each set corresponding to a regulatory element.
#' Each regionSet in the list should be named. A named list of data.tables also works. 
#' @param binNum How many bins each region should be split into for aggregation of the DNA methylation data
#' @param minReads Filter out bins with fewer than X reads before returning.
#' @param sampleNameInBSDT boolean for whether the BSDT has a sampleName column
#' @param sampleType could be case/control, tissue type, etc.
#' 
#' @return a data.table with binNum rows for each region set containing aggregated methylation data
#' Each region was split into bins; methylation was put in these bins; 
#' Output contains sum of the all corresponding bins for the regions of each region set 
#' ie for all regions in each region set: first bins summed, second bins summed, etc
#' 
#' @export
#' @example
#' R/examples/example.R
returnMIRABins = function(BSDT,GRList, binNum=11, minReads = 500, sampleNameInBSDT=TRUE,sampleType=NULL){
  
  #checking that input is in list format and converting
  if (!class(GRList) %in% c("list","GRangesList")){
    stop("GRList should be a named list/GRangesList.")
    # if (class(GRList) %in% "GRanges"){
    #   GRList=GRangesList(GRList)
    #   message("Converting...")
    # }
    # if (class(GRList) %in% "data.table"){
    #   GRList=list(GRList)
    #   message("Converting...")
    # }
  }
  
  if (is.null(names(GRList))){
    stop("GRList should be a named list/GRangesList.")
  }
  
  #checking that all objects in GRList are the same type 
  if (all(sapply(X = GRList,FUN = class) %in% "GRanges")){
    GRDTList=lapply(X = GRList,FUN = grToDt)#GRanges to data.tables
  }else if (all(sapply(X = GRList,FUN = class) %in% "data.table")){
    GRDTList=GRList #this case is okay
  }else{
    stop("GRList should be a GRangesList or a list of data.tables")
  }
  
  
  methylByBin=lapply(X = GRDTList, FUN = function(x) BSBinAggregate(BSDT = BSDT,rangeDT = x, binCount = binNum,splitFactor=NULL,minReads = minReads))
  names(methylByBin)=names(GRList)#preserving names
  #adding a feature ID column to each data.table that should identify what region set was used
  for (i in 1:length(methylByBin)){
    methylByBin[[i]][,featureID := rep(names(methylByBin)[i],nrow(methylByBin[[i]]))]
  }
  #screening out region sets that had incomplete binning
  binNumScreen=sapply(X = methylByBin,FUN = nrow)
  methylByBin=methylByBin[!(binNumScreen<binNum)]#taking out incomplete region sets
  
  bigMethylByBin=rbindlist(methylByBin)
  if (sampleNameInBSDT){
    bigMethylByBin[,sampleName := rep(BSDT[1,sampleName])] #creating new sampleName column
  }
  
  
  return(bigMethylByBin)
}

#' Function to take DNA methylation and region data and return a MIRA score;
#' a wrapper for returnMIRABins and scoreDip
#' 
#'
#' @param BSDT A single data table that has DNA methylation data on individual sites including a "chr" column with chromosome, 
#'  a "start" column with the coordinate number for the cytosine, a "methyl" column with proportion of methylation (0 to 1),
#'  a "hitCount" column with number of methylated reads for each site, and a "readCount" column with total number of reads for each site. A "sampleName" column is preferred.
#' @param GRList A GRangesList object containing region sets, each set corresponding to a regulatory element;
#' Each regionSet in the list should be named. 
#' @param binNum How many bins each region should be split into for aggregation of the DNA methylation data.
#' @param scoringMethod Method to calculate MIRA score after binning, includes "logRatio", "area". See scoreDip function.
#' @param sampleNameInBSDT boolean for whether the BSDT has a sampleName column
#' @param sampleType could be case/control, tissue type, etc.
#' 
#' @export
MIRAScore = function(BSDT,GRList, binNum=11, scoringMethod="logRatio",sampleNameInBSDT=TRUE,sampleType=NULL){
  MIRAresults=list()
  
  
  bigBin=returnMIRABins(BSDT,GRList, binNum,sampleNameInBSDT,sampleType)
  
  #using binned methylation data to calculate MIRA score
  scoreDT=bigBin[,.(score=scoreDip(methyl,binNum,method=scoringMethod)),by=.(featureID,sampleName)]
  
  return(scoreDT)
}

#' My dip scoring function - for MIRA scores;
#' That's Methylation-based Inference of Regulatory Activity
#' just calculates the ratio between the flank (shoulders) and midpoint
#'
#' @param shoulderShift The number of bins away from the center to use as the
#' shoulders. I have used 5 or 3.
#' 
#' @export
#' @example
#' R/examples/example.R
scoreDip = function(values, binCount, shoulderShift = 5,method="logRatio") {
  if (method=="logRatio"){
	centerSpot = ceiling(binCount/2)
	leftSide = centerSpot - shoulderShift  # 3
	rightSide = centerSpot + shoulderShift  # 3
	midpoint = (values[centerSpot] + values[centerSpot+1] + values[centerSpot-1] ) /3
	# log ratio...
	score=log ( ((values[leftSide] + values[rightSide])/2) / midpoint )
  }
  
	#alternate way of scoring by the area in the dip
	if (method=="area"){
	  maxMethyl=max(values)
	  score=maxMethyl*binCount-sum(values)
	}
  
  #another alternate method
  if (method=="parabola"){
    #fit2 <- lm(y~poly(x,2,raw=TRUE))
    #lines(xx, predict(fit2, data.frame(x=xx)), col="green")
    
  }
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
#' @param start 
#'
#' @return
#' A data.table, expanded to nrow= number of bins, with these id columns:
#' 		id: region ID
#' 		binID: repeating ID (this is the value to aggregate across)
#' 		ubinID: unique bin IDs
#' @export
#' @examples
#' #load example region set
#' cgIslandsDT = data.table(...)
#' binnedCGI = cgIslandsDT[, binRegion(start,end, 50)]
binRegion = function(start, end, bins, idDF=NULL) {
	#if (!is.null(idDF) & ( ! "data.frame"  %in% class(idDF))) {
	#	stop("idDF should be a data.frame")
	#}
	binSize = (end-start)/(bins)
	breaks = round(rep(start, each=(bins+1)) + (0:(bins)) * rep(binSize, each=(bins+1)))

	endpoints = (bins+1) * (1:(length(start)))
	startpoints = 1 + (bins+1)  * (0:(length(start)-1))
	#TODO: remove this code split
	if (is.null(idDF)) {
		dt = data.table(start=breaks[-endpoints], end=breaks[-startpoints], id=rep((1:length(start)), each=bins), binID= 1:bins, ubinID=1:length(breaks[-startpoints]), key="id")
	} else {
		chr = rep(idDF, each=bins)
		dt = data.table(chr, start=breaks[-endpoints], end=breaks[-startpoints], id=rep((1:length(start)), each=bins), binID= 1:bins, ubinID=1:length(breaks[-startpoints]), key="id")

	}
	return(dt)
}

#' A wrapper of BSAggregate that first bins regions and then aggregates
#' each bin across a set of regions, individually.
#'
#' Produced originally for binning Ewing RRBS data across various region sets
#'
#' @param rangeDT A data table with the sets of regions to be binned, 
#' with columns named start, end
#' @param binCount Number of bins across the region
#' @param byRegionGroup Pass along to binCount (see ?binCount)
#' @param minReads Filter out bins with fewer than X reads before returning.
#' 
#' @export
BSBinAggregate = function(BSDT, rangeDT, binCount, minReads = 500, byRegionGroup=TRUE,splitFactor="id") {
	if (! "data.table" %in% class(rangeDT)) {
		stop("rangeDT must be a data.table")
	}
	seqnamesColName = "seqnames"  # default column name
	if (! "seqnames" %in% colnames(rangeDT)) {
		if ("chr" %in% colnames(rangeDT)) {
			message("seqnames column name set to: chr")
			seqnamesColName = "chr"
		} else {
			# Got neither.
			stop("rangeDT must have a seqnames column")
		}
	}

	#if(!silent){
	  message("Binning...")
	#}
	binnedDT = rangeDT[, binRegion(start, end, binCount, get(seqnamesColName))]
	binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
	message("Aggregating...")
	binnedBSDT = BSAggregate(BSDT, regionsGRL=GRangesList(binnedGR), jCommand=buildJ(c("methyl", "readCount"), c("mean", "sum")), byRegionGroup=byRegionGroup, splitFactor=splitFactor)
	# If we aren't aggregating by bin, then don't restrict to min reads!
	if (byRegionGroup) {
		binnedBSDT = binnedBSDT[readCount >= minReads,]
	}
	
	
	return(binnedBSDT)
}

#' BSaggregate -- Aggregate a BSDT across regions or region groups,
#' for multiple samples at a time.
#' This function is as BScombineByRegion, but can handle not only multiple
#' samples in BSDT, but also simultaneously multiple region sets by passing
#' a regionsGRL (GRangesList object).
#' you can use jCommand to do other functions.

#' Given a bisulfite data table as input, with an identifier column for
#' different samples; plus a GRanges objects with regions to aggregate.
#'
#' @param BSDT The bisulfite data.table (output from one of the parsing
#' functions for methylation calls) that you wish to aggregate. It can
#' be a combined table, with individual samples identified by column passed
#' to splitFactor.
#' @param regionsGRL Regions across which you want to aggregate.
#' @param excludeGR A GenomicRanges object with regions you want to 
#' exclude from the aggregation function. These regions will be eliminated
#' from the input table and not counted.
#' @param jCommand You can pass a custom command in the j slot to data.table
#' specifying which columns to aggregate, and which functions to use. You
#' can use buildJ() to build a jCommand argument easily.
#' @param byRegionGroup You can aggregate by regionID or by regionGroupID; 
#' this reflects the regionsGRL that you pass; by default, BSAggregate will
#' aggregate each region individually -- scores will then be contiguous, and
#' the output is 1 row per region.
#' Turn on this flag to aggregate across all region groups, making the result
#' uncontiguous, and resulting in 1 row per *region group*.
#'
BSAggregate = function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, sumCols=NULL, jCommand=NULL, byRegionGroup=FALSE, keep.na=FALSE) {

	# Assert that regionsGRL is a GRL.
	# If regionsGRL is given as a GRanges, we convert to GRL
	if( "GRanges" %in% class(regionsGRL)) {
		regionsGRL = GRangesList(regionsGRL);
	} else if (! "GRangesList" %in% class(regionsGRL)) {
		stop("regionsGRL is not a GRanges or GRangesList object");
	}

	if(! is.null(excludeGR)) {
		BSDT = BSFilter(BSDT, minReads=0, excludeGR)
	}

	bsgr = BSdtToGRanges(list(BSDT));

	additionalColNames = setdiff(colnames(BSDT), c("chr","start", "end","hitCount","readCount", splitFactor));

	colModes = sapply(BSDT,mode);
	if (is.null(sumCols)) {
		sumCols = setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
		# Restrict to numeric columns.		
		sumCols = intersect(sumCols, names(colModes[which(colModes == "numeric")]))

	}
	# It's required to do a findoverlaps on each region individually,
	# Not on a GRL, because of the way overlaps with GRLs work. So,
	# we must convert the GRL to a GR, but we must keep track of which
	# regions came from which group.
	regionsGR = unlist(regionsGRL)
	
	if(is.null(regionsGRL.length)) {
		if (length(regionsGRL) > 100) {
		message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF=FALSE)
		}
		regionsGRL.length = sapply(regionsGRL, length)
		message("Done counting regionsGRL lengths.");
	}

	# Build a table to keep track of which regions belong to which group
	region2group = data.table(
		regionID=1:length(regionsGR), 
		chr=as.vector(seqnames(regionsGR)), 
		start=as.vector(start(regionsGR)), 
		end=as.vector(end(regionsGR)),
		withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
		regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
	setkey(region2group, regionID)


	message("Finding overlaps...");
	fo = findOverlaps(bsgr[[1]], regionsGR)

	setkey(BSDT, chr, start)
	# Gut check:
	# stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))

	message("Setting regionIDs...");
	BSDT = BSDT[queryHits(fo),] #restrict the table to CpGs in any region.

	if (NROW(BSDT) < 1) {
		warning("No BSDT sites in the given region list; please expand your regionsGRL")
		return(NULL)
	}

	BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
	#BSDT[queryHits(fo),regionID:=subjectHits(fo)]
	#if (!keep.na) {
	#	BSDT = BSDT[queryHits(fo),]
	#}

	if (is.null(jCommand)) {
		cols=c(sumCols, keepCols)
		funcs = c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
		jCommand = buildJ(cols, funcs)
	}
	message("jCommand: ", jCommand)
	
	# Define aggregation column. aggregate by region or by region group?
	if (byRegionGroup) {
		agCol = "regionGroupID";
	} else {
		agCol = "regionID"; # Default
	}

	# Build the by string
	if (is.null(splitFactor)) {
		byString = paste0("list(regionID)");
	} else {
		byString = paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
	}

	# Now actually do the aggregate:
	message("Combining...");
	bsCombined = BSDT[,eval(parse(text=jCommand)), by=eval(parse(text=byString))]
	setkey(bsCombined, regionID)
	# Now aggregate across groups.
	# I do this in 2 steps to avoid assigning regions to groups,
	# which takes awhile. I think this preserves memory and is faster.

	# Define aggregation column. aggregate by region or by region group?
	if (byRegionGroup) {
		# must set allow=TRUE here in case there are multiple IDs (splitCol)
		bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
		if (! is.null(splitFactor) ) { 
			byStringGroup = paste0("list(", paste("regionGroupID", paste0(splitFactor, collapse=", "), sep=", "), ")")
		} else {
			byStringGroup = "list(regionGroupID)"
		}
		bsCombined=bsCombined[,eval(parse(text=jCommand)), by=eval(parse(text=byStringGroup))]
		
		#if no strand information was given, averaging the signatures about the center...
		#...to account for unknown strand orientation, also averaging readCount about center
		if (!strand(regionsGR)@values %in% c("+","-")){
		  bsCombined[,methyl := (methyl+rev(methyl))/2]
		  bsCombined[,readCount := (readCount+rev(readCount))/2]
		}
		
		return(bsCombined);
	} else {
	  warning("Using byRegionGroup=FALSE may result in missing functionalities such as symmetrical averaging")
		e = region2group[bsCombined,]
		setkey(e, regionID);
		return(e);
	}
	# WARNING: There are now 2^2 ways to aggregate, sum vs mean
	# at each level: across regions, then across region sets. THis
	# doesn't give you a choice at this point. 
}


#' A function to plot the binned methylation of samples over a region.
#' 
#' @param binnedRegDT A datatable with specific column names containing the following:
#' bin numbers(binnedRegionDT), aggregated methylation values (methyl), name of the region set
#' (featureID), case/control column (sampleType), sample name (sampleName).
#' @param featID Region set names in a single string or vector of strings.
#' @param plotType Line or jitter (ggplot2). 
#' @export
plotMIRARegions <- function(binnedRegDT,featID=unique(binnedRegDT[,featureID]),plotType="line"){
  setkey(binnedRegDT,featureID)
  binPlot=ggplot(data=binnedRegDT[featID], mapping = aes(x=regionGroupID,y = methyl))
  if (!("sampleType" %in% names(binnedRegDT))){
    sampleType="All samples" #if no sampleType column then all lines/points will be black
    warning("sampleType column must exist if it is desired to split up sample types by color")
  }
  if (plotType=="line"){
    binPlot=binPlot+geom_line(aes(col=sampleType,group=sampleName))+facet_wrap(~featureID)
  }else if (plotType=="jitter"){
    binPlot=binPlot+geom_jitter(aes(col=sampleType))+facet_wrap(~featureID)
  }else {
    stop('The only supported values for plotType are "line" and "jitter"')
  }
  return(binPlot)
}

#' A function to plot MIRA scores and compare case/control.
#' 
#' @param scoreDT A datatable with the following columns: score, featureID (names of regions),sampleType.
#' @param featID Region set name/names in a single string or vector of strings.
#' @export
plotMIRAScores <- function(scoreDT,featID=unique(scoreDT[,featureID])){
  setkey(scoreDT,featureID)
  scorePlot=ggplot(data=scoreDT[featID], mapping = aes(x=sampleType,y=score))+
    geom_boxplot()+geom_jitter()+facet_wrap(~featureID)  
  return(scorePlot)
}


#' Adding methyl column that has proportion of reads that were methylated for each site.
#' Note: Assigns methyl column by reference with ":="
#' 
#' @param BSDTList A bisulfite datatable or list of datatables with a column for number of methylated 
#' reads (hitCount) and a column for number of total reads (readCount) for each cytosine that 
#' was measured.
#' @export
addMethCol <- function(BSDTList){
  
  #converting to a data.table list if it was a single data.table
  if ("data.table" %in% class(BSDTList)){
    BSDTList=list(BSDTList)
  }
  
  #stopping the function if the input was not data.table originally
  if (!"data.table" %in% class(BSDTList[[1]])){
    stop('Input must be a single data.table object or list of data.table objects')
  }
  
  #using anonymous function to apply operation that adds methyl column to each element of list
  #extra [] on the end is necessary for proper display/printing of the object
  BSDTList=lapply(X = BSDTList,FUN = function(x) x[, methyl := round(hitCount/readCount, 3)][])
  
  return(BSDTList)
}

#' Function to normalize case/experimental samples to the controls
#' 
#' It finds the median of the controls for each bin for each region set
#' then divides by the median (oldBinVal/medianBinVal) for each bin.
#' @param binnedDT A datatable containing bins for each region set for each sample;bins contain
#' aggregated methylation across regions for that sample; 
#' it should have a column with sample annotation so cases and controls can be split up as
#' well as annotation of the region sets (featureID column) and sampleType column (case/control).
#' 
#' @export
normalizeMIRA = function(binnedDT){
  if ("list" %in% class(binnedDT)){
    binnedDT=rbindlist(binnedDT)
  }
  if (!"data.table" %in% class(binnedDT)){
    stop("binnedDT must be a data.table (or list of data.tables)")
  }
  
  features=unique(binnedDT[,featureID]) #get a set of all features
  
  setkey(binnedDT,featureID,sampleType)
  
  #getting regionGroupIDs for the loop
  regGroupID=unique(binnedDT[,regionGroupID])
  #normalize for each feature
  for (i in 1:length(features)){ #get median methylation for each feature
    medMeth=binnedDT[.(features[i],"control"), median(methyl),by= regionGroupID]
    for (j in regGroupID){#divide each methylation values for each region group by the appropriate value
      #do operation on the rows for this feature and only one regionGroupID at a time
      binnedDT[featureID==features[i] & regionGroupID==j,methyl := methyl/medMeth[regionGroupID == j,V1] ]
    }
    
  } 
  return(binnedDT)
}

#' helper function
#' given a vector of colums, and the equally-sized vector of functions
#' to apply to those columns, constructs a j-expression for use in
#' a data.table.
#' use it in a DT[,eval(parse(text=buildJ(cols,funcs)))]
buildJ = function(cols, funcs) {
  r = paste("list(", paste(paste0(cols, "=", funcs, "(", cols, ")"), collapse=","), ")")
  return(r);
}

#Two utility functions for converting data.tables into GRanges objects
#genes = dtToGR(gModels, "chr", "txStart", "txEnd", "strand", "geneId");
dtToGrInternal = function(DT, chr, start, end=NA, strand=NA, name=NA, metaCols=NA) {
  if (is.na(end)) {
    if ("end" %in% colnames(DT)) {
      end = "end"
    } else {
      end = start;
    }
  }
  if (is.na(strand)) {
    gr=GRanges(seqnames=DT[[`chr`]], ranges=IRanges(start=DT[[`start`]], end=DT[[`end`]]), strand="*")
  } else {
    # GRanges can only handle '*' for no strand, so replace any non-accepted
    # characters with '*'
    DT[,strand:=as.character(strand)]
    DT[strand=="1", strand:="+"]
    DT[strand=="-1", strand:="-"]
    DT[[`strand`]] =  gsub("[^+-]", "*", DT[[`strand`]])
    gr=GRanges(seqnames=DT[[`chr`]], ranges=IRanges(start=DT[[`start`]], end=DT[[`end`]]), strand=DT[[`strand`]])
  }
  if (! is.na(name) ) {
    names(gr) = DT[[`name`]];
  } else {
    names(gr) = 1:length(gr);
  }
  if(! is.na(metaCols)) {
    for(x in metaCols) {
      elementMetadata(gr)[[`x`]]=DT[[`x`]]
    }
  }
  gr;
}

dtToGr = function(DT, chr="chr", start="start", end=NA, strand=NA, name=NA, splitFactor=NA, metaCols=NA) {
  if(is.na(splitFactor)) {
    return(dtToGrInternal(DT, chr, start, end, strand, name,metaCols));
  }
  if ( length(splitFactor) == 1 ) { 
    if( splitFactor %in% colnames(DT) ) {
      splitFactor = DT[, get(splitFactor)];
    }
  }
  lapply(split(1:nrow(DT), splitFactor), 
         function(x) { 
           dtToGrInternal(DT[x,], chr, start, end, strand, name,metaCols)
         }
  )
  
  
}

dtToGR = dtToGr;

#Converts a list of data.tables (From BSreadbeds) into GRanges.
BSdtToGRanges = function(dtList) {
  gList = list();
  for (i in 1:length(dtList)) {
    #dt = dtList[[i]];
    setkey(dtList[[i]], chr, start)
    #convert the data into granges object
    gList[[i]] = GRanges(seqnames=dtList[[i]]$chr, ranges=IRanges(start=dtList[[i]]$start, end=dtList[[i]]$start), strand=rep("*", nrow(dtList[[i]])), hitCount=dtList[[i]]$hitCount, readCount=dtList[[i]]$readCount)
    #I used to use end=start+1, but this targets CG instead of just a C, and it's causing edge-effects problems when I assign Cs to tiled windows using (within). Aug 2014 I'm changing to start/end at the same coordinate.
  }
  return(gList);
}

# This can run into memory problems if there are too many files...
# because of the way parallel lacks long vector support. The solution is
# to just use a single core; or to pass mc.preschedule=FALSE; This
# makes it so that each file is processed as a separate job. Much better.
#' Read in files from biseq meth caller
#' @param files	a list of filenames (use parseInputArg if necessary)
#' @param contrastList	a list of named character vectors, each with length equal to the number of items in files. These will translate into column names in the final table.
#' @param sampleNames	a vector of length length(files), name for each file. You can also just use contrastList to implement the same thing so this is really unnecessary...
#' @param cores	number of processors.
#' @export
BSreadBiSeq = function(files, contrastList=NULL, sampleNames=extractSampleName(files), cores=4, returnAsList=FALSE) {
  cores=min(length(files), cores); #not more cores than files!
  setLapplyAlias(cores);
  if (!is.null(contrastList)) {
    if( any(sapply(contrastList, length) != length(files))) {
      stop("contrastList must be a list, with each value having the same number of elements as files.");
    }
  }
  message("Reading ", length(files), " files..");
  freadList = lapplyAlias(files, fread, mc.preschedule=FALSE);
  colNames = names(contrastList)
  message("File reading finished (", length(files), " files). Parsing Biseq format...", appendLF=FALSE);
  # TODO: This parsing takes awhile, and could be done in parallel.
  freadListParsed = lapplyAlias(freadList, parseBiseq, mc.preschedule=FALSE)
  
  message("Parsing complete, building final tables and cleaning up...")
  numberOfFiles = length(freadListParsed);
  for (i in 1:numberOfFiles) {
    if (numberOfFiles > 1) {
      message(i, ": ", sampleNames[i], "; ", appendLF=FALSE)
    }
    if (numberOfFiles > 1 && i == numberOfFiles){
      message("", appendLF=TRUE)
    }
    DT = freadListParsed[[i]]; #convenience alias.
    if(!is.null(contrastList)) {
      DT[, get("colNames"):=as.list(sapply(contrastList, "[[", i))]
    }
    if (!is.null(sampleNames)) {
      DT[,sampleName:=sampleNames[i]]
    }
    freadListParsed[[i]] = DT
  }
  
  #filteredList = do.call(rbind, freadListParsed)
  #gc(); #rbind call is memory inefficient; this helps.
  # rbindlist supposedly does the same thing as do.call(rbind, list) but 
  # faster
  if (!returnAsList){ #default (returnAsList=FALSE) is to return as one combined data.table/data.frame
    filteredList = 	rbindlist(freadListParsed)
  }else{
    filteredList = freadListParsed
  }
  
  return(filteredList);
}

#' Takes a data.table from BSreadBiSeq and parses the strange x/y format
#' of methylation calls, splitting them into individual columns
#' @param DT data.table to parse
#' @return data.table with separate methylated and unmethylated columns
#' 
#' 
parseBiseq = function(DT) {
  message(".", appendLF=FALSE);
  setnames(DT, paste0("V", 1:6), c("chr", "start", "end", "meth", "rate", "strand"))
  DT[, meth:=gsub("'", "", meth)]
  #split the '12/12' format of meth calls
  ll = unlist(strsplit(DT$meth, "/", fixed=TRUE))
  idx = seq(1, length(ll), by = 2)
  DT[, `:=`(hitCount =  as.integer(ll[idx]), readCount =  as.integer(ll[idx+1]))]
  DT[,start := as.integer(start+1)] #re-index
  DT[, c("rate", "end", "meth" ):=NULL] #remove unnecessary columns
  DT[, strand:=NULL]
  DT=DT[,list(hitCount= sum(hitCount), readCount=sum(readCount)), by=list(chr, start)] #smash measurements
  setcolorder(DT,c("chr", "start", "hitCount", "readCount"));
  DT = DT[ !grep("_",chr),]; #clean Chrs
  return(DT)
}


#' convert a GenomicRanges into a data.table
#' 
#' @param a GRanges object
#' @return A data.table object.
grToDt = function(GR) {
  DF=as.data.frame(elementMetadata(GR))
  if( ncol(DF) > 0) {
    DT = data.table(chr=as.vector(seqnames(GR)), start=start(GR), end=end(GR), DF)
  } else {
    DT = data.table(chr=as.vector(seqnames(GR)), start=start(GR), end=end(GR))
  }
  return(DT)
}

#' This function is a drop-in replacement for the base list() function,
#' which automatically names your list according to the names of the 
#' variables used to construct it.
#' It seemlessly handles lists with some names and others absent,
#' not overwriting specified names while naming any unnamed parameters.
#'
#' @param ...	arguments passed to list()
#' @return A named list object.
#' @export
#' @examples
#' x=5
#' y=10
#' nlist(x,y) # returns list(x=5, y=10)
#' list(x,y) # returns unnamed list(5, 10)
nlist = function(...) {
  fcall = match.call(expand.dots=FALSE)
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
#' @param cores	Number of cpus
#' @return None
setLapplyAlias = function(cores=0) {
  if (cores < 1) {
    return(getOption("mc.cores"))
  }
  if(cores > 1) { #use multicore?
    if (requireNamespace("parallel", quietly = TRUE)) {
      options(mc.cores=cores)
    } else {
      warning("You don't have package parallel installed. Setting cores to 1.")
      options(mc.cores=1) #reset cores option.
    }
  } else {
    options(mc.cores=1) #reset cores option.
  }
}

#' Function to run lapply or mclapply, depending on the option set in
#' getOption("mc.cores"), which can be set with setLapplyAlias().
#'
#' @param ... Arguments passed lapply() or mclapply()
#' @param mc.preschedule Argument passed to mclapply
#' @return Result from lapply or parallel::mclapply
lapplyAlias = function(..., mc.preschedule=TRUE) {
  if (is.null(getOption("mc.cores"))) { setLapplyAlias(1) }
  if(getOption("mc.cores") > 1) {
    return(parallel::mclapply(..., mc.preschedule=mc.preschedule))
  } else {
    return(lapply(...))
  }
}

# extract sample names from file names as the first part of the file name (before any suffix)
extractSampleName = function(fileNames, suffixSep="\\.", pathSep="/") {
  sapply(strsplit(fileNames,pathSep),function(x) strsplit(rev(x)[1],suffixSep)[[1]][1])
}
