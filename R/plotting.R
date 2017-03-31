#Plotting functions to visualize outputs of MIRA
#Visualize MIRA signatures and MIRA scores




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
    
    if (!("sampleType" %in% names(binnedRegDT))) {
        sampleType = "All samples" 
        #if no sampleType column then all lines/points will be black
        warning("sampleType column must exist if it is desired to split up sample types by color")
    }
    if (plotType == "line") {
        binPlot = binPlot + 
            geom_line(aes(col = sampleType, group = sampleName)) + 
            facet_wrap(~featureID)
    }else if (plotType == "jitter") {
        binPlot = binPlot + geom_jitter(aes(col = sampleType)) + 
            facet_wrap(~featureID)
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
#' @return a plot of class "gg"/"ggplot" that shows MIRA scores 
#' with geom_boxplot and geom_jitter.
#' @export
#' @examples
#' data(ewingMyobigBinDT) # bigBinDT object
#' exScores = bigBinDT[, .(score=scoreDip(methyl, 
#'                                        binCount = 21, 
#'                                        shoulderShift = 5)), 
#'                     by = .(featureID, sampleName)]
#' #adding annotation
#' sampleType = rep(c("Ewing", "Muscle-related"), each = 24)
#' exScores = cbind(exScores, sampleType)
#' exScorePlot = plotMIRAScores(exScores)         
plotMIRAScores <- function(scoreDT, featID = unique(scoreDT[, featureID])){
    setkey(scoreDT, featureID)
    scorePlot = ggplot(data = scoreDT[featID], 
                       mapping = aes(x = sampleType, y = score)) + 
        geom_boxplot() + geom_jitter() + 
        facet_wrap(~featureID)  
    return(scorePlot)
}