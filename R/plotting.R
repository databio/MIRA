# Plotting functions to visualize outputs of MIRA
# Visualize MIRA signatures and MIRA scores




#' Plot summary methylation profile
#' 
#' Plot one or multiple methylation profiles. Displays each region set
#' in a different subplot.
#' If you only want to plot certain region sets,
#' subset with the `featID` parameter.
#' 
#' @param binnedRegDT A datatable with specific column names containing:
#' bin numbers(binnedRegionDT column), 
#' aggregated methylation values (methylProp column), 
#' name of the region set (featureID column), 
#' case/control column (sampleType column), 
#' sample name (sampleName column).
#' @param featID Region set names in a single string or vector of strings.
#' @param plotType Line or jitter (ggplot2). 
#' @param colBlindOption If TRUE, function will plot with a color blind
#' friendly palette which could be helpful when plotting multiple colors. 
#' 
#' @return A plot of class "gg"/ "ggplot" that shows MIRA signatures
#' @examples
#' data("exampleBins", package = "MIRA")
#' MIRAplot <- plotMIRAProfiles(binnedRegDT = exampleBins)
#' 
#' @export
plotMIRAProfiles <- function(binnedRegDT, 
                            featID = unique(binnedRegDT[, featureID]), 
                            plotType = "line",
                            colBlindOption = FALSE){
    binNum <- max(binnedRegDT[, bin])
    setkey(binnedRegDT, featureID)
    binPlot <- ggplot(data = binnedRegDT[featID], 
                     mapping = aes(x = factor(bin), 
                                   y = methylProp * 100)) +
                    theme_classic() + ylim(c(0, 100)) +
                    geom_hline(yintercept=c(0), alpha=.2) +
                    ylab("DNA Methylation (%)") + 
                    xlab("Genome Regions Surrounding Sites") +
                    scale_x_discrete(labels=xAxisForRegionPlots(binNum))
    if (colBlindOption) {
        binPlot <- binPlot + scale_color_brewer(name = "Sample Type", 
                                           palette = "Dark2")
    } else {
        binPlot <- binPlot + scale_color_brewer(name = "Sample Type", 
                                               palette = "Set1")
    }
        
    
    if (!("sampleType" %in% names(binnedRegDT))) {
        sampleType <- "All samples" 
        # if no sampleType column then all lines/points will be black
        warning(cleanws("sampleType column is required to split up 
                        sample types by color"))
    }
    if (plotType == "line") {
        binPlot <- binPlot + 
            geom_line(aes(col = sampleType, group = sampleName)) + 
            facet_wrap(~featureID)
    }else if (plotType == "jitter") {
        binPlot <- binPlot + geom_jitter(aes(col = sampleType), alpha = .4) + 
            facet_wrap(~featureID)
    }else {
        stop('The only supported values for plotType are "line" and "jitter"')
    }
    return(binPlot)
}

# A function to get right x axis numbers on the plotMIRAProfiles() plots
xAxisForRegionPlots <- function(binNum) {
    if ((binNum %% 2) == 0) { # even binNum
        tmp <- c((-1 * binNum / 2):-1, 1:(binNum / 2)) # no zero
        xAxis <- c(tmp[1], rep("", (binNum - 4) / 2), -1, 1, 
                  rep("", (binNum - 4) / 2), tmp[binNum])
    } else if ((binNum %% 2) == 1) { # odd binNum
        tmp <- (-1 * (binNum - 1) / 2):((binNum - 1) / 2)
        xAxis <- c(tmp[1], rep("", (binNum - 3) / 2), 0, 
                  rep("", (binNum - 3) / 2), tmp[binNum])
    }
    return(xAxis)
}

#' Plot MIRA scores and compare different conditions
#' 
#' Splits up samples by sample type. Displays each region set
#' in a different subplot.
#' If you only want to plot certain region sets,
#' subset with the `featID` parameter. 
#' 
#' Due to the limited number
#' of colors in the palette, a warning will be issued if 
#' there are too many (more than 9) region sets ('featureID's).
#' 
#' @param scoreDT A datatable with the following columns: 
#' score, featureID (names of regions), ideally include 'sampleType'.
#' @param featID Region set name/names in a single string or vector of strings.
#' @param colBlindOption If TRUE, function will plot with a color blind
#' friendly palette which could be helpful when plotting multiple colors. 
#' @return a plot of class "gg"/"ggplot" that shows MIRA scores 
#' with geom_boxplot and geom_jitter (or geom_violin instead
#' of boxplot if no sampleType column is given).
#' @export
#' @examples
#' data(bigBinDT2)
#' exScores <- bigBinDT2[, .(score=calcMIRAScore(methylProp,
#'                                        shoulderShift = 5)), 
#'                     by = .(featureID, sampleName)]
#' # adding annotation
#' sampleType <- rep(c("Ewing", "Muscle-related"), each = 24)
#' exScores <- cbind(exScores, sampleType)
#' exScorePlot <- plotMIRAScores(exScores)         
plotMIRAScores <- function(scoreDT, 
                           featID = unique(scoreDT[, featureID]),
                           colBlindOption = FALSE){
    # the preferred option when 'sampleType' is a column
    if ("sampleType" %in% colnames(scoreDT)) {
        sampleTypeNum <- length(unique(scoreDT[, sampleType]))
        setkey(scoreDT, featureID)
        scorePlot <- ggplot(data = scoreDT[featID], 
                           mapping = aes(x = sampleType, 
                                         y = score, 
                                         col = sampleType)) + 
            theme_classic() +
            ylab("MIRA Score") + xlab("Sample Type") +
            geom_boxplot(aes(fill = sampleType), alpha = 0.75) + 
            geom_jitter(data = scoreDT[featID], 
                        mapping = aes(x = sampleType, y = score)) + 
            scale_color_manual(guide = FALSE, values = rep("black", 
                                                           sampleTypeNum)) +
            facet_wrap(~featureID)
        if (colBlindOption) {
            scorePlot <- scorePlot + scale_fill_brewer(name = "Sample Type", 
                                                      palette="Dark2")
        } else {
            scorePlot <- scorePlot + scale_fill_brewer(name = "Sample Type", 
                                                      palette="Set1")
        }
        
    } else {
        # a less fancy plot when sampleType column is not present
        setkey(scoreDT, featureID)
        scorePlot <- ggplot(data = scoreDT[featID], 
                           mapping = aes(x = "", y = score, 
                                         col = featureID)) + 
            theme_classic() +
            theme(axis.title.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            ylab("MIRA Score") + # xlab("") +
            geom_violin(aes(fill = featureID), alpha = 0.75) +
            geom_jitter() +
            scale_color_manual(guide = FALSE, 
                               values = rep("black", length(featID))) +
            facet_wrap(~featureID)
        if (colBlindOption) {
            scorePlot <- scorePlot + scale_fill_brewer(name = "Region Set", 
                                                      palette="Dark2")
        } else {
            scorePlot <- scorePlot + scale_fill_brewer(name = "Region Set", 
                                                      palette="Set1")
        }
        
    }
        
    return(scorePlot)
}
