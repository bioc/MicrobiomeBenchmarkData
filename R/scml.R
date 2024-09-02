#' SCML: spike-in-based calibration to total microbial load
#' 
#' The \code{scml} function applies the
#' spike-in-based calibration to total microbial load (SCML) method to
#'
#' @param tse A treeSummarizedExperiment from the \code{getBenchmarkData}
#' function.
#' @param bac A character. One of the following options:
#' s = Salinibacter ruber (AF323500), r = Rhizobium radiobacter (AB247615),
#' a, = Alicyclobacillus acidiphilus (AB076660)
#' 
#' @return A TreeSummarizedExperiment with SCML data instead of counts.
#' @export
#'
#' @examples
#' tse <- getBenchmarkData("Stammler_2016_16S_spikein", dryrun = FALSE)[[1]]
#' tseSCML <- scml(tse, bac = "s")
#' 
scml <- function(tse, bac = c("s", "r", "a")) {
    bacLetter <- match.arg(bac)
    bacNames <- c(
        s = "AF323500XXXX", r = "AB247615XXXX", a = "AB076660XXXX"
    )
    bacFullNames <- c(
        AF323500XXXX = "Salinibacter ruber (AF323500)",
        AB247615XXXX = "Rhizobium radiobacter (AB247615)",
        AB076660XXXX = "Alicyclobacillus acidiphilus (AB076660)"
        
    )
    bacName <- bacNames[bacLetter]
    
    lgl <- bacName %in% rownames(tse)
    if (!lgl) {
        stop(
            "Feature", bacName, "not found.",
            "Are you sure you're using the Stammler_2016_16S_spikein",
            call. = FALSE
        )
    }
    message("Re-calibrating counts with ", bacFullNames[bacName])
    counts <- SummarizedExperiment::assay(tse, 1)
    bacAb <- counts[bacName, ]
    sizeFactor <- bacAb/mean(bacAb)
    scmlData <- counts 
    for(i in seq(ncol(scmlData))){
        scmlData[,i] <- round(scmlData[,i] / sizeFactor[i])
    }
    SummarizedExperiment::assay(tse) <- scmlData
    return(tse)
}
