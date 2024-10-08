% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scml.R
\name{scml}
\alias{scml}
\title{SCML: spike-in-based calibration to total microbial load}
\usage{
scml(tse, bac = c("s", "r", "a"))
}
\arguments{
\item{tse}{A treeSummarizedExperiment from the \code{getBenchmarkData}
function.}

\item{bac}{A character. One of the following options:
s = Salinibacter ruber (AF323500), r = Rhizobium radiobacter (AB247615),
a, = Alicyclobacillus acidiphilus (AB076660)}
}
\value{
A TreeSummarizedExperiment with SCML data instead of counts.
}
\description{
The \code{scml} function applies the
spike-in-based calibration to total microbial load (SCML) method to
}
\examples{
tse <- getBenchmarkData("Stammler_2016_16S_spikein", dryrun = FALSE)[[1]]
tseSCML <- scml(tse, bac = "s")

}
