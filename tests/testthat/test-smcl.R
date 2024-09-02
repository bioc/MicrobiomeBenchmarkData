test_that("scml works", {
    tse <- suppressWarnings(getBenchmarkData("Stammler_2016_16S_spikein", dryrun = FALSE)[[1]])
    expect_s4_class(scml(tse, bac = "s"), "TreeSummarizedExperiment")
    expect_s4_class(scml(tse, bac = "r"), "TreeSummarizedExperiment")
    expect_s4_class(scml(tse, bac = "a"), "TreeSummarizedExperiment")
})
