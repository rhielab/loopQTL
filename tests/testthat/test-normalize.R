## Tests for loopQTLnormalize() — one end-to-end smoke test on the shipped
## chr2 pseudo dataset, one saveWindows round-trip, and two error paths.

.setupPseudoObj <- function() {
    extdataDir <- system.file("extdata", package = "loopQTL")
    ss <- read.table(
        file.path(extdataDir, "example_samplesheet.tsv"),
        sep = "\t", header = TRUE, stringsAsFactors = FALSE
    )
    ss$hic   <- file.path(extdataDir, ss$hic)
    ss$loops <- file.path(extdataDir, ss$loops)
    loopQTLsetup(
        samples    = ss,
        vcf        = file.path(extdataDir, "example.vcf.gz"),
        resolution = 5000L,
        minSamples = 2L
    )
}

test_that("loopQTLnormalize() populates the phenotype matrix", {
    obj <- .setupPseudoObj()
    expect_equal(length(phenotype(obj)), 0L)

    obj <- loopQTLnormalize(
        obj, BPPARAM = BiocParallel::SerialParam()
    )
    expect_equal(nrow(phenotype(obj)), length(loops(obj)))
    expect_equal(ncol(phenotype(obj)), length(samples(obj)))
    expect_identical(colnames(phenotype(obj)), samples(obj))
    ## No infinities should leak into the normalized phenotype
    expect_false(any(is.infinite(phenotype(obj))))
})

test_that("loopQTLnormalize(saveWindows = TRUE) writes a non-empty .rds", {
    obj <- .setupPseudoObj()
    windowsOut <- tempfile(fileext = ".rds")
    on.exit(file.remove(windowsOut), add = TRUE)

    obj <- loopQTLnormalize(
        obj,
        BPPARAM     = BiocParallel::SerialParam(),
        saveWindows = TRUE,
        windowsOut  = windowsOut,
        snpWindow   = 10000L
    )
    expect_true(file.exists(windowsOut))
    w <- readRDS(windowsOut)
    expect_type(w, "list")
    expect_gt(length(w), 0L)

    ## The object should record the windows path and SNP count in metadata
    expect_equal(
        normalizePath(metadata(obj)[["windowsOut"]]),
        normalizePath(windowsOut)
    )
    expect_gt(metadata(obj)[["n_windows_saved"]], 0L)
})

test_that("loopQTLnormalize(saveWindows = TRUE) requires `windowsOut`", {
    obj <- .setupPseudoObj()
    expect_error(
        loopQTLnormalize(
            obj,
            BPPARAM     = BiocParallel::SerialParam(),
            saveWindows = TRUE
        ),
        regexp = "windowsOut"
    )
})

test_that("loopQTLnormalize() rejects non-loopQTL input", {
    expect_error(
        loopQTLnormalize(list(a = 1), BPPARAM = BiocParallel::SerialParam())
    )
})
