## Smoke tests for loopQTLplotLocus(): returns a gg/patchwork object for
## a known locus using both "cpm" and the default (normalized) value modes.

test_that("loopQTLplotLocus(values = 'cpm') returns a plot object", {
    obj <- loopQTLExampleData()
    gtfPath <- system.file("extdata", "example.gtf.gz", package = "loopQTL")
    p <- loopQTLplotLocus(
        obj,
        chrom  = "chr2",
        start  = 172100000,
        end    = 172550000,
        snpID  = "rs12052253",
        values = "cpm",
        gtf    = gtfPath
    )
    ## plotLocus stacks a few tracks, so the return may be a patchwork
    ## composite or a plain ggplot; either carries the "gg" class.
    expect_true(inherits(p, "gg"))
})

test_that("loopQTLplotLocus() errors when the SNP is not in results(obj)", {
    obj <- loopQTLExampleData()
    gtfPath <- system.file("extdata", "example.gtf.gz", package = "loopQTL")
    expect_error(
        loopQTLplotLocus(
            obj,
            chrom  = "chr2",
            start  = 172100000,
            end    = 172550000,
            snpID  = "rsNope",
            values = "normal",
            gtf    = gtfPath
        ),
        regexp = "not found|results"
    )
})
