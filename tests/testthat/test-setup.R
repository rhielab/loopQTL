## Build a samplesheet pointing at the shipped extdata, so the test is
## hermetic (no network, no user paths).
.exampleSamplesheet <- function() {
    extdataDir <- system.file("extdata", package = "loopQTL")
    s <- read.table(
        file.path(extdataDir, "example_samplesheet.tsv"),
        header = TRUE, stringsAsFactors = FALSE
    )
    s$hic <- file.path(extdataDir, s$hic)
    s$loops <- file.path(extdataDir, s$loops)
    s
}

test_that("loopQTLsetup() builds a valid loopQTL object from extdata", {
    samplesheet <- .exampleSamplesheet()
    vcfPath_ <- system.file("extdata", "example.vcf.gz", package = "loopQTL")

    obj <- loopQTLsetup(
        samples    = samplesheet,
        vcf        = vcfPath_,
        resolution = 10000,
        minSamples = 2
    )

    expect_s4_class(obj, "loopQTL")
    expect_equal(sort(samples(obj)), sort(samplesheet$Sample))
    expect_equal(resolution(obj), 10000L)
    expect_equal(vcfPath(obj), vcfPath_)
    expect_gt(length(loops(obj)), 0)
    ## Loops should be intra-chromosomal
    expect_true(all(as.character(
        GenomicRanges::seqnames(InteractionSet::anchors(loops(obj), type = "first"))
    ) == as.character(
        GenomicRanges::seqnames(InteractionSet::anchors(loops(obj), type = "second"))
    )))
})

test_that("loopQTLsetup() rejects malformed input", {
    expect_error(
        loopQTLsetup(
            samples = "not a dataframe", vcf = "x.vcf.gz",
            resolution = 10000
        ),
        regexp = "is.data.frame"
    )

    bad <- data.frame(foo = "sample1", bar = "f.hic")
    expect_error(
        loopQTLsetup(samples = bad, vcf = "x.vcf.gz", resolution = 10000),
        regexp = "Sample|hic|loops"
    )
})
