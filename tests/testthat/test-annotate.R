test_that("loopQTLannotate() adds gene/peak/Type to results", {
    ## Start from the pre-annotated example object, then strip annotation
    ## columns and re-run annotate to verify it populates them again.
    obj <- loopQTLExampleData()
    obj@results$gene <- NULL
    obj@results$peak <- NULL
    obj@results$Type <- NULL

    gtf <- system.file("extdata", "example.gtf.gz", package = "loopQTL")
    peak <- system.file("extdata", "example_annotation.bed",
        package = "loopQTL"
    )

    obj2 <- loopQTLannotate(
        obj,
        gtf             = gtf,
        peak            = peak,
        TSSwindow       = 2000,
        anchorWindow    = 5000,
        peakAnnotation  = "Enhancer",
        removeTSS       = TRUE
    )

    expect_s4_class(obj2, "loopQTL")
    expect_true(all(c("gene", "peak", "Type") %in% names(results(obj2))))
    expect_type(results(obj2)$peak, "logical")
    expect_true(any(nzchar(results(obj2)$Type)))
})

test_that("loopQTLannotate() errors on missing gtf file", {
    obj <- loopQTLExampleData()
    expect_error(
        loopQTLannotate(obj, gtf = "no-such.gtf.gz"),
        regexp = "file.exists"
    )
})
