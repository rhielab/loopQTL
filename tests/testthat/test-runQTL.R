test_that("runQTL() rejects an object with no phenotype", {
    obj <- loopQTLExampleData()
    obj@phenotype <- matrix(numeric(0), 0, 1)
    expect_error(runQTL(obj, window = 10000), regexp = "phenotype")
})

test_that("runQTL() on the example object returns a well-formed result table", {
    ## The example object was already run through runQTL(). Re-running keeps
    ## the test hermetic (no fresh normalize step) and exercises the full
    ## two-anchor merge + RINT pipeline against the shipped VCF.
    obj <- loopQTLExampleData()

    ## Give it a tiny PC matrix as covariates so the model has a fixed form.
    pcs <- matrix(
        rnorm(length(samples(obj)) * 2),
        nrow = length(samples(obj)), ncol = 2,
        dimnames = list(samples(obj), c("PC1", "PC2"))
    )

    obj2 <- runQTL(
        obj,
        window     = 10000,
        covariates = pcs,
        normal     = TRUE
    )

    r <- results(obj2)
    expect_s3_class(r, "data.frame")
    expect_gt(nrow(r), 0)
    expect_true(all(
        c("loopID", "snpID", "beta", "se", "p", "nUsed", "maf", "anchor") %in%
            names(r)
    ))
    expect_true(all(r$p >= 0 & r$p <= 1, na.rm = TRUE))
})
