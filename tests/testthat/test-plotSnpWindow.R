## Smoke test for loopQTLplotSnpWindow(): returns a ggplot for the focal
## SNP using the shipped example SNP-window cache.

test_that("loopQTLplotSnpWindow() returns a plot object for a known SNP", {
    obj <- loopQTLExampleData()
    p <- loopQTLplotSnpWindow(
        obj,
        chrom      = "chr2",
        start      = 172100000,
        end        = 172550000,
        snpID      = "rs12052253",
        normal     = TRUE,
        covariates = "auto",
        Kmax       = 2
    )
    expect_true(inherits(p, "gg"))
})
