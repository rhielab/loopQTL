test_that("loopQTLchooseK() returns K and a sample x K PC matrix (elbow)", {
    obj <- loopQTLExampleData()

    Kres <- loopQTLchooseK(obj, Kmethod = "elbow")

    expect_type(Kres, "list")
    expect_true(all(c("K", "PCs", "method") %in% names(Kres)))
    expect_equal(Kres$method, "elbow")

    expect_true(is.matrix(Kres$PCs))
    expect_equal(nrow(Kres$PCs), length(samples(obj)))
    expect_equal(ncol(Kres$PCs), Kres$K)
    expect_equal(rownames(Kres$PCs), samples(obj))

    ## K is capped at nSamples - 2 to leave df for intercept + genotype
    expect_lte(Kres$K, length(samples(obj)) - 2L)
    expect_gte(Kres$K, 1L)
})

test_that("loopQTLchooseK() with BE method returns a valid K", {
    obj <- loopQTLExampleData()

    ## Deterministic permutation runs; keep B small for speed.
    set.seed(1)
    Kres <- loopQTLchooseK(obj, Kmethod = "BE", B = 3)

    expect_equal(Kres$method, "BE")
    expect_gte(Kres$K, 1L)
    expect_equal(nrow(Kres$PCs), length(samples(obj)))
})
