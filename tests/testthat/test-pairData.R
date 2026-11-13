test_that("loopQTLpairData() returns per-sample rows for a valid pair", {
    obj <- loopQTLExampleData()

    loopID_ <- obj@results$loopID[1]
    snpID_ <- obj@results$snpID[1]

    pair <- loopQTLpairData(obj, loopID = loopID_, snpID = snpID_)

    expect_s3_class(pair, "data.frame")
    expect_true(all(
        c("sample", "loopID", "snpID", "genotype", "phenotype") %in%
            names(pair)
    ))
    expect_equal(nrow(pair), length(samples(obj)))
    expect_setequal(pair$sample, samples(obj))
    expect_true(all(pair$genotype %in% c(0L, 1L, 2L, NA)))
})

test_that("loopQTLpairData() errors on unknown loopID / snpID", {
    obj <- loopQTLExampleData()

    expect_error(
        loopQTLpairData(obj, loopID = "not_a_loop", snpID = obj@results$snpID[1]),
        regexp = "loopID"
    )
    expect_error(
        loopQTLpairData(obj, loopID = obj@results$loopID[1], snpID = "rsNope"),
        regexp = "snpID"
    )
})
