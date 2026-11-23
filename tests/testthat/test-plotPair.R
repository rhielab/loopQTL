## Smoke tests for loopQTLplotPair(): returns a ggplot for a known hit,
## fails cleanly on an unknown pair.

test_that("loopQTLplotPair() returns a ggplot for a known hit", {
    obj <- loopQTLExampleData()
    p <- loopQTLplotPair(
        obj,
        loopID = "L0000037_chr2_172207000_172457000",
        snpID  = "rs12052253"
    )
    expect_s3_class(p, "ggplot")
})

test_that("loopQTLplotPair() errors on an unknown (loop, SNP) pair", {
    obj <- loopQTLExampleData()
    expect_error(
        loopQTLplotPair(obj, loopID = "not_a_loop", snpID = "rs12052253"),
        regexp = "not found|pair"
    )
    expect_error(
        loopQTLplotPair(
            obj,
            loopID = "L0000037_chr2_172207000_172457000",
            snpID  = "rsNope"
        ),
        regexp = "not found|pair"
    )
})
