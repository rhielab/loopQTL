test_that("example results table has the documented QTL columns", {
    obj <- loopQTLExampleData()
    r <- results(obj)

    expect_gt(nrow(r), 0)
    for (col in c("beta", "se", "t", "p", "nUsed", "maf")) {
        expect_true(col %in% names(r), info = col)
        expect_type(r[[col]], "double")
    }

    ## Statistical sanity: p-values in [0, 1], MAF in [0, 0.5]
    pVals <- r$p[!is.na(r$p)]
    expect_true(all(pVals >= 0 & pVals <= 1))
    mafVals <- r$maf[!is.na(r$maf)]
    expect_true(all(mafVals >= 0 & mafVals <= 0.5))

    ## nUsed cannot exceed the sample count
    expect_true(all(r$nUsed[!is.na(r$nUsed)] <= length(samples(obj))))

    ## anchor labels are A or B
    expect_true(all(r$anchor %in% c("A", "B")))
})

test_that("annotate columns from the example object are well-formed", {
    obj <- loopQTLExampleData()
    r <- results(obj)
    ## The example object has already been annotated.
    skip_if_not(all(c("gene", "peak", "Type") %in% names(r)))

    expect_type(r$gene, "character")
    expect_type(r$peak, "logical")
    expect_type(r$Type, "character")
    expect_true(all(!is.na(r$peak)))
})
