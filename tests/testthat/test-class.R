test_that("loopQTLExampleData() returns a valid loopQTL object", {
    obj <- loopQTLExampleData()

    expect_s4_class(obj, "loopQTL")
    expect_true(methods::validObject(obj))
})

test_that("accessor methods return the expected slot types", {
    obj <- loopQTLExampleData()

    expect_type(samples(obj), "character")
    expect_true(length(samples(obj)) > 0)
    expect_false(anyDuplicated(samples(obj)) > 0)

    expect_s4_class(loops(obj), "GInteractions")
    expect_equal(length(loops(obj)), nrow(phenotype(obj)))

    expect_type(resolution(obj), "integer")
    expect_gt(resolution(obj), 0)

    expect_true(is.character(vcfPath(obj)))
    expect_true(file.exists(vcfPath(obj)))

    expect_true(is.matrix(phenotype(obj)))
    expect_equal(colnames(phenotype(obj)), samples(obj))

    expect_s3_class(results(obj), "data.frame")
    expect_true(all(
        c(
            "loopID", "snpID", "chrom", "pos", "ref", "alt",
            "beta", "se", "p", "nUsed", "maf"
        ) %in% names(results(obj))
    ))

    expect_type(metadata(obj), "list")
    expect_true("resolution" %in% names(metadata(obj)))
})

test_that("show() prints a summary without error", {
    obj <- loopQTLExampleData()
    expect_output(show(obj), "loopQTL object")
    expect_output(show(obj), "samples")
    expect_output(show(obj), "loops")
})

test_that("cis() returns a named numeric vector matching samples", {
    obj <- loopQTLExampleData()
    cisVals <- cis(obj)

    expect_type(cisVals, "double")
    if (length(cisVals) > 0) {
        expect_equal(names(cisVals), samples(obj))
        expect_true(all(cisVals > 0))
    }
})
