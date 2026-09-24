#' Per-sample boxplot for a single SNP-loop pair
#'
#' This function renders the per-sample boxplot for a chosen result SNP-loop
#' pair. It is returned as a ggplot object that can be saved with
#' `ggplot2::ggsave()`.
#'
#' @param obj A loopQTL object with `results(obj)` filled by `runQTL()`.
#' @param loopID A single loopID to plot.
#' @param snpID A single snpID to plot.
#' @param normal Plot Rank-based Inverse Normal Transformed (RINT'd) phenotype
#' values. Defaults to TRUE. Set to `FALSE` to plot pre-RINT normalized values.
#' @param jitter Add a jittered point overlay (default `TRUE`).
#'
#' @return A `ggplot` object.
#'
#' @export
#' @examples
#' ## Load the example object. In real analysis, the input
#' ## should be the output of runQTL().
#' ##Here we just use the example object. 
#' obj <- loopQTLExampleData()
#' ##Check results so that you can choose your interested loopID and snpID
#' head(results(obj))
#' ##Then you can specify the specifi SNP-loop pair data you want to plot as classical QTL box plot.
#' ##Return is a ggplot object so you can save with ggsave().
#' p <- loopQTLplotPair(
#'     obj,
#'     loopID = "L0000187_chr1_46030000_46220000",
#'     snpID  = "rs11584814")
loopQTLplotPair <- function(obj, loopID, snpID, normal = TRUE, jitter = TRUE) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(loopID), length(loopID) == 1,
        is.character(snpID),  length(snpID) == 1
    )
    if (nrow(obj@results) == 0) {
        .stopNoCall("obj@results is empty. Run runQTL(obj) first.")
    }

    ## Look up the (loop, SNP) row in obj@results for stats
    hit <- obj@results[obj@results$loopID == loopID &
        obj@results$snpID == snpID, , drop = FALSE]
    if (nrow(hit) == 0) {
        .stopNoCall("(loopID, snpID) pair not found in obj@results")
    }
    hit <- hit[1, , drop = FALSE]
    ref <- hit$ref
    alt <- hit$alt
    genotypeLabels <- c(
        paste0(ref, "/", ref),
        paste0(ref, "/", alt),
        paste0(alt, "/", alt)
    )

    ## Pull per-sample pair data
    d <- loopQTLpairData(obj, loopID, snpID, normal = normal)
    d$genotypeFactor <- factor(d$genotype, levels = c(0, 1, 2), labels = genotypeLabels)
    d <- d[!is.na(d$genotypeFactor) & !is.na(d$phenotype), , drop = FALSE]
    if (nrow(d) == 0) {
        .stopNoCall("No samples have non-NA genotype and phenotype for this pair.")
    }

    ## Per-genotype sample counts for x-tick subtitles
    nPer <- as.integer(table(d$genotypeFactor))
    xLabels <- sprintf("%s\n(n=%d)", genotypeLabels, nPer)

    plotTitle <- sprintf("%s   |   %s", snpID, loopID)
    subMain <- sprintf(
        "REF=%s   ALT=%s   |   p = %.3g   |   beta = %.3g",
        ref, alt, hit$p, hit$beta
    )
    yLabel <- if (normal) "Loop phenotype (RINT)" else "Loop phenotype (normalized)"

    p <- ggplot2::ggplot(d, ggplot2::aes(.data$genotypeFactor, .data$phenotype)) +
        ggplot2::geom_boxplot(
            fill = "lightgray", outlier.size = 0.8,
            width = 0.5
        ) +
        ggplot2::scale_x_discrete(labels = xLabels) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::labs(
            title = plotTitle,
            subtitle = subMain,
            x = "Genotype",
            y = yLabel
        ) +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 11),
            plot.subtitle = ggplot2::element_text(size = 10)
        )
    if (jitter) {
        p <- p + ggplot2::geom_jitter(
            width = 0.12, alpha = 0.6,
            size = 1.0, color = "steelblue"
        )
    }
    p
}
