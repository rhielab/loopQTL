#' Hi-C-style beta heatmap around a SNP
#'
#' This function renders a Hi-C-style heatmap centered on a chosen
#' SNP: for every bin-pair in the SNP's cached contact window (from
#' `loopQTLnormalize(saveWindows = TRUE)`), the color encodes the per-sample
#' regression coefficient (beta) of the bin-pair's contact strength on the SNP's
#' dosage. Set `normal = TRUE` and `covariates = "auto"` to match the model
#' used in `runQTL()`.
#'
#' @param obj A loopQTL object with `metadata$windowsOut` filled by
#'   `loopQTLnormalize(saveWindows = TRUE)`.
#' @param snpID A single snpID present in the saved windows file.
#' @param chrom Optional chrom for verification. If supplied, must match
#'   the SNP's chrom in the windows file.
#' @param start Optional Start genomic position in bp to crop the heatmap to.
#' @param end Optional End genomic position in bp to crop the heatmap to.
#' @param windows Either `NULL` (default, loads from `metadata(obj)$windowsOut`),
#'   a path to a windows `.rds`, or the loaded list.
#' @param mirror Logical; mirror the upper-triangle beta values to the
#'   lower triangle so the heatmap is symmetric. Default `TRUE`.
#' @param minSamplesPerCell Minimum usable (non-NA dose + complete covariates)
#'   samples required to fit the regression for a given cell. Cells with fewer samples
#'   are returned NA. Default 3.
#' @param normal Apply Rank-based Inverse Normal Transformation (RINT) per
#'   cell across samples before fitting. Matches the `normal = TRUE`
#'   default in [runQTL()]. Default `FALSE`.
#' @param covariates Covariate handling for the regression:
#'   \describe{
#'     \item{`NULL`}{Default. Fit `cell ~ dosage` per cell (no covariates).}
#'     \item{`"auto"`}{Auto-pick K phenotype PCs via PCAForQTL on
#'       `phenotype(obj)` (RINT'd if `normal = TRUE`), then fit
#'       `cell ~ dosage + PCs` per cell. Matches what [runQTL()] tested,
#'       so beta here equals the `beta` column in `results(obj)`.}
#'     \item{matrix/data.frame}{User-supplied `samples x k` covariates;
#'       partial regression as above. Rownames or matched order required.}
#'   }
#' @param Kmethod `"elbow"` (default) or `"BE"`. Passed through to
#'   [loopQTLchooseK()] when `covariates = "auto"`. Ignored otherwise.
#' @param Kmax Hard cap on the auto-chosen K when `covariates = "auto"`.
#'   Defaults to the number of samples minus 2 (same as [runQTL()]). Pass a
#'   smaller value (e.g. `floor(nSamples / 4)`) for a more conservative
#'   covariate set.
#' @param limit Symmetric color-scale limit on beta. Values beyond
#'   `+/-limit` are limited to the deepest color. Default
#'   `NULL`) in the cropped view.
#' @param palette RColorBrewer palette name for the color pattern in the heatmap. Default `"RdBu"`.
#' @param metric What to render per cell:
#'   \describe{
#'     \item{`"betaLogp"`}{Default. Use raw * -log10(p) values to plot, this consider both beta
#'     and significane of the regression analysis.}
#'     \item{`"beta"`}{Use raw beta values to plot.}
#'   }
#' @param values Which cell values to run the regression on:
#'   \describe{
#'     \item{`"normalized"`}{Default. Use the phenotype(obj) to perform the regression}
#'     \item{`"cpm"`}{Use counts per million matrix to perform the regression.}
#'     \item{`"raw"`}{Use raw straw counts to perform the regression}
#'   }
#' @param axisUnit Unit for the graph axes - `"Mb"` (default) or `"bp"`.
#'
#' @return A `ggplot` object.
#'
#' @export
#' @examples
#' ## Load the example object. In real analysis, the input
#' ## should be the output of runQTL().
#' ##Also please note that this function requires the object to have the contact matrix `rda` file. 
#' ##Check `obj@windowsOut` to see if your object contains this file.
#' ##Here we just use the example object.
#' obj <- loopQTLExampleData()
#' ##Check the results to see if you have any interested SNPs that you want to plot on.
#' head(results(obj))
#' p <- loopQTLplotSnpWindow(
#'     obj,
#'     chrom      = "chr1",
#'     start      = 45900000,
#'     end        = 46300000,
#'     snpID      = "rs11584814",
#'     normal     = TRUE,
#'     covariates = "auto")
loopQTLplotSnpWindow <- function(obj, snpID,
                                 chrom = NULL,
                                 start = NULL,
                                 end = NULL,
                                 windows = NULL,
                                 mirror = TRUE,
                                 minSamplesPerCell = 3,
                                 normal = FALSE,
                                 covariates = NULL,
                                 Kmethod = c("elbow", "BE"),
                                 Kmax = NULL,
                                 limit = NULL,
                                 palette = "RdBu",
                                 metric = c("betaLogp", "beta"),
                                 values = c("normalized", "cpm", "raw"),
                                 axisUnit = c("Mb", "bp")) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(snpID), length(snpID) == 1, nzchar(snpID),
        is.logical(mirror), length(mirror) == 1,
        is.logical(normal), length(normal) == 1, !is.na(normal),
        is.character(palette), length(palette) == 1
    )
    if (!is.null(limit)) {
        stopifnot(is.numeric(limit), length(limit) == 1, limit > 0)
    }
    Kmethod <- match.arg(Kmethod)
    metric <- match.arg(metric)
    values <- match.arg(values)
    axisUnit <- match.arg(axisUnit)

    w <- .resolveWindows(obj, windows)

    ## 2. Look up the SNP row
    i <- match(snpID, w$snpIndex$snpID)
    if (is.na(i)) {
        .stopNoCall("snpID not found in windows file: ", snpID)
    }
    row <- w$snpIndex[i, , drop = FALSE]
    cols <- row$colIndices[[1]]
    if (length(cols) == 0) {
        .stopNoCall(
            "No cells were saved for snpID = ", snpID,
            " (likely on a chrom that had an empty window slice)."
        )
    }
    if (!is.null(chrom)) {
        if (length(chrom) != 1 || !is.character(chrom)) {
            .stopNoCall("`chrom` must be a single character string")
        }
        if (!identical(chrom, row$chrom)) {
            .stopNoCall(
                "Supplied chrom = '", chrom, "' does not match the SNP's chrom ('",
                row$chrom, "') in the windows file"
            )
        }
    }

    res <- as.integer(w$resolution)
    binStart <- as.integer(row$binStart)
    binEnd <- as.integer(row$binEnd)

    ## 3. Optionally crop to user's [start, end]
    if (!is.null(start) || !is.null(end)) {
        if (is.null(start)) start <- binStart
        if (is.null(end)) end <- binEnd + res
        if (!is.numeric(start) || !is.numeric(end) ||
            length(start) != 1 || length(end) != 1 || start >= end) {
            .stopNoCall("`start` and `end` must be single numerics with start < end")
        }
        crop_lo_req <- as.integer(floor(start / res) * res)
        crop_hi_req <- as.integer(floor((end - 1) / res) * res)
        cropLo <- max(crop_lo_req, binStart)
        cropHi <- min(crop_hi_req, binEnd)
        if (cropLo > cropHi) {
            .stopNoCall(
                "Crop range [", start, ", ", end, "] does not overlap the saved ",
                "window [", binStart, ", ", binEnd + res, "] for SNP ", snpID
            )
        }
        if (cropLo != crop_lo_req || cropHi != crop_hi_req) {
            .warningNoCall(sprintf(
                paste0(
                    "Requested [%d, %d] was clipped to [%d, %d] (the saved window ",
                    "for %s covers [%d, %d]; widen `window` in loopQTLnormalize() ",
                    "to enlarge it)."
                ),
                as.integer(start), as.integer(end),
                cropLo, cropHi + res, snpID, binStart, binEnd + res
            ))
        }
    } else {
        cropLo <- binStart
        cropHi <- binEnd
    }

    ## 4. Pull genotype dosage for this SNP
    dose <- .read_snp_dosage(obj@vcfPath, row$chrom, row$pos,
        samplesVec = w$samples
    )

    ## 5. Per-cell regression slope (beta)
    ## Three modes, matching runQTL()'s semantics:
    ##   normal = FALSE, covariates = NULL          -> marginal beta on raw BNBC
    ##   normal = TRUE,  covariates = NULL          -> marginal beta on RINT'd cells
    ##   normal = TRUE,  covariates = "auto"/matrix -> partial beta matching runQTL()
    ## By Frisch-Waugh-Lovell, the partial slope equals the slope from
    ## lm(cell ~ dose + cvrt), so we residualize both cell-value and dose
    ## against [1, cvrt] then take slope = sum(YRes * dRes) / sum(dRes^2).
    ## Select the cell values matrix based on the `values` arg

    if (values %in% c("raw", "cpm") && is.null(w$cellsRaw)) {
        .stopNoCall(
            "This loopQTL object does not contain raw counts. To run ",
            "loopQTLplotSnpWindow(), please re-run loopQTLnormalize() ",
            "with saveWindows set to TRUE."
        )
    }

    cellsMatrix <- switch(values,
        normalized = w$cells,
        raw = w$cellsRaw,
        cpm = {
            if (length(obj@cis) == 0) {
                .stopNoCall(
                    "values = 'cpm' requires per-sample cis depth. ",
                    "Set the `cis` column in the samplesheet at loopQTLsetup() or ",
                    "let setup auto-compute it from .hic."
                )
            }
            libM <- obj@cis[w$samples] / 1000000
            sweep(w$cellsRaw, 1, libM, "/")
        }
    )
    cells_for_snp <- cellsMatrix[, cols, drop = FALSE]
    rownames(cells_for_snp) <- w$samples

    Y <- cells_for_snp
    if (normal) {
        Y <- apply(Y, 2, .rint)
        if (!is.matrix(Y)) Y <- as.matrix(Y)
        if (ncol(Y) != length(cols)) Y <- t(Y)
    }

    cvrtMode <- "none"
    cvrt <- NULL
    if (identical(covariates, "auto")) {
        cvrtMode <- "auto"
        if (length(obj@phenotype) == 0) {
            .stopNoCall("covariates = 'auto' needs obj@phenotype; run loopQTLnormalize() first")
        }
        Y_for_pca <- if (normal) .rintMatrix(obj@phenotype) else obj@phenotype
        Kres <- .pca_for_qtl(Y_for_pca,
            samples = obj@samples,
            Kmethod = Kmethod, Kmax = Kmax
        )
        cvrt <- Kres$PCs[w$samples, , drop = FALSE]
    } else if (!is.null(covariates)) {
        cvrtMode <- "user"
        cvDF <- as.data.frame(covariates, stringsAsFactors = FALSE)
        if (!is.null(rownames(cvDF))) {
            mm <- match(w$samples, rownames(cvDF))
            if (anyNA(mm)) {
                .stopNoCall(
                    "`covariates` is missing rows for samples: ",
                    paste(w$samples[is.na(mm)], collapse = ", ")
                )
            }
            cvDF <- cvDF[mm, , drop = FALSE]
        } else if (nrow(cvDF) != length(w$samples)) {
            .stopNoCall(
                "`covariates` has ", nrow(cvDF), " rows but ",
                length(w$samples), " samples; provide rownames or matched order"
            )
        }
        cvrt <- as.matrix(cvDF)
        rownames(cvrt) <- w$samples
    }

    usable <- !is.na(dose)
    if (!is.null(cvrt)) {
        usable <- usable & stats::complete.cases(cvrt)
    }
    nUsable <- sum(usable)

    betaVec <- rep(NA_real_, ncol(Y))
    pVec <- rep(NA_real_, ncol(Y))
    if (nUsable >= as.integer(minSamplesPerCell)) {
        Yusable <- Y[usable, , drop = FALSE]
        dU <- dose[usable]
        if (is.null(cvrt)) {
            ## Marginal: center then slope = sum(Yc * dc) / sum(dc^2)
            Y_r <- sweep(Yusable, 2, colMeans(Yusable, na.rm = TRUE), FUN = "-")
            dR <- dU - mean(dU)
            df <- nUsable - 2
        } else {
            ## Partial via FWL: residualize Y and d against [1, cvrt[usable, ]]
            X <- cbind(1, cvrt[usable, , drop = FALSE])
            qrX <- qr(X)
            Y_r <- Yusable - X %*% qr.coef(qrX, Yusable)
            dR <- as.numeric(dU - X %*% qr.coef(qrX, dU))
            df <- nUsable - ncol(X) - 1
        }
        dDenominator <- sum(dR^2)
        if (dDenominator > 0 && df > 0) {
            betaVec <- as.numeric(crossprod(Y_r, dR) / dDenominator)
            ## RSS_j = ||Y_r_j||^2 - beta_j^2 * dDenominator  (from partial regression)
            rss <- colSums(Y_r^2) - betaVec^2 * dDenominator
            rss <- pmax(rss, 0) # numerical safety
            sigma2 <- rss / df
            seVec <- sqrt(sigma2 / dDenominator)
            tVec <- betaVec / seVec
            pVec <- 2 * stats::pt(-abs(tVec), df = df)
        }
        betaVec[!is.finite(betaVec)] <- NA_real_
        pVec[!is.finite(pVec)] <- NA_real_
    }

    ## Build the display value per cell based on `metric`
    if (metric == "betaLogp") {
        ## Guard against p==0 (perfect fit); cap the -log10 to a plottable ceiling
        pSafe <- pmax(pVec, .Machine$double.xmin)
        displayV <- betaVec * (-log10(pSafe))
        displayV[!is.finite(displayV)] <- NA_real_
    } else {
        displayV <- betaVec
    }

    ## 6. Lay r values onto the bins x bins grid
    coord <- w$cellCoordinates[cols, , drop = FALSE]
    bins <- seq.int(cropLo, cropHi, by = res)
    nb <- length(bins)
    inCrop <- coord$binI >= cropLo & coord$binI <= cropHi &
        coord$binJ >= cropLo & coord$binJ <= cropHi
    ri <- match(coord$binI[inCrop], bins)
    rj <- match(coord$binJ[inCrop], bins)
    ## Shift tile centers by res/2 so each tile visually spans [binStart,
    ## binStart + res] and axis ticks at binStart align with the LEFT
    ## edge of the tile -- matching IGV / cooler / Juicebox bp conventions.
    half <- res / 2
    value <- displayV[inCrop]
    binI <- bins[ri] + half
    binJ <- bins[rj] + half
    if (mirror) {
        value <- c(value, value)
        bin_i_original <- binI
        binI <- c(binI, binJ)
        binJ <- c(binJ, bin_i_original)
        rm(bin_i_original)
    }
    long <- data.frame(
        binI = binI, binJ = binJ, value = value,
        stringsAsFactors = FALSE
    )
    long <- long[!duplicated(paste(long$binI, long$binJ, sep = "_")), ,
        drop = FALSE
    ]

    ## 7. Color scale (RdBu divergent on `value`)
    ## Value scale depends on `metric`; both are unbounded. Default limit is
    ## data-driven; pass a fixed `limit` for cross-SNP comparability.
    if (is.null(limit)) {
        vfin <- long$value[is.finite(long$value)]
        effLimit <- if (length(vfin) > 0) max(abs(vfin)) else 1
        if (!is.finite(effLimit) || effLimit <= 0) effLimit <- 1
    } else {
        effLimit <- limit
    }

    ## 8. Plot
    axisFormat <- if (axisUnit == "Mb") {
        function(x) sprintf("%.2f", x / 1000000)
    } else {
        function(x) format(x, big.mark = ",", scientific = FALSE)
    }
    axis_label_suffix <- if (axisUnit == "Mb") "Mb" else "bp"

    modeStr <- sprintf(
        "%s%s",
        if (cvrtMode == "none") "marginal beta" else "partial beta",
        if (normal) ", RINT" else ""
    )
    if (cvrtMode == "auto") {
        modeStr <- paste0(modeStr, sprintf(" | K=%d auto PCs", ncol(cvrt)))
    }
    if (cvrtMode == "user") {
        modeStr <- paste0(modeStr, sprintf(" | %d user cvrt", ncol(cvrt)))
    }
    modeStr <- paste0(modeStr, " | ", values, " values")

    metricLabel <- switch(metric,
        betaLogp = expression(beta %*% -log[10](p)),
        beta = expression(beta)
    )
    metric_title_str <- switch(metric,
        betaLogp = "beta × -log10(p)",
        beta = "beta"
    )

    title <- sprintf(
        "Per-pixel %s (contact ~ dosage) | %s",
        metric_title_str, snpID
    )
    subtitle <- sprintf(
        "%s:%d  REF=%s ALT=%s | window %s bp | %d samples (%d usable) | %s | limits ±%.3g",
        row$chrom, row$pos, row$ref, row$alt,
        format(cropHi - cropLo + res, big.mark = ","),
        length(w$samples), nUsable, modeStr, effLimit
    )

    diagonalStart <- cropLo
    diagonalEnd <- cropHi + res

    ggplot2::ggplot(
        long,
        ggplot2::aes(.data$binI, .data$binJ, fill = .data$value)
    ) +
        ggplot2::geom_raster(interpolate = FALSE) +
        ggplot2::annotate("segment",
            x = diagonalStart, y = diagonalStart,
            xend = diagonalEnd, yend = diagonalEnd,
            color = "black", linewidth = 0.5, alpha = 0.8
        ) +
        ggplot2::scale_fill_distiller(
            palette   = palette,
            direction = -1, # so + beta = red, - beta = blue
            limits    = c(-effLimit, effLimit),
            oob       = scales::squish,
            na.value  = "grey95"
        ) +
        ggplot2::scale_x_continuous(labels = axisFormat, expand = c(0, 0)) +
        ggplot2::scale_y_reverse(labels = axisFormat, expand = c(0, 0)) +
        ggplot2::coord_fixed() +
        ggplot2::labs(
            title = title, subtitle = subtitle,
            x = sprintf("binI (%s, %s)", axis_label_suffix, row$chrom),
            y = sprintf("binJ (%s, %s)", axis_label_suffix, row$chrom),
            fill = metricLabel
        ) +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(
            plot.subtitle    = ggplot2::element_text(size = 9),
            panel.grid       = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            axis.line        = ggplot2::element_line(linewidth = 0.3)
        )
}

## Internal helpers

## Resolve the `windows` arg into a loaded list, with the necessary structure.
.resolveWindows <- function(obj, windows) {
    ## If the windows argument is null, look inside the loopQTL object
    if (is.null(windows)) {
        windows <- obj@metadata$windowsOut
        if (is.null(windows)) {
            .stopNoCall(
                "The input loopQTL object does not contain a windows file ",
                "path. Please run loopQTLnormalize() with saveWindows set ",
                "to TRUE and windowsOut set to an output file path, or ",
                "specify a path to a windows file via the windows argument."
            )
        }
    }
    if (is.character(windows) && length(windows) == 1) {
        if (!file.exists(windows)) {
            .stopNoCall("Windows file does not exist: ", windows)
        }
        return(readRDS(windows))
    }
    if (is.list(windows) &&
        all(c("cells", "cellCoordinates", "snpIndex", "samples", "resolution") %in%
            names(windows))) {
        return(windows)
    }
    .stopNoCall(
        "`windows` must be NULL, a path to a windows .rds, or a loaded list ",
        "with cells/cellCoordinates/snpIndex/samples/resolution fields."
    )
}

## Targeted single-SNP dosage read from a tabix VCF. Returns a numeric vector
## of length(samplesVec) in that order; NA where missing.
.read_snp_dosage <- function(vcfPath, chrom, pos, samplesVec) {
    tab <- Rsamtools::TabixFile(vcfPath)
    vcfChroms <- Rsamtools::headerTabix(tab)$seqnames
    vcfChromNormalizerFunc <- .make_chrom_normalizer(vcfChroms)
    gr <- GenomicRanges::GRanges(
        seqnames = vcfChromNormalizerFunc(chrom),
        ranges   = IRanges::IRanges(start = as.integer(pos), end = as.integer(pos))
    )
    vcf <- VariantAnnotation::readVcf(
        vcfPath,
        genome = "unknown",
        param = VariantAnnotation::ScanVcfParam(
            which = gr, geno = "GT",
            info = NA_character_
        )
    )
    if (length(vcf) == 0) {
        .stopNoCall("VCF has no record at ", chrom, ":", pos)
    }
    ## If multi-allelic, take the first row (matching how .gt_to_dosage handles it)
    gt <- VariantAnnotation::geno(vcf)$GT
    m <- match(samplesVec, colnames(gt))
    if (anyNA(m)) {
        .stopNoCall(
            "VCF missing samples present in windows file: ",
            paste(samplesVec[is.na(m)], collapse = ", ")
        )
    }
    gt <- gt[1, m, drop = FALSE]
    dose <- .gt_to_dosage(matrix(as.character(gt), nrow = 1))[1, ]
    setNames(as.numeric(dose), samplesVec)
}
