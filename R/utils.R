## Internal function to display a warning without showing the function call
.warningNoCall <- function(...) {
    warning(..., call. = FALSE)
}

## Internal function to stop execution without showing the function call
.stopNoCall <- function(...) {
    stop(..., call. = FALSE)
}

## Internal function to display an error like stop() would, but without stopping
## execution or showing the function call
.displayError <- function(...) {
    try(stop(..., call. = FALSE))
}

## Internal function to check if an object is a single NA (since R 4.3 changed
## the behavior of if(!is.na(...))). Based on
## https://github.com/GeoBosh/gbutils/blob/16a311ce2e40813911461c6084393c1915568c5b/R/isNA.R
## (license: GPL v2 or later)
.isSingleNA <- function(x) {
    return(is.atomic(x) && length(x) == 1 && is.na(x))
}

## Bin a 1-based BEDPE-style anchor (start, end) to its analysis-bin start (0-based),
## matching strawr's binning. straw returns bin starts as the lower coordinate of
## the bin in 0-based half-open BP space, so bin = floor(midpoint / resolution) * resolution.
.binAnchor <- function(start, end, resolution) {
    mid <- (as.numeric(start) + as.numeric(end)) / 2
    as.integer(floor(mid / resolution) * resolution)
}

## Detect strawr's chromosome naming style ("1" vs "chr1") from the .hic header.
## Returns a function that maps a user-supplied chrom string (either style) to the
## style this .hic expects.
.make_chrom_normalizer <- function(hicChroms) {
    uses_chr_prefix <- any(grepl("^chr", hicChroms))
    function(x) {
        if (uses_chr_prefix) {
            ifelse(grepl("^chr", x), x, paste0("chr", x))
        } else {
            sub("^chr", "", x)
        }
    }
}

## RINT with Blom offset, NA-aware. Operates on a single vector (one loop's
## values across samples).
.rint <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    n <- sum(!is.na(r))
    if (n < 2) {
        return(rep(NA_real_, length(x)))
    }
    qnorm((r - 0.5) / n)
}

## Apply .rint() row-wise to a (loops x samples) matrix. Rows whose values are
## all-NA, fewer than 2 non-NA, or have zero variance are returned as all-NA
## (RINT on a constant produces NaN/Inf, which would poison downstream
## regressions). dimnames are preserved.
.rintMatrix <- function(mat) {
    if (length(mat) == 0) {
        return(mat)
    }
    out <- matrix(NA_real_,
        nrow = nrow(mat), ncol = ncol(mat),
        dimnames = dimnames(mat)
    )
    rowSd <- apply(mat, 1, stats::sd, na.rm = TRUE)
    ok <- is.finite(rowSd) & rowSd > 0
    if (any(ok)) {
        good <- which(ok)
        out[good, ] <- t(apply(mat[good, , drop = FALSE], 1, .rint))
    }
    out
}

## Internal: compute per-SNP window cells for saveWindows in loopQTLnormalize().
#
## For each VCF SNP that overlaps any consensus loop anchor (within `snpWindow`
## bp flank), build the SNP's "window box" - a square (binLo..binHi)^2 region
## centered on the bounding interval of the SNP position and every partner
## anchor of every loop the SNP participates in, expanded symmetrically so its
## span is at least `window` bp. Returns the per-SNP metadata plus per-chrom
## deduplicated cell-coordinate tables that the BNBC per-chrom worker will
## index after .normalizeBNBC() emits its wide table.
#
## Returns:
##   list(
##     snpsDF = data.frame(snpID, chromOrig, chromVcf, pos, ref, alt,
##                          binLo, binHi, n_loops_anchored),
##     cells_by_chrom = named list per chrom; each is data.frame(binI, binJ)
##                      of upper-triangle cells (bp), deduplicated across SNPs
##                      on that chrom,
##     snp_cells_local = list (length nrow(snpsDF)) of integer vectors of
##                      within-chrom row indices into cells_by_chrom[[chrom]]
##   )
## Returns NULL if no eligible SNPs are found (caller short-circuits).
.compute_snp_windows <- function(vcfPath, anchorA, anchorB,
                                 snpWindow, window, resolution) {
    if (!file.exists(vcfPath)) {
        .stopNoCall("VCF not found: ", vcfPath)
    }

    ## VCF chrom-naming normalizer
    tab <- Rsamtools::TabixFile(vcfPath)
    vcfChroms <- Rsamtools::headerTabix(tab)$seqnames
    if (length(vcfChroms) == 0) {
        .stopNoCall("Could not read chromosome names from VCF tabix index")
    }
    vcfChromNormalizerFunc <- .make_chrom_normalizer(vcfChroms)

    ## Build SNP target ranges = anchors padded by snpWindow
    pad <- function(gr) {
        GenomicRanges::GRanges(
            seqnames = vcfChromNormalizerFunc(as.character(GenomicRanges::seqnames(gr))),
            ranges = IRanges::IRanges(
                start = pmax(1, GenomicRanges::start(gr) - snpWindow),
                end   = GenomicRanges::end(gr) + snpWindow
            )
        )
    }
    targets <- GenomicRanges::reduce(c(pad(anchorA), pad(anchorB)))
    targets <- targets[as.character(GenomicRanges::seqnames(targets)) %in% vcfChroms]
    if (length(targets) == 0) {
        return(NULL)
    }

    message(
        "saveWindows: tabix-reading SNPs at ",
        length(targets), " anchor regions (snpWindow = ", snpWindow, " bp) ..."
    )
    param <- VariantAnnotation::ScanVcfParam(
        which = targets,
        geno = NA_character_,
        info = NA_character_
    )
    vcf <- VariantAnnotation::readVcf(vcfPath, genome = "unknown", param = param)
    if (length(vcf) == 0) {
        warning(
            "saveWindows: no SNPs in any anchor ± ",
            snpWindow, " bp window"
        )
        return(NULL)
    }
    rr <- SummarizedExperiment::rowRanges(vcf)
    snp_chrom_vcf <- as.character(GenomicRanges::seqnames(rr))
    snpPos <- as.integer(GenomicRanges::start(rr))
    snpRef <- as.character(VariantAnnotation::ref(vcf))
    altList <- VariantAnnotation::alt(vcf)
    snpAlt <- vapply(
        altList,
        function(a) paste(as.character(a), collapse = ","),
        character(1)
    )
    snpID <- rownames(vcf)
    if (is.null(snpID) || any(snpID == "" | is.na(snpID))) {
        snpID <- sprintf("%s:%d_%s/%s", snp_chrom_vcf, snpPos, snpRef, snpAlt)
    }

    ## Map SNPs to consensus loops (within snpWindow flank)  # Convert anchors to vcf-style seqnames for overlap with the SNP GRanges.
    to_vcf_gr <- function(gr) {
        GenomicRanges::GRanges(
            seqnames = vcfChromNormalizerFunc(as.character(GenomicRanges::seqnames(gr))),
            ranges = IRanges::IRanges(
                start = pmax(1, GenomicRanges::start(gr) - snpWindow),
                end   = GenomicRanges::end(gr) + snpWindow
            )
        )
    }
    padded_A_vcf <- to_vcf_gr(anchorA)
    padded_B_vcf <- to_vcf_gr(anchorB)

    snpGR <- GenomicRanges::GRanges(
        seqnames = snp_chrom_vcf,
        ranges   = IRanges::IRanges(start = snpPos, end = snpPos)
    )

    hitsA <- GenomicRanges::findOverlaps(snpGR, padded_A_vcf)
    hitsB <- GenomicRanges::findOverlaps(snpGR, padded_B_vcf)
    in_loops_by_snp <- vector("list", length(snpGR))
    if (length(hitsA) > 0) {
        spl <- split(S4Vectors::subjectHits(hitsA), S4Vectors::queryHits(hitsA))
        for (s in names(spl)) {
            in_loops_by_snp[[as.integer(s)]] <- c(in_loops_by_snp[[as.integer(s)]], spl[[s]])
        }
    }
    if (length(hitsB) > 0) {
        spl <- split(S4Vectors::subjectHits(hitsB), S4Vectors::queryHits(hitsB))
        for (s in names(spl)) {
            in_loops_by_snp[[as.integer(s)]] <- c(in_loops_by_snp[[as.integer(s)]], spl[[s]])
        }
    }
    keepSNPs <- vapply(in_loops_by_snp, function(x) length(x) > 0, logical(1))
    if (!any(keepSNPs)) {
        .warningNoCall(
            "saveWindows: no SNPs overlap any consensus loop anchor"
        )
        return(NULL)
    }
    snp_idx_keep <- which(keepSNPs)
    in_loops_by_snp <- in_loops_by_snp[snp_idx_keep]

    snpID <- snpID[snp_idx_keep]
    snp_chrom_vcf <- snp_chrom_vcf[snp_idx_keep]
    snpPos <- snpPos[snp_idx_keep]
    snpRef <- snpRef[snp_idx_keep]
    snpAlt <- snpAlt[snp_idx_keep]

    ## Native (anchor-style) chrom for output. Inverse of vcfChromNormalizerFunc is sample-specific;
    ## use anchor seqnames at the SNP's chrom directly.
    anchor_chrom_native <- as.character(GenomicRanges::seqnames(anchorA))
    anchor_chrom_vcf <- vcfChromNormalizerFunc(anchor_chrom_native)
    snp_chrom_orig <- vapply(snp_chrom_vcf, function(cVcf) {
        m <- which(anchor_chrom_vcf == cVcf)[1]
        if (is.na(m)) cVcf else anchor_chrom_native[m]
    }, character(1))

    ## Per-SNP bounding interval -> binLo/binHi
    anchorA_start_native <- GenomicRanges::start(anchorA)
    anchorA_end_native <- GenomicRanges::end(anchorA)
    anchorB_start_native <- GenomicRanges::start(anchorB)
    anchorB_end_native <- GenomicRanges::end(anchorB)

    nSnp <- length(snpID)
    binLo <- integer(nSnp)
    binHi <- integer(nSnp)
    nLoops <- integer(nSnp)

    for (i in seq_len(nSnp)) {
        li <- unique(in_loops_by_snp[[i]])
        nLoops[i] <- length(li)
        ## Bounding interval: SNP pos U all involved anchors (A and B)
        starts <- c(snpPos[i], anchorA_start_native[li], anchorB_start_native[li])
        ends <- c(snpPos[i], anchorA_end_native[li], anchorB_end_native[li])
        bLo <- min(starts)
        bHi <- max(ends)
        span <- bHi - bLo
        if (span < window) {
            center <- (bLo + bHi) / 2
            half <- window / 2
            bLo <- center - half
            bHi <- center + half
        }
        binLo[i] <- as.integer(floor(bLo / resolution) * resolution)
        binHi[i] <- as.integer(floor(bHi / resolution) * resolution)
        if (binLo[i] < 0) binLo[i] <- 0
    }

    snpsDF <- data.frame(
        snpID = snpID,
        chromOrig = snp_chrom_orig,
        chromVcf = snp_chrom_vcf,
        pos = snpPos,
        ref = snpRef,
        alt = snpAlt,
        binLo = binLo,
        binHi = binHi,
        n_loops_anchored = nLoops,
        stringsAsFactors = FALSE
    )

    ## Per-chrom dedup of cells  # cell = (binI, binJ) with binLo <= binI <= binJ <= binHi (upper tri).
    ## Group SNPs by their native (anchor-style) chrom, since BNBC emits in that style.
    cells_by_chrom <- list()
    snp_cells_local <- vector("list", nSnp)

    chromGroups <- split(seq_len(nSnp), snpsDF$chromOrig)
    for (chr in names(chromGroups)) {
        idxs <- chromGroups[[chr]]
        ## Build a master set of cells, dedupe via paste-key
        cellChunks <- vector("list", length(idxs))
        for (k in seq_along(idxs)) {
            s <- idxs[k]
            bins <- seq.int(snpsDF$binLo[s], snpsDF$binHi[s], by = resolution)
            ng <- length(bins)
            if (ng == 0) {
                cellChunks[[k]] <- data.frame(
                    binI = integer(0),
                    binJ = integer(0)
                )
                next
            }
            ## All (binI, binJ) with binI <= binJ, upper triangle inclusive
            ii <- rep(bins, times = ng:1)
            jj <- unlist(lapply(seq_len(ng), function(z) bins[z:ng]))
            cellChunks[[k]] <- data.frame(
                binI = as.integer(ii),
                binJ = as.integer(jj)
            )
        }
        combined <- do.call(rbind, cellChunks)
        keys <- paste(combined$binI, combined$binJ, sep = "_")
        uniqIdx <- !duplicated(keys)
        chromCells <- combined[uniqIdx, , drop = FALSE]
        rownames(chromCells) <- NULL
        chrom_cells_keys <- keys[uniqIdx]
        cells_by_chrom[[chr]] <- chromCells

        ## For each SNP, look up its cells' local indices
        for (k in seq_along(idxs)) {
            s <- idxs[k]
            ch <- cellChunks[[k]]
            if (nrow(ch) == 0) {
                snp_cells_local[[s]] <- integer(0)
                next
            }
            ck <- paste(ch$binI, ch$binJ, sep = "_")
            snp_cells_local[[s]] <- match(ck, chrom_cells_keys)
        }
    }

    list(
        snpsDF = snpsDF,
        cells_by_chrom = cells_by_chrom,
        snp_cells_local = snp_cells_local
    )
}

## Read a single BEDPE-like loop file. Tolerant of "#chr" headers and of
## files with or without a `loopID` column.
.read_loops_bedpe <- function(path) {
    if (!file.exists(path)) {
        .stopNoCall("Loop file not found: ", path)
    }
    hdr <- readLines(path, n = 1, warn = FALSE)
    hasHeader <- grepl("^#?chr", hdr, ignore.case = TRUE) &&
        grepl("\\bstart", hdr, ignore.case = TRUE)
    df <- utils::read.table(
        path,
        header = hasHeader, sep = "\t",
        comment.char = if (hasHeader) "" else "#",
        stringsAsFactors = FALSE, check.names = FALSE
    )
    if (!hasHeader) {
        ## assume positional: chrA startA endA chrB startB endB [loopID]
        if (ncol(df) < 6) {
            .stopNoCall("Loop file ", path, " must have at least 6 columns")
        }
        baseNames <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
        if (ncol(df) == 6) {
            names(df) <- baseNames
        } else {
            names(df) <- c(baseNames, paste0("V", seq_len(ncol(df) - 6) + 6))
        }
    } else {
        nm <- names(df)
        if (nm[1] %in% c("#chr1", "#chrA")) nm[1] <- "chrA"
        ## canonicalize common variants
        nm <- sub("^chr1$", "chrA", nm)
        nm <- sub("^x1$", "startA", nm, ignore.case = TRUE)
        nm <- sub("^y1$", "endA", nm, ignore.case = TRUE)
        nm <- sub("^chr2$", "chrB", nm)
        nm <- sub("^x2$", "startB", nm, ignore.case = TRUE)
        nm <- sub("^y2$", "endB", nm, ignore.case = TRUE)
        names(df) <- nm
    }
    need <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
    miss <- setdiff(need, names(df))
    if (length(miss) > 0) {
        .stopNoCall(
            "Loop file ", path, " missing required columns: ",
            paste(miss, collapse = ", ")
        )
    }
    df$chrA <- as.character(df$chrA)
    df$chrB <- as.character(df$chrB)
    df$startA <- as.integer(df$startA)
    df$endA <- as.integer(df$endA)
    df$startB <- as.integer(df$startB)
    df$endB <- as.integer(df$endB)
    df
}

#' Load the bundled example loopQTL object
#'
#' Convenience helper for `@examples` blocks and vignette code that need a
#' fully-populated `loopQTL` object without re-running the full pipeline.
#' Loads `example_loopQTL_obj.rda` shipped in `inst/extdata/` and rewrites
#' the file paths baked into the object so they resolve against the
#' currently installed package location.
#'
#' @return A [loopQTL-class] object with `phenotype`, `results`,
#'   and `metadata` populated.
#' @export
#' @examples
#' obj <- loopQTLExampleData()
#' obj
loopQTLExampleData <- function() {
    extdataDir <- system.file("extdata", package = "loopQTL")
    rdaPath    <- file.path(extdataDir, "example_loopQTL_obj.rda")
    e <- new.env()
    load(rdaPath, envir = e)
    obj <- e$obj
    ## Rewrite absolute paths baked into the .rda so downstream functions
    ## resolve files against the current installed location, not the machine
    ## where the .rda was built.
    obj@metadata$windowsOut <- file.path(
        extdataDir, "example_loopQTL_snp_windows.rds"
    )
    if (!is.null(obj@metadata$hicPaths)) {
        obj@metadata$hicPaths <- setNames(
            file.path(extdataDir, "example_hic",
                      paste0(names(obj@metadata$hicPaths), ".hic")),
            names(obj@metadata$hicPaths)
        )
    }
    obj@vcfPath <- file.path(extdataDir, "example.vcf.gz")
    obj
}
                                
