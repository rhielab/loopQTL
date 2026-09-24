#' Cross-normalization on Hi-C signal
#'
#' This function normalizes contact strength across samples so that
#' per-loop phenotypes are directly comparable. For each chromosome it pulls
#' the whole-genome matrix per sample, log-CPM adjusts, then does
#' band-wise quantile normalization across samples per genomic-distance band.
#' Normalized values at the consensus-loop bin-pairs are written to
#' `phenotype(obj)`.
#'
#' Please note that you need to set `saveWindows` to TRUE to additionally output
#' a per-SNP contact-window matrix as .rds file. This file is needed by the
#' `loopQTLplotSnpWindow()` and `loopQTLplotLocus()` downstream functions. 
#' `saveWindows` does not change `phenotype(obj)` or affect the `runQTL()` step.
#'
#' @param obj A loopQTL object created from `loopQTLsetup()`.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] for parallelization
#'   over chromosomes. Default [BiocParallel::SerialParam()].
#' @param saveWindows Logical. If `TRUE`, also save a local contact submatrix
#'   around each in-loop SNP for downstream functions usage.
#' @param windowsOut Required when `saveWindows = TRUE`. Path to a `.rds`
#'   file the per-SNP windows are written to.
#' @param window Minimum window extent in bp around each in-loop SNP to
#' create `.rds` file when `saveWindows = TRUE`. Default 1000000.
#' @param snpWindow Base-pair window flank added to each consensus anchor when selecting
#'   in-loop SNPs. Default 0.
#'
#' @return The loopQTL `obj` with `phenotype(obj)` filled.
#' @import data.table
#' @export
#' @examples
#' ## Load the example object. In real analysis, the input
#' ## should be the output of loopQTLsetup().
#' ##Here we just use the example object.
#' obj <- loopQTLExampleData()
#' ## Example code to run loopQTLnormalize(). By default it will NOT create
#' ## the .rds contact submatrix file that downstream functions need.
#' obj <- loopQTLnormalize(
#'     obj,
#'     BPPARAM     = BiocParallel::MulticoreParam(workers = 2),
#'     snpWindow   = 10000)
#' 
#' ## Example code to run loopQTLnormalize() and also create the `.rds`
#' ## contact submatrix file needed by downstream visualization functions.
#'
#' ## First, get a temporary directory in which to save the output file.
#' tempDirectory <- tempdir()
#' obj <- loopQTLnormalize(
#'     obj,
#'     BPPARAM     = BiocParallel::MulticoreParam(workers = 2),
#'     saveWindows = TRUE,
#'     windowsOut  = file.path(tempDirectory, "loop_snp_windows.rds"),
#'     window      = 2000000,
#'     snpWindow   = 10000)
loopQTLnormalize <- function(obj,
                             BPPARAM = BiocParallel::SerialParam(),
                             saveWindows = FALSE,
                             windowsOut = NULL,
                             window = 1000000,
                             snpWindow = 0) {
    stopifnot(methods::is(obj, "loopQTL"))

    ## saveWindows validation: path required, VCF must be set
    if (saveWindows) {
        if (is.null(windowsOut) || !is.character(windowsOut) ||
            length(windowsOut) != 1 || !nzchar(windowsOut)) {
            .stopNoCall("saveWindows = TRUE requires `windowsOut`: a single .rds path")
        }
        stopifnot(
            is.numeric(window), length(window) == 1, window > 0,
            is.numeric(snpWindow), length(snpWindow) == 1, snpWindow >= 0
        )
        window <- as.numeric(window)
        snpWindow <- as.integer(snpWindow)
        ## Best-effort writable check (touch the file)
        okTouch <- tryCatch(
            {
                con <- file(windowsOut, open = "wb")
                close(con)
                file.remove(windowsOut)
                TRUE
            },
            error = function(e) FALSE,
            warning = function(w) FALSE
        )
        if (!okTouch) {
            .stopNoCall("Cannot write to windowsOut path: ", windowsOut)
        }
        if (is.null(obj@vcfPath) || is.na(obj@vcfPath) || !nzchar(obj@vcfPath)) {
            .stopNoCall("obj@vcfPath is not set; cannot read SNPs for windows")
        }
    } else if (!is.null(windowsOut)) {
        .warningNoCall("`windowsOut` is ignored when saveWindows = FALSE")
    }

    ## BNBC's logCPM depth correction requires cis
    if (length(obj@cis) == 0) {
        .stopNoCall(
            "Normalize requires per-sample cis depth in obj@cis, but it is ",
            "empty. Options:\n",
            "  1. Include a `cis` column in the samplesheet at loopQTLsetup(), or\n",
            "  2. Re-run loopQTLsetup() (cis is auto-computed from .hic when no\n",
            "     `cis` column is provided), or\n",
            "  3. Manually attach: obj@cis <- setNames(<numeric>, obj@samples)"
        )
    }
    if (length(obj@cis) != length(obj@samples) ||
        !identical(names(obj@cis), obj@samples)) {
        .stopNoCall(
            "obj@cis must be a named numeric vector with one entry per sample, ",
            "names matching obj@samples"
        )
    }

    samplesVec <- obj@samples
    resolution <- obj@resolution
    loopsGInteractions <- obj@loops

    ## Per-sample per-chrom hic paths from metadata (we re-derive from setup invariants)
    if (is.null(obj@metadata$hicPaths)) {
        .stopNoCall("obj@metadata$hicPaths missing; loopQTLsetup must record the per-sample .hic paths")
    }
    hicPaths <- obj@metadata$hicPaths

    ## Consensus loops keyed by chrom + canonical (binA <= binB)
    anchorA <- InteractionSet::anchors(loopsGInteractions, type = "first")
    anchorB <- InteractionSet::anchors(loopsGInteractions, type = "second")
    loopChr <- as.character(GenomicRanges::seqnames(anchorA))
    ## Floor each anchor's midpoint to the bin containing it - straw returns
    ## bin-based contacts, so we look up at bin granularity even when the consensus
    ## key is exact-coord.
    loopBinA <- .binAnchor(
        GenomicRanges::start(anchorA),
        GenomicRanges::end(anchorA), resolution
    )
    loopBinB <- .binAnchor(
        GenomicRanges::start(anchorB),
        GenomicRanges::end(anchorB), resolution
    )
    swap <- loopBinA > loopBinB
    if (any(swap)) {
        tmp <- loopBinA[swap]
        loopBinA[swap] <- loopBinB[swap]
        loopBinB[swap] <- tmp
    }
    loopKey <- paste(loopChr, loopBinA, loopBinB, sep = ":")

    chroms_to_run <- unique(loopChr)
    loop_idx_by_chr <- split(seq_along(loopChr), loopChr)

    ## Make a list of functions that adjust chromosome names to match each .hic
    ## file (1 vs chr1). This allows different input HiC files to use different
    ## naming conventions.
    chrom_norm_BNBC <- lapply(hicPaths, function(p) {
        .make_chrom_normalizer(strawr::readHicChroms(p)$name)
    })

    ## saveWindows: precompute SNP -> (chrom, cells to slice)
    windowsCtx <- NULL
    if (saveWindows) {
        windowsCtx <- .compute_snp_windows(
            vcfPath = obj@vcfPath,
            anchorA = anchorA,
            anchorB = anchorB,
            snpWindow = snpWindow,
            window = window,
            resolution = resolution
        )
        if (is.null(windowsCtx)) {
            .warningNoCall(
                "saveWindows: nothing to save; ",
                "continuing without window output"
            )
            saveWindows <- FALSE
        } else {
            ## Sanity check: every SNP's chrom must be in our chroms_to_run
            winChroms <- names(windowsCtx$cells_by_chrom)
            missingChr <- setdiff(winChroms, chroms_to_run)
            if (length(missingChr) > 0) {
                .stopNoCall(
                    "saveWindows: SNP chrom(s) not in normalize's chroms_to_run: ",
                    paste(missingChr, collapse = ", ")
                )
            }
            n_total_cells <- sum(vapply(
                windowsCtx$cells_by_chrom,
                nrow, integer(1)
            ))
            estMb <- n_total_cells * length(samplesVec) * 8 / 1024^2
            message(sprintf(
                "saveWindows: %d SNP(s), %d unique cells across %d chrom(s); estimated ~%.1f MB in memory (8 bytes x cells x samples)",
                nrow(windowsCtx$snpsDF),
                n_total_cells, length(winChroms), estMb
            ))
        }
    }

    ## worker: normalize one chromosome, return (consensusRows, values matrix)
    per_chr_worker <- function(chr) {
        rows <- loop_idx_by_chr[[chr]]
        wantKeys <- loopKey[rows]
        nSamples <- length(samplesVec)

        ## 1. Pull whole-genome (i.e. whole-chrom) straw per sample
        sparseList <- lapply(seq_len(nSamples), function(si) {
            chr_in_hic <- chrom_norm_BNBC[[si]](chr)
            dat <- tryCatch(
                strawr::straw(
                    "NONE", hicPaths[si], chr_in_hic, chr_in_hic,
                    "BP", resolution
                ),
                error = function(e) NULL
            )
            if (is.null(dat) || nrow(dat) == 0) {
                return(data.frame(
                    x = integer(0), y = integer(0),
                    counts = numeric(0)
                ))
            }
            ## canonicalize x <= y
            sw <- dat$x > dat$y
            if (any(sw)) {
                t1 <- dat$x[sw]
                dat$x[sw] <- dat$y[sw]
                dat$y[sw] <- t1
            }
            dat
        })

        ## 2. BNBC normalization. Filter emission to consensus (binA, binB) pairs
        ## on this chrom; QN math is still computed over the full band, only the
        ## final triplet emission is restricted, which keeps the dcast wide-pivot
        ## from blowing past R's integer.max on large chroms with many samples.
        ## When saveWindows is on, also include the SNP window cells in the keep
        ## set so a single BNBC pass emits both consensus + window cells.
        keep_pairs_chr <- cbind(loopBinA[rows], loopBinB[rows])
        window_cells_chr <- NULL
        if (saveWindows && !is.null(windowsCtx) &&
            !is.null(windowsCtx$cells_by_chrom[[chr]])) {
            window_cells_chr <- windowsCtx$cells_by_chrom[[chr]]
            combined <- rbind(
                keep_pairs_chr,
                cbind(window_cells_chr$binI, window_cells_chr$binJ)
            )
            keep_pairs_pass <- unique(combined)
        } else {
            keep_pairs_pass <- keep_pairs_chr
        }
        cisArg <- obj@cis
        normVals <- .normalizeBNBC(sparseList, chr,
            keepPairs = keep_pairs_pass,
            cis = cisArg
        )
        ## normVals: data.table(chr, region1, region2, IF_1, ..., IF_n)

        ## 3. Index consensus loops on this chrom into the normalized table
        normKeys <- paste(normVals$chr, normVals$region1, normVals$region2, sep = ":")
        m <- match(wantKeys, normKeys)
        ifCols <- grep("^IF_", names(normVals), value = TRUE)
        chromMat <- matrix(0,
            nrow = length(rows), ncol = nSamples,
            dimnames = list(NULL, samplesVec)
        )
        for (j in seq_len(nSamples)) {
            v <- normVals[[ifCols[j]]]
            vals <- ifelse(is.na(m), 0, v[m])
            vals[is.na(vals)] <- 0
            chromMat[, j] <- vals
        }

        ## Index window cells into two (samples x n_cells) matrices:
        ##    - windowsMat:     post-BNBC normalized values (from normVals)
        ##    - windowsRawMat: raw straw counts at same cells (from sparseList)
        windowsMat <- NULL
        windowsRawMat <- NULL
        if (!is.null(window_cells_chr)) {
            winKeys <- paste(chr, window_cells_chr$binI, window_cells_chr$binJ,
                sep = ":"
            )
            mw <- match(winKeys, normKeys)
            windowsMat <- matrix(0,
                nrow = nSamples, ncol = nrow(window_cells_chr),
                dimnames = list(samplesVec, NULL)
            )
            for (j in seq_len(nSamples)) {
                v <- normVals[[ifCols[j]]]
                vals <- ifelse(is.na(mw), 0, v[mw])
                vals[is.na(vals)] <- 0
                windowsMat[j, ] <- vals
            }

            ## Raw counts: match window cells against each sample's straw sparse output
            winKeyLocal <- paste(window_cells_chr$binI, window_cells_chr$binJ,
                sep = ":"
            )
            windowsRawMat <- matrix(0,
                nrow = nSamples,
                ncol = nrow(window_cells_chr),
                dimnames = list(samplesVec, NULL)
            )
            for (s in seq_len(nSamples)) {
                dat <- sparseList[[s]]
                if (is.null(dat) || nrow(dat) == 0) next
                datKey <- paste(dat$x, dat$y, sep = ":")
                m <- match(winKeyLocal, datKey)
                windowsRawMat[s, ] <- ifelse(is.na(m), 0, dat$counts[m])
            }
        }
        list(
            rows = rows, mat = chromMat,
            windowsMat = windowsMat,
            windowsRawMat = windowsRawMat,
            window_cells_chr = window_cells_chr
        )
    }

    chromResults <- BiocParallel::bplapply(chroms_to_run, 
                                           per_chr_worker, BPPARAM = BPPARAM)
    names(chromResults) <- chroms_to_run

    ## Assemble per-chrom slices into the full (loops x samples) phenotype
    pheno <- matrix(NA_real_,
        nrow = length(loopsGInteractions), ncol = length(samplesVec),
        dimnames = list(NULL, samplesVec)
    )
    for (r in chromResults) {
        pheno[r$rows, ] <- r$mat
    }

    ## obj@phenotype stores the raw values - runQTL() stores RINT'd values
    obj@phenotype <- pheno

    ## Stitch per-chrom window outputs and write .rds
    if (saveWindows && !is.null(windowsCtx)) {
        ## Per-chrom column offsets into the genome-wide cells matrix
        chrOrder <- intersect(chroms_to_run, names(windowsCtx$cells_by_chrom))
        n_cells_per_chr <- vapply(chrOrder, function(chr) {
            r <- chromResults[[chr]]
            if (is.null(r$windowsMat)) 0 else ncol(r$windowsMat)
        }, integer(1))
        if (sum(n_cells_per_chr) == 0) {
            .warningNoCall(
                "saveWindows: every per-chrom ",
                "window slice came back empty; nothing to write"
            )
        } else {
            totalCells <- sum(n_cells_per_chr)
            cells <- matrix(0,
                nrow = length(samplesVec), ncol = totalCells,
                dimnames = list(samplesVec, NULL)
            )
            cellsRaw <- matrix(0,
                nrow = length(samplesVec), ncol = totalCells,
                dimnames = list(samplesVec, NULL)
            )
            cellCoordinates <- data.frame(
                chrom = character(totalCells),
                binI = integer(totalCells),
                binJ = integer(totalCells),
                stringsAsFactors = FALSE
            )
            offset <- 0
            chrom_col_offset <- setNames(integer(length(chrOrder)), chrOrder)
            for (chr in chrOrder) {
                r <- chromResults[[chr]]
                if (is.null(r$windowsMat) || ncol(r$windowsMat) == 0) {
                    chrom_col_offset[chr] <- NA_integer_
                    next
                }
                ncolChr <- ncol(r$windowsMat)
                rng <- seq.int(offset + 1, offset + ncolChr)
                cells[, rng] <- r$windowsMat
                cellsRaw[, rng] <- r$windowsRawMat
                cellCoordinates$chrom[rng] <- chr
                cellCoordinates$binI[rng] <- r$window_cells_chr$binI
                cellCoordinates$binJ[rng] <- r$window_cells_chr$binJ
                chrom_col_offset[chr] <- offset
                offset <- offset + ncolChr
            }

            ## Per-SNP colIndices: local within-chrom indices + that chrom's offset
            snpChrom <- windowsCtx$snpsDF$chromOrig
            nSnp <- nrow(windowsCtx$snpsDF)
            colIndices <- vector("list", nSnp)
            for (i in seq_len(nSnp)) {
                chr <- snpChrom[i]
                off <- chrom_col_offset[chr]
                loc <- windowsCtx$snp_cells_local[[i]]
                colIndices[[i]] <- if (is.na(off) || length(loc) == 0) {
                    integer(0)
                } else {
                    as.integer(loc + off)
                }
            }

            snpIndex <- windowsCtx$snpsDF
            snpIndex$chrom <- snpIndex$chromOrig
            snpIndex$chromOrig <- NULL
            snpIndex$chromVcf <- NULL
            snpIndex$colIndices <- I(colIndices)
            snpIndex <- snpIndex[, c(
                "snpID", "chrom", "pos", "ref", "alt",
                "binLo", "binHi", "n_loops_anchored",
                "colIndices"
            )]

            saveRDS(
                list(
                    cells = cells, # post-BNBC normalized values
                    cellsRaw = cellsRaw, # raw straw counts (integer)
                    cellCoordinates = cellCoordinates,
                    snpIndex = snpIndex,
                    samples = samplesVec,
                    resolution = resolution
                ),
                file = windowsOut,
                compress = "xz"
            )
            message(sprintf(
                "saveWindows: wrote %s (%.1f MB)",
                windowsOut,
                file.info(windowsOut)$size / 1024^2
            ))

            ## Store the absolute path so it can be found even if the working
            ## directory changes
            obj@metadata$windowsOut <- normalizePath(windowsOut,
                winslash = "/", mustWork = TRUE
            )
            obj@metadata$n_windows_saved <- nrow(snpIndex)
            obj@metadata$n_window_cells <- totalCells
            obj@metadata$windowBp <- as.numeric(window)
            obj@metadata$snp_window_bp <- as.integer(snpWindow)
        }
    }

    methods::validObject(obj)
    obj
}

## Internal: BNBC band-wise quantile normalization, one chrom across samples.
## Skips ComBat. Sparse triplet implementation (no dense per-sample matrices).
## QN math is delegated to preprocessCore::normalize.quantiles(), the canonical
## implementation used by BNBC and limma.
#
## Memory profile (chr1 / 50-100 samples / 10 kb): ~25-50 GB peak.
#
## Algorithm:
##   1. For each sample, store straw output as a triplet (i, j, x, d=j-i)
##      data.table in upper-triangle convention. Skip zeros entirely.
##   2. For each genomic distance band d (in bp), find the union of bin-pair
##      positions where any sample has a non-zero contact. Build a small
##      dense Ld x nSamples block with 0-fill for missing samples.
##   3. Quantile-normalize that block across samples via preprocessCore.
##   4. Emit non-zero QN'd cells as new triplets, OPTIONALLY filtered to a
##      user-supplied set of (i, j) keepPairs. QN math is always over the
##      full band; the filter only restricts what we materialize for the
##      downstream wide pivot.
##   5. Pivot the combined triplets to wide form (one IF column per sample).
#
## Why keepPairs matters: without it, the long->wide pivot (dcast) tries to
## build a grid of all genome-wide non-zero (i, j) x all samples, which can
## exceed R's integer.max for large chroms with many samples (e.g. chr3 /
## 50 samples / 10 kb gives ~6.5e9 cells in the CJ). Restricting emission
## to the consensus-loop (i, j) pairs that the caller will use shrinks this
## to a few thousand rows per chrom.
#
## Correctness: bin-pairs that are zero in every sample are absent from all
## triplets and stay absent (their QN'd value would be the row-1 reference,
## which equals mean(per-band minima); since most band positions are all-zero
## for chr-scale Hi-C, that reference is ~0). Bin-pairs non-zero in at least
## one sample are processed in full.
.normalizeBNBC <- function(sparseList, chr, keepPairs = NULL, cis = NULL) {
    nSamples <- length(sparseList)
    emptyOut <- function() {
        out <- data.table::data.table(
            chr = character(0), region1 = integer(0), region2 = integer(0)
        )
        for (k in seq_len(nSamples)) out[[paste0("IF_", k)]] <- numeric(0)
        out
    }
    if (nSamples == 0) {
        return(emptyOut())
    }

    ## Pre-compute a hash set of keepPairs for O(1) lookup per-band.
    ## keepPairs is a 2-column matrix/data.frame of (i, j) in bp.
    keepSet <- NULL
    if (!is.null(keepPairs)) {
        if (NROW(keepPairs) == 0) {
            return(emptyOut())
        }
        pairX <- as.integer(keepPairs[, 1])
        pairY <- as.integer(keepPairs[, 2])
        swap <- pairX > pairY
        if (any(swap)) {
            t1 <- pairX[swap]
            pairX[swap] <- pairY[swap]
            pairY[swap] <- t1
        }
        keepSet <- paste(pairX, pairY, sep = "_")
    }

    ## logCPM depth correction - based on bnbc/R/contactGroup_utils.R::logCPM.
    ## When cis is NULL we fall back to log2(count + 1) (no depth correction).
    if (!is.null(cis)) {
        libs <- (as.numeric(cis) + 1) / 1000000
        if (any(!is.finite(libs)) || any(libs <= 0)) {
            .stopNoCall("Invalid cis values: must be positive and finite")
        }
        cOffset <- abs(1 / min(libs) - 1 / max(libs))
        transformSample <- function(counts, s) {
            log((as.numeric(counts) + 0.5) / libs[s] + cOffset)
        }
    } else {
        transformSample <- function(counts, s) {
            log2(as.numeric(counts) + 1)
        }
    }

    ## 1. per-sample triplets, depth-corrected log, upper-triangle, keyed by band
    per_sample_dt <- vector("list", nSamples)
    for (s in seq_len(nSamples)) {
        d <- sparseList[[s]]
        if (is.null(d) || nrow(d) == 0) {
            per_sample_dt[[s]] <- data.table::data.table(
                i = integer(0), j = integer(0),
                x = numeric(0), d = integer(0)
            )
            data.table::setkey(per_sample_dt[[s]], d, i)
            next
        }
        ii <- as.integer(d$x)
        jj <- as.integer(d$y)
        swap <- ii > jj
        if (any(swap)) {
            t1 <- ii[swap]
            ii[swap] <- jj[swap]
            jj[swap] <- t1
        }
        dt <- data.table::data.table(
            i = ii, j = jj,
            x = transformSample(d$counts, s),
            d = jj - ii
        )
        data.table::setkey(dt, d, i)
        per_sample_dt[[s]] <- dt
    }

    ## Get the set of distinct bands across all samples
    allBands <- sort(unique(unlist(
        lapply(per_sample_dt, function(dt) unique(dt$d))
    )))
    if (length(allBands) == 0) {
        return(emptyOut())
    }

    ## Apply per-band quantile normalization
    outputChunks <- vector("list", length(allBands))
    for (k in seq_along(allBands)) {
        b <- allBands[k]

        ## Per-sample (i, x) for this band via keyed lookup
        sampData <- lapply(per_sample_dt, function(dt) {
            dt[.(b), .(i, x), on = "d", nomatch = 0]
        })

        ## Union of active i positions for this band
        activeI <- sort(unique(unlist(lapply(sampData, function(sd) sd$i))))
        L <- length(activeI)
        if (L == 0) next
        ## Note: L = 1 is OK - preprocessCore returns all samples equal to the mean.
        ## That degenerate cell will have zero variance downstream and be skipped.

        ## Build L x nSamples dense block. The "missing" fill is the per-sample
        ## logCPM transform of raw-count=0 (not literally 0) - because not having
        ## a triplet at (i, i+b) means the underlying raw count was 0, and
        ## log((0+0.5)/libs[s] + cOffset) is a small negative value, not 0.
        zeroFill <- vapply(
            seq_len(nSamples), function(s) transformSample(0, s),
            numeric(1)
        )
        qn <- matrix(zeroFill, nrow = L, ncol = nSamples, byrow = TRUE)
        for (s in seq_len(nSamples)) {
            sd <- sampData[[s]]
            if (nrow(sd) > 0) {
                rowIdx <- match(sd$i, activeI)
                qn[rowIdx, s] <- sd$x
            }
        }

        ## The row and column names are only retained if copy is FALSE. We
        ## don't need the raw matrix anymore anyway, and not copying it
        ## saves memory.
        preprocessCore::normalize.quantiles(qn, copy = FALSE)

        ## Filter QN output rows to keepPairs (if requested) BEFORE materializing
        ## the triplet emission. This is the fix that prevents the dcast CJ
        ## overflow on large chroms with many samples.
        if (!is.null(keepSet)) {
            bandKeys <- paste(activeI, activeI + b, sep = "_")
            keepRows <- which(bandKeys %in% keepSet)
            if (length(keepRows) == 0) next
            qn <- qn[keepRows, , drop = FALSE]
            active_i_em <- activeI[keepRows]
        } else {
            active_i_em <- activeI
        }

        ## Emit QN'd cells as triplets. With logCPM, "zero contact" maps to a
        ## per-sample baseline value (small but informative), so we emit all
        ## cells in the keep set. With log2(x+1) (no depth correction), genuine
        ## zeros stay zero and we skip them to save downstream memory.
        if (is.null(cis)) {
            nzMask <- qn != 0
            if (!any(nzMask)) next
            idxs <- which(nzMask, arr.ind = TRUE)
            vals <- qn[nzMask]
        } else {
            ## Emit every cell (all are informative under logCPM)
            idxs <- which(matrix(TRUE, nrow = nrow(qn), ncol = ncol(qn)), arr.ind = TRUE)
            vals <- as.numeric(qn)
        }
        outputChunks[[k]] <- data.table::data.table(
            s = idxs[, 2],
            i = active_i_em[idxs[, 1]],
            j = active_i_em[idxs[, 1]] + b,
            x = vals
        )
    }

    ## Free per-sample input triplets before assembling the wide output
    rm(per_sample_dt)
    invisible(gc(verbose = FALSE))

    ## Combine all output chunks into long form
    outLong <- data.table::rbindlist(outputChunks, fill = TRUE)
    rm(outputChunks)
    if (nrow(outLong) == 0) {
        return(emptyOut())
    }

    ## 5. pivot wide on (i, j) -> one IF column per sample
    wide <- data.table::dcast(outLong, i + j ~ s, value.var = "x", fill = 0)
    rm(outLong)
    invisible(gc(verbose = FALSE))
    data.table::setnames(wide, c("i", "j"), c("region1", "region2"))

    ## Ensure every sample has an IF column (some may have no data for this
    ## chromosome)
    for (s in seq_len(nSamples)) {
        nmIn <- as.character(s)
        nmOut <- paste0("IF_", s)
        if (nmIn %in% names(wide)) {
            data.table::setnames(wide, nmIn, nmOut)
        } else {
            wide[, (nmOut) := 0]
        }
    }

    wide[, chr := as.character(chr)]
    wide[, c(
        "chr", "region1", "region2",
        paste0("IF_", seq_len(nSamples))
    ), with = FALSE]
}
