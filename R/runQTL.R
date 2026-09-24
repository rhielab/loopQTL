#' Anchor-restricted association testing
#'
#' This function first reads all anchor-window SNPs from the VCF, and runs
#' MatrixEQTL twice - once treating anchor A as the gene position, once anchor 
#' B. The two result tables are outer-joined on SNP-loop level. Covariates are
#' supplied via the `covariates` argument - typically the phenotype PCs chosen
#' in the previous step. By default, rank-based inverse normal transformation
#' (RINT) is applied per loop before regression.
#'
#' @param obj A loopQTL object on which `loopQTLnormalize()` has been run.
#' @param window Base-pair window flank added to each anchor bin when matching
#' SNPs. Default 0.
#' @param covariates Optional `samples x k` matrix/data.frame. If NULL, will be 
#' auto-chosen via loopQTLchooseK().
#' @param normal Apply Rank-based Inverse Normal Transformation (RINT) per loop 
#' before regression. Default `TRUE`.
#' @param Kmethod `"elbow"` or `"BE"`; passed to `loopQTLchooseK()` when 
#' covariates is `NULL`. Default `"elbow"`.
#' @param Kmax The maximum number of PCs used. Defaults to the number of samples
#' minus two.
#' @param pvOutputThreshold MatrixEQTL p-value threshold. Default 1 (return 
#' every pair).
#'
#' @return The loopQTL object with `results(obj)` filled.
#'
#' @export
#' @examples
#' ## Load the example object. In real analysis, the input
#' ## should be the output of loopQTLnormalize().
#' ##Here we just use the example object.
#' obj <- loopQTLExampleData()
#' ##Run the function, please note that the covariates option is not specified.
#' obj <- runQTL(
#'    obj,
#'    window     = 10000)
#'
#' ## In real analysis we recommend running runQTL() with covariates
#' ## specified for more accurate regression, typically the top-K PCs.
#' ## The helper loopQTLchooseK() picks K following the PCAForQTL method
#' ## (https://github.com/heatherjzhou/PCAForQTL).
#' Kres <- loopQTLchooseK(obj)
#' obj <- runQTL(
#'     obj,
#'     window     = 10000,
#'     covariates = Kres$PCs)
runQTL <- function(obj,
                   window = 0,
                   covariates = NULL,
                   normal = TRUE,
                   Kmethod = c("elbow", "BE"),
                   Kmax = NULL,
                   pvOutputThreshold = 1) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        length(obj@phenotype) > 0,
        is.numeric(window), length(window) == 1, window >= 0,
        is.logical(normal), length(normal) == 1, !is.na(normal)
    )
    window <- as.integer(window)
    Kmethod <- match.arg(Kmethod)

    loopsGInteractions <- obj@loops
    samplesVec <- obj@samples
    loopIDs <- S4Vectors::mcols(loopsGInteractions)$loopID
    anchorA <- InteractionSet::anchors(loopsGInteractions, type = "first")
    anchorB <- InteractionSet::anchors(loopsGInteractions, type = "second")

    ## Normalize the phenotype if requested (but make a copy of it first so we
    ## don't modify the object)
    phenotype <- obj@phenotype
    if (normal) {
        phenotype <- .rintMatrix(phenotype)
    }

    n_zero_var <- sum(rowSums(!is.na(phenotype)) == 0) -
        sum(rowSums(!is.na(obj@phenotype)) == 0)
    ## n_zero_var = loops newly NA'd by RINT (had data but constant after RINT).
    ## Originally-NA loops aren't counted here.

    ## Build union of (anchor +/- window) target ranges for the VCF read
    pad <- function(gr) {
        GenomicRanges::GRanges(
            seqnames = GenomicRanges::seqnames(gr),
            ranges = IRanges::IRanges(
                start = pmax(1, GenomicRanges::start(gr) - window),
                end   = GenomicRanges::end(gr) + window
            )
        )
    }
    readTargets <- GenomicRanges::reduce(c(pad(anchorA), pad(anchorB)))

    ## tiny helper to short-circuit: write empty result table and return obj
    .returnEmpty <- function() {
        obj@results <- .empty_qtl_result()
        obj
    }

    ## Detect VCF chrom naming style and rewrite the read targets to match
    tab <- Rsamtools::TabixFile(obj@vcfPath)
    vcfChroms <- Rsamtools::headerTabix(tab)$seqnames
    if (length(vcfChroms) == 0) {
        .stopNoCall("Could not read chromosome names from VCF tabix index")
    }
    vcfChromNormalizerFunc <- .make_chrom_normalizer(vcfChroms)
    rt_chroms_vcf <- vcfChromNormalizerFunc(
            as.character(GenomicRanges::seqnames(readTargets))
    )
    keep <- rt_chroms_vcf %in% vcfChroms
    if (!any(keep)) {
        return(.returnEmpty())
    }
    readTargets <- GenomicRanges::GRanges(
        seqnames = rt_chroms_vcf[keep],
        ranges = IRanges::IRanges(
            start = GenomicRanges::start(readTargets)[keep],
            end   = GenomicRanges::end(readTargets)[keep]
        )
    )

    ## Single targeted VCF read
    param <- VariantAnnotation::ScanVcfParam(
        which = readTargets, geno = "GT", info = NA_character_
    )
    vcf <- VariantAnnotation::readVcf(
            obj@vcfPath, genome = "unknown", param = param
    )
    if (length(vcf) == 0) {
        return(.returnEmpty())
    }

    vcf_sample_names <- colnames(VariantAnnotation::geno(vcf)$GT)
    m <- match(samplesVec, vcf_sample_names)
    if (anyNA(m)) {
        .stopNoCall(
            "VCF missing samples that were validated at setup: ",
            paste(samplesVec[is.na(m)], collapse = ", ")
        )
    }

    gt <- VariantAnnotation::geno(vcf)$GT[, m, drop = FALSE]
    dosage <- .gt_to_dosage(gt)
    colnames(dosage) <- samplesVec

    rr <- SummarizedExperiment::rowRanges(vcf)
    snpChromVCF <- as.character(GenomicRanges::seqnames(rr)) 
    snpPos <- GenomicRanges::start(rr)
    snpRef <- as.character(VariantAnnotation::ref(vcf))
    altList <- VariantAnnotation::alt(vcf)
    snpAlt <- vapply(
        altList,
        function(a) paste(as.character(a), collapse = ","),
        character(1)
    )
    snp_id_raw <- rownames(vcf)
    if (is.null(snp_id_raw) || any(snp_id_raw == "" | is.na(snp_id_raw))) {
        snp_id_raw <- sprintf(
            "%s:%d_%s/%s",
            snpChromVCF, snpPos, snpRef, snpAlt
        )
    }

    ## MAF + drop SNPs with no informative variation
    af <- rowMeans(dosage, na.rm = TRUE) / 2
    mafAll <- pmin(af, 1 - af)
    snpKeep <- is.finite(mafAll) & mafAll > 0
    if (!any(snpKeep)) {
        return(.returnEmpty())
    }
    dosage <- dosage[snpKeep, , drop = FALSE]
    snpID <- snp_id_raw[snpKeep]
    snpChromVCF <- snpChromVCF[snpKeep]
    snpPos <- snpPos[snpKeep]
    snpRef <- snpRef[snpKeep]
    snpAlt <- snpAlt[snpKeep]
    mafKept <- mafAll[snpKeep]
    rownames(dosage) <- snpID

    ## Resolve covariates
    if (is.null(covariates)) {
        ## Use the already-RINT'd-or-not `phenotype` matrix to avoid double-RINT.
        Kres <- .pca_for_qtl(phenotype,
            samples = samplesVec,
            Kmethod = Kmethod, Kmax = Kmax
        )
        cvrtMat <- t(Kres$PCs) # MatrixEQTL wants k x nSamples
        obj_K_meta <- list(K = Kres$K, Kmethod = Kres$method, Ksource = "autoPCAForQTL")
    } else {
        covDF <- as.data.frame(covariates, stringsAsFactors = FALSE)
        if (!is.null(rownames(covDF))) {
            mm <- match(samplesVec, rownames(covDF))
            if (anyNA(mm)) {
                .stopNoCall(
                    "`covariates` missing rows for samples: ",
                    paste(samplesVec[is.na(mm)], collapse = ", ")
                )
            }
            covDF <- covDF[mm, , drop = FALSE]
        } else if (nrow(covDF) != length(samplesVec)) {
            .stopNoCall(
                "`covariates` has ", nrow(covDF), " rows but ",
                length(samplesVec), " samples; provide rownames or matched order"
            )
        }
        cvrtMat <- t(as.matrix(covDF))
        colnames(cvrtMat) <- samplesVec
        obj_K_meta <- list(K = nrow(cvrtMat), Kmethod = NA_character_, Ksource = "user")
    }
    obj@metadata$runQTL_K <- obj_K_meta
    obj@metadata$runQTLNormal <- normal
    obj@metadata$runQTL_n_zero_var <- max(0, n_zero_var)

    ## Build SlicedData inputs for MatrixEQTL
    snpsSd <- .asSliced(dosage)
    geneSD <- .asSliced(phenotype, rownamesDefault = loopIDs)
    cvrtSd <- .asSliced(cvrtMat)

    ## SNP positions (use original loop chrom style for output, but MatrixEQTL
    ## uses position tables only for cis-distance arithmetic, so we feed it the
    ## vcf-style chrom names and matching anchor chrom names)
    loop_chrom_vcf <- vcfChromNormalizerFunc(as.character(GenomicRanges::seqnames(anchorA)))
    snpspos <- data.frame(
        snpid = snpID,
        chr = snpChromVCF,
        pos = snpPos,
        stringsAsFactors = FALSE
    )
    geneposA <- data.frame(
        geneid = loopIDs,
        chr = loop_chrom_vcf,
        left = GenomicRanges::start(anchorA),
        right = GenomicRanges::end(anchorA),
        stringsAsFactors = FALSE
    )
    geneposB <- data.frame(
        geneid = loopIDs,
        chr = loop_chrom_vcf,
        left = GenomicRanges::start(anchorB),
        right = GenomicRanges::end(anchorB),
        stringsAsFactors = FALSE
    )

    ## Run MatrixEQTL twice
    outA <- .run_matrix_eqtl(snpsSd, geneSD, cvrtSd,
        snpspos, geneposA,
        cisDist = window,
        pvOutputThreshold = pvOutputThreshold
    )
    outB <- .run_matrix_eqtl(snpsSd, geneSD, cvrtSd,
        snpspos, geneposB,
        cisDist = window,
        pvOutputThreshold = pvOutputThreshold
    )

    ## Merge on (loopID = gene, snpID = snps); label anchor A / B / both.
    merged <- .merge_anchor_runs(outA, outB)
    if (nrow(merged) == 0) {
        obj@results <- .empty_qtl_result()
        return(obj)
    }

    ## Attach SNP metadata + MAF
    snpMeta <- data.frame(
        snpID = snpID, chrom = snpChromVCF, pos = snpPos,
        ref = snpRef, alt = snpAlt, maf = mafKept,
        stringsAsFactors = FALSE
    )
    out <- merge(merged, snpMeta, by = "snpID", all.x = TRUE, sort = FALSE)
    out$nUsed <- rowSums(!is.na(dosage[match(out$snpID, rownames(dosage)), , drop = FALSE]))

    ## Final column order
    out <- out[, c(
        "loopID", "snpID", "chrom", "pos", "ref", "alt", "anchor",
        "beta", "se", "t", "p", "nUsed", "maf"
    )]
    rownames(out) <- NULL

    ## Store the results in the object and return it
    obj@results <- out
    obj
}

## Internal helpers

.empty_qtl_result <- function() {
    data.frame(
        loopID = character(0), snpID = character(0),
        chrom = character(0), pos = integer(0),
        ref = character(0), alt = character(0),
        anchor = character(0),
        beta = numeric(0), se = numeric(0), t = numeric(0), p = numeric(0),
        nUsed = integer(0), maf = numeric(0),
        stringsAsFactors = FALSE
    )
}

## Convert a (snps x samples) GT character matrix to numeric dosage (0/1/2).
## Anything beyond biallelic 0/1 becomes NA.
.gt_to_dosage <- function(gt) {
    raw <- gsub("|", "/", gt, fixed = TRUE)
    raw[raw %in% c("", ".", "./.")] <- NA_character_
    out <- matrix(NA_real_, nrow = nrow(gt), ncol = ncol(gt))
    ok <- !is.na(raw)
    if (!any(ok)) {
        return(out)
    }
    sp <- strsplit(raw[ok], "/", fixed = TRUE)
    bad <- vapply(sp, function(x) length(x) != 2, logical(1))
    a1 <- suppressWarnings(as.integer(vapply(sp, `[`, character(1), 1)))
    a2 <- suppressWarnings(as.integer(vapply(sp, `[`, character(1), 2)))
    bad <- bad | is.na(a1) | is.na(a2) | a1 > 1 | a2 > 1
    vals <- ifelse(bad, NA_real_, a1 + a2)
    out[ok] <- vals
    out
}

## Wrap a matrix in a MatrixEQTL::SlicedData.
.asSliced <- function(mat, rownamesDefault = NULL) {
    if (!is.null(rownamesDefault) && is.null(rownames(mat))) {
        rownames(mat) <- rownamesDefault
    }
    sd <- MatrixEQTL::SlicedData$new()
    sd$CreateFromMatrix(as.matrix(mat))
    sd
}

## Run MatrixEQTL and return a simplified data.frame of cis hits
.run_matrix_eqtl <- function(snpsSd, geneSD, cvrtSd,
                             snpspos, genepos,
                             cisDist, pvOutputThreshold) {
    tmpCis <- tempfile()
    on.exit(unlink(tmpCis), add = TRUE)
    res <- MatrixEQTL::Matrix_eQTL_main(
        snps = snpsSd,
        gene = geneSD,
        cvrt = cvrtSd,
        output_file_name = NULL,
        pvOutputThreshold = 0, # disable trans output
        useModel = MatrixEQTL::modelLINEAR,
        errorCovariance = numeric(0),
        verbose = FALSE,
        output_file_name.cis = tmpCis,
        pvOutputThreshold.cis = pvOutputThreshold,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    )
    cis <- res$cis$eqtls
    if (is.null(cis) || nrow(cis) == 0) {
        return(data.frame(
            loopID = character(0), snpID = character(0),
            beta = numeric(0), se = numeric(0), t = numeric(0),
            p = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    ## MatrixEQTL columns: snps, gene, statistic, pvalue, beta (FDR ignored)
    ## Some versions also include se; if not, derive it: se = beta / statistic.
    seCol <- if ("se" %in% names(cis)) cis$se else cis$beta / cis$statistic
    data.frame(
        loopID = as.character(cis$gene),
        snpID = as.character(cis$snps),
        beta = as.numeric(cis$beta),
        se = as.numeric(seCol),
        t = as.numeric(cis$statistic),
        p = as.numeric(cis$pvalue),
        stringsAsFactors = FALSE
    )
}

## Outer-join two anchor runs on (loopID, snpID) and label `anchor`.
## When a pair appears in both runs the model is identical, so we just keep one
## row and label it `both`.
.merge_anchor_runs <- function(outA, outB) {
    if (nrow(outA) == 0 && nrow(outB) == 0) {
        return(data.frame(
            loopID = character(0), snpID = character(0),
            anchor = character(0),
            beta = numeric(0), se = numeric(0), t = numeric(0),
            p = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    key <- function(d) paste(d$loopID, d$snpID, sep = "::")
    kA <- if (nrow(outA)) key(outA) else character(0)
    kB <- if (nrow(outB)) key(outB) else character(0)

    inBoth <- intersect(kA, kB)
    onlyA <- setdiff(kA, kB)
    onlyB <- setdiff(kB, kA)

    rows <- list()
    if (length(inBoth) > 0) {
        pick <- match(inBoth, kA)
        r <- outA[pick, , drop = FALSE]
        r$anchor <- "both"
        rows[[length(rows) + 1]] <- r
    }
    if (length(onlyA) > 0) {
        pick <- match(onlyA, kA)
        r <- outA[pick, , drop = FALSE]
        r$anchor <- "A"
        rows[[length(rows) + 1]] <- r
    }
    if (length(onlyB) > 0) {
        pick <- match(onlyB, kB)
        r <- outB[pick, , drop = FALSE]
        r$anchor <- "B"
        rows[[length(rows) + 1]] <- r
    }
    do.call(rbind, rows)
}
