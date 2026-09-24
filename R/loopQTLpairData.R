#' Sample-level view of one loop × SNP pair
#'
#' This function returns the per-sample phenotype and genotype
#' (dosage 0/1/2) for a single SNP-loop pair. Useful for inspecting a
#' specific hit or exporting sample-level values for reporting.
#'
#' @param obj A loopQTL object with `results(obj)` filled.
#' @param loopID Character vector of loopIDs to extract.
#' @param snpID Character vector of snpIds to extract.
#' @param normal Apply RINT per loop before returning phenotype values. default `TRUE`.
#'
#' @return data.frame with columns `sample`, `loopID`, `snpID`, `genotype`
#'   (0 - Reference homozygou; 1 - heterozygous, 2 - Alternative homozygous), `phenotype`, `cis`.
#'
#' @export
#' @examples
#' ## Load the example object. In real analysis, the input
#' ## should be the output of runQTL().
#' ##Here we just use the example object. 
#' obj <- loopQTLExampleData()
#' ##Check results so that you can choose your interested loopID and snpID
#' head(results(obj))
#' ##Then you can specify the specifi SNP-loop pair data you want to check at samples level.
#' pair <- loopQTLpairData(
#'     obj,
#'     loopID = "L0000187_chr1_46030000_46220000",
#'     snpID  = "rs11584814")
#' head(pair)
loopQTLpairData <- function(obj, loopID, snpID, normal = TRUE) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(loopID), length(loopID) > 0,
        is.character(snpID), length(snpID) > 0,
        is.logical(normal), length(normal) == 1
    )
    if (nrow(obj@results) == 0) {
        .stopNoCall("obj@results is empty. Run runQTL(obj) first to populate it.")
    }

    loop_ids_all <- S4Vectors::mcols(obj@loops)$loopID
    missingLoops <- setdiff(loopID, loop_ids_all)
    if (length(missingLoops) > 0) {
        .stopNoCall("loopID(s) not in obj@loops: ", paste(missingLoops, collapse = ", "))
    }

    missingSNPs <- setdiff(snpID, obj@results$snpID)
    if (length(missingSNPs) > 0) {
        .stopNoCall("snpID(s) not in obj@results: ", paste(missingSNPs, collapse = ", "))
    }

    ## SNP info from obj@results (positions + ref/alt)
    snpInfo <- unique(obj@results[
        obj@results$snpID %in% snpID,
        c("snpID", "chrom", "pos", "ref", "alt"),
        drop = FALSE
    ])
    snpInfo <- snpInfo[match(snpID, snpInfo$snpID), , drop = FALSE]

    ## Genotypes via targeted tabix read  # Normalize chrom naming to whatever the VCF uses.
    tab <- Rsamtools::TabixFile(obj@vcfPath)
    vcfChroms <- Rsamtools::headerTabix(tab)$seqnames
    vcfChromNormalizerFunc <- .make_chrom_normalizer(vcfChroms)
    snp_chrom_vcf <- vcfChromNormalizerFunc(snpInfo$chrom)
    whichGR <- GenomicRanges::GRanges(
        seqnames = snp_chrom_vcf,
        ranges   = IRanges::IRanges(start = snpInfo$pos, end = snpInfo$pos)
    )
    param <- VariantAnnotation::ScanVcfParam(
        which = whichGR, geno = "GT", info = NA_character_
    )
    vcf <- VariantAnnotation::readVcf(obj@vcfPath, genome = "unknown", param = param)

    rr <- SummarizedExperiment::rowRanges(vcf)
    vcfChrom <- as.character(GenomicRanges::seqnames(rr))
    vcfPos <- GenomicRanges::start(rr)
    vcfRef <- as.character(VariantAnnotation::ref(vcf))

    ## Match by (chrom, pos, ref) to handle multi-allelic loci
    vcfKey <- paste(vcfChrom, vcfPos, vcfRef, sep = ":")
    snpKey <- paste(snp_chrom_vcf, snpInfo$pos, snpInfo$ref, sep = ":")
    snpsInVCF <- match(snpKey, vcfKey)
    if (anyNA(snpsInVCF)) {
        .stopNoCall(
            "Could not find SNP(s) in VCF at expected (chrom, pos, ref): ",
            paste(snpID[is.na(snpsInVCF)], collapse = ", ")
        )
    }

    gt <- VariantAnnotation::geno(vcf)$GT[snpsInVCF, , drop = FALSE]

    ## Reorder VCF samples to canonical
    sm <- match(obj@samples, colnames(gt))
    if (anyNA(sm)) {
        .stopNoCall(
            "VCF missing samples: ",
            paste(obj@samples[is.na(sm)], collapse = ", ")
        )
    }
    gt <- gt[, sm, drop = FALSE]
    dosage <- .gt_to_dosage(gt)
    rownames(dosage) <- snpID
    colnames(dosage) <- obj@samples

    ## Phenotype matrix for the requested loops
    loopIdx <- match(loopID, loop_ids_all)
    phenoFull <- if (normal) .rintMatrix(obj@phenotype) else obj@phenotype
    phenoSub <- phenoFull[loopIdx, , drop = FALSE]
    rownames(phenoSub) <- loopID

    ## Build long table: nSamples x nLoops x n_snps rows
    nSamples <- length(obj@samples)
    grid <- expand.grid(
        sampleIdx = seq_len(nSamples),
        loopIdx = seq_along(loopID),
        snpIdx = seq_along(snpID),
        KEEP.OUT.ATTRS = FALSE
    )
    cisVec <- if (length(obj@cis) > 0) obj@cis[obj@samples] else rep(NA_real_, nSamples)

    out <- data.frame(
        sample = obj@samples[grid$sampleIdx],
        loopID = loopID[grid$loopIdx],
        snpID = snpID[grid$snpIdx],
        genotype = dosage[cbind(grid$snpIdx, grid$sampleIdx)],
        phenotype = phenoSub[cbind(grid$loopIdx, grid$sampleIdx)],
        cis = cisVec[grid$sampleIdx],
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}
