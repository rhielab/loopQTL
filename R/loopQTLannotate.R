#' Annotate loopQTL results with target genes and regulatory peaks
#'
#' This function annotates SNP-loop QTL results with target genes, optional peak overlap,
#' and a SNP-loop Type classification.
#'
#' Adds three columns to `results(obj)`:
#' \describe{
#'   \item{`gene`}{The name of the genes found at the target anchor. Empty if no
#'     gene is found.}
#'   \item{`peak`}{`TRUE` if the SNP
#'     overlaps any peak from the supplied peak file, `FALSE` otherwise. When
#'     `removeTSS = TRUE`, peaks that overlap any TSS region are dropped
#'     before this overlap check (so promoter-like peaks don't count).}
#'   \item{`Type`}{One of:
#'     \itemize{
#'       \item `"WithinTSS"` - SNP overlaps any gene's TSS (± `TSSwindow`).
#'         Short-circuits all other annotation when `removeTSS = TRUE`.
#'       \item `"rs-gene"` - SNP not in TSS, has a target gene at the partner anchor.
#'       \item `"rs-<peakAnnotation>"` - SNP not in TSS, no gene at partner,
#'         but partner anchor overlaps a peak (only when `peak` is supplied).
#'         e.g. `"rs-Enhancer"`.
#'       \item `"rs-O"` - SNP not in TSS, no gene AND no peak at the partner anchor.
#'       \item `"rs(<peakAnnotation>)-gene"` - additionally overlaps a peak
#'         (only when `peak` is supplied), has a target gene. e.g.
#'         `"rs(Enhancer)-gene"`.
#'       \item `"rs(<peakAnnotation>)-<peakAnnotation>"` - SNP overlaps peak,
#'         partner overlaps peak, no gene at partner.
#'         e.g. `"rs(Enhancer)-Enhancer"`.
#'       \item `"rs(<peakAnnotation>)-O"` - overlaps a peak but no target gene
#'         and no peak at partner.
#'     }
#'     If both promtoer and peaks are found at the target anchor, it give priority on genes.
#'   }
#' }
#'
#' @param obj A loopQTL object with `results(obj)` filled by `runQTL()`.
#' @param gtf Path to a GTF file with `gene` entries.
#' @param peak Optional path to a peak  bed file. If supplied, an extra `peak` column
#'   is added and the `Type` labels gain the `"rs(<peakAnnotation>)-..."` annotations.
#' @param TSSwindow Base-pair window flank around each gene's TSS to define promoters. Default 2000.
#' @param anchorWindow Base-pair window flank added around each anchor when looking for
#'   overlapping promoters and peaks. Default 5000.
#' @param peakAnnotation Character label substituted into the `Type` values
#'   when `peak` is supplied. Default `"Peak"`. (e.g "Enhancer" or "H3K27ac")
#' @param removeTSS Logical; default `TRUE`. When `TRUE`, a TSS-based filtering on peaks will be
#' done to get the distal peak elements.
#'
#' @return The input `obj` with `results(obj)` augmented by `gene`, `peak`
#'   (when applicable), and `Type` columns.
#'
#' @export
#' @examples
#' ## Load the example object. In real analysis, the input
#' ## should be the output of loopQTLsetup().
#' ##Here we just use the example object.
#' obj <- loopQTLExampleData()
#' ## get the path of example annotation files
#' extdataDir <- system.file("extdata", package = "loopQTL")
#' peakPath   <- file.path(extdataDir, "example_annotation.bed")
#' gtfPath    <- file.path(extdataDir, "example.gtf.gz")
#'
#' ## Annotate SNP-loop pairs with gene, peak overlap, and Type
#' obj <- loopQTLannotate(
#'     obj,
#'     gtf            = gtfPath,
#'     peak           = peakPath,
#'     peakAnnotation = "Enhancer",
#'     TSSwindow      = 2000,
#'     removeTSS      = TRUE)
#' 
#' ##check the annotated results
#' head(results(obj))
loopQTLannotate <- function(obj,
                            gtf,
                            peak = NULL,
                            TSSwindow = 2000,
                            anchorWindow = 5000,
                            peakAnnotation = "Peak",
                            removeTSS = TRUE) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(gtf), length(gtf) == 1, file.exists(gtf),
        is.numeric(TSSwindow), length(TSSwindow) == 1, TSSwindow >= 0,
        is.numeric(anchorWindow), length(anchorWindow) == 1, anchorWindow >= 0,
        is.character(peakAnnotation), length(peakAnnotation) == 1,
        is.logical(removeTSS), length(removeTSS) == 1
    )
    if (!is.null(peak)) {
        stopifnot(is.character(peak), length(peak) == 1, file.exists(peak))
    }
    if (nrow(obj@results) == 0) {
        .stopNoCall("obj@results is empty. Run runQTL(obj) first.")
    }

    TSSwindow <- as.integer(TSSwindow)
    anchorWindow <- as.integer(anchorWindow)

    ## 1. Read GTF, extract gene entries, derive TSS
    message("reading GTF: ", gtf)
    gtfGRanges <- rtracklayer::import(gtf)
    typeCol <- if ("type" %in% names(S4Vectors::mcols(gtfGRanges))) {
        "type"
    } else if ("feature" %in% names(S4Vectors::mcols(gtfGRanges))) {
        "feature"
    } else {
        NULL
    }
    if (is.null(typeCol)) {
        .stopNoCall("GTF has neither `type` nor `feature` column; cannot find gene entries")
    }
    geneGRanges <- gtfGRanges[S4Vectors::mcols(gtfGRanges)[[typeCol]] == "gene"]
    if (length(geneGRanges) == 0) {
        .stopNoCall("No `gene` entries in GTF")
    }

    geneNames <- if ("gene_name" %in% names(S4Vectors::mcols(geneGRanges))) {
        as.character(geneGRanges$gene_name)
    } else if ("gene_id" %in% names(S4Vectors::mcols(geneGRanges))) {
        as.character(geneGRanges$gene_id)
    } else {
        .stopNoCall("GTF gene entries lack both `gene_name` and `gene_id`")
    }
    ## Some GTFs have NA gene_name for novel/predicted; fall back to gene_id.
    if ("gene_id" %in% names(S4Vectors::mcols(geneGRanges))) {
        naIdx <- which(is.na(geneNames) | geneNames == "")
        if (length(naIdx) > 0) {
            geneNames[naIdx] <- as.character(geneGRanges$gene_id[naIdx])
        }
    }

    strandChr <- as.character(GenomicRanges::strand(geneGRanges))
    tssPos <- ifelse(strandChr == "-",
        GenomicRanges::end(geneGRanges),
        GenomicRanges::start(geneGRanges)
    )
    tssRegions <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(geneGRanges),
        ranges = IRanges::IRanges(
            start = pmax(1, tssPos - TSSwindow),
            end = tssPos + TSSwindow
        )
    )
    S4Vectors::mcols(tssRegions)$gene <- geneNames
    message(length(tssRegions), " gene TSS regions (±", TSSwindow, " bp)")

    ## 2. Chrom-naming normalizer for THIS GTF  # GTF may use "chr1" (GENCODE) or "1" (Ensembl); loop anchors may use
    ## either depending on the BEDPE caller (Mustache uses "1"). Build a
    ## function that converts any chrom string to whatever the GTF uses.
    gtfChroms <- unique(as.character(GenomicRanges::seqnames(tssRegions)))
    gtf_uses_chr <- any(grepl("^chr", gtfChroms))
    normalizeChrom <- function(x) {
        x <- as.character(x)
        if (gtf_uses_chr) {
            ifelse(grepl("^chr", x), x, paste0("chr", x))
        } else {
            sub("^chr", "", x)
        }
    }

    ## 3. Optional: read peaks, filter by TSS if removeTSS
    peaksGRanges <- NULL
    if (!is.null(peak)) {
        message("reading peaks: ", peak)
        peaksGRanges <- rtracklayer::import(peak)
        ## Normalize peak chroms to GTF style so overlap actually works.
        peak_chroms_norm <- normalizeChrom(GenomicRanges::seqnames(peaksGRanges))
        peaksGRanges <- GenomicRanges::GRanges(
            seqnames = peak_chroms_norm,
            ranges   = IRanges::ranges(peaksGRanges)
        )
        if (removeTSS) {
            hits <- GenomicRanges::findOverlaps(peaksGRanges, tssRegions)
            peaks_in_tss <- unique(S4Vectors::queryHits(hits))
            nDropped <- length(peaks_in_tss)
            if (nDropped > 0) {
                peaksGRanges <- peaksGRanges[-peaks_in_tss]
                message(
                    "dropped ", nDropped,
                    " peaks overlapping TSS regions; ", length(peaksGRanges),
                    " peaks retained for SNP overlap."
                )
            }
        }
    }

    ## 4. Per-loop gene lists at anchor A and anchor B  # Critical: normalize anchor seqnames to the GTF's convention BEFORE
    ## findOverlaps. Without this, anchor seqnames like "1" silently fail to
    ## match TSS seqnames like "chr1", producing zero gene assignments. (This
    ## was the cause of the original "no rs-gene at all" bug.)
    loopsGInteractions <- obj@loops
    loop_ids_all <- S4Vectors::mcols(loopsGInteractions)$loopID
    anchorA <- InteractionSet::anchors(loopsGInteractions, type = "first")
    anchorB <- InteractionSet::anchors(loopsGInteractions, type = "second")

    padAnchor <- function(gr) {
        GenomicRanges::GRanges(
            seqnames = normalizeChrom(GenomicRanges::seqnames(gr)),
            ranges = IRanges::IRanges(
                start = pmax(1, GenomicRanges::start(gr) - anchorWindow),
                end   = GenomicRanges::end(gr) + anchorWindow
            )
        )
    }
    paddedA <- padAnchor(anchorA)
    paddedB <- padAnchor(anchorB)

    ## Anchor-side gene assignments (length = nLoops; one character vector per loop)
    build_gene_lists <- function(padded) {
        hits <- GenomicRanges::findOverlaps(padded, tssRegions)
        out <- vector("list", length(padded))
        if (length(hits) == 0) {
            return(out)
        }
        byLoop <- split(
            S4Vectors::subjectHits(hits),
            S4Vectors::queryHits(hits)
        )
        for (lp in names(byLoop)) {
            out[[as.integer(lp)]] <-
                unique(geneNames[byLoop[[lp]]])
        }
        out
    }
    genes_at_A <- build_gene_lists(paddedA)
    genes_at_B <- build_gene_lists(paddedB)

    ## Per-anchor peak overlap (logical vector, length = nLoops). NULL when no
    ## peak file was supplied -- then partner-side never gets a peak label.
    peak_at_A <- if (!is.null(peaksGRanges)) {
        GenomicRanges::countOverlaps(paddedA, peaksGRanges) > 0
    } else {
        NULL
    }
    peak_at_B <- if (!is.null(peaksGRanges)) {
        GenomicRanges::countOverlaps(paddedB, peaksGRanges) > 0
    } else {
        NULL
    }

    ## 5. Per-row annotation in obj@results
    res <- obj@results
    loopIdx <- match(res$loopID, loop_ids_all)

    ## SNP GRanges (1 bp per SNP), with chrom normalized to GTF convention.
    snpGR <- GenomicRanges::GRanges(
        seqnames = normalizeChrom(res$chrom),
        ranges   = IRanges::IRanges(start = res$pos, end = res$pos)
    )

    snp_in_tss <- if (removeTSS) {
        GenomicRanges::countOverlaps(snpGR, tssRegions) > 0
    } else {
        rep(FALSE, nrow(res))
    }
    snp_in_peak <- if (!is.null(peaksGRanges)) {
        GenomicRanges::countOverlaps(snpGR, peaksGRanges) > 0
    } else {
        rep(NA, nrow(res))
    }

    ## Resolve partner-anchor genes per row
    partnerGenes <- vector("list", nrow(res))
    partner_has_peak <- rep(FALSE, nrow(res))
    for (i in seq_len(nrow(res))) {
        li <- loopIdx[i]
        anc <- res$anchor[i]
        g <- character(0)
        if (anc == "A") {
            g <- genes_at_B[[li]]
        } else if (anc == "B") {
            g <- genes_at_A[[li]]
        } else if (anc == "both") g <- unique(c(genes_at_A[[li]], genes_at_B[[li]]))
        partnerGenes[[i]] <- if (is.null(g)) character(0) else g
        if (!is.null(peak_at_A)) {
            partner_has_peak[i] <-
                if (anc == "A") {
                    peak_at_B[li]
                } else if (anc == "B") {
                    peak_at_A[li]
                } else {
                    (peak_at_A[li] || peak_at_B[li])
                }
        }
    }

    geneCol <- vapply(
        partnerGenes,
        function(g) paste(g, collapse = ","),
        character(1)
    )
    hasGene <- nchar(geneCol) > 0

    ## Build Type column. Partner-side priority: gene > peak > O.
    rsPrefix <- if (!is.null(peaksGRanges)) {
        ifelse(snp_in_peak, paste0("rs(", peakAnnotation, ")"), "rs")
    } else {
        rep("rs", nrow(res))
    }
    partnerSuffix <- ifelse(
        hasGene, "gene",
        ifelse(partner_has_peak, peakAnnotation, "O")
    )
    typeCol <- paste0(rsPrefix, "-", partnerSuffix)
    typeCol[snp_in_tss] <- "WithinTSS"

    res$gene <- geneCol
    if (!is.null(peaksGRanges)) res$peak <- snp_in_peak
    res$Type <- typeCol

    obj@results <- res
    obj@metadata$annotate <- list(
        gtf = gtf,
        peak = peak,
        TSSwindow = TSSwindow,
        anchorWindow = anchorWindow,
        peakAnnotation = peakAnnotation,
        removeTSS = removeTSS,
        n_genes_used = length(unique(geneNames))
    )
    obj
}
