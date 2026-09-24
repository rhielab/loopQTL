#' Locus-level track view of loops, SNPs, and genes
#'
#' This function renders a locus-scale ggplot combining the consensus
#' loops, the tested SNP, and (when a GTF is supplied) gene annotations, over
#' a user-specified region. Please note that to use `values = "cpm"` or
#' `values = "raw"` to show per-loop contact strength, you need to set
#' `saveWindows = TRUE` to save the per-SNP contact-window matrices from
#' [loopQTLnormalize()].
#'
#' @param obj A loopQTL object with `results(obj)` filled by `runQTL()`.
#' @param snpID Single snpID present in `results(obj)`.
#' @param chrom Chromosome. e.g. "chr1".
#' @param start The start genomic coordinate in bp for the x-axis.
#' @param end The end genomic coordinate in bp for the x-axis.
#' @param gtf Path to a GTF file.
#' @param pThreshold Nominal p-value cutoff for "significant" loops. Default `0.05`.
#' @param values Which contact-strength matrix to color arcs by:
#'   \describe{
#'     \item{`"normalized"`}{Default. The matrix from `phenotype(obj)`.
#'       Matches what `runQTL()` tested.}
#'     \item{`"cpm"`}{CPM matrix. Requires `saveWindows` output.}
#'     \item{`"raw"`}{raw count matrix. Requires `saveWindows`.}
#'   }
#' @param windows Optional path/loaded list from `loopQTLnormalize(
#'   saveWindows = TRUE)`. Required when `values` is "cpm" or "raw".
#' @param strongColor Color for the strongest signal on the arc gradient.
#'   Default `"#FF7917"` (orange).
#' @param weakColor Color for the weakest signal. Default `"grey90"`.
#' @param limit Optional numeric cap on the color scale (values above are
#'   limited to `strongColor`).
#' @param geneBiotypes Character vector of `gene_biotype` values to keep
#'   from the GTF. Default `"protein_coding"`; pass `NULL` to keep all.
#' @param axisUnit `"Mb"` (default) or `"bp"` for x-axis tick labels.
#'
#' @return A `ggplot` (gene / SNP / arc panels
#'   stacked) object. Save with `ggplot2::ggsave()`.
#'
#' @export
#' @examples
#' ## Load the example object. In real analysis, the input
#' ## should be the output of runQTL().
#' ## If you specify `values = "cpm"` or `values = "raw"`, the object must
#' ## carry the contact matrix `.rds` (see `obj@metadata$windowsOut`);
#' ## otherwise the normalized matrix stored in the object is used.
#' ##Here we just use the example object.
#' obj <- loopQTLExampleData()
#' ##You also need a gtf file input to give the genomic annotations. Here we provide an example gtf.
#' gtfPath <- system.file("extdata", "example.gtf.gz", package = "loopQTL")
#' ## Return is a ggplot object so you can ggsave() to output.
#' p <- loopQTLplotLocus(
#'     obj,
#'     chrom = "chr1",
#'     start = 45900000,
#'     end = 46300000,
#'     snpID = "rs11584814",
#'     values = "cpm",
#'     gtf = gtfPath)
loopQTLplotLocus <- function(obj, snpID,
                             chrom, start, end,
                             gtf,
                             pThreshold = 0.05,
                             values = c("normalized", "cpm", "raw"),
                             windows = NULL,
                             strongColor = "#FF7917",
                             weakColor = "grey90",
                             limit = NULL,
                             geneBiotypes = "protein_coding",
                             axisUnit = c("Mb", "bp")) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(snpID), length(snpID) == 1, nzchar(snpID),
        is.character(chrom), length(chrom) == 1,
        is.numeric(start), length(start) == 1,
        is.numeric(end), length(end) == 1, start < end,
        is.character(gtf), length(gtf) == 1, file.exists(gtf),
        is.numeric(pThreshold), length(pThreshold) == 1,
        pThreshold > 0, pThreshold <= 1
    )
    values <- match.arg(values)
    axisUnit <- match.arg(axisUnit)
    if (nrow(obj@results) == 0) {
        .stopNoCall("obj@results is empty. Run runQTL(obj) first.")
    }

    start <- as.integer(start)
    end <- as.integer(end)

    ## 1. Filter obj@results to (this SNP, p <= threshold) and to loops in the region
    resultsForThisSNP <- obj@results[obj@results$snpID == snpID, , drop = FALSE]
    if (nrow(resultsForThisSNP) == 0) {
        .stopNoCall("SNP ", snpID, " not found in obj@results")
    }
    if (!identical(unique(as.character(resultsForThisSNP$chrom)), chrom)) {
        .stopNoCall(
            "Supplied chrom = '", chrom, "' does not match the SNP's chrom ('",
            resultsForThisSNP$chrom[1], "') in obj@results"
        )
    }
    snpPos <- resultsForThisSNP$pos[1]
    snpRef <- resultsForThisSNP$ref[1]
    snpAlt <- resultsForThisSNP$alt[1]

    significantResults <- resultsForThisSNP[resultsForThisSNP$p <= pThreshold & is.finite(resultsForThisSNP$p), ,
        drop = FALSE
    ]
    if (nrow(significantResults) == 0) {
        .stopNoCall("No significant (p <= ", pThreshold, ") loops for SNP ", snpID)
    }

    ## Look up anchor coords for each significant loop; filter to region
    loop_ids_all <- S4Vectors::mcols(obj@loops)$loopID
    ancAAll <- InteractionSet::anchors(obj@loops, "first")
    ancBAll <- InteractionSet::anchors(obj@loops, "second")
    li <- match(significantResults$loopID, loop_ids_all)
    loopsDF <- data.frame(
        loopID = significantResults$loopID,
        chrom = as.character(GenomicRanges::seqnames(ancAAll[li])),
        ancAStart = GenomicRanges::start(ancAAll[li]),
        ancAEnd = GenomicRanges::end(ancAAll[li]),
        ancBStart = GenomicRanges::start(ancBAll[li]),
        ancBEnd = GenomicRanges::end(ancBAll[li]),
        beta = significantResults$beta,
        p = significantResults$p,
        rowI = li,
        stringsAsFactors = FALSE
    )
    inRegion <- loopsDF$chrom == chrom &
        pmin(loopsDF$ancAStart, loopsDF$ancBStart) >= start &
        pmax(loopsDF$ancAEnd, loopsDF$ancBEnd) <= end
    loopsDF <- loopsDF[inRegion, , drop = FALSE]
    if (nrow(loopsDF) == 0) {
        .stopNoCall(
            "No significant loops fully inside [", start, ", ", end,
            "] for SNP ", snpID
        )
    }
    loopsDF$ancAMid <- (loopsDF$ancAStart + loopsDF$ancAEnd) / 2
    loopsDF$ancBMid <- (loopsDF$ancBStart + loopsDF$ancBEnd) / 2

    ## 2. Per-sample dosage for this SNP
    dose <- .read_snp_dosage(obj@vcfPath, chrom, snpPos,
        samplesVec = obj@samples
    )
    ## Split samples into genotype classes 0/1/2 (drop NA)
    genotypeGroups <- split(
        names(dose)[!is.na(dose)],
        factor(dose[!is.na(dose)], levels = c(0, 1, 2))
    )

    ## 3. For each loop, get per-sample contact values, then mean per genotype
    valuesMat <- .locus_values_matrix(obj, loopsDF, values, windows)

    mean_per_geno <- vapply(names(genotypeGroups), function(g) {
        ss <- genotypeGroups[[g]]
        if (length(ss) == 0) {
            return(rep(NA_real_, nrow(loopsDF)))
        }
        rowMeans(valuesMat[, ss, drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(loopsDF)))

    ## as.matrix() is needed in case there is only one loop - if so, the result
    ## would not be a matrix, and colnames() would fail.
    if (length(dim(mean_per_geno)) < 2) {
        mean_per_geno <- t(as.matrix(mean_per_geno))
    }

    colnames(mean_per_geno) <- names(genotypeGroups)

    nPer <- vapply(genotypeGroups, length, integer(1))
    genotypeLabels <- c(
        paste0(snpRef, "/", snpRef),
        paste0(snpRef, "/", snpAlt),
        paste0(snpAlt, "/", snpAlt)
    )
    panelTitles <- sprintf("%s (n=%d)", genotypeLabels, nPer)

    ## 4. Color scale: use ACTUAL data range so small differences show up.
    ## BNBC values can be small (or negative); a fixed 0..max scale often
    ## squashes real differences into the pale end of the gradient.
    vfin <- mean_per_geno[is.finite(mean_per_geno)]
    if (!is.null(limit)) {
        effLo <- 0
        effHi <- limit
    } else if (length(vfin) > 0) {
        effLo <- min(vfin)
        effHi <- max(vfin)
        if (effHi - effLo < .Machine$double.eps) {
            effLo <- effLo - 0.5
            effHi <- effHi + 0.5
        }
    } else {
        effLo <- 0
        effHi <- 1
    }

    geneDF <- .load_genes_in_region(gtf, chrom, start, end, geneBiotypes)
    exonDF <- .load_exons_in_region(gtf, chrom, start, end, geneDF$gene)

    axisFormat <- if (axisUnit == "Mb") {
        function(x) sprintf("%.2f", x / 1000000)
    } else {
        function(x) format(x, big.mark = ",", scientific = FALSE)
    }

    ## 7. Build the four ggplot panels
    pGene <- .build_gene_track(geneDF, exonDF, start, end, chrom, axisFormat)
    pSnp <- .build_snp_track(snpPos, snpID, start, end, axisFormat)

    arcPanels <- lapply(seq_along(genotypeGroups), function(gI) {
        genotypeLabel <- names(genotypeGroups)[gI]
        arcDF <- data.frame(
            loopID = loopsDF$loopID,
            x1 = loopsDF$ancAMid,
            x2 = loopsDF$ancBMid,
            value = mean_per_geno[, genotypeLabel]
        )
        .build_arc_panel(arcDF, start, end, axisFormat,
            effLo, effHi, strongColor, weakColor,
            panelTitles[gI],
            showX = (gI == length(genotypeGroups))
        )
    })

    ## Stack the plots with patchwork.
    ## guides = "collect" consolidates all three arc legends into one.
    heights <- c(1.2, 0.4, rep(1.2, length(arcPanels)))
    wrap <- patchwork::wrap_plots(
        c(list(pGene, pSnp), arcPanels),
        ncol = 1,
        heights = heights
    ) + patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(
            title = sprintf(
                "%s | %s:%d-%d",
                snpID, chrom, start, end
            ),
            subtitle = sprintf(
                "REF=%s  ALT=%s  |  %d significant loop(s) (p <= %.3g)  |  %s values",
                snpRef, snpAlt, nrow(loopsDF),
                pThreshold, values
            )
        )
    wrap
}


## Internal helpers

## Pull per-sample contact values for the given loops based on `values`.
## Returns a (nLoops x nSamples) matrix in obj@samples order.
.locus_values_matrix <- function(obj, loopsDF, values, windows) {
    loop_ids_all <- S4Vectors::mcols(obj@loops)$loopID
    li <- match(loopsDF$loopID, loop_ids_all)
    if (values == "normalized") {
        if (length(obj@phenotype) == 0) {
            .stopNoCall(
                "loopQTLnormalize() must be run before using the ",
                "values = \"normalized\" option."
            )
        }
        return(obj@phenotype[li, obj@samples, drop = FALSE])
    }
    ## cpm / raw modes need windows data
    w <- .resolveWindows(obj, windows)
    if (is.null(w$cellsRaw)) {
        .stopNoCall(
            "loopQTLnormalize() must be run with saveWindows set to TRUE to ",
            "use values = \"", values, "\"."
        )
    }
    res <- obj@resolution
    nLoops <- nrow(loopsDF)
    out <- matrix(0,
        nrow = nLoops, ncol = length(obj@samples),
        dimnames = list(loopsDF$loopID, obj@samples)
    )
    for (i in seq_len(nLoops)) {
        binA <- .binAnchor(loopsDF$ancAStart[i], loopsDF$ancAEnd[i], res)
        binB <- .binAnchor(loopsDF$ancBStart[i], loopsDF$ancBEnd[i], res)
        lo <- min(binA, binB)
        hi <- max(binA, binB)
        hit <- which(w$cellCoordinates$chrom == loopsDF$chrom[i] &
            w$cellCoordinates$binI == lo &
            w$cellCoordinates$binJ == hi)
        if (length(hit) == 0) next
        rawVals <- w$cellsRaw[obj@samples, hit[1]]
        if (values == "cpm") {
            if (length(obj@cis) == 0) {
                .stopNoCall("values = 'cpm' needs obj@cis (per-sample cis depth).")
            }
            rawVals <- rawVals / (obj@cis[obj@samples] / 1000000)
        }
        out[i, ] <- rawVals
    }
    out
}

## Read genes overlapping [start, end] on chrom from a GTF. Handles chr vs no-chr.
.load_genes_in_region <- function(gtf, chrom, start, end, biotypes) {
    gr <- rtracklayer::import(gtf)
    typeCol <- if ("type" %in% names(S4Vectors::mcols(gr))) {
        "type"
    } else if ("feature" %in% names(S4Vectors::mcols(gr))) {
        "feature"
    } else {
        .stopNoCall("GTF has neither `type` nor `feature` column")
    }
    genes <- gr[S4Vectors::mcols(gr)[[typeCol]] == "gene"]
    if (length(genes) == 0) .stopNoCall("No `gene` entries in GTF")

    gtfChroms <- unique(as.character(GenomicRanges::seqnames(genes)))
    ## Adjust the chromosome names to match the GTF
    chromQ <- if (any(grepl("^chr", gtfChroms))) {
        (if (grepl("^chr", chrom)) chrom else paste0("chr", chrom))
    } else {
        sub("^chr", "", chrom)
    }

    keep <- as.character(GenomicRanges::seqnames(genes)) == chromQ &
        GenomicRanges::end(genes) >= start &
        GenomicRanges::start(genes) <= end
    genes <- genes[keep]
    if (!is.null(biotypes) && "gene_biotype" %in% names(S4Vectors::mcols(genes))) {
        genes <- genes[as.character(S4Vectors::mcols(genes)$gene_biotype) %in% biotypes]
    }

    nameCol <- if ("gene_name" %in% names(S4Vectors::mcols(genes))) {
        "gene_name"
    } else if ("gene_id" %in% names(S4Vectors::mcols(genes))) {
        "gene_id"
    } else {
        NULL
    }
    data.frame(
        gene = if (is.null(nameCol)) {
            NA_character_
        } else {
            as.character(S4Vectors::mcols(genes)[[nameCol]])
        },
        start = pmax(start, GenomicRanges::start(genes)),
        end = pmin(end, GenomicRanges::end(genes)),
        strand = as.character(GenomicRanges::strand(genes)),
        stringsAsFactors = FALSE
    )
}

## Read exons of the given gene names in the region.
.load_exons_in_region <- function(gtf, chrom, start, end, geneNames) {
    if (length(geneNames) == 0) {
        return(data.frame(
            gene = character(0), start = integer(0),
            end = integer(0), strand = character(0),
            stringsAsFactors = FALSE
        ))
    }
    gr <- rtracklayer::import(gtf)
    typeCol <- if ("type" %in% names(S4Vectors::mcols(gr))) {
        "type"
    } else if ("feature" %in% names(S4Vectors::mcols(gr))) {
        "feature"
    } else {
        return(data.frame(
            gene = character(0), start = integer(0),
            end = integer(0), strand = character(0),
            stringsAsFactors = FALSE
        ))
    }
    exons <- gr[S4Vectors::mcols(gr)[[typeCol]] == "exon"]
    if (length(exons) == 0) {
        return(data.frame(
            gene = character(0), start = integer(0),
            end = integer(0), strand = character(0),
            stringsAsFactors = FALSE
        ))
    }
    gtfChroms <- unique(as.character(GenomicRanges::seqnames(exons)))
    chromQ <- if (any(grepl("^chr", gtfChroms))) {
        (if (grepl("^chr", chrom)) chrom else paste0("chr", chrom))
    } else {
        sub("^chr", "", chrom)
    }
    nameCol <- if ("gene_name" %in% names(S4Vectors::mcols(exons))) {
        "gene_name"
    } else if ("gene_id" %in% names(S4Vectors::mcols(exons))) {
        "gene_id"
    } else {
        NULL
    }
    if (is.null(nameCol)) {
        return(data.frame(
            gene = character(0), start = integer(0),
            end = integer(0), strand = character(0),
            stringsAsFactors = FALSE
        ))
    }
    keep <- as.character(GenomicRanges::seqnames(exons)) == chromQ &
        GenomicRanges::end(exons) >= start &
        GenomicRanges::start(exons) <= end &
        as.character(S4Vectors::mcols(exons)[[nameCol]]) %in% geneNames
    exons <- exons[keep]
    if (length(exons) == 0) {
        return(data.frame(
            gene = character(0), start = integer(0),
            end = integer(0), strand = character(0),
            stringsAsFactors = FALSE
        ))
    }
    data.frame(
        gene = as.character(S4Vectors::mcols(exons)[[nameCol]]),
        start = pmax(start, GenomicRanges::start(exons)),
        end = pmin(end, GenomicRanges::end(exons)),
        strand = as.character(GenomicRanges::strand(exons)),
        stringsAsFactors = FALSE
    )
}

## Non-overlapping row assignment: greedy first-fit into rows so text labels don't collide.
.assign_gene_rows <- function(gDF) {
    n <- nrow(gDF)
    if (n == 0) {
        return(integer(0))
    }
    ord <- order(gDF$start)
    rows <- integer(n)
    rowEnds <- numeric(0) # last end per row
    for (k in ord) {
        placed <- FALSE
        for (r in seq_along(rowEnds)) {
            if (gDF$start[k] > rowEnds[r]) {
                rows[k] <- r
                rowEnds[r] <- gDF$end[k]
                placed <- TRUE
                break
            }
        }
        if (!placed) {
            rows[k] <- length(rowEnds) + 1
            rowEnds <- c(rowEnds, gDF$end[k])
        }
    }
    rows
}

.build_gene_track <- function(geneDF, exonDF, start, end, chrom, axisFormat) {
    if (nrow(geneDF) == 0) {
        return(
            ggplot2::ggplot() +
                ggplot2::geom_blank() +
                ggplot2::coord_cartesian(xlim = c(start, end)) +
                ggplot2::scale_x_continuous(labels = axisFormat, expand = c(0, 0)) +
                ggplot2::labs(y = "Genes") +
                ggplot2::theme_classic(base_size = 10) +
                ggplot2::theme(
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(),
                    axis.title.x = ggplot2::element_blank()
                )
        )
    }
    geneDF$row <- .assign_gene_rows(geneDF)

    ## Attach row assignment to exons via gene name
    if (nrow(exonDF) > 0) {
        exonDF$row <- geneDF$row[match(exonDF$gene, geneDF$gene)]
        exonDF <- exonDF[!is.na(exonDF$row), , drop = FALSE]
    }

    ## Strand direction chevrons along the gene body (thin line)
    ## ">" for + strand, "<" for - strand, placed every ~5% of region width
    span <- end - start
    step <- max(as.integer(span / 30), 1)
    chev <- do.call(rbind, lapply(seq_len(nrow(geneDF)), function(i) {
        lo <- geneDF$start[i] + step
        hi <- geneDF$end[i] - step
        if (hi <= lo) {
            return(NULL)
        } # gene too short for chevrons
        xs <- seq(lo, hi, by = step)
        if (length(xs) == 0) {
            return(NULL)
        }
        data.frame(
            x = xs, y = geneDF$row[i],
            lab = if (geneDF$strand[i] == "-") "<" else ">",
            stringsAsFactors = FALSE
        )
    }))

    gg <- ggplot2::ggplot() +
        ## thin line for the whole gene body (intron backbone)
        ggplot2::geom_segment(
            data = geneDF,
            ggplot2::aes(
                x = .data$start, xend = .data$end,
                y = .data$row, yend = .data$row
            ),
            color = "steelblue4", linewidth = 0.3
        )
    ## strand chevrons on the intron line
    if (!is.null(chev) && nrow(chev) > 0) {
        gg <- gg + ggplot2::geom_text(
            data = chev,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
            color = "steelblue4", size = 2.5
        )
    }
    ## Thick exon rectangles
    if (nrow(exonDF) > 0) {
        gg <- gg + ggplot2::geom_rect(
            data = exonDF,
            ggplot2::aes(
                xmin = .data$start, xmax = .data$end,
                ymin = .data$row - 0.28, ymax = .data$row + 0.28
            ),
            fill = "steelblue4", color = NA
        )
    }
    ## Gene name labels above each gene row
    topRow <- max(geneDF$row)
    gg <- gg + ggplot2::geom_text(
        data = geneDF,
        ggplot2::aes(
            x = (.data$start + .data$end) / 2,
            y = .data$row + 0.55,
            label = .data$gene
        ),
        size = 3, fontface = "italic"
    )

    gg +
        ggplot2::coord_cartesian(
            xlim = c(start, end),
            ylim = c(0.3, topRow + 0.9)
        ) +
        ## Axis at TOP (IGV / UCSC ruler convention)
        ggplot2::scale_x_continuous(
            labels = axisFormat, expand = c(0, 0),
            position = "top"
        ) +
        ## IGV-style stacked left-margin labels: bold `chr1` tag on top,
        ## `Genes` axis label below.
        ggplot2::labs(y = "Genes", tag = chrom) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(
                size = 10,
                angle = 0, vjust = 0.5,
                margin = ggplot2::margin(r = 4)
            ),
            axis.title.x = ggplot2::element_blank(),
            ## bottom edge fully hidden (axis is on top now)
            axis.line.x.bottom = ggplot2::element_blank(),
            axis.text.x.bottom = ggplot2::element_blank(),
            axis.ticks.x.bottom = ggplot2::element_blank(),
            ## top-left chrom tag (same size as the Genes axis label, unbold)
            plot.tag = ggplot2::element_text(size = 10, hjust = 0),
            plot.tag.position = c(0, 0.9)
        )
}

.build_snp_track <- function(snpPos, snpID, start, end, axisFormat) {
    ggplot2::ggplot() +
        ggplot2::geom_vline(
            xintercept = snpPos,
            color = "#FF7917", linewidth = 0.6
        ) +
        ggplot2::annotate("text",
            x = snpPos, y = 0.5,
            label = snpID, hjust = -0.05,
            size = 3, fontface = "bold"
        ) +
        ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(0, 1)) +
        ggplot2::scale_x_continuous(labels = axisFormat, expand = c(0, 0)) +
        ggplot2::labs(y = "SNP") +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank()
        )
}

## Generate parabolic arcs (t in [0, 1]; y = 4*t*(1-t) scaled by loop size).
.arcPath <- function(x1, x2) {
    t <- seq(0, 1, length.out = 60)
    loopSize <- abs(x2 - x1)

    return(data.frame(
        x = x1 + t * (x2 - x1),
        y = -(4 * t * (1 - t)) * loopSize
    ))
}

.build_arc_panel <- function(arcDF, start, end, axisFormat,
                             effLo, effHi, strongColor, weakColor,
                             panelTitle, showX = FALSE) {
    ## Expand each arc into a path of points
    paths <- do.call(rbind, lapply(seq_len(nrow(arcDF)), function(i) {
        p <- .arcPath(arcDF$x1[i], arcDF$x2[i])
        p$loopID <- arcDF$loopID[i]
        p$value <- arcDF$value[i]
        p
    }))
    yMin <- min(paths$y, na.rm = TRUE)
    if (!is.finite(yMin) || yMin >= 0) yMin <- -1
    gg <- ggplot2::ggplot(paths, ggplot2::aes(
        x = .data$x, y = .data$y,
        group = .data$loopID,
        color = .data$value
    )) +
        ggplot2::geom_path(linewidth = 0.9, na.rm = TRUE) +
        ggplot2::scale_color_gradient(
            low = weakColor, high = strongColor,
            limits = c(effLo, effHi),
            oob = scales::squish,
            na.value = weakColor
        ) +
        ggplot2::coord_cartesian(
            xlim = c(start, end),
            ylim = c(yMin * 1.1, 0)
        ) +
        ggplot2::scale_x_continuous(labels = axisFormat, expand = c(0, 0)) +
        ggplot2::labs(y = panelTitle, color = "Mean strength") +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(size = 9),
            panel.grid = ggplot2::element_blank()
        )
    ## X-axis ruler is on top of the gene track (IGV-style); always hide it here.
    gg + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
    )
}
