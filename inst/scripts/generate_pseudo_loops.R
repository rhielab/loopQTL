#!/usr/bin/env Rscript
## generate_pseudo_loops.R
##
## Create N pseudo per-sample loop BEDPE files from a single source
## BEDPE by, for each sample:
##   (1) randomly dropping a Uniform(5%, 10%) fraction of the loops
##   (2) then, for a random 20% of the loops that remain, shifting each
##       anchor independently by +/- one bin (5 kb), while keeping the
##       count and q-value fields unchanged.
##
## The source BEDPE is assumed to be the FitHiC output with 8 columns
##   chrA startA endA chrB startB endB count q-value
## and the loops are filtered to the target chromosome before any
## perturbation, so that every loop references bins that exist in the
## pseudo .hic files (which are chr2-only at 5 kb).
##
## Usage:
##   generate_pseudo_loops.R <source.bedpe> <target_chrom> <outdir> \
##                           <nsamples> <bin_bp> [seed]
##
## Example:
##   generate_pseudo_loops.R \
##     C42B..._FitHiC_q0.05_30KB-2Mb.bedpe chr2 chr2_loops/ 10 5000

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
    stop(
        "Usage: generate_pseudo_loops.R <source.bedpe> <target_chrom> ",
        "<outdir> <nsamples> <bin_bp> [seed]",
        call. = FALSE
    )
}
src_bedpe <- args[1]
target_chrom <- args[2]
outdir <- args[3]
N <- as.integer(args[4])
bin_bp <- as.integer(args[5])
seed <- if (length(args) >= 6L) as.integer(args[6]) else 42L

stopifnot(
    file.exists(src_bedpe),
    is.finite(N), N > 0L,
    is.finite(bin_bp), bin_bp > 0L
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

cat(sprintf("[%s] Reading %s ...\n",
            format(Sys.time(), "%H:%M:%S"), src_bedpe))
dt <- fread(src_bedpe, header = FALSE, sep = "\t",
            col.names = c("chrA", "startA", "endA",
                          "chrB", "startB", "endB",
                          "count", "qvalue"))
cat(sprintf("[%s] Loaded %s loops across %d chromosomes\n",
            format(Sys.time(), "%H:%M:%S"),
            format(nrow(dt), big.mark = ","),
            length(unique(dt$chrA))))

## Keep only loops on the target chromosome (both anchors)
dt <- dt[chrA == target_chrom & chrB == target_chrom]
if (nrow(dt) == 0L) {
    stop("No loops found on target chromosome: ", target_chrom, call. = FALSE)
}
cat(sprintf("[%s] %s loops remain on %s\n",
            format(Sys.time(), "%H:%M:%S"),
            format(nrow(dt), big.mark = ","), target_chrom))

n_src <- nrow(dt)

for (i in seq_len(N)) {
    ## (1) drop a Uniform(5%, 10%) fraction
    drop_frac <- runif(1L, min = 0.05, max = 0.10)
    n_drop <- round(n_src * drop_frac)
    drop_idx <- sample.int(n_src, n_drop, replace = FALSE)
    kept <- dt[-drop_idx]
    n_kept <- nrow(kept)

    ## (2) shift 20% of the remaining loops by +/- one bin per anchor,
    ##     independently for each anchor.
    n_shift <- round(n_kept * 0.20)
    shift_idx <- sample.int(n_kept, n_shift, replace = FALSE)

    dA <- sample(c(-1L, 1L), n_shift, replace = TRUE) * bin_bp
    dB <- sample(c(-1L, 1L), n_shift, replace = TRUE) * bin_bp

    kept[shift_idx, `:=`(
        startA = pmax(startA + dA, 0L),
        endA   = pmax(endA   + dA, 0L),
        startB = pmax(startB + dB, 0L),
        endB   = pmax(endB   + dB, 0L)
    )]

    ## Keep A <= B canonicalization in case a shift reversed the order
    swap <- kept$startA > kept$startB
    if (any(swap)) {
        tmp_s <- kept$startA[swap]; tmp_e <- kept$endA[swap]
        kept[swap, `:=`(
            startA = startB, endA = endB,
            startB = tmp_s,  endB = tmp_e
        )]
    }

    out_bedpe <- file.path(outdir, sprintf("sample%d_loops.bedpe", i))
    cat(sprintf(
        "[%s] sample%d: dropped %s (%.1f%%), shifted %s (20%%); writing %s\n",
        format(Sys.time(), "%H:%M:%S"), i,
        format(n_drop, big.mark = ","), drop_frac * 100,
        format(n_shift, big.mark = ","), out_bedpe
    ))

    fwrite(kept, file = out_bedpe, sep = "\t",
           col.names = FALSE, quote = FALSE)
}

cat(sprintf("[%s] done\n", format(Sys.time(), "%H:%M:%S")))
