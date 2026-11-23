#!/usr/bin/env Rscript
## generate_pseudo_hic.R
##
## Create N pseudo Hi-C sample pixel tables from a single source BG2 by
## independently perturbing each non-zero pixel's raw count. The
## perturbation for each pixel is sign * magnitude with:
##   sign      ~ Uniform({-1, +1})
##   magnitude ~ Uniform(1..10)
## so delta is drawn uniformly from {-10,...,-1, 1,...,10} (never 0).
## Perturbed counts are clamped at 0; zero-count pixels are dropped.
##
## While writing the per-sample BG2, chromosome names are rewritten
## from "<N>" style to "chr<N>" style so the resulting .hic files match
## GRCh38 / hg38 conventions used in GTF / VCF annotations.
##
## Usage:
##   generate_pseudo_hic.R <source_pixels.bg2> <outdir> <nsamples> [seed]
##
## Inputs
##   source_pixels.bg2  Output of `hictk dump --join` on the source .hic
##                      restricted to a single chromosome at one resolution.
##                      Columns: chrA startA endA chrB startB endB count.
##   outdir             Directory to write sample1.bg2 .. sampleN.bg2
##   nsamples           Positive integer
##   seed               Optional integer for reproducibility (default: 42)

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
    stop("Usage: generate_pseudo_hic.R <source.bg2> <outdir> <nsamples> [seed]",
         call. = FALSE)
}

src_bg2 <- args[1]
outdir  <- args[2]
N       <- as.integer(args[3])
seed    <- if (length(args) >= 4L) as.integer(args[4]) else 42L

stopifnot(
    file.exists(src_bg2),
    is.finite(N), N > 0L
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

cat(sprintf("[%s] Reading %s ...\n",
            format(Sys.time(), "%H:%M:%S"), src_bg2))
dt <- fread(src_bg2, header = FALSE, sep = "\t",
            col.names = c("chrA","startA","endA",
                          "chrB","startB","endB","count"))
cat(sprintf("[%s] Loaded %s pixels (unique chroms: %s)\n",
            format(Sys.time(), "%H:%M:%S"),
            format(nrow(dt), big.mark = ","),
            paste(unique(dt$chrA), collapse = ",")))

## Rewrite "N" -> "chrN" so the output .hic uses standard hg38 chrom names
chrA_out <- paste0("chr", dt$chrA)
chrB_out <- paste0("chr", dt$chrB)

n_pix <- nrow(dt)
for (i in seq_len(N)) {
    ## Per-pixel perturbation: sign in {-1, +1}, magnitude in 1..10
    sign_i <- sample(c(-1L, 1L), n_pix, replace = TRUE)
    mag_i  <- sample.int(10L, n_pix, replace = TRUE)
    delta  <- sign_i * mag_i

    perturbed <- pmax(dt$count + delta, 0L)
    keep <- perturbed > 0L

    out_bg2 <- file.path(outdir, sprintf("sample%d.bg2", i))
    cat(sprintf("[%s] sample%d: writing %s non-zero pixels to %s\n",
                format(Sys.time(), "%H:%M:%S"), i,
                format(sum(keep), big.mark = ","), out_bg2))

    fwrite(
        data.table(
            chrA_out[keep], dt$startA[keep], dt$endA[keep],
            chrB_out[keep], dt$startB[keep], dt$endB[keep],
            perturbed[keep]
        ),
        file = out_bg2, sep = "\t", col.names = FALSE, quote = FALSE
    )
}

cat(sprintf("[%s] done\n", format(Sys.time(), "%H:%M:%S")))
