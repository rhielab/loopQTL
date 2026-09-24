# loopQTL 0.99.0

Initial submission to Bioconductor.

## New features

`loopQTL` maps chromatin loop QTLs from Hi-C: which SNPs at a loop's
anchors associate with the loop's contact strength across samples. The
package tests only SNPs that fall inside a loop's **anchor bins** rather
than every SNP in a wide cis-window, so signal is concentrated on the
variants biologically relevant to whether that specific loop forms.

The full workflow is documented in the package vignette with below functions included:

- `loopQTLsetup()` -- Validate inputs and build loopQTL object.
- `loopQTLnormalize()` -- Cross-normalization on Hi-C signal.
- `loopQTLchooseK()` -- Auto-select the number of PCs for covariates.
- `runQTL()` -- Anchor-restricted association testing.
- `loopQTLpairData()` -- Sample-level view of one loop × SNP pair.
- `loopQTLannotate()` -- Annotate loopQTL results with target genes and
  regulatory peaks.
- `loopQTLplotPair()` -- Per-sample boxplot for a single SNP-loop pair.
- `loopQTLplotLocus()` -- Locus-level track view of loops, SNPs, and genes.
- `loopQTLplotSnpWindow()` -- Hi-C-style beta heatmap around a SNP.

## Example data

A 10-sample Hi-C test dataset is bundled under `inst/extdata/`. So the vignette, and the example
function documentation run end-to-end without external downloads.
