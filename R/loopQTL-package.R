#' loopQTL: Anchor-restricted loop QTL mapping from Hi-C
#'
#' loopQTL identifies SNPs associated with chromatin interaction signals
#' across samples. It implements a band-wise joint normalization approach
#' (referenced from Brant et al.) that systematically tests the association
#' between variants and their anchored loops, together with downstream
#' functions for annotation and visualization for the user's convenience.
#'
#' The typical workflow is:
#' \enumerate{
#'   \item [loopQTLsetup()] - build a consensus loop set and container
#'     from per-sample `.hic` files and a VCF.
#'   \item [loopQTLnormalize()] - compute per-loop phenotype values via
#'     band-wise quantile normalization.
#'   \item [runQTL()] - run anchor-restricted QTL tests.
#'   \item [loopQTLannotate()] - annotate SNP-loop pairs against a GTF
#'     (and optional peak set).
#'   \item [loopQTLplotLocus()], [loopQTLplotPair()],
#'     [loopQTLplotSnpWindow()] - visualize results.
#' }
#'
#' @author
#' Zexun Wu,
#' Andrew Vu,
#' Matthew Salcedo,
#' Ethan Nelson-Moore,
#' Suhn K. Rhie*
#'
#' @keywords internal
#' @name loopQTL-package
#' @aliases loopQTL loopQTL-package
#' @importFrom methods new is slot slotNames validObject
#' @importMethodsFrom methods show
#' @importFrom stats prcomp qnorm setNames
#' @importFrom utils head
"_PACKAGE"

## Suppress R CMD check NOTEs about undefined globals used by data.table NSE
## and by ggplot2 aes() via the .data pronoun.
utils::globalVariables(c(".data", "i", "x", "."))
