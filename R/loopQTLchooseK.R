#' Auto-select the number of PCs for covariates
#'
#' This function runs PCA on the normalized phenotype matrix and picks a
#' number of principal components K to use as covariates in the downstream
#' regression. This function returns a list with the chosen K and the
#' sample × K matrix of PCs.
#'
#' @param obj A loopQTL object with `phenotype` matrix filled by `loopQTLnormalize()`.
#' @param Kmethod `"elbow"` (default) or `"BE"`.
#' @param normal Apply RINT per loop before PCA. Default `TRUE`.
#' @param Kmax The maximum number of PCs returned. Default
#'   `length(samples(obj)) - 2`,
#' @param ... Forwarded to PCAForQTL's `runElbow` or `runBE` function.
#'
#' @return A list:
#'   \describe{
#'     \item{`K`}{Chosen number of PCs (integer).}
#'     \item{`PCs`}{`samples x K` matrix of principal components.}
#'     \item{`method`}{Which K-selection method was used.}
#'   }
#'
#' @export
#' @examples
#' obj <- loopQTLExampleData()
#' Kres <- loopQTLchooseK(obj, Kmethod = "elbow")
#' Kres$K
loopQTLchooseK <- function(obj,
                           Kmethod = c("elbow", "BE"),
                           normal = TRUE,
                           Kmax = NULL,
                           ...) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        length(obj@phenotype) > 0
    )
    Kmethod <- match.arg(Kmethod)

    if (is.null(Kmax)) Kmax <- max(1, length(obj@samples) - 2)
    Kmax <- as.integer(Kmax)

    Y <- obj@phenotype
    if (normal) Y <- .rintMatrix(Y)

    .pca_for_qtl(Y,
        samples = obj@samples, Kmethod = Kmethod,
        Kmax = Kmax, ...
    )
}

## Internal: PCA + K selection on a (loops x samples) matrix. Shared by
## loopQTLchooseK() and runQTL() (so we don't RINT twice when runQTL
## already RINT'd the matrix it's passing).
.pca_for_qtl <- function(Y,
                         samples,
                         Kmethod = c("elbow", "BE"),
                         Kmax = NULL,
                         ...) {
    Kmethod <- match.arg(Kmethod)
    if (is.null(Kmax)) Kmax <- max(1, length(samples) - 2)
    Kmax <- as.integer(Kmax)

    keepRows <- stats::complete.cases(Y)
    if (sum(keepRows) < 2) {
        .stopNoCall("Not enough non-NA loops to run PCA")
    }
    Yt <- t(Y[keepRows, , drop = FALSE])
    pcaResults <- stats::prcomp(Yt, center = TRUE, scale. = FALSE)

    Kraw <- switch(Kmethod,
        elbow = .PCAForQTLrunElbow(prcompResult = pcaResults, ...),
        BE    = .PCAForQTLrunBE(prcompResult = pcaResults, ...)
    )
    Kraw <- as.integer(Kraw)
    if (is.na(Kraw) || Kraw < 1) Kraw <- 1
    Kraw <- min(Kraw, ncol(pcaResults$x))

    K <- min(Kraw, Kmax)
    if (K < Kraw) {
        message(sprintf(
            "PCAForQTL %s chose K = %d, capped to Kmax = %d (nSamples - 2 leaves df for intercept + genotype).",
            Kmethod, Kraw, K
        ))
    }

    PCs <- pcaResults$x[, seq_len(K), drop = FALSE]
    rownames(PCs) <- samples
    colnames(PCs) <- paste0("PC", seq_len(K))

    list(K = K, PCs = PCs, method = Kmethod, Kraw = Kraw, Kmax = Kmax)
}
