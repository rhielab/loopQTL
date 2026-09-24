## Modified based on https://github.com/heatherjzhou/PCAForQTL/blob/b8f7d09e77cb5b50dd712ad99f3fbccab899b73d/R/22.01.04_main1.2_runBE.R (retrieved 2026/07/31) - license: GPLv3

## Implements the BE algorithm, a permutation-based approach for choosing K in
## PCA. Intuitively, the BE algorithm retains PCs that explain more variance in
## the data than by random chance and discards those that do not.

## Given X, the data matrix (must be observation by feature),
## B, the number of permutations (default is 20),
## and alpha, the significance level (default is 0.05),
## run BE and return p-values, alpha, and numOfPCsChosen.

## FIXME: For reproducibility, make sure to change the RNG type (necessary
## unless mc.cores is 1) and set the seed before using this function.

.PCAForQTLrunBE <- function(X, B = 20, alpha = 0.05,
                            mc.cores = max(1, min(B, BiocParallel::bpnworkers(BiocParallel::MulticoreParam()) - 1))) {
    if (alpha < 0 || alpha > 1) {
        .stopNoCall("alpha must be between 0 and 1.")
    }

    n <- nrow(X) # Number of observations.
    p <- ncol(X) # Number of features.
    d <- min(n, p) # Total number of PCs.

    results <- BiocParallel::bplapply(seq_len(B), FUN = function(b) {
        ## Permute each column of X. That is, permute the observations in each
        ## feature.
        XPermuted <- matrix(data = NA, nrow = n, ncol = p)
        for (j in seq_len(p)) {
            XPermuted[, j] <- sample(x = X[, j], size = n, replace = FALSE)
        }

        prcompResultPerm <- prcomp(x = XPermuted, center = TRUE, scale. = TRUE)
        importanceTablePerm <- summary(prcompResultPerm)$importance
        PVEsPerm <- importanceTablePerm[2, ]
        return(PVEsPerm)
    }, BPPARAM = BiocParallel::MulticoreParam(workers = mc.cores))
    ## results is a list of vectors
    temp <- unlist(results)
    testStatsPerm <- matrix(data = temp, nrow = d, byrow = FALSE)

    prcompResult <- prcomp(x = X, center = TRUE, scale. = TRUE)
    importanceTable <- summary(prcompResult)$importance
    PVEs <- importanceTable[2, ]
    pValues <- (rowSums(testStatsPerm >= PVEs) + 1) / (B + 1)

    ## The p-value for the jth PC is calculated as, roughly speaking, the
    ## proportion of permutations where the PVE of the jth PC is greater than or
    ## equal to PVE_j.
    for (j in 2:d) {
        ## Enforce monotone increase of the p-values.
        if (pValues[j] < pValues[j - 1]) {
            pValues[j] <- pValues[j - 1]
        }
    }

    numOfPCsChosen <- sum(pValues <= alpha)
    toReturn <- list(pValues = pValues, alpha = alpha, numOfPCsChosen = numOfPCsChosen)
    return(toReturn)
}
