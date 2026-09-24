## Modified based on https://github.com/heatherjzhou/PCAForQTL/blob/b8f7d09e77cb5b50dd712ad99f3fbccab899b73d/R/22.01.04_main1.1_runElbow.R
## Retrieved 2026/07/31 - license: GPLv3

## Implements the elbow method for choosing K in PCA

## Given the output from running prcomp() on X, select the number of PCs by
## maximizing the distance to the diagonal line (see below for details, or the
## paper: Zhou HJ, Li L, Li Y, Li W, Li JJ. PCA outperforms popular hidden
## variable inference methods for molecular QTL mapping. Genome Biol. 2022 Oct
## 11;23(1):210. doi: 10.1186/s13059-022-02761-4. PMID: 36221136; PMCID:
## PMC9552461.).
.PCAForQTLrunElbow <- function(prcompResult) {
    if (!isa(prcompResult, "prcomp")) {
        .stopNoCall("prcompResult must be a prcomp object returned by the function prcomp().")
    }

    importanceTable <- summary(prcompResult)$importance
    x <- seq_len(ncol(importanceTable)) # PC indices
    y <- importanceTable[2, ] # PVEs

    ## Given x and y, calculate the distance between each point and the diagonal
    ## line (the line connecting the first and last points).
    diagonalLineStartX <- x[1]
    diagonalLineStartY <- y[1]
    diagonalLineEndX <- x[length(x)]
    diagonalLineEndY <- y[length(y)]

    diagonalLineWidth <- diagonalLineEndX - diagonalLineStartX
    diagonalLineHeight <- diagonalLineEndY - diagonalLineStartY

    diagonalLineLength <- sqrt(diagonalLineWidth^2 + diagonalLineHeight^2)
    distancesNumerator <- abs(diagonalLineWidth * (diagonalLineStartY - y) - (diagonalLineStartX - x) * diagonalLineHeight)
    distances <- distancesNumerator / diagonalLineLength

    numOfPCsChosen <- unname(which.max(distances))
    return(numOfPCsChosen)
}
