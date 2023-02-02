plotLigerHeatmap <- function(
        object,
        slot = c("raw.data", "norm.data", "scale.data", "H", "H.norm", "W", "V"),
        cellIdx = NULL,
        geneIdx = NULL,
        factorIdx = NULL,
        cellAnnotation = NULL,
        geneAnnotation = NULL,
        factorAnnotation = NULL,
        cellSplitBy = NULL,
        geneSplitBy = NULL,
        factorSplitBy = NULL,
        ...
) {
    if (slot %in% c("raw.data", "norm.data", "scale.data")) featureIdx = geneIdx
    else featureIdx = factorIdx
    hmData <- .retrieveCellFeature(object, slot = slot, feature = featureIdx)

}
