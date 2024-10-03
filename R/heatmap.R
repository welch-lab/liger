#' Plot Heatmap of Gene Expression or Factor Loading
#' @param object A \linkS4class{liger} object, with data to be plot available.
#' @param features,factors Character vector of genes of interests or numeric
#' index of factor to be involved. \code{features} is required, while
#' \code{factors} is by default all the factors (reads object recorded k value
#' in \code{uns} slot).
#' @param cellIdx Valid index to subscribe cells to be included. See
#' \code{\link{subsetLiger}}. Default \code{NULL} use all cells.
#' @param slot Use the chosen matrix for heatmap. For \code{plotGeneHeatmap},
#' default \code{"normData"}, alternatively \code{"rawData"},
#' \code{"scaleData"} or \code{"scaleUnsharedData"}. For
#' \code{plotFactorHeatmap}, default \code{"H.norm"}, alternatively \code{"H"}.
#' @param useCellMeta Character vector of available variable names in
#' \code{cellMeta}, variables will be added as annotation to the heatmap.
#' Default \code{NULL}.
#' @param cellAnnotation data.frame object for using external annotation, with
#' each column a variable and each row is a cell. Row names of this data.frame
#' will be used for matching cells involved in heatmap. For cells not found in
#' this data.frame, \code{NA}s will be added with warning. Default \code{NULL}.
#' @param featureAnnotation,factorAnnotation Similar as \code{cellAnnotation},
#' while each row would be a gene or factor, respectively. Default \code{NULL}.
#' @param cellSplitBy Character vector of variable names available in annotation
#' given by \code{useCellMeta} and \code{cellAnnotation}. This slices the
#' heatmap by specified variables. Default \code{NULL}.
#' @param featureSplitBy,factorSplitBy Similar as \code{cellSplitBy}. Default
#' \code{NULL}
#' @param trim Numeric vector of two numbers. Higher value limits the maximum
#' value and lower value limits the minimum value. Default \code{c(0, 0.03)}.
#' @param viridisOption See \code{option} argument of
#' \code{\link[viridisLite]{viridis}}. Default \code{"C"} (plasma) for
#' \code{plotGeneHeatmap} and \code{"D"} (viridis) for \code{plotFactorHeatmap}.
#' @inheritDotParams .plotHeatmap transpose showCellLabel showCellLegend showFeatureLabel showFeatureLegend cellAnnColList featureAnnColList scale baseSize cellTextSize featureTextSize cellTitleSize featureTitleSize legendTextSize legendTitleSize viridisDirection RColorBrewerOption
#' @return \code{\link[ComplexHeatmap]{HeatmapList-class}} object
#' @export
#' @rdname plotHeatmap
#' @examples
#' \donttest{
#' plotGeneHeatmap(pbmcPlot, varFeatures(pbmcPlot))
#' plotGeneHeatmap(pbmcPlot, varFeatures(pbmcPlot),
#'                 useCellMeta = c("leiden_cluster", "dataset"),
#'                 cellSplitBy = "leiden_cluster")
#'
#' plotFactorHeatmap(pbmcPlot)
#' plotFactorHeatmap(pbmcPlot, cellIdx = pbmcPlot$leiden_cluster %in% 1:3,
#'                   useCellMeta = c("leiden_cluster", "dataset"),
#'                   cellSplitBy = "leiden_cluster")
#' }
plotGeneHeatmap <- function(
        object,
        features,
        cellIdx = NULL,
        slot = c("normData", "rawData", "scaleData", "scaleUnsharedData"),
        useCellMeta = NULL,
        cellAnnotation = NULL,
        featureAnnotation = NULL,
        cellSplitBy = NULL,
        featureSplitBy = NULL,
        viridisOption = "C",
        ...
) {
    slot <- match.arg(slot)
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    hmData <- retrieveCellFeature(object, slot = slot, cellIdx = cellIdx,
                                  feature = features)
    if (slot == "normData") dataScaleFunc <- function(x) log2(10000*x + 1)
    else dataScaleFunc <- NULL
    # Organize annotation

    cellAnn <- .formatAnn(object, charIdx = colnames(object)[cellIdx],
                          useMeta = useCellMeta,
                          annDF = cellAnnotation,
                          splitBy = cellSplitBy)
    featureAnn <- .formatAnn(object, charIdx = features,
                             useMeta = NULL, annDF = featureAnnotation,
                             splitBy = featureSplitBy)

    .plotHeatmap(dataMatrix = t(hmData), dataName = "Gene\nExpression",
                 cellDF = cellAnn$ann,
                 featureDF = featureAnn$ann,
                 cellSplitVar = cellAnn$split,
                 featureSplitVar = featureAnn$split,
                 dataScaleFunc = dataScaleFunc,
                 viridisOption = viridisOption,
                 ...)
}

#' @export
#' @rdname plotHeatmap
plotFactorHeatmap <- function(
        object,
        factors = NULL,
        cellIdx = NULL,
        slot = c("H.norm", "H"),
        useCellMeta = NULL,
        cellAnnotation = NULL,
        factorAnnotation = NULL,
        cellSplitBy = NULL,
        factorSplitBy = NULL,
        trim = c(0, 0.03),
        viridisOption = "D",
        ...
) {
    slot <- match.arg(slot)
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    if (is.null(factors)) factors <- seq(object@uns$factorization$k)
    hmData <- retrieveCellFeature(object, slot = slot, cellIdx = cellIdx,
                                  feature = factors)
    cellAnn <- .formatAnn(object, charIdx = colnames(object)[cellIdx],
                          useMeta = useCellMeta,
                          annDF = cellAnnotation,
                          splitBy = cellSplitBy)
    featureAnn <- .formatAnn(object, charIdx = factors,
                             useMeta = NULL, annDF = factorAnnotation,
                             splitBy = factorSplitBy)

    .plotHeatmap(dataMatrix = t(hmData), dataName = "Factor\nLoading",
                 cellDF = cellAnn$ann,
                 featureDF = featureAnn$ann,
                 cellSplitVar = cellAnn$split,
                 featureSplitVar = featureAnn$split,
                 viridisOption = viridisOption,
                 trim = trim,
                 ...)
}


#' General heatmap plotting with prepared matrix and data.frames
#' @description This is not an exported function. This documentation just
#' serves for a manual of extra arguments that users can use when generating
#' heatmaps with \code{\link{plotGeneHeatmap}} or
#' \code{\link{plotFactorHeatmap}}.
#'
#' Note that the following arguments are pre-occupied by upstream wrappers so
#' users should not include them in a function call: \code{dataMatrix},
#' \code{dataName}, \code{cellDF}, \code{featureDF}, \code{cellSplitVar},
#' \code{featureSplitVar}.
#'
#' The following arguments of \code{\link[ComplexHeatmap]{Heatmap}} is occupied
#' by this function, so users should include them in a function call as well:
#' \code{matrix}, \code{name}, \code{col}, \code{heatmap_legend_param},
#' \code{top_annotation}, \code{column_title_gp}, \code{column_names_gp},
#' \code{show_column_names}, \code{column_split}, \code{column_gap},
#' \code{left_annotation}, \code{row_title_gp}, \code{row_names_gp},
#' \code{show_row_names}, \code{row_split}, \code{row_gap}.
#' @param dataMatrix Matrix object with features/factors as rows and cells as
#' columns.
#' @param dataName Text for heatmap color bar title. Default \code{Value}.
#' @param cellDF data.frame object. Number of rows must match with number of
#' columns of \code{dataMatrix}.
#' @param featureDF data.frame object. Number of columns must match with number
#' of rows of \code{dataMatrix}.
#' @param transpose Logical, whether to "rotate" the heatmap by 90 degrees so
#' that cell information is displayed by row. Default \code{FALSE}.
#' @param cellSplitVar,featureSplitVar Subset columns of \code{cellDF} or
#' \code{featureDF}, respectively.
#' @param dataScaleFunc A function object, applied to \code{dataMatrix}.
#' @param showCellLabel,showFeatureLabel Logical, whether to show cell barcodes,
#' gene symbols or factor names. Default \code{TRUE} for gene/factors but
#' \code{FALSE} for cells.
#' @param showCellLegend,showFeatureLegend Logical, whether to show cell or
#' feature legends. Default \code{TRUE}. Can be a scalar for overall control
#' or a vector matching with each given annotation variable.
#' @param cellAnnColList,featureAnnColList List object, with each element a
#' named vector of R-interpretable color code. The names of the list elements
#' are used for matching the annotation variable names. The names of the colors
#' in the vectors are used for matching the levels of a variable (factor object,
#' categorical). Default \code{NULL} generates ggplot-flavor categorical colors.
#' @param scale Logical, whether to take z-score to scale and center gene
#' expression. Applied after \code{dataScaleFunc}. Default \code{FALSE}.
#' @param trim Numeric vector of two values. Limit the z-score value into this
#' range when \code{scale = TRUE}. Default \code{c(-2, 2)}.
#' @param baseSize One-parameter control of all text sizes. Individual text
#' element sizes can be controlled by other size arguments. "Title" sizes are
#' 2 points larger than "text" sizes when being controlled by this.
#' @param cellTextSize,featureTextSize,legendTextSize Size of cell barcode
#' labels, gene/factor labels, or legend values. Default \code{NULL}.
#' @param cellTitleSize,featureTitleSize,legendTitleSize Size of titles of the
#' cell slices, gene/factor slices, or the legends. Default \code{NULL}.
#' @param viridisOption,viridisDirection See argument \code{option} and
#' \code{direction} of \code{\link[viridisLite]{viridis}}. Default \code{"A"}
#' and \code{-1}.
#' @param RColorBrewerOption When \code{scale = TRUE}, heatmap color will be
#' mapped with \code{\link[RColorBrewer]{brewer.pal}}. This is passed to
#' \code{name}. Default \code{"RdBu"}.
#' @param ... Additional arguments to be passed to
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#' @return \code{\link[ComplexHeatmap]{HeatmapList-class}} object
.plotHeatmap <- function(
        dataMatrix,
        dataName = "Value",
        cellDF = NULL,
        featureDF = NULL,
        transpose = FALSE,
        cellSplitVar = NULL,
        featureSplitVar = NULL,
        dataScaleFunc = NULL,
        showCellLabel = FALSE,
        showCellLegend = TRUE,
        showFeatureLabel = TRUE,
        showFeatureLegend = TRUE,
        cellAnnColList = NULL,
        featureAnnColList = NULL,
        scale = FALSE,
        trim = c(-2, 2),
        baseSize = 8,
        cellTextSize = NULL,
        featureTextSize = NULL,
        cellTitleSize = NULL,
        featureTitleSize = NULL,
        legendTextSize = NULL,
        legendTitleSize = NULL,
        viridisOption = "A",
        viridisDirection = -1,
        RColorBrewerOption = "RdBu",
        ...
) {
    # Final data processing
    if (!is.null(dataScaleFunc)) dataMatrix <- dataScaleFunc(dataMatrix)
    if (isTRUE(scale)) {
        dataMatrix <- .zScore(dataMatrix, trim = trim)
    }
    # Viridis option checks
    viridisAvail <- c(
        "magma", "A", "inferno", "B", "plasma", "C", "viridis", "D",
        "cividis", "E", "rocket", "F", "mako", "G", "turbo", "H"
    )
    if (length(viridisOption) != 1 ||
        !viridisOption %in% viridisAvail)
        cli::cli_abort(
            c("{.var viridisOption} has to be one value from the available choices: ",
              "{.val {viridisAvail}}")
        )

    ## Font-size specification
    # Broadcast one-param setting to each
    cellText <- featureText <- legendText <- baseSize
    cellTitle <- featureTitle <- legendTitle <- baseSize + 2
    # And set specific ones if specified
    if (!is.null(cellTextSize)) cellText <- cellTextSize
    if (!is.null(cellTitleSize)) cellTitle <- cellTitleSize
    if (!is.null(featureTextSize)) featureText <- featureTextSize
    if (!is.null(featureTitleSize)) featureTitle <- featureTitleSize
    if (!is.null(legendTextSize)) legendText <- legendTextSize
    if (!is.null(legendTitleSize)) legendTitle <- legendTitleSize

    ## Customized color mapping requires a function object returned by
    ## colorRamp2
    if (!isTRUE(scale)) {
        # If not scaled (min = 0)
        col_fun <- circlize::colorRamp2(
            breaks = c(0, max(dataMatrix) / 2, max(dataMatrix)),
            colors = viridis::viridis(
                n = 20, option = viridisOption,
                direction = viridisDirection
            )[c(1, 10, 20)]
        )
    } else {
        # If scaling/centering the matrix (z-score), use blue-white-red color
        # palette and let center value (zero) be white
        col_fun <- circlize::colorRamp2(
            breaks = c(min(dataMatrix), 0, max(dataMatrix)),
            colors = RColorBrewer::brewer.pal(
                9, RColorBrewerOption
            )[c(8, 5, 2)]
        )
    }
    ## Construct HeatmapAnnotation
    cellHA <- .constructHA(cellDF, legendTitleSize = legendTitle,
                           legendTextSize = legendText,
                           which = ifelse(transpose, "row", "column"),
                           showLegend = showCellLegend,
                           colList = cellAnnColList)
    featureHA <- .constructHA(featureDF, legendTitleSize = legendTitle,
                              legendTextSize = legendText,
                              which = ifelse(transpose, "column", "row"),
                              showLegend = showFeatureLegend,
                              colList = featureAnnColList)

    if (!isTRUE(transpose)) {
        hm <- ComplexHeatmap::Heatmap(
            # General settings
            matrix = dataMatrix, name = dataName,
            col = col_fun,
            heatmap_legend_param = list(
                title_gp = grid::gpar(fontsize = legendTitle,
                                      fontface = "bold"),
                labels_gp = grid::gpar(fontsize = legendText)
            ),

            # Column settings
            top_annotation = cellHA,
            column_title_gp = grid::gpar(fontsize = cellTitle),
            column_names_gp = grid::gpar(fontsize = cellText),
            show_column_names = showCellLabel,
            column_split = cellSplitVar,
            column_gap = grid::unit(0, "mm"),

            # Row settings
            left_annotation = featureHA,
            row_names_gp = grid::gpar(fontsize = featureText),
            row_title_gp = grid::gpar(fontsize = featureTitle),
            show_row_names = showFeatureLabel,
            row_split = featureSplitVar,
            row_gap = grid::unit(0, "mm"),
            ...
        )
    } else {
        hm <- ComplexHeatmap::Heatmap(
            # General settings
            matrix = t(dataMatrix), name = dataName,
            col = col_fun,
            heatmap_legend_param = list(
                title_gp = grid::gpar(fontsize = legendTitle,
                                      fontface = "bold"),
                labels_gp = grid::gpar(fontsize = legendText)
            ),

            # Column settings
            top_annotation = featureHA,
            column_title_gp = grid::gpar(fontsize = featureTitle),
            column_names_gp = grid::gpar(fontsize = featureText),
            show_column_names = showFeatureLabel,
            column_split = featureSplitVar,
            column_gap = grid::unit(0, "mm"),

            # Row settings
            left_annotation = cellHA,
            row_title_gp = grid::gpar(fontsize = cellTitle),
            row_names_gp = grid::gpar(fontsize = cellText),
            show_row_names = showCellLabel,
            row_split = cellSplitVar,
            row_gap = grid::unit(0, "mm"),
            ...
        )
    }
    grDevices::pdf(nullfile())
    hml <- ComplexHeatmap::draw(hm, merge_legend = TRUE)
    grDevices::dev.off()
    return(hml)
}

.formatAnn <- function(
        object,
        charIdx,
        useMeta = NULL,
        annDF = NULL,
        splitBy = NULL
) {
    ### Check and format the information in cellMeta
    # TODO: Only have cellMeta for now but not gonna use featureMeta
    # See if also allow featureMeta in the future
    AnnDF <- cellMeta(object, columns = useMeta, cellIdx = charIdx,
                       as.data.frame = TRUE, drop = FALSE)

    ### Check and append customized annotation
    if (inherits(annDF, c("data.frame", "DFrame"))) {
        notFound <- !(charIdx %in% rownames(annDF))
        if (any(notFound))
            warning(sum(notFound), " selected could not be found in ",
                    "given annotation.")
        # Convert to data.frame first so missing value can be filled with NA
        if (!is.data.frame(annDF))
            annDF <- as.data.frame(annDF)
        annDF <- annDF[charIdx, , drop = FALSE]
        if (is.null(AnnDF)) AnnDF <- annDF
        else AnnDF <- cbind(AnnDF, annDF)
    } else if (!is.null(annDF)) {
        warning("Annotation of class ", class(annDF),
                " is not supported yet.")
    }

    if (!is.null(splitBy)) {
        notFound <- !splitBy %in% colnames(AnnDF)
        if (any(notFound))
            warning("Variables in `cell/featureSplitBy` not detected in specified ",
                    "annotation: ",
                    paste(splitBy[notFound], collapse = ", "))
        splitBy <- splitBy[!notFound]
    }
    if (length(splitBy) > 0) cellSplitVar <- AnnDF[,splitBy]
    else cellSplitVar <- NULL

    return(list(ann = AnnDF, split = cellSplitVar))
}

# HA - HeatmapAnnotation()
.constructHA <- function(df, legendTitleSize, legendTextSize,
                         which = c("row", "column"), showLegend = TRUE,
                         colList = NULL) {
    which <- match.arg(which)
    if (!is.null(df) && ncol(df) > 0) {
        annCol <- list()
        for (var in colnames(df)) {
            if (is.factor(df[[var]])) {
                if (var %in% names(colList)) {
                    df[[var]] <- droplevels(df[[var]])
                    if (any(!levels(df[[var]]) %in% names(colList[[var]]))) {
                        cli::cli_abort(
                            "Given customized annotation color must have names matching to all available levels in the annotation."
                        )
                    }
                    annCol[[var]] <- colList[[var]][levels(df[[var]])]
                } else {
                    # Automatic generate with ggplot2 strategy,
                    # with level awareness
                    if (nlevels(df[[var]]) > length(scPalette)) {
                        annCol[[var]] <- scales::hue_pal()(nlevels(df[[var]]))
                    } else {
                        annCol[[var]] <- scPalette[1:nlevels(df[[var]])]
                    }
                    names(annCol[[var]]) <- levels(df[[var]])
                    df[[var]] <- droplevels(df[[var]])
                }
            }
        }
        ha <- ComplexHeatmap::HeatmapAnnotation(
            df = df, which = which, col = annCol,
            show_legend = showLegend,
            annotation_legend_param = list(
                title_gp = grid::gpar(fontsize = legendTitleSize,
                                      fontface = "bold"),
                labels_gp = grid::gpar(fontsize = legendTextSize)
            )
        )
    }
    else ha <- NULL
}

.zScore <- function(x, trim = NULL) {
    if (inherits(x, "dgCMatrix")) {
        m <- Matrix::rowMeans(x)
        v <- rowVars_sparse_rcpp(x, m)
    } else {
        m <- rowMeans(x)
        v <- rowVarsDense(x, m)
    }
    x <- (x - m) / sqrt(v)
    if (!is.null(trim)) {
        if (!is.numeric(trim) || length(trim) != 2)
            warning("`trim` must be a numeric vector of two values")
        else {
            x[x > max(trim)] <- max(trim)
            x[x < min(trim)] <- min(trim)
        }
    }
    return(x)
}


