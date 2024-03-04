#' Make dot plot of gene expression in cell groups
#' @description This function produces dot plots. Each column represent a group
#' of cells specified by \code{groupBy}, each row is a gene specified by
#' \code{features}. The color of dots reflects mean of normalized expression of
#' specified genes in each cell group and sizes reflects the percentage of cells
#' expressing each gene in a group. We utilize
#' \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html}{ComplexHeatmap}
#' for simplified management of adding annotation and slicing subplots. This was
#' inspired by the implementation in
#' \href{https://samuel-marsh.github.io/scCustomize/reference/Clustered_DotPlot.html}{scCustomize}.
#' @details For \code{...}, please notice that arguments \code{colorMat},
#' \code{sizeMat}, \code{featureAnnDF}, \code{cellSplitVar}, \code{cellLabels}
#' and \code{viridisOption} from \code{\link{.complexHeatmapDotPlot}} are
#' already occupied by this function internally. A lot of arguments from
#' \code{\link[ComplexHeatmap]{Heatmap}} have also been occupied: \code{matrix,
#' name, heatmap_legend_param, rect_gp, col, layer_fun, km, border, border_gp,
#' column_gap, row_gap, cluster_row_slices, cluster_rows, row_title_gp,
#' row_names_gp, row_split, row_labels, cluster_column_slices, cluster_columns,
#' column_split, column_title_gp, column_title, column_labels, column_names_gp,
#' top_annotation}.
#' @param object A \linkS4class{liger} object
#' @param features Use a character vector of gene names to make plain dot plot
#' like a heatmap. Use a data.frame where the first column is gene names and
#' second column is a grouping variable (e.g. subset \code{runMarkerDEG} output)
#' @param groupBy The names of the columns in \code{cellMeta} slot storing
#' categorical variables. Expression data would be aggregated basing on these,
#' together with \code{splitBy}. Default uses default clusters.
#' @param splitBy The names of the columns in \code{cellMeta} slot storing
#' categorical variables. Dotplot panel splitting would be based on these.
#' Default \code{NULL}.
#' @param featureScaleFunc A function object applied to normalized data for
#' scaling the value for better visualization. Default \code{function(x)
#' log2(10000*x + 1)}
#' @param cellIdx Valid cell subscription. See \code{\link{subsetLiger}}.
#' Default \code{NULL} for using all cells.
#' @param legendColorTitle Title for colorbar legend. Default
#' \code{"Mean\\nExpression"}.
#' @param legendSizeTitle Title for size legend. Default
#' \code{"Percent\\nExpressed"}
#' @param viridisOption Name of available viridis palette. See
#' \code{\link[viridis]{viridis}}. Default \code{"magma"}.
#' @param verbose Logical. Whether to show progress information. Mainly when
#' subsetting data. Default \code{FALSE}.
#' @param ... Additional theme setting arguments passed to
#' \code{\link{.complexHeatmapDotPlot}} and heatmap setting arguments passed to
#' \code{\link[ComplexHeatmap]{Heatmap}}. See Details.
#' @return \code{\link[ComplexHeatmap]{HeatmapList}} object.
#' @export
#' @examples
#' # Use character vector of genes
#' features <- varFeatures(pbmcPlot)[1:10]
#' plotClusterGeneDot(pbmcPlot, features = features)
#'
#' # Use data.frame with grouping information, with more tweak on plot
#' features <- data.frame(features, rep(letters[1:5], 2))
#' plotClusterGeneDot(pbmcPlot, features = features,
#'                    clusterFeature = TRUE, clusterCell = TRUE, maxDotSize = 6)
plotClusterGeneDot <- function(
        object,
        features,
        groupBy = NULL,
        splitBy = NULL,
        featureScaleFunc = function(x) log2(10000*x + 1),
        cellIdx = NULL,
        legendColorTitle = "Mean\nExpression",
        legendSizeTitle = "Percent\nExpressed",
        viridisOption = "magma",
        verbose = FALSE,
        ...
) {
    groupBy <- groupBy %||% object@uns$defaultCluster
    allVars <- c(groupBy, splitBy)
    grouping <- .fetchCellMetaVar(object, variables = c(groupBy, splitBy),
                                  checkCategorical = TRUE, cellIdx = cellIdx,
                                  returnList = TRUE)

    # Get data that need to plot
    # Retrieved a matrix with features as columns, cells as rows
    if (is.data.frame(features)) {
        # tibble object passes is.data.frame, but does not work with many
        # trivial data.frame operation
        features <- as.data.frame(features)
        features <- features[,c(1,2)]
    } else if (is.character(features)) {
        features <- data.frame(feature = features)
    }
    # Retrieved expression matrix would always be cell x feature, designed for
    # ggplot
    mat <- retrieveCellFeature(object, feature = unique(features[,1]),
                               slot = "normData", cellIdx = cellIdx,
                               verbose = verbose)
    if (any(duplicated(features[,1]))) {
        mat <- mat[,features[,1]]
    }
    # In case specified features not found
    features <- features[features[,1] %in% colnames(mat), , drop = FALSE]
    allFeatures <- make.unique(features[,1])
    # Make sure everything consistent
    colnames(mat) <- allFeatures
    rownames(features) <- allFeatures
    # Calculate values as the aggregated representation,
    # format the final input for ComplexHeatmap

    if (!is.null(featureScaleFunc)) mat <- featureScaleFunc(mat)
    expMat <- stats::aggregate(mat, grouping, FUN = mean)
    # cellAnn can also be generated by coercing allVars to data.frame
    cellAnn <- expMat[, allVars, drop = FALSE]
    rownames(expMat) <- apply(cellAnn, 1, paste, collapse = ".")
    cellLabels <- apply(cellAnn[,!allVars %in% splitBy, drop = FALSE], 1,
                        paste, collapse = ".")
    expMat <- t(expMat[,-seq_along(allVars)])
    percMat <- stats::aggregate(mat, grouping, FUN = function(x) mean(x > 0))
    percMat <- t(percMat[,-seq_along(allVars)])

    cellSplitVar <- cellAnn[, splitBy, drop = FALSE]
    if (ncol(cellSplitVar) == 0) cellSplitVar <- NULL

    .complexHeatmapDotPlot(
        colorMat = expMat, sizeMat = percMat, featureAnnDF = features,
        legendColorTitle = legendColorTitle, legendSizeTitle = legendSizeTitle,
        cellSplitVar = cellSplitVar, cellLabels = cellLabels, ...
    )
}

#' Make dot plot of factor loading in cell groups
#' @description This function produces dot plots. Each column represent a group
#' of cells specified by \code{groupBy}, each row is a factor specified by
#' \code{useDims}. The color of dots reflects mean of factor loading of
#' specified factors in each cell group and sizes reflects the percentage of
#' cells that have loadings of a factor in a group. We utilize
#' \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html}{ComplexHeatmap}
#' for simplified management of adding annotation and slicing subplots. This was
#' inspired by the implementation in
#' \href{https://samuel-marsh.github.io/scCustomize/reference/Clustered_DotPlot.html}{scCustomize}.
#' @details For \code{...}, please notice that arguments \code{colorMat},
#' \code{sizeMat}, \code{featureAnnDF}, \code{cellSplitVar}, \code{cellLabels}
#' and \code{viridisOption} from \code{\link{.complexHeatmapDotPlot}} are
#' already occupied by this function internally. A lot of arguments from
#' \code{\link[ComplexHeatmap]{Heatmap}} have also been occupied: \code{matrix,
#' name, heatmap_legend_param, rect_gp, col, layer_fun, km, border, border_gp,
#' column_gap, row_gap, cluster_row_slices, cluster_rows, row_title_gp,
#' row_names_gp, row_split, row_labels, cluster_column_slices, cluster_columns,
#' column_split, column_title_gp, column_title, column_labels, column_names_gp,
#' top_annotation}.
#' @param object A \linkS4class{liger} object
#' @param groupBy The names of the columns in \code{cellMeta} slot storing
#' categorical variables. Loading data would be aggregated basing on these,
#' together with \code{splitBy}. Default uses default clusters.
#' @param useDims A Numeric vector to specify exact factors of interests.
#' Default \code{NULL} uses all available factors.
#' @param useRaw Whether to use un-aligned cell factor loadings (\eqn{H}
#' matrices). Default \code{FALSE}.
#' @param splitBy The names of the columns in \code{cellMeta} slot storing
#' categorical variables. Dotplot panel splitting would be based on these.
#' Default \code{NULL}.
#' @param factorScaleFunc A function object applied to factor loading matrix for
#' scaling the value for better visualization. Default \code{NULL}.
#' @param cellIdx Valid cell subscription. See \code{\link{subsetLiger}}.
#' Default \code{NULL} for using all cells.
#' @param legendColorTitle Title for colorbar legend. Default
#' \code{"Mean Factor\nLoading"}.
#' @param legendSizeTitle Title for size legend. Default
#' \code{"Percent\nLoaded"}
#' @param viridisOption Name of available viridis palette. See
#' \code{\link[viridis]{viridis}}. Default \code{"viridis"}.
#' @param verbose Logical. Whether to show progress information. Mainly when
#' subsetting data. Default \code{FALSE}.
#' @param ... Additional theme setting arguments passed to
#' \code{\link{.complexHeatmapDotPlot}} and heatmap setting arguments passed to
#' \code{\link[ComplexHeatmap]{Heatmap}}. See Details.
#' @return \code{\link[ComplexHeatmap]{HeatmapList}} object.
#' @export
#' @examples
#' plotClusterFactorDot(pbmcPlot)
plotClusterFactorDot <- function(
        object,
        groupBy = NULL,
        useDims = NULL,
        useRaw = FALSE,
        splitBy = NULL,
        factorScaleFunc = NULL,
        cellIdx = NULL,
        legendColorTitle = "Mean Factor\nLoading",
        legendSizeTitle = "Percent\nLoaded",
        viridisOption = "viridis",
        verbose = FALSE,
        ...
) {
    groupBy <- groupBy %||% object@uns$defaultCluster
    allVars <- c(groupBy, splitBy)
    grouping <- .fetchCellMetaVar(object, variables = c(groupBy, splitBy),
                                  checkCategorical = TRUE, cellIdx = cellIdx,
                                  returnList = TRUE)
    # Retrieved expression matrix would always be cell x feature, designed for
    # ggplot
    if (is.null(useDims)) useDims <- seq(object@uns$factorization$k)
    mat <- retrieveCellFeature(object, feature = useDims,
                               slot = ifelse(useRaw, "H", "H.norm"),
                               cellIdx = cellIdx, verbose = verbose)

    # Calculate values as the aggregated representation,
    # format the final input for ComplexHeatmap

    if (!is.null(factorScaleFunc)) mat <- factorScaleFunc(mat)
    expMat <- stats::aggregate(mat, grouping, FUN = mean)
    # cellAnn can also be generated by coercing allVars to data.frame
    cellAnn <- expMat[, allVars, drop = FALSE]
    rownames(expMat) <- apply(cellAnn, 1, paste, collapse = ".")
    cellLabels <- apply(cellAnn[,!allVars %in% splitBy, drop = FALSE], 1,
                        paste, collapse = ".")
    expMat <- t(expMat[,-seq_along(allVars)])
    percMat <- stats::aggregate(mat, grouping, FUN = function(x) mean(x > 0))
    percMat <- t(percMat[,-seq_along(allVars)])

    cellSplitVar <- cellAnn[, splitBy, drop = FALSE]
    if (ncol(cellSplitVar) == 0) cellSplitVar <- NULL

    .complexHeatmapDotPlot(
        colorMat = expMat, sizeMat = percMat,
        featureAnnDF = data.frame(colnames(mat)),
        legendColorTitle = legendColorTitle, legendSizeTitle = legendSizeTitle,
        cellSplitVar = cellSplitVar, cellLabels = cellLabels,
        viridisOption = viridisOption, ...
    )
}

#' Generate dot plot from input matrix with ComplexHeatmap
#' @param colorMat,sizeMat Matrix of the same size. Values in \code{colorMat}
#' will be visualized with color while values in \code{sizeMat} will be
#' reflected by dot size.
#' @param featureAnnDF Data frame of features containing feature names and
#' grouping labels.
#' @param cellSplitVar Split the cell orientation (default columns) by this
#' variable.
#' @param cellLabels Label to be shown on cell orientation.
#' @param maxDotSize The maximum dot size. Default \code{4}.
#' @param clusterFeature,clusterCell Whether the feature/cell orientation
#' (default rows/column, respectively) should be clustered. Default
#' \code{FALSE}.
#' @param legendColorTitle,legendSizeTitle The title for color bar and dot size
#' legends, repectively. Default see \code{"Matrix Value"} and \code{"Fraction
#' Value"}.
#' @param transpose Logical, whether to rotate the dot plot orientation. i.e.
#' rows as cell aggregation and columns as features. Default \code{FALSE}.
#' @param baseSize One-parameter control of all text sizes. Individual text
#' element sizes can be controlled by other size arguments. "Title" sizes are
#' 2 points larger than "text" sizes when being controlled by this. Default
#' \code{8}.
#' @param cellTextSize,featureTextSize,legendTextSize Size of cell labels,
#' feature label and legend text. Default \code{NULL} controls by
#' \code{baseSize}.
#' @param cellTitleSize,featureTitleSize,legendTitleSize Size of titles on
#' cell and feature orientation and legend title. Default \code{NULL} controls
#' by \code{baseSize + 2}.
#' @param featureGrpRot Number of degree to rotate the feature grouping label.
#' Default \code{0}.
#' @param viridisOption,viridisDirection See argument \code{option} and
#' \code{direction} of \code{\link[viridisLite]{viridis}}. Default \code{"A"}
#' and \code{-1}.
#' @param ... Additional arguments passed to
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#' @return A \code{\link[ComplexHeatmap]{HeatmapList}} object.
.complexHeatmapDotPlot <- function(
        colorMat,
        sizeMat,
        featureAnnDF = NULL,
        cellSplitVar = NULL,
        cellLabels = NULL,
        maxDotSize = 4,
        clusterFeature = FALSE,
        clusterCell = FALSE,
        legendColorTitle = "Matrix Value",
        legendSizeTitle = "Fraction Value",
        transpose = FALSE,
        baseSize = 8,
        cellTextSize = NULL,
        featureTextSize = NULL,
        cellTitleSize = NULL,
        featureTitleSize = NULL,
        legendTextSize = NULL,
        legendTitleSize = NULL,
        featureGrpRot = 0,
        viridisOption = "C",
        viridisDirection = -1,
        ...
) {
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

    if (isTRUE(transpose)) {
        colorMat <- t(colorMat)
        sizeMat <- t(sizeMat)
    }
    ## Customized color mapping requires a function object returned by
    ## colorRamp2
    col_fun <- circlize::colorRamp2(
        breaks = c(0, max(colorMat)/2, max(colorMat)),
        colors = viridis::viridis(n = 20, option = viridisOption,
                                  direction = viridisDirection)[c(1, 10, 20)]
    )

    ## Draw the sized colored dot, with original heatmap rectangles disabled.
    layer_fun = function(j, i, x, y, w, h, fill) {
        grid::grid.rect(x = x, y = y, width = w, height = h,
                        gp = grid::gpar(col = NA, fill = NA))
        grid::grid.circle(
            x = x, y = y,
            r = ComplexHeatmap::pindex(sizeMat, i, j) *
                grid::unit(maxDotSize, "pt"),
            gp = grid::gpar(
                fill = col_fun(ComplexHeatmap::pindex(colorMat, i, j)),
                col = NA
            )
        )
    }

    ## Hand made size legend
    sizeLgd <- list(ComplexHeatmap::Legend(
        labels = c(0, 0.25, 0.5, 0.75),
        labels_gp = grid::gpar(fontsize = legendText),
        title = legendSizeTitle,
        title_gp = grid::gpar(fontsize = legendTitle),
        graphics = list(
            function(x, y, w, h)
                grid::grid.circle(
                    x = x, y = y, r = 0 * grid::unit(maxDotSize, "pt"),
                    gp = grid::gpar(fill = "black")),
            function(x, y, w, h)
                grid::grid.circle(
                    x = x, y = y, r = 0.25 * grid::unit(maxDotSize, "pt"),
                    gp = grid::gpar(fill = "black")),
            function(x, y, w, h)
                grid::grid.circle(
                    x = x, y = y, r = 0.5 * grid::unit(maxDotSize, "pt"),
                    gp = grid::gpar(fill = "black")),
            function(x, y, w, h)
                grid::grid.circle(
                    x = x, y = y, r = 0.75 * grid::unit(maxDotSize, "pt"),
                    gp = grid::gpar(fill = "black"))
        )
    ))


    featureHA <- NULL
    if (isFALSE(colorMat)) featureLabels <- rownames(colorMat)
    else featureLabels <- colnames(colorMat)
    featureSplitVar <- NULL
    if (!is.null(featureAnnDF)) {
        if (is.data.frame(featureAnnDF) && ncol(featureAnnDF) > 1) {
            sliceLabel <- featureAnnDF[,2]
            if (!is.factor(sliceLabel)) sliceLabel <- factor(sliceLabel)
            blockAnn <- ComplexHeatmap::anno_block(
                gp = grid::gpar(fill = 0, col = "black"),
                labels = levels(sliceLabel),
                labels_gp = grid::gpar(fontsize = featureTitle),
                which = ifelse(transpose, "column", "row"), labels_rot = featureGrpRot
            )
            if (isFALSE(transpose))
                featureHA <- ComplexHeatmap::rowAnnotation(group = blockAnn)
            else if (isTRUE(transpose))
                featureHA <- ComplexHeatmap::columnAnnotation(group = blockAnn)
            featureLabels <- featureAnnDF[,1]
            featureSplitVar <- sliceLabel
        } else {
            # TODO whether to have better error messages on unknown classes?
            featureLabels <- featureAnnDF[,1]
        }
    }
    if (isFALSE(transpose)) {
        # If features are passed in a grouped structure, have the plot splitted by
        # specified grouping, and have group label added.
        hm <- ComplexHeatmap::Heatmap(
            # General settings
            matrix = colorMat,
            name = legendColorTitle,
            heatmap_legend_param = list(
                title_gp = grid::gpar(fontsize = legendTitle),
                labels_gp = grid::gpar(fontsize = legendText)
            ),
            ## This removes the heatmap small rectangles
            rect_gp = grid::gpar(type = "none"),
            col = col_fun,
            ## And then replace with the dots
            layer_fun = layer_fun,
            km = NULL,
            border = "grey",
            border_gp = grid::gpar(lwd = 0.2),
            # Column settings
            cluster_columns = clusterCell,
            cluster_column_slices = clusterCell,
            column_title_gp = grid::gpar(fontsize = cellTitle),
            column_names_gp = grid::gpar(fontsize = cellText),
            column_gap = grid::unit(0, "mm"),
            column_split = cellSplitVar,
            column_labels = cellLabels,
            # Row settings
            cluster_row_slices = clusterFeature,
            cluster_rows = clusterFeature,
            row_split = featureSplitVar,
            row_title_gp = grid::gpar(type = "none"),
            row_title = NULL,
            row_labels = featureLabels,
            row_names_gp = grid::gpar(fontsize = featureText),
            row_gap = grid::unit(0, "mm"),
            left_annotation = featureHA,
            ...
        )
    } else if (isTRUE(transpose)) {
        # It is pretty annoying that ComplexHeatmap does not provide an easy
        # method to flip things
        hm <- ComplexHeatmap::Heatmap(
            # General settings
            matrix = colorMat,
            name = legendColorTitle,
            heatmap_legend_param = list(
                title_gp = grid::gpar(fontsize = legendTitle),
                labels_gp = grid::gpar(fontsize = legendText)
            ),
            ## This removes the heatmap small rectangles
            rect_gp = grid::gpar(type = "none"),
            col = col_fun,
            ## And then replace with the dots
            layer_fun = layer_fun,
            km = NULL,
            border = "grey",
            border_gp = grid::gpar(lwd = 0.2),
            column_gap = grid::unit(0, "mm"),
            row_gap = grid::unit(0, "mm"),
            # row settings for cells
            cluster_row_slices = clusterCell,
            cluster_rows = clusterCell,
            row_title_gp = grid::gpar(fontsize = cellTitle),
            row_names_gp = grid::gpar(fontsize = cellText),
            row_split = cellSplitVar,
            row_labels = cellLabels,
            # Column settings for features
            cluster_column_slices = clusterFeature,
            cluster_columns = clusterFeature,
            column_split = featureSplitVar,
            column_title_gp = grid::gpar(type = "none"),
            column_title = NULL,
            column_labels = featureLabels,
            column_names_gp = grid::gpar(fontsize = featureText),
            top_annotation = featureHA,
            ...
        )
    }
    grDevices::pdf(nullfile())
    dp <- ComplexHeatmap::draw(hm, annotation_legend_list = sizeLgd,
                               merge_legend = TRUE)
    grDevices::dev.off()
    return(dp)
}

