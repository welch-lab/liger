#' @importFrom Matrix colSums rowSums t
#' @importFrom grDevices dev.off pdf
NULL

#' Plot t-SNE coordinates of cells across datasets
#'
#' Generates two plots of all cells across datasets, one colored by dataset and one colored by
#' cluster. These are useful for visually examining the alignment and cluster distributions,
#' respectively. If clusters have not been set yet (quantileAlignSNF not called), will plot by
#' single color for second plot. It is also possible to pass in another clustering (as long as
#' names match those of cells).
#'
#' @param object \code{liger} object. Should call runTSNE or runUMAP before calling.
#' @param clusters Another clustering to use for coloring second plot (must have same names as
#'   clusters slot) (default NULL).
#' @param title Plot titles (list or vector of length 2) (default NULL).
#' @param pt.size Controls size of points representing cells (default 0.3).
#' @param text.size Controls size of plot text (cluster center labels) (default 3).
#' @param do.shuffle Randomly shuffle points so that points from same dataset are not plotted
#'   one after the other (default TRUE).
#' @param rand.seed Random seed for reproducibility of point shuffling (default 1).
#' @param axis.labels Vector of two strings to use as x and y labels respectively.
#' @param do.legend Display legend on plots (default TRUE).
#' @param legend.size Size of legend on plots (default 5).
#' @param reorder.idents logical whether to reorder the datasets from default order before plotting (default FALSE).
#' @param new.order new dataset factor order for plotting.  must set reorder.idents = TRUE.
#' @param return.plots Return ggplot plot objects instead of printing directly (default FALSE).
#' @param legend.fonts.size Controls the font size of the legend.
#' @param raster Rasterization of points (default NULL). Automatically convert to raster format if
#'   there are over 100,000 cells to plot.
#'
#' @return List of ggplot plot objects (only if return.plots TRUE, otherwise prints plots to
#'   console).
#'
#' @importFrom ggplot2 ggplot geom_point geom_text ggtitle guides guide_legend aes theme xlab ylab
#' @importFrom dplyr %>% group_by summarize
#' @importFrom scattermore geom_scattermore
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' # get tsne.coords for normalized data
#' ligerex <- runTSNE(ligerex)
#' # plot to console
#' plotByDatasetAndCluster(ligerex)
#' # return list of plots
#' plots <- plotByDatasetAndCluster(ligerex, return.plots = TRUE)
#' }

plotByDatasetAndCluster <- function(object, clusters = NULL, title = NULL, pt.size = 0.3,
                                    text.size = 3, do.shuffle = TRUE, rand.seed = 1,
                                    axis.labels = NULL, do.legend = TRUE, legend.size = 5,
                                    reorder.idents = FALSE, new.order = NULL,
                                    return.plots = FALSE, legend.fonts.size = 12, raster = NULL) {
    # check raster and set by number of cells total if NULL
    if (is.null(x = raster)) {
        if (nrow(x = object@cell.data) > 1e5) {
            raster <- TRUE
            message("NOTE: Points are rasterized as number of cells/nuclei plotted exceeds 100,000.
              \n To plot in vector form set `raster = FALSE`.")
        } else {
            raster <- FALSE
        }
    }

    tsne_df <- data.frame(object@tsne.coords)
    colnames(tsne_df) <- c("Dim1", "Dim2")
    tsne_df[['Dataset']] <- unlist(lapply(1:length(object@H), function(x) {
        rep(names(object@H)[x], nrow(object@H[[x]]))
    }))
    if (reorder.idents == TRUE){
        tsne_df$Dataset <- factor(tsne_df$Dataset, levels = new.order)
    }
    c_names <- names(object@clusters)
    if (is.null(clusters)) {
        # if clusters have not been set yet
        if (length(object@clusters) == 0) {
            clusters <- rep(1, nrow(object@tsne.coords))
            names(clusters) <- c_names <- rownames(object@tsne.coords)
        } else {
            clusters <- object@clusters
            c_names <- names(object@clusters)
        }
    }
    tsne_df[['Cluster']] <- clusters[c_names]
    if (do.shuffle) {
        set.seed(rand.seed)
        idx <- sample(1:nrow(tsne_df))
        tsne_df <- tsne_df[idx, ]
    }


    if (isTRUE(x = raster)) {
        p1 <- ggplot(tsne_df, aes_string(x = 'Dim1', y = 'Dim2', color = 'Dataset')) + theme_bw() +
            theme_cowplot(legend.fonts.size) + geom_scattermore(pointsize = pt.size) +
            guides(color = guide_legend(override.aes = list(size = legend.size)))

        centers <- tsne_df %>% group_by(.data[['Cluster']]) %>% summarize(
            Dim1 = median(x = .data[['Dim1']]),
            Dim2 = median(x = .data[['Dim2']])
        )

        p2 <- ggplot(tsne_df, aes_string(x = 'Dim1', y = 'Dim2', color = 'Cluster')) +
            theme_cowplot(legend.fonts.size) + geom_scattermore(pointsize = pt.size) +
            geom_text(data = centers, mapping = aes_string(label = 'Cluster'), colour = "black", size = text.size) +
            guides(color = guide_legend(override.aes = list(size = legend.size)))
    } else {
        p1 <- ggplot(tsne_df, aes_string(x = 'Dim1', y = 'Dim2', color = 'Dataset')) + theme_bw() +
            theme_cowplot(legend.fonts.size) + geom_point(size = pt.size, stroke = 0.2) +
            guides(color = guide_legend(override.aes = list(size = legend.size)))

        centers <- tsne_df %>% group_by(.data[['Cluster']]) %>% summarize(
            Dim1 = median(x = .data[['Dim1']]),
            Dim2 = median(x = .data[['Dim2']])
        )

        p2 <- ggplot(tsne_df, aes_string(x = 'Dim1', y = 'Dim2', color = 'Cluster')) +
            theme_cowplot(legend.fonts.size) + geom_point(size = pt.size, stroke = 0.2) +
            geom_text(data = centers, mapping = aes_string(label = 'Cluster'), colour = "black", size = text.size) +
            guides(color = guide_legend(override.aes = list(size = legend.size)))
    }


    if (!is.null(title)) {
        p1 <- p1 + ggtitle(title[1])
        p2 <- p2 + ggtitle(title[2])
    }
    if (!is.null(axis.labels)) {
        p1 <- p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
        p2 <- p2 + xlab(axis.labels[1]) + ylab(axis.labels[2])
    }
    p1 <- p1 + theme_cowplot(12)
    p2 <- p2 + theme_cowplot(12)
    if (!do.legend) {
        p1 <- p1 + theme(legend.position = "none")
        p2 <- p2 + theme(legend.position = "none")
    }
    if (return.plots) {
        return(list(p1, p2))
    } else {
        print(p1)
        print(p2)
    }
}

#' Plot specific feature on t-SNE coordinates
#'
#' Generates one plot for each dataset, colored by chosen feature (column) from cell.data slot.
#' Feature can be categorical (factor) or continuous.
#' Can also plot all datasets combined with by.dataset = FALSE.
#'
#' @param object \code{liger} object. Should call runTSNE or runUMAP before calling.
#' @param feature Feature to plot (should be column from cell.data slot).
#' @param by.dataset Whether to generate separate plot for each dataset (default TRUE).
#' @param discrete Whether to treat feature as discrete; if left NULL will infer from column class
#'   in cell.data (if factor, treated like discrete) (default NULL).
#' @param title Plot title (default NULL).
#' @param pt.size Controls size of points representing cells (default 0.3).
#' @param text.size Controls size of plot text (cluster center labels) (default 3).
#' @param do.shuffle Randomly shuffle points so that points from same dataset are not plotted
#'   one after the other (default TRUE).
#' @param rand.seed Random seed for reproducibility of point shuffling (default 1).
#' @param do.labels Print centroid labels for categorical features (default FALSE).
#' @param axis.labels Vector of two strings to use as x and y labels respectively.
#' @param do.legend Display legend on plots (default TRUE).
#' @param legend.size Size of legend spots for discrete data (default 5).
#' @param option Colormap option to use for ggplot2's scale_color_viridis (default 'plasma').
#' @param cols.use Vector of colors to form gradient over instead of viridis colormap (low to high).
#'   Only applies to continuous features (default NULL).
#' @param zero.color Color to use for zero values (no expression) (default '#F5F5F5').
#' @param return.plots Return ggplot plot objects instead of printing directly (default FALSE).
#'
#' @return List of ggplot plot objects (only if return.plots TRUE, otherwise prints plots to
#'   console).
#'
#' @importFrom ggplot2 ggplot geom_point geom_text ggtitle aes guides guide_legend labs
#' scale_color_viridis_c scale_color_gradientn theme xlab ylab
#' @importFrom dplyr %>% group_by summarize
#' @importFrom stats median
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' # get tsne.coords for normalized data
#' ligerex <- runTSNE(ligerex)
#' # plot nUMI to console
#' plotFeature(ligerex, feature = 'nUMI')
#' }

plotFeature <- function(object, feature, by.dataset = TRUE, discrete = NULL, title = NULL,
                        pt.size = 0.3, text.size = 3, do.shuffle = TRUE, rand.seed = 1, do.labels = FALSE,
                        axis.labels = NULL, do.legend = TRUE, legend.size = 5, option = 'plasma',
                        cols.use = NULL, zero.color = '#F5F5F5', return.plots = FALSE) {
    dr_df <- data.frame(object@tsne.coords)
    colnames(dr_df) <- c("dr1", "dr2")
    if (!(feature %in% colnames(object@cell.data))) {
        stop('Please select existing feature in cell.data, or add it before calling.')
    }
    dr_df$feature <- object@cell.data[, feature]
    if (is.null(discrete)) {
        if (class(dr_df$feature) != "factor") {
            discrete <- FALSE
        } else {
            discrete <- TRUE
        }
    }
    if (!discrete){
        dr_df$feature[dr_df$feature == 0] <- NA
    }
    if (by.dataset) {
        dr_df$dataset <- object@cell.data$dataset
    } else {
        dr_df$dataset <- factor("single")
    }
    if (do.shuffle) {
        set.seed(rand.seed)
        idx <- sample(1:nrow(dr_df))
        dr_df <- dr_df[idx, ]
    }
    p_list <- list()
    for (sub_df in split(dr_df, f = dr_df$dataset)) {
        ggp <- ggplot(sub_df, aes_string(x = 'dr1', y = 'dr2', color = 'feature')) + geom_point(size = pt.size)

        # if data is discrete
        if (discrete) {
            ggp <- ggp + guides(color = guide_legend(override.aes = list(size = legend.size))) +
                labs(col = feature)
            if (do.labels) {
                centers <- sub_df %>% group_by(feature) %>% summarize(
                    dr1 = median(x = sub_df[['dr1']]),
                    dr2 = median(x = sub_df[['dr2']])
                )
                ggp <- ggp + geom_text(data = centers, mapping = aes(label = feature),
                                       colour = "black", size = text.size)
            }
        } else {
            if (is.null(cols.use)) {
                ggp <- ggp + scale_color_viridis_c(option = option,
                                                   direction = -1,
                                                   na.value = zero.color) + labs(col = feature)
            } else {
                ggp <- ggp + scale_color_gradientn(colors = cols.use,
                                                   na.value = zero.color) + labs(col = feature)
            }

        }
        if (by.dataset) {
            base <- as.character(sub_df$dataset[1])
        } else {
            base <- ""
        }
        if (!is.null(title)) {
            base <- paste(title, base)
        }
        ggp <- ggp + ggtitle(base)
        if (!is.null(axis.labels)) {
            ggp <- ggp + xlab(axis.labels[1]) + ylab(axis.labels[2])
        }
        if (!do.legend) {
            ggp <- ggp + theme(legend.position = "none")
        }
        p_list[[as.character(sub_df$dataset[1])]] <- ggp
    }
    if (by.dataset) {
        p_list <- p_list[names(object@raw.data)]
    }

    if (return.plots){
        if (length(p_list) == 1) {
            return(p_list[[1]])
        } else {
            return(p_list)
        }
    } else {
        for (plot in p_list) {
            print(plot)
        }
    }
}

#' Plot scatter plots of unaligned and aligned factor loadings
#'
#' @description
#' Generates scatter plots of factor loadings vs cells for both unaligned and aligned
#' (normalized) factor loadings. This allows for easier visualization of the changes made to the
#' factor loadings during the alignment step. Lists a subset of highly loading genes for each factor.
#' Also provides an option to plot t-SNE coordinates of the cells colored by aligned factor loadings.
#'
#' It is recommended to call this function into a PDF due to the large number of
#' plots produced.
#'
#' @param object \code{liger} object. Should call quantileAlignSNF before calling.
#' @param num.genes Number of genes to display for each factor (default 10).
#' @param cells.highlight Names of specific cells to highlight in plot (black) (default NULL).
#' @param plot.tsne Plot t-SNE coordinates for each factor (default FALSE).
#' @param verbose Print messages (TRUE by default)
#'
#' @return Plots to console (1-2 pages per factor)
#'
#' @importFrom graphics legend par plot
#' @importFrom grDevices rainbow
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' # get tsne.coords for normalized data
#' ligerex <- runTSNE(ligerex)
#' # factor plots into pdf file
#' # pdf("plot_factors.pdf")
#' plotFactors(ligerex)
#' # dev.off()
#' }

plotFactors <- function(object, num.genes = 10, cells.highlight = NULL, plot.tsne = FALSE, verbose = TRUE) {
    k <- ncol(object@H.norm)
    if (verbose) {
        pb <- txtProgressBar(min = 0, max = k, style = 3)
    }
    W <- t(object@W)
    rownames(W) <- colnames(object@H[[1]])
    Hs_norm <- object@H.norm
    # restore default settings when the current function exits
    init_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(init_par))
    for (i in 1:k) {
        graphics::par(mfrow = c(2, 1))
        top_genes.W <- rownames(W)[order(W[, i], decreasing = TRUE)[1:num.genes]]
        top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
        factor_textstring <- paste0("Factor", i)

        plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
        cols <- rep("gray", times = nrow(Hs_norm))
        names(cols) <- rownames(Hs_norm)
        cols.use <- grDevices::rainbow(length(object@H))

        for (cl in 1:length(object@H)) {
            cols[rownames(object@H[[cl]])] <- rep(cols.use[cl], times = nrow(object@H[[cl]]))
        }
        if (!is.null(cells.highlight)) {
            cols[cells.highlight] <- rep("black", times = length(cells.highlight))
        }
        graphics::plot(1:nrow(Hs_norm), do.call(rbind, object@H)[, i],
                       cex = 0.2, pch = 20,
                       col = cols, main = plot_title1, xlab = "Cell", ylab = "Raw H Score"
        )
        graphics::legend("top", names(object@H), pch = 20, col = cols.use, horiz = TRUE, cex = 0.75)
        graphics::plot(1:nrow(Hs_norm), object@H.norm[, i],
                       pch = 20, cex = 0.2,
                       col = cols, xlab = "Cell", ylab = "H_norm Score"
        )
        if (plot.tsne) {
            graphics::par(mfrow = c(1, 1))
            fplot(object@tsne.coords, object@H.norm[, i], title = paste0("Factor ", i))
        }
        if (verbose) {
            setTxtProgressBar(pb, i)
        }
    }
}

#' Generate word clouds and t-SNE plots
#'
#' @description
#' Plots t-SNE coordinates of all cells by their loadings on each factor. Underneath it displays the
#' most highly loading shared and dataset-specific genes, with the size of the marker indicating
#' the magnitude of the loading.
#'
#' It is recommended to call this function into a PDF due to the large number of
#' plots produced.
#'
#' @param object \code{liger} object. Should call runTSNE before calling.
#' @param dataset1 Name of first dataset (by default takes first two datasets for dataset1 and 2)
#' @param dataset2 Name of second dataset
#' @param num.genes Number of genes to show in word clouds (default 30).
#' @param min.size Size of smallest gene symbol in word cloud (default 1).
#' @param max.size Size of largest gene symbol in word cloud (default 4).
#' @param factor.share.thresh Use only factors with a dataset specificity less than or equalt to
#'   threshold (default 10).
#' @param log.fc.thresh Lower log-fold change threshold for differential expression in markers
#'   (default 1).
#' @param pval.thresh Upper p-value threshold for Wilcoxon rank test for gene expression
#'   (default 0.05).
#' @param do.spec.plot Include dataset specificity plot in printout (default TRUE).
#' @param return.plots Return ggplot objects instead of printing directly (default FALSE).
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return List of ggplot plot objects (only if return.plots TRUE, otherwise prints plots to
#'   console).
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes aes_string geom_point ggtitle scale_color_gradient scale_size
#' scale_x_continuous scale_y_continuous coord_fixed labs
#' @importFrom grid roundrectGrob
#' @importFrom grid gpar
#' @importFrom cowplot draw_grob
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object based on in-memory datasets), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' ligerex <- runTSNE(ligerex)
#' # pdf('word_clouds.pdf')
#' plotWordClouds(ligerex, num.genes = 20)
#' # dev.off()
#' # ligerex (liger object based on datasets in HDF5 format), factorization complete input
#' ligerex <- readSubset(ligerex, slot.use = "norm.data", max.cells = 5000)
#' plotWordClouds(ligerex, num.genes = 20)
#' }

plotWordClouds <- function(object, dataset1 = NULL, dataset2 = NULL, num.genes = 30, min.size = 1,
                           max.size = 4, factor.share.thresh = 10, log.fc.thresh = 1, pval.thresh = 0.05,
                           do.spec.plot = TRUE, return.plots = FALSE, verbose = TRUE) {
    if (is.null(dataset1) | is.null(dataset2)) {
        dataset1 <- names(object@H)[1]
        dataset2 <- names(object@H)[2]
    }

    if(class(object@raw.data[[1]])[1] == "H5File"){
        sample.idx = unlist(lapply(object@sample.data, colnames))
        H_aligned = object@H.norm[sample.idx, ]
        tsne_coords <- object@tsne.coords[sample.idx, ]
    } else {
        H_aligned <- object@H.norm
        tsne_coords <- object@tsne.coords
    }

    W <- t(object@W)
    V1 <- t(object@V[[dataset1]])
    V2 <- t(object@V[[dataset2]])
    W <- pmin(W + V1, W + V2)

    dataset.specificity <- calcDatasetSpecificity(object, dataset1 = dataset1,
                                                  dataset2 = dataset2, do.plot = do.spec.plot)
    factors.use <- which(abs(dataset.specificity[[3]]) <= factor.share.thresh)

    markers <- getFactorMarkers(object, dataset1 = dataset1, dataset2 = dataset2,
                                factor.share.thresh = factor.share.thresh,
                                num.genes = num.genes, log.fc.thresh = log.fc.thresh,
                                pval.thresh = pval.thresh,
                                dataset.specificity = dataset.specificity,
                                verbose = verbose
    )

    rownames(W) <- rownames(V1) <- rownames(V2) <- object@var.genes
    loadings_list <- list(V1, W, V2)
    names_list <- list(dataset1, "Shared", dataset2)
    if (verbose) {
        pb <- txtProgressBar(min = 0, max = length(factors.use), style = 3)
    }
    return_plots <- list()
    for (i in factors.use) {
        tsne_df <- data.frame(H_aligned[, i], tsne_coords)
        factorlab <- paste("Factor", i, sep = "")
        colnames(tsne_df) <- c(factorlab, "Dim1", "Dim2")
        factor_ds <- paste("Factor", i, "Dataset Specificity:", dataset.specificity[[3]][i])
        p1 <- ggplot(tsne_df, aes_string(x = "Dim1", y = "Dim2", color = factorlab)) + geom_point() +
            scale_color_gradient(low = "yellow", high = "red") + ggtitle(label = factor_ds)

        top_genes_V1 <- markers[[1]]$gene[markers[[1]]$factor_num == i]
        top_genes_W <- markers[[2]]$gene[markers[[2]]$factor_num == i]
        top_genes_V2 <- markers[[3]]$gene[markers[[3]]$factor_num == i]

        top_genes_list <- list(top_genes_V1, top_genes_W, top_genes_V2)
        plot_list <- lapply(seq_along(top_genes_list), function(x) {
            top_genes <- top_genes_list[[x]]
            gene_df <- data.frame(
                genes = top_genes,
                loadings = loadings_list[[x]][top_genes, i]
            )
            if (length(top_genes) == 0) {
                gene_df <- data.frame(genes = c("no genes"), loadings = c(1))
            }
            out_plot <- ggplot(gene_df, aes(x = 1, y = 1, size = loadings, label = .data[['genes']])) +
                geom_text_repel(force = 100, segment.color = NA) +
                scale_size(range = c(min.size, max.size), guide = FALSE) +
                scale_y_continuous(breaks = NULL) +
                scale_x_continuous(breaks = NULL) +
                labs(x = "", y = "") + ggtitle(label = names_list[[x]]) + coord_fixed() + ggplot2::theme_void()
            return(out_plot)
        })

        p2 <- (plot_grid(plotlist = plot_list, align = "hv", nrow = 1)
               + draw_grob(roundrectGrob(
                   x = 0.33, y = 0.5, width = 0.67, height = 0.70,
                   gp = gpar(fill = "khaki1", col = "Black", alpha = 0.5, lwd = 2)
               ))
               + draw_grob(roundrectGrob(
                   x = 0.67, y = 0.5, width = 0.67, height = 0.70,
                   gp = gpar(fill = "indianred1", col = "Black", alpha = 0.5, lwd = 2)
               )))
        return_plots[[i]] <- plot_grid(p1, p2, nrow = 2, align = "h")
        if (!return.plots) {
            print(return_plots[[i]])
        }
        if (verbose) {
            setTxtProgressBar(pb, i)
        }
    }
    if (return.plots) {
        return(return_plots)
    }
}

#' Generate t-SNE plots and gene loading plots
#'
#' @description
#' Plots t-SNE coordinates of all cells by their loadings on each factor. Underneath it displays the
#' most highly loading shared and dataset-specific genes, along with the overall gene loadings
#' for each dataset.
#'
#' It is recommended to call this function into a PDF due to the large number of
#' plots produced.
#'
#' @param object \code{liger} object. Should call runTSNE before calling.
#' @param dataset1 Name of first dataset (by default takes first two datasets for dataset1 and 2)
#' @param dataset2 Name of second dataset
#' @param num.genes Number of genes to show in word clouds (default 30).
#' @param num.genes.show Number of genes displayed as y-axis labels in the gene loading plots at
#' the bottom (default 12)
#' @param mark.top.genes Plot points corresponding to top loading genes in different color (default
#'   TRUE).
#' @param factor.share.thresh Use only factors with a dataset specificity less than or equal to
#'   threshold (default 10).
#' @param log.fc.thresh Lower log-fold change threshold for differential expression in markers
#'   (default 1).
#' @param umi.thresh Lower UMI threshold for markers (default 30).
#' @param frac.thresh Lower threshold for fraction of cells expressing marker (default 0).
#' @param pval.thresh Upper p-value threshold for Wilcoxon rank test for gene expression
#'   (default 0.05).
#' @param do.spec.plot Include dataset specificity plot in printout (default TRUE).
#' @param max.val Value between 0 and 1 at which color gradient should saturate to max color. Set to
#'   NULL to revert to default gradient scaling. (default 0.1)
#' @param pt.size Point size for plots (default 0.4).
#' @inheritParams plotGene
#' @param return.plots Return ggplot objects instead of printing directly (default FALSE).
#' @param axis.labels Vector of two strings to use as x and y labels respectively (default NULL).
#' @param do.title Include top title with cluster and Dataset Specificity (default FALSE).
#' @param verbose Print progress bar/messages (TRUE by default)
#' @param raster Rasterization of points (default NULL). Automatically convert to raster format if
#'   there are over 100,000 cells to plot.
#'
#' @return List of ggplot plot objects (only if return.plots TRUE, otherwise prints plots to
#'   console).
#'
#' @importFrom ggplot2 aes aes_string annotate coord_cartesian element_blank ggplot geom_point
#' ggtitle scale_color_viridis_c theme
#' theme_bw
#' @importFrom grid gpar unit
#' @import patchwork
#' @importFrom stats loadings
#' @importFrom cowplot theme_cowplot
#' @importFrom scattermore geom_scattermore
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object based on in-memory datasets), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' ligerex <- runUMAP(ligerex)
#' # pdf("gene_loadings.pdf")
#' plotGeneLoadings(ligerex, num.genes = 20)
#' # dev.off()
#' # ligerex (liger object based on datasets in HDF5 format), factorization complete input
#' ligerex <- readSubset(ligerex, slot.use = "norm.data", max.cells = 5000)
#' plotGeneLoadings(ligerex, num.genes = 20)
#' }
#'
plotGeneLoadings <- function(object, dataset1 = NULL, dataset2 = NULL, num.genes.show = 12,
                             num.genes = 30, mark.top.genes = TRUE, factor.share.thresh = 10,
                             log.fc.thresh = 1, umi.thresh = 30, frac.thresh = 0,
                             pval.thresh = 0.05, do.spec.plot = TRUE, max.val = 0.1, pt.size = 0.4,
                             option = "plasma", zero.color = "#F5F5F5", return.plots = FALSE,
                             axis.labels = NULL, do.title = FALSE, verbose = TRUE, raster = NULL) {
    # check raster and set by number of cells total if NULL
    if (is.null(x = raster)) {
        if (nrow(x = object@cell.data) > 1e5) {
            raster <- TRUE
            message("NOTE: Points are rasterized as number of cells/nuclei plotted exceeds 100,000.
              \n To plot in vector form set `raster = FALSE`.")
        } else {
            raster <- FALSE
        }
    }

    if (is.null(dataset1) | is.null(dataset2)) {
        dataset1 <- names(object@H)[1]
        dataset2 <- names(object@H)[2]
    }

    if(class(object@raw.data[[1]])[1] == "H5File"){
        sample.idx = unlist(lapply(object@sample.data, colnames))
        H_aligned = object@H.norm[sample.idx, ]
        tsne_coords <- object@tsne.coords[sample.idx, ]
    } else {
        H_aligned <- object@H.norm
        tsne_coords <- object@tsne.coords
    }

    W_orig <- t(object@W)
    V1 <- t(object@V[[dataset1]])
    V2 <- t(object@V[[dataset2]])
    W <- pmin(W_orig + V1, W_orig + V2)

    dataset.specificity <- calcDatasetSpecificity(object,
                                                  dataset1 = dataset1,
                                                  dataset2 = dataset2, do.plot = do.spec.plot
    )

    factors.use <- which(abs(dataset.specificity[[3]]) <= factor.share.thresh)


    markers <- getFactorMarkers(object,
                                dataset1 = dataset1, dataset2 = dataset2,
                                factor.share.thresh = factor.share.thresh,
                                num.genes = num.genes, log.fc.thresh = log.fc.thresh,
                                pval.thresh = pval.thresh,
                                dataset.specificity = dataset.specificity,
                                verbose = verbose
    )

    rownames(W) <- rownames(V1) <- rownames(V2) <- rownames(W_orig) <- object@var.genes
    loadings_list <- list(V1, W, V2)
    names_list <- list(dataset1, "Shared", dataset2)
    if (verbose) {
        pb <- txtProgressBar(min = 0, max = length(factors.use), style = 3)
    }
    return_plots <- list()
    for (i in factors.use) {
        tsne_df <- data.frame(H_aligned[, i], tsne_coords)
        factorlab <- paste("Factor", i, sep = "")
        colnames(tsne_df) <- c(factorlab, "Dim1", "Dim2")
        tsne_df[[factorlab]][tsne_df[[factorlab]] == 0] <- NA
        factor_ds <- paste("Factor", i, "Dataset Specificity:", dataset.specificity[[3]][i])
        data.max <- max(object@H.norm[, i])
        # plot TSNE
        if (!is.null(max.val)) {
            values <- c(0, max.val, 1)
        } else {
            values <- NULL
        }

        if (isTRUE(x = raster)) {
            p1 <- ggplot(tsne_df, aes_string(x = "Dim1", y = "Dim2", color = factorlab)) +
                geom_scattermore(pointsize = pt.size) +
                scale_color_viridis_c(
                    option = option,
                    direction = -1,
                    na.value = zero.color, values = values
                ) +
                theme_cowplot(12)
        } else {
            p1 <- ggplot(tsne_df, aes_string(x = "Dim1", y = "Dim2", color = factorlab)) +
                geom_point(size = pt.size) +
                scale_color_viridis_c(
                    option = option,
                    direction = -1,
                    na.value = zero.color, values = values
                ) +
                theme_cowplot(12)
        }


        if (!is.null(axis.labels)) {
            p1 <- p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
        }
        if (do.title) {
            p1 <- p1 + ggtitle(label = factor_ds)
        }

        # subset to specific factor and sort by p-value
        top_genes_V1 <- markers[[1]][markers[[1]]$factor_num == i, ]
        top_genes_V1 <- top_genes_V1[order(top_genes_V1$p_value), ]$gene
        # don't sort for W
        top_genes_W <- markers[[2]][markers[[2]]$factor_num == i, ]$gene
        top_genes_V2 <- markers[[3]][markers[[3]]$factor_num == i, ]
        top_genes_V2 <- top_genes_V2[order(top_genes_V2$p_value), ]$gene

        top_genes_list <- list(top_genes_V1, top_genes_W, top_genes_V2)
        # subset down to those which will be shown if sorting by p-val

        top_genes_list <- lapply(top_genes_list, function(x) {
            if (length(x) > num.genes.show) {
                # to avoid subset warning
                x <- x[1:num.genes.show]
            }
            x
        })

        plot_list <- lapply(seq_along(top_genes_list), function(x) {
            top_genes <- top_genes_list[[x]]
            # make dataframe for cum gene loadings plot
            sorted <- sort(loadings_list[[x]][, i])
            # sort by loadings instead - still only showing num.genes.show
            # look through top num.genes in loadings
            top_loaded <- names(rev(sorted[(length(sorted) - num.genes + 1):length(sorted)]))
            top_genes <- top_loaded[which(top_loaded %in% top_genes)]
            if (length(top_genes) == 0) {
                top_genes <- c("no genes")
            }

            gene_df <- data.frame(
                loadings = sorted,
                xpos = seq(0, 1, length.out = length(sorted)),
                top_k = names(sorted) %in% top_genes
            )
            y_lim_text <- max(gene_df$loadings)
            # plot and annotate with top genes

            out_plot <- ggplot(gene_df, aes_string(x = 'xpos', y = 'loadings')) +
                geom_point(size = pt.size) +
                theme_bw() +
                theme(
                    axis.ticks.x = element_blank(),
                    axis.line.x = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank()
                ) +
                ggtitle(label = names_list[[x]]) +
                annotate("text",
                         x = 1.1,
                         y = seq(y_lim_text, 0, length.out = num.genes.show)[1:length(top_genes)],
                         label = top_genes, hjust = 0, col = "#8227A0"
                ) +
                coord_cartesian(
                    xlim = c(0, 1), # This focuses the x-axis on the range of interest
                    clip = "off"
                ) +
                theme(plot.margin = unit(c(1, 4, 1, 1), "lines"))

            if (mark.top.genes) {
                out_plot <- out_plot + geom_point(
                    data = subset(gene_df, gene_df[['top_k']] == TRUE),
                    aes_string('xpos', 'loadings'),
                    col = "#8227A0", size = 0.5
                )
            }
            return(out_plot)
        })

        # p2 <- plot_grid(plotlist = plot_list, nrow = 1)

        return_plots[[i]] <- p1 / (plot_list[[1]] | plot_list[[2]] | plot_list[[3]])
        # if can figure out how to make cowplot work, might bring this back
        # return_plots[[i]] <- plot_grid(p1, p2, nrow = 2, align = "h")
        if (!return.plots) {
            print(return_plots[[i]])
        }
        if (verbose) {
            setTxtProgressBar(pb, i)
        }
    }
    if (return.plots) {
        return(return_plots)
    }
}

#' Plot violin plots for gene expression
#'
#' Generates violin plots of expression of specified gene for each dataset.
#'
#' @param object \code{liger} object.
#' @param gene Gene for which to plot relative expression.
#' @param methylation.indices Indices of datasets in object with methylation data (this data is not
#'   magnified and put on log scale).
#' @param by.dataset Plots gene expression for each dataset separately (default TRUE).
#' @param return.plots Return ggplot objects instead of printing directly to console (default
#'   FALSE).
#'
#' @return List of ggplot plot objects (only if return.plots TRUE, otherwise prints plots to
#'   console).
#'
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string ggplot geom_point geom_boxplot geom_violin ggtitle labs
#' scale_color_gradient2 theme
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object based on in-memory datasets), factorization complete
#' # plot expression for CD4 and return plots
#' violin_plots <- plotGeneViolin(ligerex, "CD4", return.plots = TRUE)
#' # ligerex (liger object based on datasets in HDF5 format), factorization complete input
#' ligerex <- readSubset(ligerex, slot.use = "norm.data", max.cells = 5000)
#' violin_plots <- plotGeneViolin(ligerex, "CD4", return.plots = TRUE)
#' }

plotGeneViolin <- function(object, gene, methylation.indices = NULL,
                           by.dataset = TRUE, return.plots = FALSE) {
    if (class(object@raw.data[[1]])[1] == "H5File"){
        if (object@h5file.info[[1]][["sample.data.type"]] != "norm.data"){
            stop("norm.data should be sampled for making violin plots.")
        }
    }

    gene_vals <- c()
    gene_df <- data.frame(object@tsne.coords)
    rownames(gene_df) <- names(object@clusters)

    for (i in 1:length(object@raw.data)) {
        if (class(object@raw.data[[i]])[1] == "H5File"){
            if (i %in% methylation.indices) {
                gene_vals <- c(gene_vals, object@sample.data[[i]][gene, ])
            } else {
                if (gene %in% rownames(object@sample.data[[i]])) {
                    gene_vals_int <- log2(10000 * object@sample.data[[i]][gene, ] + 1)
                }
                else {
                    gene_vals_int <- rep(list(0), ncol(object@sample.data[[i]]))
                    names(gene_vals_int) <- colnames(object@sample.data[[i]])
                }
                gene_vals <- c(gene_vals, gene_vals_int)
            }
        } else {
            if (i %in% methylation.indices) {
                gene_vals <- c(gene_vals, object@norm.data[[i]][gene, ])
            } else {
                if (gene %in% rownames(object@norm.data[[i]])) {
                    gene_vals_int <- log2(10000 * object@norm.data[[i]][gene, ] + 1)
                }
                else {
                    gene_vals_int <- rep(list(0), ncol(object@norm.data[[i]]))
                    names(gene_vals_int) <- colnames(object@norm.data[[i]])
                }
                gene_vals <- c(gene_vals, gene_vals_int)
            }
        }
    }

    gene_df$Gene <- as.numeric(gene_vals[rownames(gene_df)])
    colnames(gene_df) <- c("Dim1", "Dim2", "gene")
    gene_plots <- list()
    for (i in 1:length(object@scale.data)) {
        if (by.dataset) {
            gene_df.sub <- gene_df[rownames(object@H[[i]]), ]
            gene_df.sub$Cluster <- object@clusters[rownames(object@H[[i]])]
            title <- names(object@scale.data)[i]
        } else {
            gene_df.sub <- gene_df
            gene_df.sub$Cluster <- object@clusters
            title <- "All Datasets"
        }
        max_v <- max(gene_df.sub["gene"], na.rm = TRUE)
        min_v <- min(gene_df.sub["gene"], na.rm = TRUE)
        midpoint <- (max_v - min_v) / 2
        plot_i <- ggplot(gene_df.sub, aes_string(x = "Cluster", y = "gene", fill = "Cluster")) +
            geom_boxplot(position = "dodge", width = 0.4, outlier.shape = NA, alpha = 0.7) +
            geom_violin(position = "dodge", alpha = 0.7) +
            ggtitle(title)
        gene_plots[[i]] <- plot_i + theme(legend.position = "none") + labs(y = gene)
        if (i == 1 & !by.dataset) {
            break
        }
    }
    if (return.plots) {
        return(gene_plots)
    } else {
        for (i in 1:length(gene_plots)) {
            print(gene_plots[[i]])
        }
    }
}

#' Plot gene expression on dimensional reduction (t-SNE) coordinates
#'
#' Generates plot of dimensional reduction coordinates (default t-SNE) colored by expression of
#' specified gene. Data can be scaled by dataset or selected feature column from cell.data (or across
#' all cells). Data plots can be split by feature.
#'
#' @param object \code{liger} object. Should call runTSNE before calling.
#' @param gene Gene for which to plot expression.
#' @param use.raw Plot raw UMI values instead of normalized, log-transformed data (default FALSE).
#' @param use.scaled Plot values scaled across specified groups of cells (with log transformation)
#'   (default FALSE).
#' @param scale.by Grouping of cells by which to scale gene (can be any factor column in cell.data
#'   or 'none' for scaling across all cells) (default 'dataset').
#' @param log2scale Whether to show log2 transformed values or original normalized, raw, or scaled
#'   values (as stored in object). Default value is FALSE if use.raw = TRUE, otherwise TRUE.
#' @param methylation.indices Indices of datasets in object with methylation data (this data is not
#'   log transformed and must use normalized values). (default NULL)
#' @param plot.by How to group cells for plotting (can be any factor column in cell.data or 'none'
#'   for plotting all cells in a single plot). Note that this can result in large number of plots.
#'   Users are encouraged to use same value as for scale.by (default 'dataset').
#' @param set.dr.lims Whether to keep dimensional reduction coordinates consistent when multiple
#'   plots created (default FALSE).
#' @param pt.size Point size for plots (default 0.1).
#' @param min.clip Minimum value for expression values plotted. Can pass in quantile (0-1) or
#'   absolute cutoff (set clip.absolute = TRUE). Can also pass in vector if expecting multiple plots;
#'   users are encouraged to pass in named vector (from levels of desired feature) to avoid
#'   mismatches in order (default NULL).
#' @param max.clip Maximum value for expression values plotted. Can pass in quantile (0-1) or
#'   absolute cutoff (set clip.absolute = TRUE). Can also pass in vector if expecting multiple plots;
#'   users are encouraged to pass in named vector (from levels of desired feature) to avoid
#'   mismatches in order (default NULL).
#' @param clip.absolute Whether to treat clip values as absolute cutoffs instead of quantiles
#'   (default FALSE).
#' @param points.only Remove axes, background, and legend when plotting coordinates (default FALSE).
#' @param option Colormap option to use for ggplot2's scale_color_viridis (default 'plasma').
#' @param cols.use Vector of colors to form gradient over instead of viridis colormap (low to high).
#'   (default NULL).
#' @param zero.color Color to use for zero values (no expression) (default '#F5F5F5').
#' @param axis.labels Vector of two strings to use as x and y labels respectively. (default NULL)
#' @param do.legend Display legend on plots (default TRUE).
#' @param return.plots Return ggplot objects instead of printing directly (default FALSE).
#' @param keep.scale Maintain min/max color scale across all plots when using plot.by (default FALSE)
#' @param raster Rasterization of points (default NULL). Automatically convert to raster format if
#'   there are over 100,000 cells to plot.
#'
#' @return If returning single plot, returns ggplot object; if returning multiple plots; returns
#'   list of ggplot objects.
#'
#' @importFrom dplyr %>% group_by mutate_at vars group_cols
#' @importFrom ggplot2 ggplot geom_point aes_string element_blank ggtitle labs xlim ylim
#' scale_color_viridis_c scale_color_gradientn theme
#' @importFrom stats quantile
#' @importFrom scattermore geom_scattermore
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object based on in-memory datasets), factorization complete
#' ligerex
#' ligerex <- runTSNE(ligerex)
#' # plot expression for CD4 and return plots
#' gene_plots <- plotGene(ligerex, "CD4", return.plots = TRUE)
#' # ligerex (liger object based on datasets in HDF5 format), factorization complete input
#' ligerex <- readSubset(ligerex, slot.use = "norm.data", max.cells = 5000)
#' gene_plots <- plotGene(ligerex, "CD4", return.plots = TRUE)
#' }

plotGene <- function(object, gene, use.raw = FALSE, use.scaled = FALSE, scale.by = 'dataset',
                     log2scale = NULL, methylation.indices = NULL, plot.by = 'dataset',
                     set.dr.lims = FALSE, pt.size = 0.1, min.clip = NULL, max.clip = NULL,
                     clip.absolute = FALSE, points.only = FALSE, option = 'plasma', cols.use = NULL,
                     zero.color = '#F5F5F5', axis.labels = NULL, do.legend = TRUE, return.plots = FALSE,
                     keep.scale = FALSE, raster = NULL) {
    if ((plot.by != scale.by) & (use.scaled)) {
        warning("Provided values for plot.by and scale.by do not match; results may not be very
            interpretable.")
    }

    # check raster and set by number of cells total if NULL
    if (is.null(x = raster)) {
        if (nrow(x = object@cell.data) > 1e5) {
            raster <- TRUE
            message("NOTE: Points are rasterized as number of cells/nuclei plotted exceeds 100,000.
              \n To plot in vector form set `raster = FALSE`.")
        } else {
            raster <- FALSE
        }
    }


    if (use.raw) {
        if (is.null(log2scale)) {
            log2scale <- FALSE
        }
        # drop only outer level names
        if (class(object@raw.data[[1]])[1] == "H5File") {
            if (object@h5file.info[[1]][["sample.data.type"]] != "raw.data"){
                stop("raw.data should be sampled for this plot.")
            }
            gene_vals <- getGeneValues(object@sample.data, gene, log2scale = log2scale)
        } else {
            gene_vals <- getGeneValues(object@raw.data, gene, log2scale = log2scale)
        }
    } else {
        if (is.null(log2scale)) {
            log2scale <- TRUE
        }
        # rescale in case requested gene not highly variable
        if (use.scaled) {
            # check for feature
            if (!(scale.by %in% colnames(object@cell.data)) & scale.by != 'none') {
                stop("Please select existing feature in cell.data to scale.by, or add it before calling.")
            }
            if (class(object@raw.data[[1]])[1] == "H5File") {
                if (object@h5file.info[[1]][["sample.data.type"]] != "norm.data"){
                    stop("norm.data should be sampled for this plot.")
                }
                gene_vals <- getGeneValues(object@sample.data, gene)
                cells <- unlist(lapply(object@sample.data, colnames))
            } else {
                gene_vals <- getGeneValues(object@norm.data, gene)
                cells <- unlist(lapply(object@norm.data, colnames))
            }
            cellnames <- names(gene_vals)
            # set up dataframe with groups
            gene_df <- data.frame(gene = gene_vals)
            if (scale.by == 'none') {
                gene_df[['scaleby']] = 'none'
            } else {
                gene_df[['scaleby']] = factor(object@cell.data[cells,][[scale.by]])
            }
            gene_df1 <- gene_df %>%
                group_by(.data[['scaleby']]) %>%
                # scale by selected feature
                mutate_at(vars(-group_cols()), function(x) { scale(x, center = FALSE)})
            gene_vals <- gene_df1$gene
            if (log2scale) {
                gene_vals <- log2(10000 * gene_vals + 1)
            }
            names(gene_vals) <- cellnames
        } else {
            # using normalized data
            # indicate methylation indices here
            if (class(object@raw.data[[1]])[1] == "H5File") {
                if (object@h5file.info[[1]][["sample.data.type"]] != "norm.data"){
                    stop("norm.data should be sampled for this plot.")
                }
                gene_vals <- getGeneValues(object@sample.data, gene, methylation.indices = methylation.indices,
                                           log2scale = log2scale)
            } else {
                gene_vals <- getGeneValues(object@norm.data, gene, methylation.indices = methylation.indices,
                                           log2scale = log2scale)
            }
        }
    }
    gene_vals[gene_vals == 0] <- NA
    # Extract min and max expression values for plot scaling if keep.scale = TRUE
    if (keep.scale){
        max_exp_val <- max(gene_vals, na.rm = TRUE)
        min_exp_val <- min(gene_vals, na.rm = TRUE)
    }

    if (class(object@raw.data[[1]])[1] == "H5File") {
        cells <- unlist(lapply(object@sample.data, colnames))
        dr_df <- data.frame(object@tsne.coords[cells,])
    } else {
        dr_df <- data.frame(object@tsne.coords)
        rownames(dr_df) <- rownames(object@cell.data)
    }
    dr_df$gene <- as.numeric(gene_vals[rownames(dr_df)])
    colnames(dr_df) <- c("dr1", "dr2", "gene")
    # get dr limits for later
    lim1 <- c(min(dr_df$dr1), max(dr_df$dr1))
    lim2 <- c(min(dr_df$dr2), max(dr_df$dr2))

    if (plot.by != 'none') {
        if (!(plot.by %in% colnames(object@cell.data))) {
            stop("Please select existing feature in cell.data to plot.by, or add it before calling.")
        }
        dr_df$plotby <- factor(object@cell.data[rownames(dr_df),][[plot.by]])
    } else {
        dr_df$plotby <- factor("none")
    }
    # expand clip values if only single provided
    num_levels <- length(levels(dr_df$plotby))
    if (length(min.clip) == 1) {
        min.clip <- rep(min.clip, num_levels)
        names(min.clip) <- levels(dr_df$plotby)
    }
    if (length(max.clip) == 1) {
        max.clip <- rep(max.clip, num_levels)
        names(max.clip) <- levels(dr_df$plotby)
    }
    if (!is.null(min.clip) & is.null(names(min.clip))) {
        if (num_levels > 1) {
            message("Adding names to min.clip according to levels in plot.by group; order may not be
              preserved as intended if multiple clip values passed in. Pass in named vector to
              prevent this.")
        }
        names(min.clip) <- levels(dr_df$plotby)
    }
    if (!is.null(max.clip) & is.null(names(max.clip))) {
        if (num_levels > 1) {
            message("Adding names to max.clip according to levels in plot.by group; order may not be
              preserved as intended if multiple clip values passed in. Pass in named vector to
              prevent this.")
        }
        names(max.clip) <- levels(dr_df$plotby)
    }
    p_list <- list()
    for (sub_df in split(dr_df, f = dr_df$plotby)) {
        # maybe do quantile cutoff here
        group_name <- as.character(sub_df$plotby[1])
        if (!clip.absolute) {
            max_v <- quantile(sub_df$gene, probs = max.clip[group_name], na.rm = TRUE)
            min_v <- quantile(sub_df$gene, probs = min.clip[group_name], na.rm = TRUE)
        } else {
            max_v <- max.clip[group_name]
            min_v <- min.clip[group_name]
        }
        sub_df$gene[sub_df$gene < min_v & !is.na(sub_df$gene)] <- min_v
        sub_df$gene[sub_df$gene > max_v & !is.na(sub_df$gene)] <- max_v

        if (isTRUE(x = raster)) {
            ggp <- ggplot(sub_df, aes_string(x = 'dr1', y = 'dr2', color = 'gene')) + geom_scattermore(pointsize = pt.size) +
                labs(col = gene)
        } else {
            ggp <- ggplot(sub_df, aes_string(x = 'dr1', y = 'dr2', color = 'gene')) + geom_point(size = pt.size) +
                labs(col = gene)
        }

        if (!is.null(cols.use)) {
            if (keep.scale) {
                ggp <- ggp + scale_color_gradientn(colors = cols.use,
                                                   na.value = zero.color,
                                                   limits = c(min_exp_val, max_exp_val))
            } else {
                ggp <- ggp + scale_color_gradientn(colors = cols.use,
                                                   na.value = zero.color)
            }
        } else {
            if (keep.scale) {
                ggp <- ggp + scale_color_viridis_c(option = option,
                                                   direction = -1,
                                                   na.value = zero.color,
                                                   limits = c(min_exp_val, max_exp_val))
            } else {
                ggp <- ggp + scale_color_viridis_c(option = option,
                                                   direction = -1,
                                                   na.value = zero.color)
            }
        }
        if (set.dr.lims) {
            ggp <- ggp + xlim(lim1) + ylim(lim2)
        }

        if (plot.by != 'none') {
            base <- as.character(sub_df$plotby[1])
        } else {
            base <- ""
        }
        ggp <- ggp + ggtitle(base)

        if (!is.null(axis.labels)) {
            ggp <- ggp + xlab(axis.labels[1]) + ylab(axis.labels[2])
        }
        if (!do.legend) {
            ggp <- ggp + theme(legend.position = "none")
        }
        if (points.only) {
            ggp <- ggp + theme(
                axis.line = element_blank(), axis.text.x = element_blank(),
                axis.text.y = element_blank(), axis.ticks = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(), legend.position = "none",
                panel.background = element_blank(), panel.border = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                plot.background = element_blank(), plot.title = element_blank()
            )
        }
        p_list[[as.character(sub_df$plotby[1])]] <- ggp + theme_cowplot(12)
    }
    if (plot.by == 'dataset') {
        p_list <- p_list[names(object@raw.data)]
    }

    if (return.plots){
        if (length(p_list) == 1) {
            return(p_list[[1]])
        } else {
            return(p_list)
        }
    } else {
        for (plot in p_list) {
            print(plot)
        }
    }
}

#' Plot expression of multiple genes
#'
#' Uses plotGene to plot each gene (and dataset) on a separate page. It is recommended to call this
#' function into a PDF due to the large number of plots produced.
#'
#' @param object \code{liger} object. Should call runTSNE before calling.
#' @param genes Vector of gene names.
#' @param ... arguments passed from \code{\link[rliger]{plotGene}}
#'
#' @return If returning single plot, returns ggplot object; if returning multiple plots; returns
#'   list of ggplot objects.
#'
#' @importFrom ggplot2 ggplot geom_point aes_string scale_color_gradient2 ggtitle
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete input
#' ligerex <- runTSNE(ligerex)
#' # plot expression for CD4 and FCGR3A
#' # pdf("gene_plots.pdf")
#' plotGenes(ligerex, c("CD4", "FCGR3A"))
#' # dev.off()
#' }

plotGenes <- function(object, genes, ...) {
    for (i in 1:length(genes)) {
        print(genes[i])
        plotGene(object, genes[i], ...)
    }
}

#' Generate a river (Sankey) plot
#'
#' Creates a riverplot to show how separate cluster assignments from two datasets map onto a
#' joint clustering. The joint clustering is by default the object clustering, but an external one
#' can also be passed in. Uses the riverplot package to construct riverplot object and then plot.
#'
#' @param object \code{liger} object. Should run quantileAlignSNF before calling.
#' @param cluster1 Cluster assignments for dataset 1. Note that cluster names should be distinct
#'   across datasets.
#' @param cluster2 Cluster assignments for dataset 2. Note that cluster names should be distinct
#'   across datasets.
#' @param cluster_consensus Optional external consensus clustering (to use instead of object
#'   clusters)
#' @param min.frac Minimum fraction of cluster for edge to be shown (default 0.05).
#' @param min.cells Minumum number of cells for edge to be shown (default 10).
#' @param river.yscale y-scale to pass to riverplot -- scales the edge with values by this factor,
#'   can be used to squeeze vertically (default 1).
#' @param river.lty Line style to pass to riverplot (default 0).
#' @param river.node_margin Node_margin to pass to riverplot -- how much vertical space to keep
#'   between the nodes (default 0.1).
#' @param label.cex Size of text labels (default 1).
#' @param label.col Color of text labels (defualt "black").
#' @param lab.srt Angle of text labels (default 0).
#' @param river.usr Coordinates at which to draw the plot in form (x0, x1, y0, y1).
#' @param node.order Order of clusters in each set (list with three vectors of ordinal numbers).
#'   By default will try to automatically order them appropriately.
#'
#' @return A riverplot object
#'
#' @importFrom plyr mapvalues
#' @importFrom riverplot makeRiver
#' @importFrom riverplot riverplot
#' @importFrom grDevices hcl
#' @importFrom utils capture.output
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete input
#' # toy clusters
#' cluster1 <- sample(c('type1', 'type2', 'type3'), ncol(ligerex@raw.data[[1]]), replace = TRUE)
#' names(cluster1) <- colnames(ligerex@raw.data[[1]])
#' cluster2 <- sample(c('type4', 'type5', 'type6'), ncol(ligerex@raw.data[[2]]), replace = TRUE)
#' names(cluster2) <- colnames(ligerex@raw.data[[2]])
#' # create riverplot
#' makeRiverplot(ligerex, cluster1, cluster2)
#' }

makeRiverplot <- function(object, cluster1, cluster2, cluster_consensus = NULL, min.frac = 0.05,
                          min.cells = 10, river.yscale = 1, river.lty = 0, river.node_margin = 0.1,
                          label.cex = 1, label.col = "black", lab.srt = 0, river.usr = NULL,
                          node.order = "auto") {
    cluster1 <- droplevels(cluster1)
    cluster2 <- droplevels(cluster2)
    if (is.null(cluster_consensus)) {
        cluster_consensus <- droplevels(object@clusters)
    }
    # Make cluster names unique if necessary
    if (length(intersect(levels(cluster1), levels(cluster2))) > 0 |
        length(intersect(levels(cluster1), levels(cluster_consensus))) > 0 |
        length(intersect(levels(cluster2), levels(cluster_consensus))) > 0) {
        message("Duplicate cluster names detected. Adding 1- and 2- to make unique names.")
        cluster1 <- mapvalues(cluster1, from = levels(cluster1),
                              to = paste("1", levels(cluster1), sep = "-"))
        cluster2 <- mapvalues(cluster2, from = levels(cluster2),
                              to = paste("2", levels(cluster2), sep = "-"))
    }
    cluster1 <- cluster1[intersect(names(cluster1), names(cluster_consensus))]
    cluster2 <- cluster2[intersect(names(cluster2), names(cluster_consensus))]

    # set node order
    if (identical(node.order, "auto")) {
        tab.1 <- table(cluster1, cluster_consensus[names(cluster1)])
        tab.1 <- sweep(tab.1, 1, rowSums(tab.1), "/")
        tab.2 <- table(cluster2, cluster_consensus[names(cluster2)])
        tab.2 <- sweep(tab.2, 1, rowSums(tab.2), "/")
        whichmax.1 <- apply(tab.1, 1, which.max)
        whichmax.2 <- apply(tab.2, 1, which.max)
        ord.1 <- order(whichmax.1)
        ord.2 <- order(whichmax.2)
        cluster1 <- factor(cluster1, levels = levels(cluster1)[ord.1])
        cluster2 <- factor(cluster2, levels = levels(cluster2)[ord.2])
    } else {
        if (is.list(node.order)) {
            cluster1 <- factor(cluster1, levels = levels(cluster1)[node.order[[1]]])
            cluster_consensus <- factor(cluster_consensus,
                                        levels = levels(cluster_consensus)[node.order[[2]]])
            cluster2 <- factor(cluster2, levels = levels(cluster2)[node.order[[3]]])
        }
    }
    cluster1 <- cluster1[!is.na(cluster1)]
    cluster2 <- cluster2[!is.na(cluster2)]
    nodes1 <- levels(cluster1)[table(cluster1) > 0]
    nodes2 <- levels(cluster2)[table(cluster2) > 0]
    nodes_middle <- levels(cluster_consensus)[table(cluster_consensus) > 0]
    node_Xs <- c(
        rep(1, length(nodes1)), rep(2, length(nodes_middle)),
        rep(3, length(nodes2))
    )

    # first set of edges
    edge_list <- list()
    for (i in 1:length(nodes1)) {
        temp <- list()
        i_cells <- names(cluster1)[cluster1 == nodes1[i]]
        for (j in 1:length(nodes_middle)) {
            if (length(which(cluster_consensus[i_cells] == nodes_middle[j])) / length(i_cells) > min.frac &
                length(which(cluster_consensus[i_cells] == nodes_middle[j])) > min.cells) {
                temp[[nodes_middle[j]]] <- sum(cluster_consensus[i_cells] ==
                                                   nodes_middle[j]) / length(cluster1)
            }
        }
        edge_list[[nodes1[i]]] <- temp
    }
    # second set of edges
    cluster3 <- cluster_consensus[names(cluster2)]
    for (i in 1:length(nodes_middle)) {
        temp <- list()
        i_cells <- names(cluster3)[cluster3 == nodes_middle[i]]
        for (j in 1:length(nodes2)) {
            j_cells <- names(cluster2)[cluster2 == nodes2[j]]
            if (length(which(cluster_consensus[j_cells] == nodes_middle[i])) / length(j_cells) > min.frac &
                length(which(cluster_consensus[j_cells] == nodes_middle[i])) > min.cells) {
                if (!is.na(sum(cluster2[i_cells] == nodes2[j]))) {
                    temp[[nodes2[j]]] <- sum(cluster2[i_cells] ==
                                                 nodes2[j]) / length(cluster2)
                }
            }
        }
        edge_list[[nodes_middle[i]]] <- temp
    }
    # set cluster colors
    node_cols <- list()
    ggplotColors <- function(g) {
        d <- 360 / g
        h <- cumsum(c(15, rep(d, g - 1)))
        grDevices::hcl(h = h, c = 100, l = 65)
    }
    pal <- ggplotColors(length(nodes1))
    for (i in 1:length(nodes1)) {
        node_cols[[nodes1[i]]] <- list(col = pal[i], textcex = label.cex,
                                       textcol = label.col, srt = lab.srt)
    }
    pal <- ggplotColors(length(nodes_middle))
    for (i in 1:length(nodes_middle)) {
        node_cols[[nodes_middle[i]]] <- list(col = pal[i], textcex = label.cex,
                                             textcol = label.col, srt = lab.srt)
    }
    pal <- ggplotColors(length(nodes2))
    for (i in 1:length(nodes2)) {
        node_cols[[nodes2[i]]] <- list(col = pal[i], textcex = label.cex,
                                       textcol = label.col, srt = lab.srt)
    }
    # create nodes and riverplot object
    nodes <- list(nodes1, nodes_middle, nodes2)
    node.limit <- max(unlist(lapply(nodes, length)))

    node_Ys <- lapply(1:length(nodes), function(i) {
        seq(1, node.limit, by = node.limit / length(nodes[[i]]))
    })
    rp <- makeRiver(c(nodes1, nodes_middle, nodes2), edge_list,
                    node_xpos = node_Xs, node_ypos = unlist(node_Ys), node_styles = node_cols
    )
    # prevent normal riverplot output being printed to console
    invisible(capture.output(riverplot(rp,
                                       yscale = river.yscale, lty = river.lty,
                                       node_margin = river.node_margin, usr = river.usr
    )))
}

#' Plot cluster proportions by dataset
#'
#' Generates plot of clusters sized by the proportion of total cells
#'
#' @param object \code{liger} object. Should call quantileAlignSNF before calling.
#' @param return.plot Return ggplot object (default FALSE)
#'
#' @return print plot to console (return.plot = FALSE); ggplot object (return.plot = TRUE)
#'   list of ggplot objects.
#'
#' @importFrom grid unit
#' @importFrom ggplot2 ggplot aes coord_fixed element_blank geom_point guides guide_legend
#' scale_size scale_y_discrete theme
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete input
#' ligerex <- quantile_norm(ligerex)
#' # plot cluster proportions
#' plotClusterProportions(ligerex)
#' }

plotClusterProportions <- function(object, return.plot = FALSE) {

    sample_names <- unlist(lapply(seq_along(object@H), function(i) {
        rep(names(object@H)[i], nrow(object@H[[i]]))
    }))
    freq_table <- data.frame(rep(object@clusters, length(object@scale.data)),
                             sample_names)
    freq_table <- table(freq_table[,1], freq_table[,2])
    for (i in 1:ncol(freq_table)) {
        freq_table[, i] <- freq_table[, i] / sum(freq_table[, i])
    }
    freq_table <- as.data.frame(freq_table)
    colnames(freq_table) <- c("Cluster", "Sample", "Proportion")
    p1 <- ggplot(freq_table, aes_string(x = "Cluster", y = "Sample")) +
        geom_point(aes_string(size = 'Proportion', fill = 'Cluster', color = 'Cluster')) +
        scale_size(guide = "none") + theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            legend.position = 'bottom',
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            legend.justification = "center"
        ) + scale_y_discrete(position = "right") +
        guides(fill = guide_legend(ncol = 6, override.aes = list(size = 4))) +
        coord_fixed(ratio = 0.5)
    if (return.plot) {
        return(p1)
    }
    else {
        print(p1)
    }
}

#' Plot heatmap of cluster/factor correspondence
#'
#' Generates matrix of cluster/factor correspondence, using sum of row-normalized factor loadings
#' for every cell in each cluster. Plots heatmap of matrix, with red representing high total
#' loadings for a factor, black low. Optionally can also include dendrograms and sorting for
#' factors and clusters.
#'
#' @param object \code{liger} object.
#' @param use.aligned Use quantile normalized factor loadings to generate matrix (default FALSE).
#' @param Rowv Determines if and how the row dendrogram should be computed and reordered. Either a
#'   dendrogram or a vector of values used to reorder the row dendrogram or NA to suppress any row
#'   dendrogram (and reordering) (default NA for no dendrogram).
#' @param Colv Determines if and how the column dendrogram should be reordered. Has the same options
#'   as the Rowv argument (default 'Rowv' to match Rowv).
#' @param col Color map to use (defaults to red and black)
#' @param return.data Return matrix of total factor loadings for each cluster (default FALSE).
#' @param ... Additional parameters to pass on to heatmap()
#'
#' @return If requested, matrix of size num_cluster x num_factor
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats heatmap
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete input
#' # plot expression for CD4 and return plots
#' loading.matrix <- plotClusterFactors(ligerex, return.data = TRUE)
#' }

plotClusterFactors <- function(object, use.aligned = FALSE, Rowv = NA, Colv = "Rowv", col = NULL,
                               return.data = FALSE, ...) {
    if (use.aligned) {
        data.mat <- object@H.norm
    } else {
        scaled <- lapply(object@H, function(i) {
            scale(i, center = FALSE, scale = TRUE)
        })
        data.mat <- Reduce(rbind, scaled)
    }
    row.scaled <- t(apply(data.mat, 1, function(x) {
        x / sum(x)
    }))
    cluster.bars <- list()
    for (cluster in levels(object@clusters)) {
        cluster.bars[[cluster]] <- colSums(row.scaled[names(object@clusters)
                                                      [which(object@clusters == cluster)], ])

    }
    cluster.bars <- Reduce(rbind, cluster.bars)
    if (is.null(col)) {
        colfunc <- grDevices::colorRampPalette(c("black", "red"))
        col <- colfunc(15)
    }
    rownames(cluster.bars) <- levels(object@clusters)
    colnames(cluster.bars) <- 1:ncol(cluster.bars)
    title <- ifelse(use.aligned, "H.norm", "raw H")
    stats::heatmap(cluster.bars,
                   Rowv = Rowv, Colv = Rowv, col = col, xlab = "Factor", ylab = "Cluster",
                   main = title, ...
    )
    if (return.data) {
        return(cluster.bars)
    }
}
