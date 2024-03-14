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
#' ligerex <- readSubset(ligerex, slot.use = "normData", max.cells = 5000)
#' plotWordClouds(ligerex, num.genes = 20)
#' }

plotWordClouds <- function(object, dataset1 = NULL, dataset2 = NULL, num.genes = 30, min.size = 1,
                           max.size = 4, factor.share.thresh = 10, log.fc.thresh = 1, pval.thresh = 0.05,
                           do.spec.plot = TRUE, return.plots = FALSE, verbose = TRUE) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package \"ggrepel\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages(\"ggrepel\")",
             call. = FALSE)
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
        pb <- utils::txtProgressBar(min = 0, max = length(factors.use), style = 3)
    }
    return_plots <- list()
    for (i in factors.use) {
        tsne_df <- data.frame(H_aligned[, i], tsne_coords)
        factorlab <- paste("Factor", i, sep = "")
        colnames(tsne_df) <- c(factorlab, "Dim1", "Dim2")
        factor_ds <- paste("Factor", i, "Dataset Specificity:", dataset.specificity[[3]][i])
        p1 <- ggplot2::ggplot(tsne_df, ggplot2::aes_string(x = "Dim1", y = "Dim2", color = factorlab)) + ggplot2::geom_point() +
            ggplot2::scale_color_gradient(low = "yellow", high = "red") + ggplot2::ggtitle(label = factor_ds)

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
            out_plot <- ggplot2::ggplot(gene_df, ggplot2::aes(x = 1, y = 1, size = .data[["loadings"]], label = .data[['genes']])) +
                ggrepel::geom_text_repel(force = 100, segment.color = NA) +
                ggplot2::scale_size(range = c(min.size, max.size), guide = FALSE) +
                ggplot2::scale_y_continuous(breaks = NULL) +
                ggplot2::scale_x_continuous(breaks = NULL) +
                ggplot2::labs(x = "", y = "") + ggplot2::ggtitle(label = names_list[[x]]) + ggplot2::coord_fixed() + ggplot2::theme_void()
            return(out_plot)
        })

        p2 <- (cowplot::plot_grid(plotlist = plot_list, align = "hv", nrow = 1)
               + cowplot::draw_grob(grid::roundrectGrob(
                   x = 0.33, y = 0.5, width = 0.67, height = 0.70,
                   gp = grid::gpar(fill = "khaki1", col = "Black", alpha = 0.5, lwd = 2)
               ))
               + cowplot::draw_grob(grid::roundrectGrob(
                   x = 0.67, y = 0.5, width = 0.67, height = 0.70,
                   gp = grid::gpar(fill = "indianred1", col = "Black", alpha = 0.5, lwd = 2)
               )))
        return_plots[[i]] <- cowplot::plot_grid(p1, p2, nrow = 2, align = "h")
        if (!return.plots) {
            print(return_plots[[i]])
        }
        if (verbose) {
            utils::setTxtProgressBar(pb, i)
        }
    }
    if (return.plots) {
        return(return_plots)
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
    if (!requireNamespace("riverplot", quietly = TRUE))
        stop("Package \"riverplot\" needed for this function to work. ",
             "Please install it by command:\n",
             "BiocManager::install('riverplot')",
             call. = FALSE)
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
    rp <- riverplot::makeRiver(c(nodes1, nodes_middle, nodes2), edge_list,
                    node_xpos = node_Xs, node_ypos = unlist(node_Ys), node_styles = node_cols
    )
    # prevent normal riverplot output being printed to console
    invisible(utils::capture.output(riverplot::riverplot(rp,
                                       yscale = river.yscale, lty = river.lty,
                                       node_margin = river.node_margin, usr = river.usr
    )))
}
