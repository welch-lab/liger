#' `r lifecycle::badge("experimental")` Suggest optimal K value for the factorization
#' @export
#' @description
#' This function sweeps through a series of k values (number of ranks the
#' datasets are factorized into). For each k value, it repeats the factorization
#' for a number of random starts and obtains the objective errors from each run.
#' The optimal k value is recommended to be the one with the lowest variance.
#'
#' \bold{We are currently actively testing the methodology and the function is
#' subject to change. Please report any issues you encounter.}
#'
#' Currently we have identified that a wider step of k values (e.g. 5, 10, 15,
#' ...) shows a more stable variance than a narrower step (e.g. 5, 6, 7, ...).
#'
#' Note that this function is supposed to take a long time when a larger number
#' of random starts is requested (e.g. 50) for a robust suggestion. It is safe
#' to interupt the progress (e.g. Ctrl+C) and the function will still return
#' the recorded objective errors already completed.
#' @param object A \linkS4class{liger} object.
#' @param kTest A numeric vector of k values to be tested. Default 5, 10, 15,
#' ..., 50.
#' @param nRandomStart Number of random starts for each k value. Default
#' \code{10}.
#' @param lambda Regularization parameter. Default \code{5}.
#' @param nIteration Number of iterations for each run. Default \code{30}.
#' @param nCores Number of cores to use for each run. Default \code{1L}.
#' @param verbose Whether to print progress messages. Default \code{TRUE}.
#' @return A list containing:
#' \item{stats}{A data frame containing the k values, objective errors, and
#' random starts.}
#' \item{figure}{A ggplot2 object showing the objective errors and variance
#' for each k value.}
#' @examples
#' pbmcPlot <- scaleNotCenter(pbmcPlot)
#' suggestK(pbmcPlot, kTest = c(10, 15, 20), nRandomStart = 5)
suggestK <- function(
        object,
        kTest = seq(5, 50, 5),
        nRandomStart = 10,
        lambda = 5,
        nIteration = 30,
        nCores = 1L,
        verbose = getOption("ligerVerbose", TRUE)
) {
    if (isTRUE(verbose)) {
        cli::cli_alert_info("The progress might take long. Completed result will still be returned even if interrupted.")
    }
    scaledList <- scaleData(object)
    resultDF <- data.frame(
        # Should be like 5, 5, 5, 10, 10, 10, ..
        k = rep(kTest, each = nRandomStart),
        objErr = rep(NA, length(kTest*nRandomStart)),
        randomStart = rep(seq(nRandomStart), length(kTest))
    )
    disp <- numeric(length(kTest))
    on.exit({
        return(list(stats = resultDF, figure = .plotSuggestK(resultDF), disp = disp))
    })
    for (k in kTest) {
        cli::cli_progress_bar(
            name = sprintf("Working on k = %d", k),
            total = nRandomStart,
            type = 'tasks'
        )
        conn_mat_trip_accum <- list()
        # conn_mat <- NULL
        for (i in seq(nRandomStart)) {
            res <- RcppPlanc::inmf(
                scaledList,
                k = k,
                lambda = lambda,
                niter = nIteration,
                Hinit = NULL,
                Vinit = NULL,
                Winit = NULL,
                verbose = FALSE,
                nCores = nCores
            )
            # if (is.null(conn_mat)) {
            #     conn_mat <- .H_to_conn_mat(H = res$H)
            # } else {
            #     conn_mat <- conn_mat + .H_to_conn_mat(H = res$H)
            # }
            resultDF[resultDF$k == k & resultDF$randomStart == i, "objErr"] <- res$objErr
            cli::cli_progress_update()
        }
        # conn_mat@x <- conn_mat@x / nRandomStart
        # disp_k <- sum(conn_mat - conn_mat * conn_mat)/2/ncol(object)/(ncol(object) - 1)*8
        # disp[kTest == k] <- disp_k
        cli::cli_progress_done()
    }
    return(list(stats = resultDF, figure = .plotSuggestK(resultDF)))#, disp = disp))
}

.plotSuggestK <- function(stats, ...) {
    bandDF <- stats %>%
        dplyr::group_by(.data[['k']]) %>%
        dplyr::summarise(
            min = min(.data[['objErr']], na.rm = TRUE),
            max = max(.data[['objErr']], na.rm = TRUE)
        )
    varDF <- stats %>%
        dplyr::group_by(.data[['k']]) %>%
        dplyr::summarise(
            variance = stats::var(.data[['objErr']], na.rm = TRUE)
        )
    band_y_top <- max(stats$objErr, na.rm = TRUE)
    band_y_bottom <- min(stats$objErr, na.rm = TRUE)
    bar_y_top <- max(varDF$variance, na.rm = TRUE)
    bar_y_min <- min(varDF$variance, na.rm = TRUE)
    bar_y_bottom <- bar_y_min - 0.1*(bar_y_top - bar_y_min)
    star_y <- bar_y_min - 0.5*(bar_y_min - bar_y_bottom)
    rescale_bar <- function(y2) {
        (y2 - bar_y_bottom) / (bar_y_top - bar_y_bottom)*(band_y_top - band_y_bottom) + band_y_bottom
    }
    best_k <- varDF$k[which.min(varDF$variance)]

    p <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(
            mapping = ggplot2::aes(
                x = .data[['k']],
                ymin = .data[['min']],
                ymax = .data[['max']]
            ),
            data = bandDF,
            fill = 'grey'
        ) +
        ggplot2::geom_point(
            mapping = ggplot2::aes(
                x = .data[['k']],
                y = .data[['objErr']]
            ),
            data = stats
        ) +
        ggplot2::geom_line(
            mapping = ggplot2::aes(
                x = .data[['k']],
                y = rescale_bar(.data[['variance']])
            ),
            data = varDF,
            color = '#54B0E4',
            linewidth = 2,
            alpha = 0.8
        ) +
        ggplot2::geom_point(
            mapping = ggplot2::aes(
                x = best_k,
                y = rescale_bar(star_y)
            ),
            data = data.frame(k = best_k, star_y = star_y),
            color = '#FF0000',
            size = 3,
            shape = 8
        ) +
        # ggplot2::scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
        ggplot2::scale_x_continuous(
            breaks = sort(unique(stats$k))
        ) +
        ggplot2::scale_y_continuous(
            # Create the second y axis for variance
            name = "Objective errors (dot and band)",
            # limits = c(band_y_bottom/2, band_y_top*1.5),
            sec.axis = ggplot2::sec_axis(
                transform = ~ (. - band_y_bottom)/(band_y_top - band_y_bottom)*(bar_y_top - bar_y_bottom) + bar_y_bottom,
                name = "Variance (blue)"
            )#,
            # expand = c(0, 0.05 * max(stats$objErr))
        )
    .ggplotLigerTheme(p, ...)
}

# Make sure H is cell x factor
# Experimental method to examine the dispersion of results
.H_to_conn_mat <- function(H) {
    # Get max factor loading of all cells and concatenate
    H <- Reduce(rbind, H)
    # 1-based clustering assignment returned from CPP code
    H_clust <- max_factor_rcpp(H, dims_use = seq_len(ncol(H)), center = TRUE)
    # Build connectivity matrix from clustering assignment:
    # N x N symmetric upper-triangle sparse matrix that is 1 if two cells are
    # in the same cluster
    accumulate <- list()
    for (i in seq_len(ncol(H))) {
        # For each factor/cluster, find the cells assigned this label
        cellIdx <- which(H_clust == i)
        if (length(cellIdx) < 2) {
            next
        }
        triplets <- t(utils::combn(cellIdx, 2))
        triplets <- cbind(triplets, 1)
        accumulate <- c(accumulate, list(triplets))
    }
    all_trips <- Reduce(rbind, accumulate)
    rm(accumulate)
    conn_mat <- Matrix::sparseMatrix(
        i = all_trips[, 1],
        j = all_trips[, 2],
        x = all_trips[, 3],
        symmetric = TRUE,
        repr = "C"
    )
    return(conn_mat)
}
