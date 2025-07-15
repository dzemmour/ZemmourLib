#' Plot Differential Abundance (DA) Results as Beeswarm Plot
#'
#' This function generates a quasi-random (beeswarm) plot for differential abundance
#' results from Milo, visualizing log fold change and significance.
#'
#' @param da.res A data frame containing differential abundance results (e.g., from `miloR::DA_results`).
#' @param group.by Character string, column in `da.res` to group the beeswarm points by. Default is NULL.
#' @param alpha Numeric, significance threshold for `SpatialFDR`. Default is 0.1.
#' @param subset.nhoods Numeric or character vector, subset of neighborhood IDs to plot. Default is NULL.
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom dplyr mutate arrange
#' @importFrom rlang .data
#' @return A ggplot object.
#' @export
MyplotDAbeeswarm = function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL)
{
    requireNamespace("ggbeeswarm", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("rlang", quietly = TRUE)

    if (!is.null(group.by)) {
        if (!group.by %in% colnames(da.res)) {
            stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ",
                 group.by, ")?")
        }
        # if (is.numeric(da.res[, group.by])) { # Original logic, but it just assigns if numeric.
        # }
        da.res <- dplyr::mutate(da.res, group_by = da.res[, group.by])
    }
    else {
        da.res <- dplyr::mutate(da.res, group_by = "g1")
    }
    if (!is.factor(da.res[, "group_by"])) {
        message("Converting group_by to factor...")
        da.res <- dplyr::mutate(da.res, group_by = factor(.data$group_by,
                                                          levels = unique(.data$group_by)))
    }
    if (!is.null(subset.nhoods)) {
        da.res <- da.res[subset.nhoods, ]
    }
    beeswarm_pos <- ggplot_build(da.res %>% dplyr::mutate(is_signif = ifelse(.data$PValue <
                                                                                 alpha, 1, 0)) %>% dplyr::arrange(.data$group_by) %>% ggplot(aes(.data$group_by,
                                                                                                                                                 .data$logFC)) + ggbeeswarm::geom_quasirandom())
    pos_x <- beeswarm_pos$data[[1]]$x
    pos_y <- beeswarm_pos$data[[1]]$y
    n_groups <- unique(da.res$group_by) %>% length()
    da.res %>% dplyr::mutate(is_signif = ifelse(.data$SpatialFDR < alpha,
                                                1, 0)) %>% dplyr::mutate(logFC_color = ifelse(.data$is_signif == 1,
                                                                                              .data$logFC, NA)) %>% dplyr::arrange(.data$group_by) %>% dplyr::mutate(Nhood = factor(.data$Nhood,
                                                                                                                                                                                    levels = unique(.data$Nhood))) %>% dplyr::mutate(pos_x = pos_x, pos_y = pos_y) %>%
        ggplot(aes(.data$pos_x, .data$pos_y, color = .data$logFC, size = .data$PValue <
                       alpha))  + scale_color_gradient2() + xlab(group.by) + ylab("Log Fold Change") + #guides(color = "none")
        scale_x_continuous(breaks = seq(1, n_groups), labels = setNames(levels(da.res$group_by),
                                                                        seq(1, n_groups))) + geom_point() + coord_flip() +
        theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0))
}
