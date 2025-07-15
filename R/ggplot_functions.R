#' Theme Element: No Grid Lines
#'
#' This helper function returns a `ggplot2` theme element that removes both
#' major and minor grid lines from a plot.
#'
#' @param ... Additional arguments passed to `theme()`.
#' @import ggplot2
#' @return A `theme` object for `ggplot2`.
#' @export
NoGrid <- function(...) {
    requireNamespace("ggplot2", quietly = TRUE)
    no.grid.theme <- theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        validate = TRUE,
        ...
    )
    return(no.grid.theme)
}

#' Theme Element: No Legend
#'
#' This helper function returns a `ggplot2` theme element that removes the legend from a plot.
#' Note: A similar function `NoLegend` exists in Seurat. This provides a standalone version.
#'
#' @param ... Additional arguments passed to `theme()`.
#' @import ggplot2
#' @return A `theme` object for `ggplot2`.
#' @export
NoLegend <- function (...)
{
    requireNamespace("ggplot2", quietly = TRUE)
    no.legend.theme <- theme(legend.position = "none", validate = TRUE,
                             ...)
    return(no.legend.theme)
}
