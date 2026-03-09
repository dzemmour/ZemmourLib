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

#' Create a Raster Cache for Cell Embeddings
#'
#' Generates a high-resolution PNG of cell embeddings (e.g., UMAP) to be used
#' as a fast background layer in subsequent plots.
#'
#' @param so A Seurat object.
#' @param reduction The name of the reduction to use (e.g., "umap").
#' @param cells A vector of cell IDs to include. Defaults to all cells.
#' @param outfile_prefix Prefix for the generated .png and .rds files.
#' @param width Width of the output image in pixels.
#' @param height Height of the output image in pixels.
#' @param point_size Size of the points.
#' @param point_alpha Transparency of the points.
#' @param point_color Color of the points.
#' @return Invisibly returns a list containing metadata about the cache.
#' @importFrom Seurat Embeddings
#' @importFrom ggplot2 ggplot aes geom_point coord_cartesian theme_void theme margin ggsave
#' @export
MakeEmbeddingRasterCache <- function(
        so,
        reduction,
        cells = NULL,
        outfile_prefix,
        width = 2048,
        height = 2048,
        point_size = 1,
        point_alpha = 1,
        point_color = "grey"
) {

    emb <- Seurat::Embeddings(so[[reduction]])
    if (!is.null(cells)) emb <- emb[cells, , drop = FALSE]

    df <- data.frame(x = emb[,1], y = emb[,2])
    xlim <- range(df$x)
    ylim <- range(df$y)

    pngfile <- paste0(outfile_prefix, ".png")
    metafile <- paste0(outfile_prefix, ".rds")

    # Generate background using ggplot to ensure coordinate alignment
    bg_plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_point(color = point_color, alpha = point_alpha, size = point_size, stroke = 0) +
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0))

    # Save with ggsave to maintain alignment
    ggplot2::ggsave(pngfile, plot = bg_plot, width = width/300, height = height/300, dpi = 300, bg = "transparent")

    meta <- list(
        png = normalizePath(pngfile),
        xlim = xlim,
        ylim = ylim,
        reduction = reduction,
        n_cells = nrow(df)
    )
    saveRDS(meta, metafile)
    invisible(meta)
}

#' Add a Library Raster Background to a ggplot
#'
#' Reads a pre-cached background from the ZemmourLib internal data (inst/extdata)
#' and adds it as a layer using \code{annotation_raster}.
#'
#' @param lineage Character string of the cell type, options: allT, CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP
#'   Must correspond to an .rds file in the library's extdata folder.
#' @param interpolate Logical; whether to interpolate the raster image.
#' @param alpha Transparency of the background layer.
#' @param color Optional; a color to tint the background points (e.g., "grey").
#' @param resolution Optional; a numeric vector of length 2 to resample the image.
#' @return A ggplot2 \code{annotation_raster} layer.
#' @importFrom png readPNG
#' @importFrom grDevices col2rgb
#' @importFrom ggplot2 annotation_raster
#' @importFrom magick image_read image_resize
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # 1. Create your standard DimPlot
#' p <- Seurat::DimPlot(
#'   so_query,
#'   reduction = "mde_incremental",
#'   group.by = "level2",
#'   pt.size = 2,
#'   raster = FALSE
#' )
#'
#' # 2. Generate the background layer from the library
#' bg <- AddEmbeddingRasterBackground(
#'   lineage = "CD4",
#'   interpolate = FALSE,
#'   color = "grey",
#'   alpha = 0.5
#' )
#'
#' # 3. Insert the background layer BEHIND the DimPlot points
#' p$layers <- c(list(bg), p$layers)
#'
#' # 4. Final adjustments
#' p + ggplot2::coord_fixed()
#' }
AddEmbeddingRasterBackground <- function(
        lineage,
        interpolate = TRUE,
        alpha = 1,
        color = NULL,
        resolution = NULL
) {

    # 1. Locate the internal RDS file
    meta_path <- system.file("extdata", sprintf("immgenT_%s_MDE_background.Rds", lineage), package = "ZemmourLib")

    if (meta_path == "") {
        stop(paste0("Could not find background metadata for '", lineage,
                    "' in ZemmourLib. Ensure the .rds exists in inst/extdata."))
    }

    # 2. Load metadata and read the associated PNG
    meta <- readRDS(meta_path)

    # Ensure the PNG path is correctly handled (it might have been absolute during creation)
    # We look for the PNG in the same folder as the RDS
    png_path <- system.file("extdata", sprintf("immgenT_%s_MDE_background.png", lineage), package = "ZemmourLib")

    if (png_path == "") {
        stop(paste0("Could not find background image for '", lineage, "'."))
    }

    img <- png::readPNG(png_path)

    # 3. Recolor if requested
    if (!is.null(color)) {
        rgb_col <- grDevices::col2rgb(color) / 255
        img[,,1] <- rgb_col[1]
        img[,,2] <- rgb_col[2]
        img[,,3] <- rgb_col[3]
    }

    # 4. Adjust transparency
    img[,,4] <- img[,,4] * alpha

    # 5. Resample image if resolution requested
    if (!is.null(resolution)) {
        stopifnot(length(resolution) == 2)
        img_magick <- magick::image_read(img)
        img_resized <- magick::image_resize(img_magick, paste0(resolution[1], "x", resolution[2], "!"))
        img <- as.raster(img_resized)
    }

    # 6. Return the raster layer
    ggplot2::annotation_raster(
        raster = img,
        xmin = meta$xlim[1],
        xmax = meta$xlim[2],
        ymin = meta$ylim[1],
        ymax = meta$ylim[2],
        interpolate = interpolate
    )
}
