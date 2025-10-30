#' Plot Gene Expression on a Subset of Cells with Full Dataset Background
#'
#' Creates a FeaturePlot-like visualization by plotting gene expression for a
#' subset of cells over a grey background of all cells in the dataset.
#'
#' @param seurat_object A Seurat object containing the expression data and UMAP embeddings.
#' @param feature A single character string of the gene name (e.g., "Foxp3").
#' @param cells_to_highlight A character vector of cell IDs to use for the foreground plot.
#' @param umap_to_plot A character string specifying the name of the DimReduc object (e.g., "umap", "mde2_totalvi_20241006").
#' @param assay The assay slot to pull expression data from (e.g., "RNA", "SCT").
#' @param colors A vector of colors for the continuous gradient scale (e.g., output of viridis(10)).
#' @param title A character string for the plot title.
#' @param pixels A numeric vector of length 2 defining the resolution for the background scattermore plot (e.g., c(512, 512)).
#' @param highlight_size A numeric value for the size of the points in the foreground (highlighted cells).
#' @param highlight_alpha A numeric value for the transparency of the points in the foreground (0 to 1).
#' @param xlim A numeric vector of length 2 for the x-axis limits. If NULL, limits are determined automatically.
#' @param ylim A numeric vector of length 2 for the y-axis limits. If NULL, limits are determined automatically.
#' @param mycols DEPRECATED: Use \code{colors} instead. A vector of colors for the continuous gradient scale.
#' @param background_alpha A numeric value for the transparency of the background points (0 to 1).
#' @return A ggplot object.
#' @importFrom Seurat GetAssayData
#' @importFrom scattermore geom_scattermore
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn theme_minimal labs coord_cartesian
#' @export
MyFeaturePlotHighlight <- function(seurat_object, feature, cells_to_highlight, umap_to_plot = "umap", assay = "RNA",
                          colors = viridis::viridis(10), title = "", pixels = c(512, 512),
                          highlight_size = 1, highlight_alpha = 1, xlim = NULL, ylim = NULL,
                          background_alpha = 0.1) {

    # 1. Prepare data for the full background
    dim_data <- seurat_object[[umap_to_plot]]@cell.embeddings
    dim1 <- dim_data[, 1]
    dim2 <- dim_data[, 2]
    tmp_all <- data.frame(seurat_object@meta.data, dim1 = dim1, dim2 = dim2)

    # 2. Create the background plot (all cells in grey)
    p_base <- ggplot2::ggplot(tmp_all) +
        scattermore::geom_scattermore(ggplot2::aes(x = dim1, y = dim2),
                                      color = "grey",
                                      pointsize = 1, # pointsize for scattermore often remains 1
                                      alpha = background_alpha, # NEW: Added background_alpha
                                      pixels = pixels)          # NEW: Added pixels

    # 3. Prepare data for the foreground (highlighted cells with feature expression)
    feature_expr <- Seurat::GetAssayData(seurat_object, assay = assay, slot = "data")[feature, ]

    # Subset the metadata/UMAP data to ONLY the cells being highlighted
    tmp_highlight <- tmp_all[cells_to_highlight, ]
    tmp_highlight$feature_value <- feature_expr[cells_to_highlight]

    # 4. Add the foreground layer
    p_final <- p_base +
        ggplot2::geom_point(data = tmp_highlight,
                            ggplot2::aes(x = dim1, y = dim2, color = feature_value),
                            size = highlight_size, # NEW: Used highlight_size
                            alpha = highlight_alpha) + # NEW: Used highlight_alpha
        ggplot2::scale_color_gradientn(colours = colors, name = feature) +
        ggplot2::theme_void() +
        NoGrid() + # Assuming NoGrid is defined in your package
        ggplot2::labs(title = title, # NEW: Used title
                      x = paste0(umap_to_plot, "_1"),
                      y = paste0(umap_to_plot, "_2"))

    # NEW: Apply axis limits if provided
    if (!is.null(xlim) || !is.null(ylim)) {
        p_final <- p_final + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
    }

    return(p_final)
}


#' Custom DotPlot Function for Seurat Objects
#'
#' Creates a dot plot to visualize gene expression by identity and feature. This is a modified version of Seurat's DotPlot with custom alpha and shape controls.
#'
#' @param object A Seurat object.
#' @param features A vector of features to plot.
#' @param assay The assay to pull the data from.
#' @param cols A vector of colors to use.
#' @param col.min Minimum value for color scaling.
#' @param col.max Maximum value for color scaling.
#' @param dot.min Minimum dot size.
#' @param dot.scale Scale factor for dot sizes.
#' @param idents Identity classes to include in the plot.
#' @param group.by A grouping variable.
#' @param split.by A variable to split the plot by.
#' @param cluster.idents Whether to cluster identities.
#' @param scale Whether to scale the expression.
#' @param scale.by Method for scaling dot size.
#' @param scale.min Minimum value for dot size.
#' @param scale.max Maximum value for dot size.
#' @param point.shape Shape of the points.
#' @param alpha.range A vector of length 2 for alpha range.
#' @param show.alpha.legend Whether to show the alpha legend.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_point scale_radius scale_alpha theme_minimal theme element_blank labs facet_grid scale_fill_identity scale_fill_distiller scale_fill_gradient scale_fill_gradientn guides guide_colorbar guide_legend
#' @importFrom Seurat `%||%` DefaultAssay CellsByIdentities FetchData MinMax
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom grid unit
#' @export
#' @examples
#' \dontrun{
#' # This example demonstrates the use of ImmgenTFeaturePlots.
#' MyDotPlot(object = so, top_genes, assay = "RNA", group.by = "annotation_level2", dot.scale = 6, scale = T, scale.by = "radius",point.shape = 21, alpha.range = c(0, 1), show.alpha.legend = F)  +
#' labs(y = topic, x = NULL) + scale_fill_gradientn(colours =  colorRampPalette(c("grey",'white','red'))(10)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
#' } # End \dontrun

MyDotPlot <- function (
        object, features, assay = NULL,
        cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5,
        dot.min = 0, dot.scale = 6, idents = NULL, group.by = NULL,
        split.by = NULL, cluster.idents = FALSE, scale = TRUE,
        scale.by = "radius", scale.min = NA, scale.max = NA,
        point.shape = 21,
        alpha.range = c(0.2, 1.0),
        show.alpha.legend = FALSE
) {
    # This code has been placed directly into your package.
    # The original 'if (!requireNamespace(...))' checks are no longer needed
    # as the package dependencies are handled by the Imports field in DESCRIPTION.

    `%||%` <- Seurat::`%||%`
    assay <- assay %||% Seurat::DefaultAssay(object)
    Seurat::DefaultAssay(object) <- assay

    split.colors <- !is.null(split.by) && !any(cols %in% rownames(RColorBrewer::brewer.pal.info))
    scale.func <- switch(
        scale.by,
        size   = ggplot2::scale_size,
        radius = ggplot2::scale_radius,
        stop("'scale.by' must be either 'size' or 'radius'")
    )

    feature.groups <- NULL
    if (is.list(features) || any(!is.na(names(features)))) {
        feature.groups <- unlist(lapply(seq_along(features), function(i) rep(names(features)[i], length(features[[i]]))))
        if (any(is.na(feature.groups))) warning("Some feature groups are unnamed.", call. = FALSE, immediate. = TRUE)
        features <- unlist(features)
        names(feature.groups) <- features
    }

    cells <- unlist(Seurat::CellsByIdentities(object, cells = colnames(object[[assay]]), idents = idents))
    data.features <- Seurat::FetchData(object, vars = features, cells = cells)

    data.features$id <- if (is.null(group.by)) Seurat::Idents(object)[cells, drop = TRUE] else object[[group.by, drop = TRUE]][cells, drop = TRUE]
    if (!is.factor(data.features$id)) data.features$id <- factor(data.features$id)
    id.levels <- levels(data.features$id)
    data.features$id <- as.vector(data.features$id)

    if (!is.null(split.by)) {
        splits <- Seurat::FetchData(object, vars = split.by)[cells, split.by]
        if (split.colors) {
            if (length(unique(splits)) > length(cols)) stop(sprintf("Need at least %d colors in 'cols'", length(unique(splits))))
            cols <- cols[seq_along(unique(splits))]
            names(cols) <- unique(splits)
        }
        data.features$id <- paste(data.features$id, splits, sep = "_")
        us <- unique(splits)
        id.levels <- paste0(rep(id.levels, each = length(us)), "_", rep(us, times = length(id.levels)))
    }

    data.plot <- lapply(unique(data.features$id), function(ident) {
        df <- data.features[data.features$id == ident, 1:(ncol(data.features) - 1), drop = FALSE]
        avg.exp <- apply(df, 2, function(x) mean(expm1(x)))
        pct.exp <- apply(df, 2, Seurat:::PercentAbove, threshold = 0)
        list(avg.exp = avg.exp, pct.exp = pct.exp)
    })
    names(data.plot) <- unique(data.features$id)

    if (cluster.idents) {
        mat <- do.call(rbind, lapply(data.plot, unlist))
        mat <- scale(mat)
        id.levels <- id.levels[hclust(dist(mat))$order]
    }

    data.plot <- do.call(rbind, lapply(names(data.plot), function(x) {
        d <- as.data.frame(data.plot[[x]])
        d$features.plot <- rownames(d); d$id <- x; d
    }))
    if (!is.null(id.levels)) data.plot$id <- factor(data.plot$id, levels = id.levels)

    ngroup <- length(levels(data.plot$id))
    if (ngroup == 1) { scale <- FALSE; warning("Only one identity present; expression values will not be scaled", call. = FALSE, immediate. = TRUE) }
    else if (ngroup < 5 && scale) { warning("Scaling with few groups can be misleading", call. = FALSE, immediate. = TRUE) }

    # color-like variable (kept for fill mapping)
    avg.exp.scaled <- sapply(unique(data.plot$features.plot), function(x) {
        v <- data.plot[data.plot$features.plot == x, "avg.exp"]
        if (scale) Seurat:::MinMax(scale(log1p(v)), min = col.min, max = col.max) else log1p(v)
    })
    avg.exp.scaled <- as.vector(t(avg.exp.scaled))
    if (split.colors) avg.exp.scaled <- as.numeric(cut(avg.exp.scaled, breaks = 20))
    data.plot$avg.exp.scaled <- avg.exp.scaled

    # alpha driver: per-feature log1p → [0,1]
    avg.exp.alpha <- sapply(unique(data.plot$features.plot), function(x) {
        v <- data.plot[data.plot$features.plot == x, "avg.exp"]; v <- log1p(v); Seurat:::MinMax(v, min = 0, max = 1)
    })
    data.plot$alpha <- as.vector(t(avg.exp.alpha))

    data.plot$features.plot <- factor(data.plot$features.plot, levels = features)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100

    if (split.colors) {
        splits.use <- unlist(lapply(data.plot$id, function(x)
            sub(paste0(".*_(", paste(sort(unique(splits), decreasing = TRUE), collapse = "|"), ")$"), "\\1", x)))
        data.plot$colors <- mapply(function(color, value) colorRampPalette(c("grey", color))(20)[value],
                                   color = cols[splits.use], value = avg.exp.scaled)
    }

    fill.by <- if (split.colors) "colors" else "avg.exp.scaled"

    if (!is.na(scale.min)) data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    if (!is.na(scale.max)) data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    if (!is.null(feature.groups)) data.plot$feature.groups <- factor(feature.groups[data.plot$features.plot], levels = unique(feature.groups))

    plot <- ggplot2::ggplot(data.plot, ggplot2::aes(x = .data$features.plot, y = .data$id)) +
        ggplot2::geom_point(
            ggplot2::aes(size = .data$pct.exp,
                         fill = .data[[fill.by]],
                         alpha = .data$alpha),
            shape = point.shape, colour = "black", stroke = 0.25
        ) +
        scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
        ggplot2::scale_alpha(range = alpha.range, limits = c(0, 1)) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())

    if (!is.null(feature.groups)) {
        plot <- plot +
            ggplot2::facet_grid(~feature.groups, scales = "free_x", space = "free_x", switch = "y") +
            ggplot2::theme(panel.spacing = grid::unit(1, "lines"), strip.background = ggplot2::element_blank())
    }

    if (split.colors) {
        plot <- plot + ggplot2::scale_fill_identity()
    } else if (length(cols) == 1) {
        plot <- plot + ggplot2::scale_fill_distiller(palette = cols)
    } else if (length(cols) == 2) {
        plot <- plot + ggplot2::scale_fill_gradient(low = cols[1], high = cols[2])
    } else {
        plot <- plot + ggplot2::scale_fill_gradientn(colours = cols)
    }

    if (!split.colors) {
        plot <- plot +
            NoGrid() +
            ggplot2::guides(
                fill  = ggplot2::guide_colorbar(title = "Average Expression"),
                size  = ggplot2::guide_legend(title = "Percent Expressed"),
                alpha = if (show.alpha.legend) ggplot2::guide_legend(title = "Avg exp (α)") else "none"
            )
    }
    return(plot)
}

#' Generate Feature Plots for ImmgenT Single-Cell Data
#'
#' This function creates a series of UMAP and Feature plots tailored for ImmgenT
#' single-cell data, visualizing overall expression and expression within
#' different annotation levels (e.g., CD4, CD8 subsets). It assumes data is
#' already normalized.
#'
#' @param so_orig A Seurat object containing the overall dataset.
#' @param so_list A named list of Seurat objects, typically created by
#'   `SplitObject(so_orig, split.by = "annotation_level1")`. Defaults to this split.
#' @param gene Character string, the gene to plot expression for.
#' @param slot Character string, the data slot to pull gene expression from (e.g., "data", "counts"). Defaults to "data".
#' @param reduction_overall Character string, the name of the dimensionality reduction
#'   to use for the overall (unsplit) plots. Defaults to "mde2_totalvi_20241006".
#' @param reduction_level2 Character string, the name of the dimensionality reduction
#'   to use for the level 2 (split) plots. Defaults to "mde_incremental".
#' @param mypal_level1 A character vector of colors for `annotation_level1`.
#'   Defaults to `ZemmourLib::immgent_colors[["level1"]]`.
#' @param mypal_level2 A character vector of colors for `annotation_level2`.
#'   Defaults to `ZemmourLib::immgent_colors[["level2"]]`.
#' @importFrom Seurat DimPlot FeaturePlot SplitObject NormalizeData
#' @importFrom ggplot2 scale_color_manual ggtitle theme_void
#' @return Prints several ggplot objects and arranges them using `grid.arrange`.
#' @export
#' @examples
#' \dontrun{
#' # This example demonstrates the use of ImmgenTFeaturePlots.
#' so_orig = NormalizeData(so_orig,assay = "RNA", normalization.method = "LogNormalize",verbose = T)
#' so_list = SplitObject(so_orig, split.by = "annotation_level1") # Split the seurat object
#' so_list[["thymocyte"]] = NULL # Remove thymocytes from the lists
#' mypal_level2 = immgent_colors$level2
#' mygene = "Cd4" # The gene to plot
#' pdf(sprintf("FeaturePlot_%s.pdf",mygene), width = 20, height = 20, useDingbats = F)
#' ImmgenTFeaturePlots(so_orig, so_list, gene = mygene, slot = "data")
#' dev.off()
#' } # End \dontrun
ImmgenTFeaturePlots = function(so_orig, so_list = Seurat::SplitObject(so_orig, split.by = "annotation_level1"), gene = "Cd4", slot = "data", reduction_overall = "mde2_totalvi_20241006", reduction_level2 = "mde_incremental",mypal_level1 = ZemmourLib::immgent_colors[["level1"]], mypal_level2 = ZemmourLib::immgent_colors[["level2"]]) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for this function but is not installed. Please install it using install.packages('gridExtra').", call. = FALSE)
    }
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for this function but is not installed. Please install it using install.packages('gridExtra').", call. = FALSE)
    }

    level1_annotations = c("CD4", "CD8", "Treg", "gdT", "CD8aa", "nonconv", "DN", "DP")

    ps1_overall = Seurat::DimPlot(so_orig, reduction = reduction_overall, group.by = "annotation_level1", raster = T, raster.dpi = c(1024,1024), label = T) +
        ggplot2::scale_color_manual(values = mypal_level1) +
        ggplot2::ggtitle("") +
        ggplot2::theme_void() +
        ZemmourLib::NoLegend() # Using your package's NoLegend function
    print(ps1_overall)

    ps2_overall = Seurat::FeaturePlot(so_orig, reduction = reduction_overall, slot = slot, features = gene, order = FALSE, raster = T, raster.dpi = c(1024,1024)) +
        ggplot2::ggtitle(sprintf("%s", gene)) +
        ggplot2::theme_void() +
        ZemmourLib::NoLegend()
    print(ps2_overall)

    ps1_level2_list = list()
    ps2_level2_list = list()

    for (l in level1_annotations) {
        if (!is.null(so_list[[l]])) { # Check if the split object for this level exists
            print(l)
            ps1_level2_list[[l]] = Seurat::DimPlot(so_list[[l]], reduction = reduction_level2, group.by = "annotation_level2", raster = T, raster.dpi = c(1024,1024), label = T) +
                ggplot2::scale_color_manual(values = mypal_level2) +
                ggplot2::ggtitle(sprintf("MDE %s", l)) +
                ggplot2::theme_void() +
                ZemmourLib::NoLegend()
            ps2_level2_list[[l]] = Seurat::FeaturePlot(so_list[[l]], reduction = reduction_level2, slot = slot, features = gene, order = FALSE, raster = T, raster.dpi = c(1024,1024)) +
                ggplot2::ggtitle(sprintf("MDE %s, gene %s", l, gene)) +
                ggplot2::theme_void() +
                ZemmourLib::NoLegend()
        } else {
            message(paste("Skipping", l, "as corresponding Seurat object not found in so_list."))
        }
    }
    gridExtra::grid.arrange(grobs = ps1_level2_list, ncol = 3, nrow = ceiling(length(ps1_level2_list)/3)) # Adjusted nrow
    gridExtra::grid.arrange(grobs = ps2_level2_list, ncol = 3, nrow = ceiling(length(ps2_level2_list)/3)) # Adjusted nrow

}


#' Generate Feature Scatter Plots with Highlighted Cell Density
#'
#' This function creates scatter plots (e.g., for ADT features) highlighting a
#' specific population by displaying their density, similar to flow cytometry plots.
#'
#' @param so A Seurat object. Default is `so`.
#' @param assay Character string, name of the assay to pull features from (e.g., "ADT", "RNA"). Default is "ADT".
#' @param slot Character string, data slot to use (e.g., "data", "counts", "scale.data"). Default is "data".
#' @param feature1 Character string, name of the first feature (e.g., gene or ADT).
#' @param feature2 Character string, name of the second feature.
#' @param group.by Character string, metadata column to define the highlighted group (should be logical TRUE/FALSE). Default is `sprintf("is_%s", cl)`.
#' @param split.by Character string, metadata column to split plots by using `facet_wrap`. Default is NULL.
#' @param raster Logical, whether to use `geom_scattermore` for highlighted points (TRUE) or `geom_point` (FALSE). Default is TRUE.
#' @param cols A character vector of colors to define the density color ramp. Default is `rev(rainbow(10, end = 4/6))`.
#' @param highlight_size Numeric, point size for highlighted cells. Default is 1.
#' @param highlight_alpha Numeric, alpha transparency for highlighted cells. Default is 1.
#' @import ggplot2
#' @importFrom scattermore geom_scattermore
#' @importFrom grDevices densCols colorRampPalette
#' @importFrom grDevices rainbow
#' @return A ggplot object.
#' @export
MyFeatureScatter = function(so = so,
                            assay = "ADT",
                            slot = "data",
                            feature1 = NULL,
                            feature2 = NULL,
                            group.by = sprintf("is_%s", cl),
                            split.by = NULL,
                            raster = T,
                            cols = rev(rainbow(10, end = 4/6)),
                            highlight_size = 1,
                            highlight_alpha = 1) {
    requireNamespace("scattermore", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)

    if (is.null(feature1) || is.null(feature2)) {
        stop("Please provide feature1 and feature2 for MyFeatureScatter.")
    }
    # Ensure features exist in the specified assay and slot
    if (!feature1 %in% rownames(so[[assay]][slot])) {
        stop(paste0("Feature1 '", feature1, "' not found in assay '", assay, "', slot '", slot, "'."))
    }
    if (!feature2 %in% rownames(so[[assay]][slot])) {
        stop(paste0("Feature2 '", feature2, "' not found in assay '", assay, "', slot '", slot, "'."))
    }

    df = data.frame(feature1 = so[[assay]][slot][feature1,],feature2 = so[[assay]][slot][feature2,], group.by = so@meta.data[,group.by])
    if (!is.null(split.by)) {
        df$split.by = so@meta.data[,split.by]
    }

    df2 = df[df$group.by == T,]
    df2$density = densCols(df2$feature1, df2$feature2, colramp = colorRampPalette(cols), nbin = 500)

    p1 = ggplot(df) + scattermore::geom_scattermore(aes(feature1, feature2), color = "grey", pixels = c(512, 512))
    if (raster == T ){
        p2 = scattermore::geom_scattermore(data = df2, aes(feature1, feature2, color = density), pointsize = highlight_size, pixels = c(216, 216))
    } else {
        p2 = geom_point(data = df2, aes(feature1, feature2, color = density), size = highlight_size, alpha = highlight_alpha)
    }
    if(is.null(split.by)) {
        p3 = p1 + p2 + xlab(feature1) + ylab(feature2) +scale_colour_identity()+  theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    } else {
        p3 = p1 + p2 + xlab(feature1) + ylab(feature2) +scale_colour_identity()+  theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~split.by)
    }
    return(p3)
}

#' Generate UMAP Plots for Seurat Object
#'
#' This function generates multiple UMAP plots for a Seurat object,
#' including coloring by a metadata column, splitting by metadata columns,
#' and visualizing gene expression.
#'
#' @param seurat_object A Seurat object. Default is `so`.
#' @param dim1 Numeric vector, UMAP dimension 1 embeddings. Default is `so[["umap_unintegrated"]]@cell.embeddings[,1]`.
#' @param dim2 Numeric vector, UMAP dimension 2 embeddings. Default is `so[["umap_unintegrated"]]@cell.embeddings[,2]`.
#' @param color_by Character string, metadata column to color points by. Default is "spleen_standard".
#' @param split_by1 Character string, metadata column to split plots by (rows in facet_grid). Default is "IGT". Can be NULL.
#' @param split_by2 Character string, metadata column to split plots by (columns in facet_grid). Default is NULL.
#' @param genes Character vector, genes to plot expression for. Default is `c("Foxp3", "Il2ra")`.
#' @param cluster_key Character string, key for cluster information in metadata. Default is "ClusterSCVI_Res".
#' @param mypal A color palette to use for discrete colors. Default is `glasbey()`.
#' @import ggplot2
#' @importFrom ggrastr geom_point_rast
#' @return Prints multiple ggplot objects.
#' @export
MyPlots = function (seurat_object = so, dim1 = seurat_object[["umap_unintegrated"]]@cell.embeddings[,1], dim2 = seurat_object[["umap_unintegrated"]]@cell.embeddings[,2], color_by = "spleen_standard", split_by1 = "IGT", split_by2 =  NULL, genes = c("Foxp3", "Il2ra"), cluster_key = "ClusterSCVI_Res", mypal = NULL) {

    so = seurat_object
    #so@meta.data[,"split_by"] = so@meta.data[,split_by] # This line modifies the Seurat object, generally not recommended in plot functions
    so@meta.data[,"color_by"] = factor(so@meta.data[,color_by])
    so@meta.data[,"split_by1"] = so@meta.data[,split_by1]
    so@meta.data[,"split_by2"] = so@meta.data[,split_by2]

    alpha = 0.5

    message("Plot 1: UMAP")
    plot1 = ggplot(data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)) +
        ggrastr::geom_point_rast(aes(dim1, dim2, color = color_by), alpha = I(alpha), raster.dpi = 100) +
        theme_bw()
    if (!is.null(mypal)) {
        plot1 = plot1 + scale_color_manual(values = mypal)
    }
    print(plot1)

    message("Plot 2: UMAP split")
    tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)
    bkgrd = data.frame(dim1 = dim1, dim2 = dim2)

    p = ggplot(bkgrd) + ggrastr::geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, alpha = 0.2, raster.dpi = 50)
    p2 = ggrastr::geom_point_rast(data = tmp, aes(dim1, dim2, color = color_by), size = 1,  alpha = alpha)

    if (is.null(split_by1)) {
        plot2 = p + p2
        if (!is.null(mypal)) {
            plot2 = plot2 + scale_color_manual(values = mypal)
        }
        plot2 = plot2 + theme_bw() + ggtitle(label = sprintf("color: %s", color_by))
    } else if (is.null(split_by2)) {
        plot2 = p + p2
        if (!is.null(mypal)) {
            plot2 = plot2 + scale_color_manual(values = mypal)
        }
        plot2 = plot2 + theme_bw() + facet_wrap(facets = vars(so@meta.data[,"split_by1"])) + ggtitle(label = sprintf("color: %s, grid: %s", color_by, split_by1))
    } else {
        plot2 = p + p2
        if (!is.null(mypal)) {
            plot2 = plot2 + scale_color_manual(values = mypal)
        }
        plot2 = plot2 + theme_bw() + facet_grid(rows = vars(so@meta.data[,"split_by1"]), cols = vars(so@meta.data[,"split_by2"])) + ggtitle(label = sprintf("color: %s, grid: %s x %s", color_by, split_by1, split_by2))
    }
    print(plot2)

    message("Plot 3: UMAP genes")

    for (g in genes) {
        print(g)
        tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2, size = so@assays$RNA$counts[rownames(so@assays$RNA$counts) %in% g,])
        bkgrd = data.frame(dim1 = dim1, dim2 = dim2)

        p = ggplot(bkgrd) + ggrastr::geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, raster.dpi = 50)
        p2 = geom_point(data = tmp, aes(dim1, dim2, color = size > 0, size = size,  alpha = size > 0))


        if (is.null(split_by1)) {
            plot3 = p + p2  +
                scale_color_manual(values = c("black", "red")) +
                scale_alpha_manual(values = c(0.5,1))  +
                theme_bw()  +
                ggtitle(label = sprintf("gene: %s", g))
        } else if (is.null(split_by2)) {
            plot3 = p + p2  +
                scale_color_manual(values = c("black", "red")) +
                scale_alpha_manual(values = c(0.5,1))  +
                theme_bw()  +
                facet_wrap(facets = vars(so@meta.data[,"split_by1"])) +
                ggtitle(label = sprintf("gene: %s, grid: %s", g, split_by1))
        } else {
            plot3 = p + p2  +
                scale_color_manual(values = c("black", "red")) +
                scale_alpha_manual(values = c(0.5,1))  +
                theme_bw() +
                facet_grid(rows = vars(so@meta.data[,"split_by1"]), cols = vars(so@meta.data[,"split_by2"])) +
                ggtitle(label = sprintf("gene: %s, grid: %s x %s", g, split_by1, split_by2)) # Corrected placeholder
        }
        print(plot3)

    }

}



#' Generate UMAP Plot with Highlighted Cells
#'
#' This function creates UMAP plots where a subset of cells can be highlighted
#' and optionally labeled. It returns a list of ggplot objects.
#'
#' @param seurat_object A Seurat object. Default is `so`.
#' @param umap_to_plot Character string, name of the UMAP embedding to use. Default is "mde_totalvi_20241201_gdT_rmIGTsample".
#' @param cells_to_highlight Character vector, names of cells to highlight. Default highlights cells where `so$nonconv_tcr_recog == TRUE`.
#' @param highlight_column_name Character string, name of the metadata column used for highlighting/labeling. Default is "nonconv_tcr_recog".
#' @param title Character string, plot title. Default is "".
#' @param labelclusters Logical, whether to add labels for highlighted clusters. Default is TRUE.
#' @param pixels Numeric vector of length 2, pixel dimensions for `geom_scattermore`. Default is `c(512, 512)`.
#' @param background_color Logical, whether to color the backrgound cells with same highlight_column_name. Default is FALSE (grey)
#' @param background_alpha Numeric, alpha transparency for background cells (only when background_color = T). Default is 0.1.
#' @param highlight_size Numeric, point size for highlighted cells. Default is 1.
#' @param highlight_alpha Numeric, alpha transparency for highlighted cells. Default is 1.
#' @param xlim Numeric vector of length 2, x-axis limits. Default is NULL (auto).
#' @param ylim Numeric vector of length 2, y-axis limits. Default is NULL (auto).
#' @param mycols A named character vector of colors for the `highlight_column_name`. Default is `mypal`.
#' @param print_plot1 Logical, whether to include the first plot (with legend) in the returned list. Default is TRUE.
#' @param print_plot2 Logical, whether to include the second plot (no legend) in the returned list. Default is TRUE.
#' @import ggplot2
#' @importFrom scattermore geom_scattermore
#' @importFrom rlang sym
#' @importFrom dplyr filter group_by summarise
#' @importFrom Seurat NoLegend
#' @return A list containing ggplot objects: `plot1` (with legend), `plot2` (no legend), and optionally `plot3` (with labels).
#' @export
MyDimPlotHighlight = function (seurat_object = so, umap_to_plot = "mde_totalvi_20241201_gdT_rmIGTsample",
                                cells_to_highlight = names(which(seurat_object$nonconv_tcr_recog ==
                                                                     T)), highlight_column_name = "nonconv_tcr_recog", title = "",
                                labelclusters = T, pixels = c(512, 512), highlight_size = 1,
                                highlight_alpha = 1, xlim = NULL, ylim = NULL, mycols = NULL,
                                print_plot1 = TRUE, print_plot2 = TRUE, background_color = F, background_alpha = 0.1)
{
    requireNamespace("scattermore", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("rlang", quietly = TRUE)
    requireNamespace("Seurat", quietly = TRUE)
    if (is.null(mycols)) {
        mycols <- c(`TRUE` = "red", `FALSE` = "grey")
    }
    missing_names = setdiff(unique(seurat_object@meta.data[[highlight_column_name]]),
                            names(mycols))
    if (length(missing_names) > 0) {
        if (any(is.na(names(mycols)))) {
            new_colors = mycols[which(is.na(names(mycols)))][seq_len(length(missing_names))]
        }
        else {
            new_colors = mycols[seq_len(length(missing_names))]
        }
        names(new_colors) <- missing_names
        color_mapping <- c(mycols[!is.na(names(mycols))], new_colors)
    }
    else {
        color_mapping = mycols
    }
    dim1 <- seurat_object[[umap_to_plot]]@cell.embeddings[, 1]
    dim2 <- seurat_object[[umap_to_plot]]@cell.embeddings[, 2]
    tmp <- data.frame(seurat_object@meta.data, dim1 = dim1, dim2 = dim2)
    tmp2 <- tmp[cells_to_highlight, ]
    highlight_values <- unique(tmp2[[highlight_column_name]])
    if (background_color) {
        bkrg <- ggplot(tmp) + scattermore::geom_scattermore(aes(dim1,
                                                                dim2, color = !!rlang::sym(highlight_column_name)), pointsize = 0.5, alpha = background_alpha,
                                                            pixels = pixels)
    } else {
        bkrg <- ggplot(tmp) + scattermore::geom_scattermore(aes(dim1,
                                                                dim2), color = "grey50", pointsize = 0.5, alpha = 0.5,
                                                            pixels = pixels)
    }
    p2 <- geom_point(data = tmp2, aes(dim1, dim2, color = !!rlang::sym(highlight_column_name)),
                     size = highlight_size, alpha = highlight_alpha)

    plot1 <- bkrg + p2 + scale_color_manual(values = color_mapping) +
        scale_x_continuous(limits = xlim) + scale_y_continuous(limits = ylim) +
        ggtitle(title) + theme_minimal() + ZemmourLib::NoGrid()
    plot2 <- bkrg + p2 + scale_color_manual(values = color_mapping) +
        scale_x_continuous(limits = xlim) + scale_y_continuous(limits = ylim) +
        theme_void() + Seurat::NoLegend() + ZemmourLib::NoGrid()
    plot3 <- NULL
    if (labelclusters) {
        tmp_labels <- tmp2 %>% dplyr::filter(!!rlang::sym(highlight_column_name) %in%
                                                 highlight_values) %>% dplyr::group_by(!!rlang::sym(highlight_column_name)) %>%
            dplyr::summarise(dim1 = median(dim1), dim2 = median(dim2))
        plot3 <- plot2 + geom_text(data = tmp_labels, aes(x = dim1,
                                                          y = dim2, label = !!rlang::sym(highlight_column_name),
                                                          color = !!rlang::sym(highlight_column_name)), size = 4,
                                   show.legend = FALSE) + scale_color_manual(values = color_mapping) +
            scale_x_continuous(limits = xlim) + scale_y_continuous(limits = ylim) +
            ggtitle(title) + ZemmourLib::NoGrid()
    }
    plot_list = list(plot1 = plot1, plot2 = plot2, plot3 = plot3)
    return(plot_list)
}

# MyDimPlotHighlight <- function(seurat_object = so,
#                                umap_to_plot = "mde_totalvi_20241201_gdT_rmIGTsample",
#                                cells_to_highlight = names(which(seurat_object$nonconv_tcr_recog == T)),
#                                highlight_column_name = "nonconv_tcr_recog",
#                                title = "",
#                                labelclusters = T,
#                                pixels = c(512, 512),
#                                highlight_size = 1,
#                                highlight_alpha = 1,
#                                xlim = NULL,
#                                ylim = NULL,
#                                mycols = NULL,
#                                print_plot1 = TRUE,
#                                print_plot2 = TRUE) {
#
#     requireNamespace("scattermore", quietly = TRUE) # Ensures scattermore is available
#     requireNamespace("ggplot2", quietly = TRUE) # Ensures ggplot2 is available
#     requireNamespace("dplyr", quietly = TRUE)
#     requireNamespace("rlang", quietly = TRUE)
#     requireNamespace("Seurat", quietly = TRUE)
#
#     # Handle default mycols if not provided
#     if (is.null(mycols)) {
#         # Using a simple default if mypal is not defined globally or passed
#         mycols <- c("TRUE" = "red", "FALSE" = "grey")
#     }
#
#     missing_names = setdiff(unique(seurat_object@meta.data[[highlight_column_name]]), names(mycols))
#
#     if (length(missing_names) > 0) {
#         if (any(is.na(names(mycols)))) {
#             new_colors = mycols[which(is.na(names(mycols)))][seq_len(length(missing_names))]
#         } else {
#             # Fallback to a simple palette if no NA names, just pick from existing
#             new_colors = mycols[seq_len(length(missing_names))]
#         }
#         names(new_colors) <- missing_names
#         color_mapping <- c(mycols[!is.na(names(mycols))], new_colors)
#     } else {
#         color_mapping = mycols
#     }
#
#     # UMAP embeddings
#     dim1 <- seurat_object[[umap_to_plot]]@cell.embeddings[, 1]
#     dim2 <- seurat_object[[umap_to_plot]]@cell.embeddings[, 2]
#     tmp <- data.frame(seurat_object@meta.data, dim1 = dim1, dim2 = dim2)
#
#     # Subset for highlighted cells
#     tmp2 <- tmp[cells_to_highlight, ]
#
#     # Filter unique cluster labels present in cells_to_highlight
#     highlight_values <- unique(tmp2[[highlight_column_name]])
#
#     # Background plot
#     bkrg <- ggplot(tmp) +
#         scattermore::geom_scattermore(aes(dim1, dim2), color = "grey50", pointsize = 0.5, alpha = 0.5, pixels = pixels)
#
#     # Highlighted points
#     p2 <- geom_point(data = tmp2,
#                      aes(dim1, dim2, color = !!rlang::sym(highlight_column_name)),
#                      size = highlight_size,
#                      alpha = highlight_alpha)
#
#     # Plots with consistent colors
#     plot1 <- bkrg + p2 +
#         scale_color_manual(values = color_mapping) +
#         scale_x_continuous(limits = xlim) +
#         scale_y_continuous(limits = ylim) +
#         ggtitle(title) +
#         theme_minimal()
#
#     plot2 <- bkrg + p2 +
#         scale_color_manual(values = color_mapping) +
#         scale_x_continuous(limits = xlim) +
#         scale_y_continuous(limits = ylim) +
#         theme_void() + Seurat::NoLegend()
#
#     # Add cluster labels dynamically for only the present values
#     plot3 <- NULL
#     if (labelclusters) {
#         tmp_labels <- tmp2 %>%
#             dplyr::filter(!!rlang::sym(highlight_column_name) %in% highlight_values) %>%
#             dplyr::group_by(!!rlang::sym(highlight_column_name)) %>%
#             dplyr::summarise(dim1 = median(dim1), dim2 = median(dim2))
#
#         plot3 <- plot2 +
#             geom_text(data = tmp_labels,
#                       aes(x = dim1, y = dim2, label = !!rlang::sym(highlight_column_name),
#                           color = !!rlang::sym(highlight_column_name)),
#                       size = 4, show.legend = FALSE) +
#             scale_color_manual(values = color_mapping) +
#             scale_x_continuous(limits = xlim) +
#             scale_y_continuous(limits = ylim) +
#             ggtitle(title)
#     }
#
#     plot_list = list(plot1 = plot1, plot2 = plot2, plot3 = plot3)
#     return(plot_list)
# }

#' Generate UMAP Plot with Highlighted Cell Density
#'
#' This function creates a UMAP plot, highlighting a specific group of cells
#' by showing their density on top of a grey background of all cells.
#'
#' @param seurat_object A Seurat object. Default is `so`.
#' @param umap_to_plot Character string, name of the UMAP embedding to use. Default is "mde_incremental".
#' @param group.by Character string, metadata column to define the highlighted group (should be logical TRUE/FALSE). Default is `sprintf("is_%s", cl)`.
#' @param split.by Character string, metadata column to split plots by using `facet_wrap`. Default is NULL.
#' @param raster Logical, whether to use `geom_scattermore` for highlighted points (TRUE) or `geom_point` (FALSE). Default is TRUE.
#' @param highlight_size Numeric, point size for highlighted cells if `raster` is FALSE. Default is 0.5.
#' @param highlight_alpha Numeric, alpha transparency for highlighted cells if `raster` is FALSE. Default is 0.5.
#' @param cols A character vector of colors to define the density color ramp. Default is `rev(rainbow(10, end = 4/6))`.
#' @import ggplot2
#' @importFrom scattermore geom_scattermore
#' @importFrom grDevices densCols colorRampPalette
#' @importFrom grDevices rainbow
#' @return A ggplot object.
#' @export
MyDimPlotHighlightDensity = function(seurat_object = so, umap_to_plot = "mde_incremental",group.by = sprintf("is_%s", cl), split.by = NULL, raster = T, highlight_size = 0.5, highlight_alpha = 0.5, cols = rev(rainbow(10, end = 4/6))) {
    requireNamespace("scattermore", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)

    df = data.frame(feature1 = seurat_object[[umap_to_plot]]@cell.embeddings[, 1], feature2 = seurat_object[[umap_to_plot]]@cell.embeddings[, 2], group.by = seurat_object@meta.data[,group.by])
    if (!is.null(split.by)) {
        df$split.by = seurat_object@meta.data[,split.by]
    }

    df2 = df[df$group.by == T,]
    df2$density = densCols(df2$feature1, df2$feature2, colramp = colorRampPalette(cols), nbin = 500)

    p1 = ggplot(df) + scattermore::geom_scattermore(aes(feature1, feature2), color = "grey", pixels = c(512, 512))
    if (raster == T ){
        p2 = scattermore::geom_scattermore(data = df2, aes(feature1, feature2, color = density),  pixels = c(216, 216))
    } else {
        p2 = geom_point(data = df2, aes(feature1, feature2, color = density), size = highlight_size, alpha = highlight_alpha)
    }
    if(is.null(split.by)) {
        p3 = p1 + p2  +scale_colour_identity()+  theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    } else {
        p3 = p1 + p2 +scale_colour_identity()+  theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~split.by)
    }
    return(p3)
}


#' Convert Seurat V5 Object to Seurat V3 Object Structure
#'
#' This function takes a Seurat V5 object and converts it into a structure
#' resembling a Seurat V3 object, extracting specified assays (e.g., RNA and ADT).
#'
#' @param so A Seurat V5 object.
#' @param assay1 Character string, name of the first assay to extract (e.g., "RNA"). Default is "RNA".
#' @param assay2 Character string, name of the second assay to extract (e.g., "ADT"). Default is "ADT".
#' @import Seurat
#' @return A Seurat V3-like object.
#' @export
ConvertS5toS3 = function(so, assay1 = "RNA", assay2 = "ADT") {
    requireNamespace("Seurat", quietly = TRUE)
    so.v3 = Seurat::CreateAssayObject(counts = so[[assay1]]$counts)
    so.v3 = Seurat::CreateSeuratObject(so.v3)
    so.v3[[assay2]] = Seurat::CreateAssayObject(counts = so[[assay2]]$counts) #samples@assays$ADT$counts
    #print(all(colnames(so.v3@assays$RNA$counts) == colnames(so.v3@assays$ADT$counts)))
    #all(colnames(so.v3) == colnames(so))
    so.v3@meta.data = so@meta.data
    return(so.v3)
}

#' Generate Plots After Seurat Integration
#'
#' This function generates a series of UMAP plots often used after Seurat integration,
#' including plots colored by a split variable, split by a variable, gene expression, and clusters.
#'
#' @param seurat_object A Seurat object. Default is `so`.
#' @param dim1 Numeric vector, UMAP dimension 1 embeddings. Default is `so@reductions$umap@cell.embeddings[,1]`.
#' @param dim2 Numeric vector, UMAP dimension 2 embeddings. Default is `so@reductions$umap@cell.embeddings[,2]`.
#' @param split_by Character string, metadata column to split plots by (faceting). Default is "IGT".
#' @param genes Character vector, genes to plot expression for. Default is `c("Foxp3", "Il2ra")`.
#' @param cluster_key Character string, key for cluster information in metadata. Default is "ClusterSCVI_Res".
#' @param mypal A color palette to use. Default is `mypal`.
#' @import ggplot2
#' @importFrom ggrastr geom_point_rast
#' @return Prints multiple ggplot objects.
#' @export
PlotsAfterIntegration = function (seurat_object = so, dim1 = seurat_object@reductions$umap@cell.embeddings[,1], dim2 = seurat_object@reductions$umap@cell.embeddings[,2], split_by = "IGT", genes = c("Foxp3", "Il2ra"), cluster_key = "ClusterSCVI_Res", mypal = NULL) {

    so = seurat_object
    so@meta.data[,"split_by"] = so@meta.data[,split_by]

    alpha = 0.5
    # Assuming `mypal` and `sample_name` are available or passed
    # This part might need adjustment based on how mypal and sample_name_colors are handled globally
    # For now, making mypal an argument.
    if (is.null(mypal)) {
        mypal = RColorBrewer::brewer.pal(length(unique(so@meta.data[,split_by])), "Set3") # A default palette
    }
    sample_name_colors = mypal[seq_len(length(unique(so@meta.data[,split_by])))]
    names(sample_name_colors) = levels(so$sample_name) # This assumes 'sample_name' exists and is a factor.

    message("Plot 1: UMAP")
    p1 = ggplot(data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)) + ggrastr::geom_point_rast(aes(dim1, dim2, color = .data[[split_by]]), alpha = I(alpha), raster.dpi = 100) + theme_bw() + scale_color_manual(values = mypal)
    print(p1)

    message("Plot 2: UMAP split")
    tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)
    bkgrd = data.frame(dim1 = dim1, dim2 = dim2)

    p = ggplot(bkgrd) + ggrastr::geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, alpha = 0.2, raster.dpi = 50)
    p2 = geom_point(data = tmp, aes(dim1, dim2, color = .data[[split_by]]), size = 1,  alpha = alpha) # Using .data for NSE

    q = p + p2 + scale_color_manual(values = sample_name_colors)  + theme_bw() + facet_wrap(~ .data[[split_by]])
    print(q)

    message("Plot 3: UMAP genes")

    for (g in genes) {
        print(g)
        # Ensure gene exists before trying to subset
        if (g %in% rownames(so@assays$RNA$counts)) {
            tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2, size = so@assays$RNA$counts[g,])
        } else {
            message(paste("Gene", g, "not found in RNA counts. Skipping."))
            next
        }
        bkgrd = data.frame(dim1 = dim1, dim2 = dim2)

        p = ggplot(bkgrd) + ggrastr::geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, raster.dpi = 50)
        p2 = geom_point(data = tmp, aes(dim1, dim2, color = size > 0, size = size,  alpha = size > 0))

        p3 = p + p2  + scale_color_manual(values = c("black", "red")) + scale_alpha_manual(values = c(0.5,1))  + theme_bw() + facet_wrap(~.data[[split_by]]) + ggtitle(g)
        print(p3)
    }


    message("Plot 4: UMAP Clusters")
    # This loop needs careful handling of `cluster_key` as it's a single string,
    # but the original code iterates over columns grepping "cluster_key".
    # Assuming 'cluster_key' is the exact column name or a pattern to find related cluster columns.
    # I'll modify to iterate over actual metadata columns that match the cluster_key pattern.
    cluster_cols = colnames(so@meta.data)[grepl(cluster_key, colnames(so@meta.data))]
    if (length(cluster_cols) == 0) {
        message(paste("No cluster columns found matching pattern:", cluster_key))
    }
    for (i in cluster_cols) {
        print(i)
        p = ggplot(data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)) +
            ggrastr::geom_point_rast(aes(dim1, dim2), color = so@meta.data[,i], alpha = I(alpha), raster.dpi = 100) +
            ggtitle(i) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        if (!is.null(mypal)) {
            p = p + scale_color_manual(values = mypal)
        }
        print(p + facet_wrap(~.data[[split_by]])) # Using .data for NSE
    }

}


#' Add Latent Data (e.g., UMAP, MDE, PCA) to a Seurat Object
#'
#' This function reads latent space embeddings (e.g., from totalVI, MDE, PCA)
#' from a CSV file and adds them as a `DimReducObject` to a Seurat object.
#' Optionally, it can also calculate UMAP on these new embeddings.
#'
#' @param so A Seurat object.
#' @param latent_file Character string, path to the CSV file containing the latent embeddings.
#' @param prefix Character string, prefix for the new `DimReduc` object name (e.g., "totalvi_mde").
#' @param calculate_umap Logical, whether to run UMAP on the newly added latent space. Default is FALSE.
#' @import Seurat
#' @return The modified Seurat object with the new `DimReduc` object.
#' @export
AddLatentData = function(so, latent_file, prefix, calculate_umap = F) {
    requireNamespace("Seurat", quietly = TRUE)

    totalvi = read.csv(file = latent_file, header = T, sep = ",", row.names = 1)
    colnames(totalvi) = as.character(seq(1:ncol(totalvi)))
    totalvi = as.matrix(totalvi)

    # Check if cell names in latent_file match Seurat object colnames
    # This pattern check seems specific to certain naming conventions, consider making it more robust
    # or relying purely on matching colnames(so) and rownames(totalvi)
    # pattern = "^IGT\\d{1,2}\\.[A-Z]{16}$"
    # matches = grepl(pattern, rownames(totalvi))
    # if (!all(matches)) {
    #     print("updating cell names in csv file, assuming same order as seurat object" )
    #     rownames(totalvi) = colnames(so) # This assumes perfect order, which is risky
    # }

    # Robust cell matching: intersect and reorder
    common_cells <- intersect(colnames(so), rownames(totalvi))
    if (length(common_cells) == 0) {
        warning("No common cells found between Seurat object and latent file. Returning original Seurat object.")
        return(so)
    }

    # Ensure ordering is consistent and only common cells are used
    totalvi <- totalvi[common_cells, ]
    so <- so[, common_cells] # Subset Seurat object to common cells

    if (all(colnames(so) == rownames(totalvi))) {
        print("All cells in Seurat object found and matched in the latent_file.")
        so[[prefix]] = Seurat::CreateDimReducObject(embeddings = totalvi, key = prefix, assay = "RNA")
        if (calculate_umap == T) {
            so = Seurat::RunUMAP(so, dims = 1:ncol(totalvi), reduction = prefix, n.components = 2, reduction.name = sprintf("umap_%s", prefix))
        }
        return(so)
    } else {
        # This branch should ideally not be reached if common_cells logic is applied,
        # but kept for safety/debugging.
        print("Mismatch in cell names after processing. Returning original Seurat object.")
        return(so)
    }
}

