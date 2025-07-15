#' Extract Numeric Values from a String
#'
#' This function takes a string and extracts all numeric characters,
#' then converts them to a numeric type.
#'
#' @param x A character string.
#' @return A numeric value extracted from the string.
#' @examples
#' extract_numeric("cell_123")
#' extract_numeric("ABC-456DEF")
extract_numeric = function(x) {
    as.numeric(gsub("\\D", "", x))
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
MyDimPlotHighlight <- function(seurat_object = so,
                               umap_to_plot = "mde_totalvi_20241201_gdT_rmIGTsample",
                               cells_to_highlight = names(which(seurat_object$nonconv_tcr_recog == T)),
                               highlight_column_name = "nonconv_tcr_recog",
                               title = "",
                               labelclusters = T,
                               pixels = c(512, 512),
                               highlight_size = 1,
                               highlight_alpha = 1,
                               xlim = NULL,
                               ylim = NULL,
                               mycols = NULL,
                               print_plot1 = TRUE,
                               print_plot2 = TRUE) {

    requireNamespace("scattermore", quietly = TRUE) # Ensures scattermore is available
    requireNamespace("ggplot2", quietly = TRUE) # Ensures ggplot2 is available
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("rlang", quietly = TRUE)
    requireNamespace("Seurat", quietly = TRUE)

    # Handle default mycols if not provided
    if (is.null(mycols)) {
        # Using a simple default if mypal is not defined globally or passed
        mycols <- c("TRUE" = "red", "FALSE" = "grey")
    }

    missing_names = setdiff(unique(seurat_object@meta.data[[highlight_column_name]]), names(mycols))

    if (length(missing_names) > 0) {
        if (any(is.na(names(mycols)))) {
            new_colors = mycols[which(is.na(names(mycols)))][seq_len(length(missing_names))]
        } else {
            # Fallback to a simple palette if no NA names, just pick from existing
            new_colors = mycols[seq_len(length(missing_names))]
        }
        names(new_colors) <- missing_names
        color_mapping <- c(mycols[!is.na(names(mycols))], new_colors)
    } else {
        color_mapping = mycols
    }

    # UMAP embeddings
    dim1 <- seurat_object[[umap_to_plot]]@cell.embeddings[, 1]
    dim2 <- seurat_object[[umap_to_plot]]@cell.embeddings[, 2]
    tmp <- data.frame(seurat_object@meta.data, dim1 = dim1, dim2 = dim2)

    # Subset for highlighted cells
    tmp2 <- tmp[cells_to_highlight, ]

    # Filter unique cluster labels present in cells_to_highlight
    highlight_values <- unique(tmp2[[highlight_column_name]])

    # Background plot
    bkrg <- ggplot(tmp) +
        scattermore::geom_scattermore(aes(dim1, dim2), color = "grey50", pointsize = 0.5, alpha = 0.5, pixels = pixels)

    # Highlighted points
    p2 <- geom_point(data = tmp2,
                     aes(dim1, dim2, color = !!rlang::sym(highlight_column_name)),
                     size = highlight_size,
                     alpha = highlight_alpha)

    # Plots with consistent colors
    plot1 <- bkrg + p2 +
        scale_color_manual(values = color_mapping) +
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(limits = ylim) +
        ggtitle(title) +
        theme_minimal()

    plot2 <- bkrg + p2 +
        scale_color_manual(values = color_mapping) +
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(limits = ylim) +
        theme_void() + Seurat::NoLegend()

    # Add cluster labels dynamically for only the present values
    plot3 <- NULL
    if (labelclusters) {
        tmp_labels <- tmp2 %>%
            dplyr::filter(!!rlang::sym(highlight_column_name) %in% highlight_values) %>%
            dplyr::group_by(!!rlang::sym(highlight_column_name)) %>%
            dplyr::summarise(dim1 = median(dim1), dim2 = median(dim2))

        plot3 <- plot2 +
            geom_text(data = tmp_labels,
                      aes(x = dim1, y = dim2, label = !!rlang::sym(highlight_column_name),
                          color = !!rlang::sym(highlight_column_name)),
                      size = 4, show.legend = FALSE) +
            scale_color_manual(values = color_mapping) +
            scale_x_continuous(limits = xlim) +
            scale_y_continuous(limits = ylim) +
            ggtitle(title)
    }

    plot_list = list(plot1 = plot1, plot2 = plot2, plot3 = plot3)
    return(plot_list)
}

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
#' @importFrom graphics rainbow
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

#' Ensure Directory Exists
#'
#' This helper function checks if a directory exists and creates it if it doesn't.
#'
#' @param path Character string, the path to the directory.
#' @return No explicit return value, invisibly returns TRUE if directory exists or was created, FALSE otherwise.
#' @export
ensure_directory <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
        message(paste("Directory created:", path))
    } else {
        message(paste("Directory already exists:", path))
    }
    invisible(dir.exists(path))
}

#' Run Limma-Trend Differential Gene Expression with Confounding Variables
#'
#' This function performs differential gene expression analysis using `limma-trend`
#' given TMM-normalized expression data, grouping variables, confounding variables,
#' and specified contrasts.
#'
#' @param tmm An `EList` object or a matrix of TMM-normalized expression data (log-transformed).
#' @param group A factor or vector indicating the primary grouping for differential expression.
#' @param confoundings A list of vectors or factors representing confounding variables.
#' @param formula.mod.matrix Character string, the R formula for the design matrix (e.g., "~ group + confoundings[[1]]").
#' @param contrasts A named list of character strings defining the contrasts to be made.
#' @param gene_symbol A vector of gene symbols corresponding to the rows of `tmm`.
#' @import limma
#' @import edgeR
#' @return A named list of `topTable` results for each contrast.
#' @export
run_limmatrend_contrasts_counfoundings = function(tmm = tmm, group = group, confoundings = confoundings, formula.mod.matrix = formula.mod.matrix, contrasts =contrasts, gene_symbol= gene_symbol) { #count = count, dge = dge,
    suppressPackageStartupMessages(library(limma))
    suppressPackageStartupMessages(library(edgeR))
    message("limmatrend")

    design = model.matrix(formula(formula.mod.matrix))
    colnames(design) = gsub(paste(c("group", sprintf("confoundings\\[\\[%s\\]\\]", seq_along(confoundings))), collapse = "|"), "", colnames(design)) # Adjusted to seq_along(confoundings)
    cont.matrix = makeContrasts( contrasts = contrasts, levels = design)
    colnames(cont.matrix) = names(contrasts)

    message("lmFit")
    fit = lmFit(tmm, design = design)
    fit2 = contrasts.fit(fit, cont.matrix)
    fit2 = eBayes(fit2, trend = TRUE, robust = TRUE)

    tt_list = list()
    for (i in colnames(cont.matrix)) {
        print(i)
        tt_list[[i]] = topTable(fit2,coef = i ,n = Inf, adjust.method = "BH", sort.by = "none")
        tt_list[[i]]$SYMBOL = gene_symbol
    }

    return(tt_list)
}

#' Get Groups for Differential Expression and Check Balance
#'
#' This function identifies cell IDs for two specified groups based on filter
#' expressions and performs chi-squared tests to check for imbalances in other
#' metadata variables between these groups.
#'
#' @param metadata Data frame, containing cell metadata.
#' @param group1_filter An unquoted expression to filter for group 1 (e.g., `condition == "control"`).
#' @param group2_filter An unquoted expression to filter for group 2 (e.g., `condition == "treated"`).
#' @param id_column Character string, the name of the column containing unique identifiers (e.g., "cell_id").
#' @import dplyr
#' @import rlang
#' @return A list containing `group1_ids` (character vector), `group2_ids` (character vector),
#'         and `results` (named list of p-values from chi-squared tests).
#' @export
GetGroups <- function(metadata, group1_filter, group2_filter, id_column) {
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("rlang", quietly = TRUE)

    # Evaluate the filter expressions within the context of metadata
    group1 <- metadata %>% dplyr::filter(!! rlang::enquo(group1_filter)) # Use enquo and !!
    group2 <- metadata %>% dplyr::filter(!! rlang::enquo(group2_filter)) # Use enquo and !!

    # Pull the specified identifier column
    group1_ids <- group1 %>% dplyr::pull(!! rlang::sym(id_column))
    group2_ids <- group2 %>% dplyr::pull(!! rlang::sym(id_column))

    # Combine groups and add group identifiers
    group1$group <- "group1"
    group2$group <- "group2"
    combined_groups <- rbind(group1, group2)

    # Function to perform Chi-squared test
    test_balanced <- function(var) {
        obs <- table(combined_groups$group, combined_groups[[var]])
        test <- tryCatch({
            chisq.test(obs)
        }, error = function(e) NA)  # Return NA on error
        if (!is.na(test$p.value)) {
            return(test$p.value)
        } else {
            return(1)  # Return non-significant p-value if test could not be performed
        }
    }

    # Apply test to each relevant column and suppress warnings
    results <- suppressWarnings(
        lapply(colnames(metadata)[!colnames(metadata) %in% c("sample_id", id_column)], test_balanced)
    )
    names(results) <- colnames(metadata)[!colnames(metadata) %in% c("sample_id", id_column)]

    # Check for significant results
    if (any(na.omit(unlist(results) < 0.05))) {
        message("WARNING: IMBALANCED GROUPS FOR:")
        print(names(results)[which(unlist(results) < 0.05)])
    }

    return(list(group1_ids = group1_ids, group2_ids = group2_ids, results = results))
}

#' Create Comparison Name from Filter Expressions
#'
#' This function generates a descriptive comparison name from two filter expressions,
#' typically used for naming differential expression contrasts. It extracts quoted strings from the expressions.
#'
#' @param group1_filter An unquoted expression for group 1 (e.g., `condition == "control"`).
#' @param group2_filter An unquoted expression for group 2 (e.g., `condition == "treated"`).
#' @import rlang
#' @import stringr
#' @return A character string representing the comparison name (e.g., "control_vs_treated").
#' @export
CreateComparisonName <- function(group1_filter, group2_filter) {
    requireNamespace("rlang", quietly = TRUE)
    requireNamespace("stringr", quietly = TRUE)
    # Convert the expressions to strings
    group1_name <- rlang::expr_text(rlang::enquo(group1_filter)) # Use enquo
    group2_name <- rlang::expr_text(rlang::enquo(group2_filter)) # Use enquo

    # Function to extract quoted substrings
    extract_quoted_parts <- function(name) {
        stringr::str_match_all(name, '"([^"]+)"')[[1]][,2]
    }

    # Extract quoted parts for each group
    group1_parts <- extract_quoted_parts(group1_name)
    group2_parts <- extract_quoted_parts(group2_name)

    # Concatenate the extracted parts
    group1_clean <- paste0(group1_parts, collapse = "")
    group2_clean <- paste0(group2_parts, collapse = "")

    # Create the comparison name
    comparison_name <- paste0(group1_clean, "_vs_", group2_clean)

    return(comparison_name)
}

#' Collapse Limma-Trend Differential Expression Results
#'
#' This function takes a list of `topTable` results from `limma-trend`
#' and collapses them into a single data frame, merging log fold changes,
#' average expression, p-values, and adjusted p-values.
#'
#' @param l A named list of data frames, where each data frame is a `topTable` result.
#' @return A data frame with collapsed differential expression results, including columns
#'         for logFC, AveExpr, P.Value, and adj.P.Val for each contrast.
#' @export
CollapseDiff_limmatrend = function(l) {
    l = lapply(l, function(x) x = data.frame(x, rownames = rownames(x)))

    fc = data.frame(l[[1]][,c("rownames", "logFC")])
    for (i in 2:length(l)) { fc = merge(x = fc, y = l[[i]][,c("rownames", "logFC")], by.x = "rownames", by.y = "rownames", all = T) }

    mean = data.frame(l[[1]][,c("rownames", "AveExpr")])
    for (i in 2:length(l)) { mean = merge(x = mean, y = l[[i]][,c("rownames", "AveExpr")], by.x = "rownames", by.y = "rownames", all = T) }

    pv = data.frame(l[[1]][,c("rownames", "P.Value")])
    for (i in 2:length(l)) { pv = merge(x = pv, y = l[[i]][,c("rownames", "P.Value")], by.x = "rownames", by.y = "rownames", all = T) }

    qcl = data.frame(l[[1]][,c("rownames", "adj.P.Val")])
    for (i in 2:length(l)) { qcl = merge(x = qcl, y = l[[i]][,c("rownames", "adj.P.Val")], by.x = "rownames", by.y = "rownames", all = T) }

    rownames(fc) = fc[,1]
    rownames(mean) = fc[,1]
    rownames(pv) = pv[,1]
    rownames(qcl) = qcl[,1]
    fc = fc[,-1]
    mean = mean[,-1]
    pv = pv[,-1]
    qcl = qcl[,-1]
    colnames(fc) = paste("logFC", names(l), sep = "_")
    colnames(mean) = paste("AveExpr", names(l), sep = "_")
    colnames(pv) = paste("P.Value", names(l), sep = "_")
    colnames(qcl) = paste("adj.P.Val", names(l), sep = "_")
    # All checks: all(rownames(fc) == rownames(pv)), all(rownames(fc) == rownames(mean)), all(rownames(fc) == rownames(qcl))
    fc[is.na(fc)] = 0
    pv[is.na(pv)] = 1
    qcl[is.na(qcl)] = 1
    diff = data.frame(fc, mean, pv, qcl)
    #which(is.na(diff)) # No need to print this.

    return(diff)

}

#' Create a V-Plot (Volcano Plot-like with Fold Change on X)
#'
#' This function generates a V-plot to visualize gene expression changes,
#' with fold change on the x-axis (log-transformed) and p-value on the y-axis (log-transformed).
#'
#' @param vplot Data frame with columns `fc` (fold change) and `pval` (p-value).
#' @param xlab Character string, label for the x-axis. Default is "FC".
#' @param xlimits Numeric vector of length 2, x-axis limits (e.g., `c(0.1, 10)`).
#' @param ylimits Numeric vector of length 2, y-axis limits (e.g., `c(10^-300,1)`).
#' @import ggplot2
#' @importFrom ggrastr geom_point_rast
#' @importFrom scales log_trans trans_breaks trans_format math_format
#' @return A ggplot object.
#' @export
Vplot = function(vplot, xlab = "FC", xlimits = c(0.1, 10), ylimits = c(10^-300,1)) {
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("ggrastr", quietly = TRUE)
    requireNamespace("scales", quietly = TRUE)

    p = ggplot(data = vplot) + ggrastr::geom_point_rast(aes(x = .data$fc, y = .data$pval), colour = "grey", alpha = I(1), size = I(1), raster.dpi = 100) +
        scale_x_continuous(trans = scales::log_trans(10), limits = xlimits) +
        scale_y_continuous(trans = scales::log_trans(10), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = ylimits) +
        annotation_logticks(sides = "b") +
        xlab(xlab) +
        theme_bw() +
        ylab("p value") +
        theme(axis.text.x  = element_text(size=20,angle = 0, hjust = 0.5), axis.text.y  = element_text(size=20), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20))
    return(p)
}

#' Add Significance Information to a V-Plot
#'
#' This function takes a V-plot object and adds colored points for significant
#' up-regulated and down-regulated genes. It also performs and displays
#' chi-squared test results for the distribution of significant genes.
#'
#' @param p A ggplot object generated by `Vplot`.
#' @param vplot Data frame with columns `fc` (fold change), `pval` (p-value), and `sig` (significance category, e.g., "up", "down").
#' @param y_text Numeric, y-coordinate for the annotation text. Default is `10^-100`.
#' @param sigtoplot Character string, title for the plot, likely related to the significance filter applied.
#' @import ggplot2
#' @importFrom ggrastr geom_point_rast
#' @return A ggplot object with significance highlights and annotations.
#' @export
VplotAddSig = function(p, vplot, y_text = 10^-100, sigtoplot = "Significance") { # Added sigtoplot argument
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("ggrastr", quietly = TRUE)

    a = length(which(vplot$fc < 1))
    b = length(which(vplot$fc > 1))
    c_up = length(which(vplot$fc[vplot$sig == "up"] < 1))
    d_up = length(which(vplot$fc[vplot$sig == "up"] > 1))
    cont_table_up = rbind(total = c(round(a/(a+b)*(c_up+d_up)), round(b/(a+b)*(c_up+d_up))), geneset = c(c_up, d_up)) # Renamed to avoid conflict
    colnames(cont_table_up) = c("FC<1", "FC>1")
    c_up_test = chisq.test(cont_table_up) # Renamed to avoid conflict
    pvalsigup = c_up_test$p.value

    a = length(which(vplot$fc < 1))
    b = length(which(vplot$fc > 1))
    c_down = length(which(vplot$fc[vplot$sig == "down"] < 1))
    d_down = length(which(vplot$fc[vplot$sig == "down"] > 1))
    cont_table_down = rbind(total = c(round(a/(a+b)*(c_down+d_down)), round(b/(a+b)*(c_down+d_down))), geneset = c(c_down, d_down)) # Renamed
    colnames(cont_table_down) = c("FC<1", "FC>1")
    c_down_test = chisq.test(cont_table_down) # Renamed
    pvalsigdown = c_down_test$p.value

    q = p + ggrastr::geom_point_rast(aes(.data$fc,.data$pval, color = .data$sig), data = vplot[vplot$sig %in% c("up", "down"), ], raster.dpi = 100) +
        scale_color_manual(values = c("blue", "red")) +
        ggtitle(label = sigtoplot) + # Uses the passed sigtoplot
        annotate("text", label = sprintf("total = %s - %s \n sigup: p = %s, observed = %s - %s \n sigdown: p = %s, observed = %s - %s",
                                         a,b,round(pvalsigup, 6), c_up, d_up, round(pvalsigdown, 6), c_down, d_down), x = 1, y = y_text)
    return(q)

}

#' Create a Fold Change vs. Fold Change Plot (FCFC Plot)
#'
#' This function generates a scatter plot comparing fold changes from two different
#' conditions or analyses, with options to label specific genes.
#'
#' @param vplot Data frame with columns `fc1`, `fc2` (fold changes for two conditions), and optionally `SYMBOL` (gene symbol).
#' @param xlab Character string, label for the x-axis (first fold change). Default is "fc1".
#' @param ylab Character string, label for the y-axis (second fold change). Default is "fc2".
#' @param main Character string, plot title. Default is "".
#' @param xlimits Numeric vector of length 2, x-axis limits. Default is `c(0.1, 10)`.
#' @param ylimits Numeric vector of length 2, y-axis limits. Default is `c(0.1, 10)`.
#' @param genes_to_label Character vector, gene symbols to label on the plot. Default is NULL.
#' @import ggplot2
#' @importFrom scales log_trans
#' @importFrom ggrepel geom_text_repel
#' @return A ggplot object.
#' @export
FCFCplot = function(vplot, xlab = "fc1", ylab = "fc2", main = "",
                    xlimits = c(0.1, 10), ylimits = c(0.1, 10), genes_to_label = NULL) {
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("scales", quietly = TRUE)
    requireNamespace("ggrepel", quietly = TRUE)

    # Create the base plot
    p = ggplot(data = vplot) +
        geom_point(aes(x = .data$fc1, y = .data$fc2), colour = "black", alpha = I(1), size = I(0.5)) +
        scale_x_continuous(trans = scales::log_trans(10), breaks = c(0.125,0.5, 1, 2,5, 10),
                           labels = c(0.125,0.5, 1, 2,5, 10), limits = xlimits) +
        scale_y_continuous(trans = scales::log_trans(10), breaks = c(0.125,0.5, 1, 2,5, 10),
                           labels = c(0.125,0.5, 1, 2,5, 10), limits = ylimits) +
        geom_hline(aes(yintercept = 1), linetype="dashed", color = "brown") +
        geom_vline(aes(xintercept = 1), linetype="dashed", color = "brown") +
        annotation_logticks(sides = "bl") +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(main) +
        theme_bw() +
        theme(axis.text.x  = element_text(size=15,angle = 0, hjust = 1),
              axis.text.y  = element_text(size=15),
              legend.text=element_text(size=20),
              axis.title.x = element_text(size=20) ,
              axis.title.y = element_text(size=20))

    # Add labels for specific genes if requested
    if (!is.null(genes_to_label)) {
        vplot2 = vplot[vplot$SYMBOL %in% genes_to_label,]
        p = p + ggrepel::geom_text_repel(data = vplot2, aes(x = .data$fc1, y = .data$fc2, label = .data$SYMBOL),
                                         size = 3, color = "red", box.padding = 0.35,
                                         point.padding = 0.5, segment.color = 'grey50', max.overlaps = 20)
    }

    return(p)
}

#' Generate a Color Palette with Split Alpha Values
#'
#' This function takes an existing color palette and creates a new one
#' where each color is split into `n` shades of increasing alpha transparency.
#' This is useful for visualizing nested clustering or categories.
#'
#' @param pal Character vector of colors, the base palette. Default is `mypal[1:length(levels(cds$clustersCCA9))]`.
#' @param splitvector Numeric vector, indicating how many sub-shades each color in `pal` should have.
#'        Default calculates this based on `cds$clustersCCA9` and `cds$clustersCCA9Split`.
#' @return A character vector of new colors with varying alpha values.
#' @export
SplitColors = function(pal = NULL, splitvector = NULL) {
    if (is.null(pal)) {
        stop("Please provide a base color palette 'pal'.")
    }
    if (is.null(splitvector)) {
        stop("Please provide 'splitvector' indicating the number of sub-shades for each color.")
    }
    newpal = c()
    for (c_idx in seq_along(splitvector)){ # Changed 'c' to 'c_idx' to avoid conflict with base R 'c'
        n = splitvector[c_idx]
        if (n == 0) next # Skip if no shades needed for this color
        newcol_rgb = col2rgb(pal[c_idx], alpha = TRUE) # Get RGBA
        newcol_rgba_matrix = matrix(rep(newcol_rgb, n), ncol = n)
        # newcol_rgba_matrix[4,] = seq(50, 255, length.out = n) # This line might lead to non-monotonic alpha if n=1
        if (n > 1) {
            newcol_rgba_matrix[4,] = seq(50, 255, length.out = n)
        } else {
            newcol_rgba_matrix[4,] = 255 # Full opacity for single shades
        }

        newpal = c(newpal, rgb(red = newcol_rgba_matrix[1,]/255, blue = newcol_rgba_matrix[2,]/255, green = newcol_rgba_matrix[3,]/255, alpha = newcol_rgba_matrix[4,]/255))
    }
    return(newpal)
}

#' @rdname SplitColors
#' @export
sp <- SplitColors


#' Remove Duplicate Columns from a Data Frame
#'
#' This function identifies and removes columns from a data frame that are
#' exact duplicates of other columns.
#'
#' @param df A data frame.
#' @return A data frame with only unique columns. Prints the names of removed duplicate columns.
#' @export
removeDuplicateColumns <- function(df) {
    # Determine unique columns by comparing all pairs
    # Note: `duplicated(as.list(df))` works by comparing elements of the list,
    # which are the columns of the dataframe.
    uniqueCols <- !duplicated(as.list(df))
    # Subset dataframe to keep only unique columns
    dfUnique <- df[, uniqueCols]
    print(colnames(df)[!uniqueCols])
    return(dfUnique)
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
#' @import ggplot2
#' @importFrom scattermore geom_scattermore
#' @importFrom grDevices densCols colorRampPalette
#' @importFrom graphics rainbow
#' @return A ggplot object.
#' @export
MyFeatureScatter = function(so = so, assay = "ADT", slot = "data", feature1 = NULL, feature2 = NULL, group.by = sprintf("is_%s", cl), split.by = NULL, raster = T, cols = rev(rainbow(10, end = 4/6))) {
    requireNamespace("scattermore", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)

    if (is.null(feature1) || is.null(feature2)) {
        stop("Please provide feature1 and feature2 for MyFeatureScatter.")
    }
    # Ensure features exist in the specified assay and slot
    if (!feature1 %in% rownames(so[[assay]][[slot]])) {
        stop(paste0("Feature1 '", feature1, "' not found in assay '", assay, "', slot '", slot, "'."))
    }
    if (!feature2 %in% rownames(so[[assay]][[slot]])) {
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
        p2 = scattermore::geom_scattermore(data = df2, aes(feature1, feature2, color = density),  pixels = c(216, 216))
    } else {
        p2 = geom_point(data = df2, aes(feature1, feature2, color = density))
    }
    if(is.null(split.by)) {
        p3 = p1 + p2 + xlab(feature1) + ylab(feature2) +scale_colour_identity()+  theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    } else {
        p3 = p1 + p2 + xlab(feature1) + ylab(feature2) +scale_colour_identity()+  theme_bw() + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~split.by)
    }
    return(p3)
}


#' Calculate True Positives, True Negatives, False Positives, and False Negatives
#'
#' This helper function calculates the counts of TP, TN, FP, and FN given two
#' binary vectors (e.g., one from a gating strategy and one from a true cluster).
#'
#' @param x A logical vector (e.g., `colnames(so) %in% gate_list[["cl2"]]`).
#' @param y A logical vector (e.g., `so$is_cl2`).
#' @return A named numeric vector with counts for tp, tn, fp, fn.
#' @examples
#' tp_tn_fp_fn_example <- function() {
#'   set.seed(123)
#'   true_status <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.3, 0.7))
#'   gate_result <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.25, 0.75))
#'   TpTnFpFn(gate_result, true_status)
#' }
#' tp_tn_fp_fn_example()
#' @export
TpTnFpFn = function(x, y) {
    tmp = table(x,y)
    # Ensure all combinations are present, filling with 0 if not
    tmp_matrix <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
    tmp_matrix[rownames(tmp), colnames(tmp)] <- tmp

    tp = tmp_matrix["TRUE","TRUE"]
    tn = tmp_matrix["FALSE","FALSE"]
    fp = tmp_matrix["TRUE","FALSE"]
    fn = tmp_matrix["FALSE","TRUE"]
    return(c(tp = tp, tn = tn, fp=fp, fn=fn))
}

#' Calculate Sensitivity, Specificity, PPV, and NPV
#'
#' This function calculates sensitivity (Se), specificity (Sp), positive predictive value (PPV),
#' and negative predictive value (NPV) from true positive, true negative, false positive, and false negative counts.
#'
#' @param tp Numeric, true positives count.
#' @param tn Numeric, true negatives count.
#' @param fp Numeric, false positives count.
#' @param fn Numeric, false negatives count.
#' @return A named numeric vector with se, sp, ppv, npv.
#' @examples
#' counts <- TpTnFpFn(x = c(T,T,F,F,T), y = c(T,F,T,F,T))
#' SeSpPPVNPV(tp = counts["tp"], tn = counts["tn"], fp = counts["fp"], fn = counts["fn"])
#' @export
SeSpPPVNPV = function(tp = NULL, tn = NULL, fp = NULL, fn = NULL) {
    if (is.null(tp) || is.null(tn) || is.null(fp) || is.null(fn)) {
        stop("All arguments (tp, tn, fp, fn) must be provided.")
    }
    # Handle division by zero gracefully, return NA if denominator is zero
    se = ifelse((tp + fn) > 0, tp / (tp+fn), NA)
    sp = ifelse((fp + tn) > 0, tn / (fp+tn), NA)
    ppv = ifelse((tp + fp) > 0, tp / (tp+fp), NA)
    npv = ifelse((tn + fn) > 0, tn / (tn+fn), NA)
    res = c(se, sp, ppv, npv)
    names(res) =  c("se", "sp", "ppv", "npv")
    return(res)
}

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
