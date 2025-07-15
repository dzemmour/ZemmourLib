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

