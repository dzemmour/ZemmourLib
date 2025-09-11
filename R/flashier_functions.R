# flashier_functions.R

#' Creates a structured bar plot for topic expression across cell groups.
#'
#' @param loadings_matrix A matrix or data frame where rows are cells and columns are topics. Rownames must be cell IDs.
#' @param grouping_vector A vector with group assignments for each cell, in the same order as the rows in loadings_matrix.
#' @param topics A character vector of the topic column names to include in the plot.
#' @param colors A named vector of colors for the topics.
#' @param y_limit A numeric vector of length 2 for the y-axis limits, e.g., c(0, 4).
#' @param bar_width The width of the bars in the plot.
#' @return A ggplot object.
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate
#' @importFrom reshape2 melt
#' @import ggplot2
MyStructurePlot <- function(loadings_matrix, grouping_vector, topics, colors = mypal_topics, y_limit = c(0, 4), bar_width = 1) {

    # --- 1. Combine inputs into a single data frame ---
    plot_data <- tibble::rownames_to_column(as.data.frame(loadings_matrix), var = "cellID") %>%
        dplyr::mutate(grouping = grouping_vector)

    # --- 2. Reshape and prepare data for plotting ---
    plot_data <- reshape2::melt(plot_data, id.vars = c("cellID", "grouping"))
    plot_data <- plot_data[plot_data$variable %in% topics, ]

    # --- 3. Create the Plot ---
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cellID, y = value, fill = variable)) +
        ggplot2::geom_col(width = bar_width) +
        ggplot2::facet_grid(. ~ grouping,
                            scales = "free_x",
                            space = "free_x") +
        ggplot2::scale_fill_manual(values = colors) +
        ggplot2::coord_cartesian(ylim = y_limit) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            strip.placement = "bottom",
            strip.background = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(face = "bold", size = 10),
            panel.spacing.x = ggplot2::unit(0.5, "cm"),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(size = 10, face = "bold")
        ) +
        NoGrid() +
        ggplot2::labs(
            title = "",
            y = "Topic Expression",
            fill = "Topic"
        )

    return(p)
}


#' Create a tile plot of top gene expression for specified topics.
#'
#' @param gene_factor_matrix A matrix with genes as rows and topics as columns. Rownames must be gene names.
#' @param topics A character vector of the topic column names to plot.
#' @param n_genes_per_topic The number of top genes to select from each topic.
#' @param title A character vector of the title of the plot.
#' @param alpha_range A numeric vector of the range of values of expression to apply alpha. Default NULL is determined automatically.
#' @return A ggplot object.
#' @importFrom reshape2 melt
#' @import ggplot2
MyGeneTilePlot <- function(gene_factor_matrix,
                           topics = c("F68", "F27"),
                           n_genes_per_topic = 15,
                           alpha_range = NULL,
                           title = "Dot Plot (Squares) of Gene Expression") {

    # Safety checks
    if (is.null(colnames(gene_factor_matrix))) stop("gene_factor_matrix must have column names (topic names).")
    if (is.null(rownames(gene_factor_matrix))) stop("gene_factor_matrix must have rownames (gene names).")
    missing_topics <- setdiff(topics, colnames(gene_factor_matrix))
    if (length(missing_topics)) stop(sprintf("These topics are not in gene_factor_matrix: %s", paste(missing_topics, collapse = ", ")))

    # Get top genes per topic (indices), then unique gene names
    idx_mat <- apply(gene_factor_matrix[, topics, drop = FALSE], 2, function(col)
        head(order(col, decreasing = TRUE), n_genes_per_topic)
    )
    idx_vec <- unique(as.vector(idx_mat))
    top_genes <- rownames(gene_factor_matrix)[idx_vec]

    # Build plotting data
    F_sub <- gene_factor_matrix[top_genes, topics, drop = FALSE]
    plot_data <- reshape2::melt(F_sub, varnames = c("Gene", "Factor"), value.name = "Expression")
    plot_data$Factor <- factor(plot_data$Factor, levels = topics)

    # Ensure genes appear ordered by first topic (highest to lowest)
    ordering_topic <- topics[1]
    gene_order <- names(sort(gene_factor_matrix[top_genes, ordering_topic], decreasing = TRUE))
    plot_data$Gene <- factor(plot_data$Gene, levels = rev(gene_order))

    # Plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Factor, y = Gene, fill = Expression, alpha = Expression)) +
        ggplot2::geom_tile(width = 0.9, height = 0.9) +
        ggplot2::scale_fill_gradient(low = "white", high = "brown", limits = c(0, NA), na.value = "white")
    if (!is.null(alpha_range)) {
        p <- p + ggplot2::scale_alpha_continuous(limits = alpha_range)
    }

    p <- p +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank()
        ) +
        ggplot2::coord_fixed(ratio = 1) +
        NoGrid() +
        ggplot2::labs(title = title, x = "", y = "", fill = "Expression (Color)", alpha = "Expression (alpha)")

    return(p)
}
