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
#' @export
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
#' @export
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

#' Perform differential gene expression (DGE) analysis with flashier results.
#'
#' This function calculates and visualizes differential expression for both gene programs (factors) and individual genes
#' between two specified groups of cells. It returns two MA plots and data tables of the results.
#'
#' @param F1 A gene-by-factor matrix from a flashier fit.
#' @param L1 A cell-by-factor matrix from a flashier fit.
#' @param group1 A character vector of cell IDs belonging to the first group.
#' @param group2 A character vector of cell IDs belonging to the second group.
#' @param title_plot A character string for the title of the plots.
#' @return A list containing:
#'   \itemize{
#'     \item \code{p1}: A ggplot object of the MA plot for factors.
#'     \item \code{p2}: A ggplot object of the MA plot for genes.
#'     \item \code{diff_factors}: A data frame of the differential expression results for factors.
#'     \item \code{diff_genes}: A data frame of the differential expression results for genes.
#'   }
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes geom_point xlim theme_minimal labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom scattermore geom_scattermore
#' @export
#'
#' @examples
#' \dontrun{
#' # Note: The following example assumes you have a 'fit.Rds' file
#' # and a 'metadata_cells.txt' file in your working directory.
#'
#' # Get F1 and L1: Recommended for large datasets: Import from a flashier fit object
#' fit = readRDS("fit.Rds")
#' res = ldf(fit, type = "i")
#' L1 = with(res, L %*% diag(D) )
#' L1 = L1[mdata_orig$cellID,]
#' colnames(L1) = sprintf("F%s", 1:ncol(L1))
#' F1 = with(ldf(fit, type = "i"), F)
#' colnames(F1) = paste("F", 1:ncol(F1), sep = "")
#'
#' # We'll create some dummy data for a runnable example.
#' # In your actual use case, replace this with your real data loading code.
#' F1 <- matrix(runif(1000 * 20), nrow = 1000, ncol = 20)
#' L1 <- matrix(runif(500 * 20), nrow = 500, ncol = 20)
#' rownames(L1) <- mdata_orig$cellID[1:500]
#' rownames(F1) <- paste0("gene", 1:1000)
#' colnames(L1) <- paste0("F", 1:20)
#' colnames(F1) <- paste0("F", 1:20)
#'
#' # Define your groups based on metadata
#' mdata_orig <- read.table("metadata_cells.txt", header = TRUE, sep = "\t")
#' # Assuming mdata_orig has columns "cellID", "annotation_level1", etc.
#' group1 <- mdata_orig %>%
#'   dplyr::filter(annotation_level1 %in% c("Treg") &
#'                 annotation_level2 %in% c("Treg_cl8") &
#'                 organ_simplified == "colon" &
#'                 condition_broad == "healthy") %>%
#'   dplyr::pull(cellID)
#'
#' group2 <- mdata_orig %>%
#'   dplyr::filter(annotation_level1 %in% c("Treg") &
#'                 annotation_level2 %in% c("Treg_cl2") &
#'                 organ_simplified == "colon" &
#'                 condition_broad == "healthy") %>%
#'   dplyr::pull(cellID)
#'
#' # Call the function
#' results <- FlashierDGE(F1, L1, group1, group2,
#'                        title_plot = "Treg_cl8 vs Treg_cl2 in colon healthy")
#'
#' # Access and display the plots and data
#' results$p1    # MA plot of factors
#' results$p2    # MA plot of genes
#'
#' results$diff_factors # table for fold changes of gene programs
#'
#' # Access and sort the gene fold changes
#' results$diff_genes %>%
#'   dplyr::arrange(dplyr::desc(log2FC))
#' }

FlashierDGE = function(F1, L1, group1, group2, title_plot = "") {

    # Calculate the loadings
    loadings_group1 = colMeans(L1[group1,])
    loadings_group2 = colMeans(L1[group2,])
    loadings_groups = colMeans(L1[c(group1, group2),])

    # Mean genes for each group
    mean_genes_group1 = F1 %*% loadings_group1
    mean_genes_group2 = F1 %*% loadings_group2
    mean_genes = F1 %*% colMeans(L1[c(group1, group2),])

    # Fold Change between the two groups
    fc_loadings = loadings_group1 - loadings_group2
    fc_genes = F1 %*% fc_loadings %>% as.data.frame()

    # Create vplot for fold changes (MA plot of factors)
    vplot = data.frame(SYMBOL = names(fc_loadings), log2FC = fc_loadings / log(2), AveExpr = loadings_groups)

    # Maximum FC for symmetrical x-axis range
    max_fc = ceiling(max(abs(vplot$log2FC)))
    top_genes = vplot %>%
        dplyr::arrange(dplyr::desc(abs(log2FC))) %>%
        utils::head(50)

    # MA Plot (Factors)
    p1 = ggplot2::ggplot(data = vplot) +
        ggplot2::geom_point(ggplot2::aes(x = log2FC, y = AveExpr), colour = "black", alpha = I(1), size = I(1)) +
        ggplot2::xlim(-max_fc, max_fc) +
        ggrepel::geom_text_repel(data = top_genes, ggplot2::aes(x = log2FC, y = AveExpr, label = SYMBOL),
                                 size = 3, color = "red", box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps  = 20) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Fold Change (log2)",
            y = "Average Expression",
            title = title_plot
        )

    # Create vplot for gene expression (MA plot of genes)
    vplot_genes = data.frame(SYMBOL = rownames(fc_genes), log2FC = fc_genes[,1] / log(2), AveExpr = mean_genes[,1])

    # Maximum FC for symmetrical x-axis range (for genes)
    max_fc_genes = ceiling(max(abs(vplot_genes$log2FC)))
    top_genes_genes <- vplot_genes %>%
        dplyr::arrange(dplyr::desc(abs(log2FC))) %>%
        utils::head(50)

    # MA Plot (Genes)
    p2 <- ggplot2::ggplot(data = vplot_genes) +
        scattermore::geom_scattermore(ggplot2::aes(x = log2FC, y = AveExpr), colour = "black", alpha = I(1), size = I(1), pixels = c(512, 512)) +
        ggplot2::xlim(-max_fc_genes, max_fc_genes) +
        ggrepel::geom_text_repel(data = top_genes_genes, ggplot2::aes(x = log2FC, y = AveExpr, label = SYMBOL),
                                 size = 3, color = "red", box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps  = 20) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Fold Change (log2)",
            y = "Average Expression",
            title = title_plot
        )

    # Return the plots and vplot tables
    list(p1 = p1, p2 = p2, diff_factors = vplot, diff_genes = vplot_genes)
}
