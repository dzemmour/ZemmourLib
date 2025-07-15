



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
