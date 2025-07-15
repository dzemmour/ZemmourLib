# ZemmourLib <img src="man/figures/logo.png" align="center" height="138" />

A collection of custom R functions primarily designed for analyzing single-cell RNA-seq, especially useful with [Seurat](https://satijalab.org/seurat/) objects. It includes utilities for UMAP visualization, differential expression analysis with `limma-trend`, and other data processing tasks.

## Installation

You can install the development version of `ZemmourLib` from GitHub with `devtools`:

```R
# install.packages("devtools")
devtools::install_github("dzemmour/ZemmourLib")
````

## Core Features

`ZemmourLib` provides a suite of functions categorized for ease of use:

### 1\. Seurat Object Utilities & Visualization

  * **`MyPlots()`**: Generate various UMAP plots for Seurat objects, including coloring by metadata, splitting by conditions, and visualizing gene expression.
  * **`MyDimPlotHighlight()`**: Create UMAP plots with specific cell subsets highlighted and optionally labeled.
  * **`MyDimPlotHighlightDensity()`**: Visualize UMAPs with the density of highlighted cell populations.
  * **`MyFeatureScatter()`**: Generate scatter plots (e.g., for ADT features) with highlighted population density, akin to flow cytometry plots.
  * **`ConvertS5toS3()`**: Convert Seurat V5 objects to a Seurat V3-like structure.
  * **`PlotsAfterIntegration()`**: A specialized plotting function for visualizing Seurat objects post-integration.
  * **`AddLatentData()`**: Add external latent space embeddings (e.g., from totalVI, MDE, PCA) into a Seurat object.

### 2\. Differential Expression Analysis (DEA) with `limma-trend`

  * **`ensure_directory()`**: A helper to ensure output directories exist.
  * **`run_limmatrend_contrasts_counfoundings()`**: Perform differential gene expression analysis using `limma-trend` with support for confounding variables and multiple contrasts.
  * **`GetGroups()`**: Identify cell IDs for defined groups and check for metadata balance between them, crucial for robust DEA.
  * **`CreateComparisonName()`**: Generate consistent comparison names from filter expressions.
  * **`CollapseDiff_limmatrend()`**: Consolidate `limma-trend` `topTable` results from multiple contrasts into a single, merged data frame.
  * **`Vplot()`**: Create a "V-plot" (volcano plot variant with log Fold Change on X-axis).
  * **`VplotAddSig()`**: Add significance highlighting and chi-squared test statistics to a V-plot.
  * **`FCFCplot()`**: Generate a Fold Change vs. Fold Change plot to compare results across two conditions or analyses, with optional gene labeling.

### 3\. Milo Differential Abundance Analysis

  * **`MyplotDAbeeswarm()`**: Visualize differential abundance results from `miloR` as a beeswarm plot.

### 4\. Utility Functions

  * **`extract_numeric()`**: Extract numeric characters from strings.
  * **`removeDuplicateColumns()`**: Remove columns that are exact duplicates in a data frame.
  * **`SplitColors()` / `sp()`**: Generate a color palette with varying alpha transparencies for sub-categories.
  * **`TpTnFpFn()`**: Calculate True Positives, True Negatives, False Positives, and False Negatives.
  * **`SeSpPPVNPV()`**: Calculate Sensitivity, Specificity, Positive Predictive Value (PPV), and Negative Predictive Value (NPV).
  * **`NoGrid()`**: A `ggplot2` theme helper to remove grid lines.
  * **`NoLegend()`**: A `ggplot2` theme helper to remove the legend.

## Usage Examples

Here are some brief examples to get you started. For more detailed usage and parameters, please refer to the function documentation (e.g., `?MyPlots` after installing the package).

```r
library(ZemmourLib)
library(Seurat)
library(ggplot2)
libary(dplyr)

# --- Example 1: Using a Plotting Function ---
# Assuming 'so' is a loaded Seurat object with UMAP embeddings and metadata
# For illustrative purposes, let's create a dummy Seurat object if you don't have one ready
if (!exists("so")) {
  # Create a small dummy Seurat object for demonstration
  data <- matrix(rpois(100*50, lambda = 5), nrow = 50, ncol = 100)
  rownames(data) <- paste0("gene", 1:50)
  colnames(data) <- paste0("cell", 1:100)
  so <- CreateSeuratObject(counts = data)
  so[["umap_unintegrated"]] <- CreateDimReducObject(
    embeddings = matrix(rnorm(100*2), ncol = 2, dimnames = list(colnames(so), c("UMAP_1", "UMAP_2"))),
    key = "UMAP_unintegrated_", assay = "RNA"
  )
  so$spleen_standard <- sample(c("WT", "KO"), 100, replace = TRUE)
  so$IGT <- sample(c("GroupA", "GroupB", "GroupC"), 100, replace = TRUE)
  so$nonconv_tcr_recog <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.2, 0.8))
  so$cell_type <- sample(LETTERS[1:5], 100, replace = TRUE)
}

# Example of MyPlots
MyPlots(so, color_by = "spleen_standard", split_by1 = "IGT", genes = c("gene1", "gene2"))

# Example of MyDimPlotHighlight
highlight_colors <- c("TRUE" = "purple", "FALSE" = "lightgrey")
plot_list <- MyDimPlotHighlight(
  seurat_object = so,
  umap_to_plot = "umap_unintegrated",
  cells_to_highlight = names(which(so$nonconv_tcr_recog == TRUE)),
  highlight_column_name = "nonconv_tcr_recog",
  title = "Cells with Non-Conventional TCR Recognition",
  mycols = highlight_colors
)
print(plot_list$plot1)


# --- Example 2: Differential Expression Analysis (Conceptual) ---
# This requires actual expression data (tmm) and metadata
# For example, assume 'tmm_data' is your TMM-normalized expression matrix
# and 'sample_metadata' is a data frame with 'sample_id', 'condition', 'batch' columns.

# Example: GetGroups
# result_groups <- GetGroups(
#   metadata = sample_metadata,
#   group1_filter = condition == "Control",
#   group2_filter = condition == "Treated",
#   id_column = "sample_id"
# )
# print(result_groups$results)

# Example: CreateComparisonName
# comp_name <- CreateComparisonName(condition == "Control", condition == "Treated")
# print(comp_name) # Output: "Control_vs_Treated"

# Example: Vplot (conceptual data)
vplot_data <- data.frame(
  fc = runif(1000, 0.1, 10),
  pval = runif(1000, 1e-10, 1),
  sig = sample(c("up", "down", "none"), 1000, replace = TRUE, prob = c(0.05, 0.05, 0.9))
)
p_vplot <- Vplot(vplot_data)
p_vplot_sig <- VplotAddSig(p_vplot, vplot_data)
print(p_vplot_sig)
```

## Contributing

Contributions to `ZemmourLib` are welcome\! If you have suggestions for improvements, bug reports, or new functionalities, please open an issue or submit a pull request on the [GitHub repository](https://www.google.com/search?q=https://github.com/dzemmour/ZemmourLib).

## License

This package is licensed under the MIT License. See the `LICENSE` file for more details.

