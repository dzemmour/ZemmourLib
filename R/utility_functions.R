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

#test
