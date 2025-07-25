% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_functions.R
\name{MyPlots}
\alias{MyPlots}
\title{Generate UMAP Plots for Seurat Object}
\usage{
MyPlots(
  seurat_object = so,
  dim1 = seurat_object[["umap_unintegrated"]]@cell.embeddings[, 1],
  dim2 = seurat_object[["umap_unintegrated"]]@cell.embeddings[, 2],
  color_by = "spleen_standard",
  split_by1 = "IGT",
  split_by2 = NULL,
  genes = c("Foxp3", "Il2ra"),
  cluster_key = "ClusterSCVI_Res",
  mypal = NULL
)
}
\arguments{
\item{seurat_object}{A Seurat object. Default is \code{so}.}

\item{dim1}{Numeric vector, UMAP dimension 1 embeddings. Default is \code{so[["umap_unintegrated"]]@cell.embeddings[,1]}.}

\item{dim2}{Numeric vector, UMAP dimension 2 embeddings. Default is \code{so[["umap_unintegrated"]]@cell.embeddings[,2]}.}

\item{color_by}{Character string, metadata column to color points by. Default is "spleen_standard".}

\item{split_by1}{Character string, metadata column to split plots by (rows in facet_grid). Default is "IGT". Can be NULL.}

\item{split_by2}{Character string, metadata column to split plots by (columns in facet_grid). Default is NULL.}

\item{genes}{Character vector, genes to plot expression for. Default is \code{c("Foxp3", "Il2ra")}.}

\item{cluster_key}{Character string, key for cluster information in metadata. Default is "ClusterSCVI_Res".}

\item{mypal}{A color palette to use for discrete colors. Default is \code{glasbey()}.}
}
\value{
Prints multiple ggplot objects.
}
\description{
This function generates multiple UMAP plots for a Seurat object,
including coloring by a metadata column, splitting by metadata columns,
and visualizing gene expression.
}
