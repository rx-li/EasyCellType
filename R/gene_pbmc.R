#' Differential expressed marker genes in 9 clusters.
#'
#' A data frame containing marker genes, clusters as well as the average of log 2 fold changes.
#' The original data set is from 10X genomics, and we followed the standard workflow 
#' provided by Seurat package to process data, and then format to get the data frame.
#' 
#' @usage data(gene_pbmc)
#' 
#' @format A data frame with 727 rows and 3 variables:
#' \describe{
#'   \item{gene}{Entrez IDs of the marker genes}
#'   \item{cluster}{Cluster}
#'   \item{score}{Average of log 2 fold changes getting from the process procedure}
#'   
#' }
#' @source \url{https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
"gene_pbmc"