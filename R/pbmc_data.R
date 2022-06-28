#' Peripheral Blood Mononuclear Cells (PBMC) data.
#'
#' Count matrix of Peripheral Blood Mononuclear Cells (PBMC). 
#' The original data set is from 10X genomics.
#' 
#' @format A large dgCMatrix: 32378 * 2700
#' \describe{
#'   \item{i}{Row index of the non-zero values}
#'   \item{p}{A vector to refer the column index of the non-zero values}
#'   \item{Dim}{Dimension of the matrix}
#'   \item{Dimnames}{A list of length 2 containing the row names and column names of the matrix}
#'   \item{x}{Vector containing all the non-zero values}
#'   
#' }
#' @source \url{https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
"pbmc_data"