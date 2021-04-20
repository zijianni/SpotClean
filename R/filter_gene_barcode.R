#' @title Filter out low count genes and barcodes from count matrix (deprecated)
#'
#' @description This function is used for filtering out low count genes and barcodes
#' from count matrix based on total gene expression counts (row sums) and
#' barcode expression counts (column sums).
#'
#' @param count_mat (matrix of num) Input count matrix to be filtered.
#'
#' @param g_threshold (nonnegative num) Filtering
#' threshold for genes. Any gene whose total expression count is less or
#' equal to \code{g_threshold} will be filtered out. Default: \code{0}.
#'
#' @param b_threshold (nonnegative num) Filtering
#' threshold for barcodes. Any barcode whose total count is less or equal
#' to \code{b_threshold} will be filtered out. Default: \code{0}.
#'
#' @return A filtered matrix with the same format as input matrix.

#' @examples
#' data(MbrainSmall)
#' dim(MbrainSmall)
#' MbrainSmall_f <- FilterGB(MbrainSmall)
#' dim(MbrainSmall_f)
#'
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#'

# FilterGB <- function(count_mat,
#                     g_threshold = 0,
#                     b_threshold = 0) {
#     #filter barcodes and genes
#     bc <- colSums(count_mat)
#     count_mat <- count_mat[, bc > b_threshold]
#     gc <- rowSums(count_mat)
#     count_mat <- count_mat[gc > g_threshold, ]
#     return(count_mat)
# }
