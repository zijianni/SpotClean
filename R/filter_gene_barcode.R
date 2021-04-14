#' Filter out low count genes and barcodes from count matrix
#' 
#' This function is used for filtering out low count genes and barcodes 
#' from count matrix based on total gene expression count (row sums) and 
#' barcode expression count (column sums). \code{CB2FindCell}
#' has already integrated this function into it with \code{g_threshold = 0} 
#' and \code{b_threshold = 0}. If users plan to customize their filtering 
#' threshold, this function can be applied to the raw expression
#' count matrix prior to running \code{CB2FindCell}.
#' 
#' @param dat Input count matrix to be filtered.
#' 
#' @param g_threshold Nonnegative integer. Default: \code{0}. Filtering 
#' threshold for genes. Any gene whose total expression count is less or 
#' equal to \code{g_threshold} will be filtered out.
#' 
#' @param b_threshold Nonnegative integer. Default: \code{0}. Filtering 
#' threshold for barcodes. Any barcode whose total count is less or equal 
#' to \code{b_threshold} will be filtered out.
#' 
#' @return A filtered matrix with the same format as input matrix.

#' @examples 
#' data(mbrainSub)
#' dim(mbrainSub)
#' mbrainSub_f <- FilterGB(mbrainSub)
#' dim(mbrainSub_f)
#' 
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#'
#' @export
FilterGB <- function(dat,
                    g_threshold = 0,
                    b_threshold = 0) {
    #filter barcodes and genes
    bc <- colSums(dat)
    dat <- dat[, bc > b_threshold]
    gc <- rowSums(dat)
    dat <- dat[gc > g_threshold, ]
    return(dat)
}
