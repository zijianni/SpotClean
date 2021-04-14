#' @title Read 10x Space Ranger output data
#'
#' @description \code{Read10xRaw} is a one-line handy function for reading
#' the raw expression data from 10x Space Ranger outputs and producing a count
#' matrix as an R object.
#'
#' @param count_dir The directory of 10x output matrix data. The directory should include
#' three files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
#'
#' @return If \code{meta = TRUE}, \code{Read10xRaw()} or \code{Read10xRawH5()}
#'
#' @examples
#'
#' @importFrom utils read.delim
#' @importFrom Matrix readMM
#' @importFrom Matrix sparseMatrix
#'
#'
#' @export



CreateSlide <- function(count_mat, slide){



}
