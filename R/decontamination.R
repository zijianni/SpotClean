#' @title Decontaminate spot swapping effect in spatial transcriptomics data
#'
#' @description This is the main function for decontaminating spot swapping
#' effect in spatial transcriptomics data.
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
