#' Subset of 1k Brain Cells from an E18 Mouse
#'
#' 1k Brain Cells from an E18 Mouse is a public dataset from 10X Genomics.
#' This subset is the first 50,000 barcodes of original data. 
#'
#' @docType data
#'
#' @usage data(mbrainSub)
#'
#' @format An object of class \code{"dgCMatrix"}.
#'
#' @keywords datasets
#'
#' @source \href{http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900}{1k Brain Cells from an E18 Mouse}
#'
#' @examples
#' data(mbrainSub)
#' str(mbrainSub)
#' 
#' @importClassesFrom Matrix dgCMatrix
"mbrainSub"