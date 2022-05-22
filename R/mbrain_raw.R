#' @title Example 10x Visium spatial data: raw count matrix
#'
#' @description  This dataset contains one sparse count matrix
#' \code{mbrain_raw} of gene expressions. The original dataset can be found at
#' https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain.
#' For simplicity, we only keep the top 100 highest expressed genes in this
#' example data.
#'
#' If users read raw data using \code{read10xRaw} (or \code{read10xRawH5}),
#' they should get exactly the same format as the object in this dataset.
#'
#' @docType data
#'
#' @usage data(mbrain_raw)
#'
#' @format An object of class \code{"dgCMatrix"}.
#'
#' @keywords datasets
#'
#' @source \href{https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain}{V1_Adult_Mouse_Brain}
#'
#' @examples
#' data(mbrain_raw)
#' str(mbrain_raw)
#'
#' @importClassesFrom Matrix dgCMatrix
"mbrain_raw"
