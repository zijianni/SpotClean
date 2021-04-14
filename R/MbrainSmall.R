#' Example 10x Visium spatial data
#'
#' This dataset contains one sparse count matrix \code{mbrain_raw} of
#' gene expressions, one data.frame \code{mbrain_slide} of the spatial
#' information for spots, and one grob object \code{mbrain_grob} of the tissue
#' image. The original dataset can be found at
#' https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain.
#' For simplicity, we only keep the top 100 highest expressed genes in this
#' example data.
#'
#' @docType data
#'
#' @usage data(MbrainSmall)
#'
#' @format An object of class \code{"dgCMatrix"}, an object of class
#' \code{"data.frame"} and an object of class \code{"grob"}.
#'
#'
#' @keywords datasets
#'
#' @source \href{https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain}
#'
#' @examples
#' data(MbrainSmall)
#' str(mbrain_raw)
#' str(mbrain_slide)
#' str(mbrain_grob)
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom grid grob rastergrob gDesc
"MbrainSmall"
