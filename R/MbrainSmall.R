#' Example 10x Visium spatial data
#'
#' This dataset contains one sparse count matrix \code{mbrain_raw} of
#' gene expressions and one list \code{mbrain_slide_info}. The first slot
#' \code{slide} of \code{mbrain_slide_info} is a data.frame of the spatial
#' information for spots. The second slot \code{grob} is a Grob object of
#' the tissue image. The original dataset can be found at
#' https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain.
#' For simplicity, we only keep the top 100 highest expressed genes in this
#' example data.
#'
#' If users read raw data using \code{Read10xRaw} (or \code{Read10xRawH5}) and
#' \code{Read10xSlide}, they should get exactly the same format as the objects
#' in this dataset.
#'
#' @docType data
#'
#' @usage data(MbrainSmall)
#'
#' @format An object of class \code{"dgCMatrix"}, a list of slide information
#' containing an object of class
#' \code{"data.frame"} and an object of class \code{"grob"}.
#'
#' @keywords datasets
#'
#' @source \href{https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain}
#'
#' @examples
#' data(MbrainSmall)
#' str(mbrain_raw)
#' str(mbrain_slide_info)
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom grid grob rastergrob gDesc
"MbrainSmall"
