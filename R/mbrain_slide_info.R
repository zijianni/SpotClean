#' @title Example 10x Visium spatial data: slide information
#'
#' @description  This dataset contains one list \code{mbrain_slide_info}.
#' The first slot \code{slide} of \code{mbrain_slide_info} is a
#' data.frame of the spatial information for spots. The second slot \code{grob}
#' is a Grob object of the tissue image. The original dataset can be found at
#' https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain.
#'
#' If users read raw data using \code{Read10xSlide}, they should get
#' exactly the same format as the objects in this dataset.
#'
#' @docType data
#'
#' @usage data(mbrain_slide_info)
#'
#' @format A list of slide information containing an object of class
#' \code{"data.frame"} and an object of class \code{"grob"}.
#'
#' @keywords datasets
#'
#' @source \href{https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Adult_Mouse_Brain}{V1_Adult_Mouse_Brain}
#'
#' @examples
#' data(mbrain_slide_info)
#' str(mbrain_slide_info)
#'
#' @importFrom grid grob
"mbrain_slide_info"
