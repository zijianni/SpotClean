#' @title Create a new slide object
#'
#' @description This function takes input of the count matrix
#' (from \code{Read10xRaw} or \code{Read10xRawH5}) and slide information
#' (from \code{Read10xSlide}) and outputs a \code{SummarizedExperiment} object
#' as our slide object for downstream decontamination and visualization.
#'
#'
#' @param count_mat (matrix of num) The raw gene-by-barcode count matrix.
#' Can be either standard matrix format or sparse matrix format.
#'
#' @param slide_info (list) The slide information from \code{Read10xSlide()}.
#'
#' @return A \code{SummarizedExperiment} object containing
#' gene expression and spot metadata.
#'
#' @examples
#'
#' data(MbrainSmall)
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' mbrain_obj


#' @import SummarizedExperiment
#'
#' @export

CreateSlide <- function(count_mat, slide_info){

    if(!identical(sort(colnames(count_mat)),
                 sort(slide_info$slide$barcode)
                 )){
        stop("Barcodes in count matrix do not match ",
             "barcodes in slide information.")
    }

    # rearrange barcodes
    count_mat <- count_mat[,slide_info$slide$barcode]

    total_counts <- colSums(count_mat)
    slide_info$slide$total_counts <- total_counts

    obj <- SummarizedExperiment(assays=list(raw=count_mat),
                                metadata=slide_info)
    return(obj)
}
