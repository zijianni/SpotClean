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
#' @param gene_cutoff (num) Filter out genes with average expressions
#' among all spots below this cutoff.
#' Default: 0.1.
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
#' @import Matrix
#'
#' @export

CreateSlide <- function(count_mat, slide_info, gene_cutoff=0.1){

    if(!identical(sort(colnames(count_mat)),
                 sort(slide_info$slide$barcode)
                 )){
        stop("Barcodes in count matrix do not match ",
             "barcodes in slide information.")
    }

    # rearrange barcodes
    count_mat <- count_mat[,slide_info$slide$barcode]

    # filter genes
    good_gene <- rowMeans(count_mat)>=gene_cutoff
    message("Filtered out ",sum(!good_gene)
            ," genes with average expressions below ",gene_cutoff, ".")
    count_mat <- count_mat[good_gene,]

    obj <- SummarizedExperiment(assays=list(raw=count_mat),
                                metadata=slide_info)
    return(obj)
}
