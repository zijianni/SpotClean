#' @title Create a new slide object
#'
#' @description This function takes input of the count matrix
#' (from \code{Read10xRaw} or \code{Read10xRawH5}) and slide information
#' (from \code{Read10xSlide} or manually specified data frame) and outputs
#' a \code{SummarizedExperiment} object
#' as our slide object for downstream decontamination and visualization.
#'
#'
#' @param count_mat (matrix of num) The raw gene-by-barcode count matrix.
#' Can be either standard matrix format or sparse matrix format.
#'
#' @param slide_info (list or data.frame) A list of slide information from
#' \code{Read10xSlide()}, or a data frame only containing spot information
#' like barcode, tissue, imagerow, imagecol, etc.
#'
#' @param gene_cutoff (num) Filter out genes with average expressions
#' among all spots below or equal to this cutoff.
#' Default: 0.1.
#'
#' @param verbose (logical) Whether print progress information.
#' Default: \code{TRUE}
#'
#' @return A \code{SummarizedExperiment} object containing
#' gene expression and spot metadata.
#'
#' @examples
#'
#' data(mbrain_raw)
#' data(mbrain_slide_info)
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' mbrain_obj


#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Matrix rowMeans
#'
#'
#' @export

CreateSlide <- function(count_mat, slide_info, gene_cutoff=0.1, verbose=TRUE){

    if(is.data.frame(slide_info) &
       all(c("barcode","tissue","imagerow","imagecol"
             )%in%colnames(slide_info))){
        warning("Input slide info does not contain image.\n")
        slide_info <- list(slide=slide_info)
    }

    if(!identical(sort(colnames(count_mat)),
                 sort(slide_info$slide$barcode)
                 )){
        stop("Barcodes in count matrix do not match ",
             "barcodes in slide information.")
    }

    # rearrange barcodes
    count_mat <- count_mat[,slide_info$slide$barcode]

    # filter genes
    gene_cutoff <- max(gene_cutoff,0)
    good_gene <- rowMeans(count_mat)>gene_cutoff

    if(verbose){
        message("Filtered out ",sum(!good_gene)
                ," genes with average expressions below or equal to ",
                gene_cutoff, ".")
    }

    count_mat <- count_mat[good_gene,]

    obj <- SummarizedExperiment(assays=list(raw=count_mat),
                                metadata=slide_info)
    return(obj)
}
