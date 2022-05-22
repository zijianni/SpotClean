#' @title Create a new slide object
#'
#' @description This function takes input of the count matrix
#' (from \code{read10xRaw} or \code{read10xRawH5}) and slide information
#' (from \code{read10xSlide} or manually specified data frame) and outputs
#' a \code{SummarizedExperiment} object
#' as our slide object for downstream decontamination and visualization.
#'
#'
#' @param count_mat (matrix of num) The raw gene-by-barcode count matrix.
#' Can be either standard matrix format or sparse matrix format.
#'
#' @param slide_info (list or data.frame) A list of slide information from
#' \code{read10xSlide()}, or a data frame only containing spot information
#' like barcode, tissue, imagerow, imagecol, etc.
#'
#' @param gene_cutoff (num) Filter out genes with average expressions
#' among tissue spots below or equal to this cutoff.
#' Default: 0.1
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
#' mbrain_obj <- createSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' mbrain_obj


#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Matrix rowMeans
#'
#'
#' @export

createSlide <- function(count_mat, slide_info, gene_cutoff=0.1, verbose=TRUE){

    # validate arguments
    if(!is.numeric(gene_cutoff)){
        stop("invalid argument (gene_cutoff)")
    }
    if(gene_cutoff<0) {
        stop("gene_cutoff should be non-negative")
    }
    if(!is.logical(verbose)){
        stop("invalid argument (verbose)")
    }
    if(is.data.frame(slide_info) &
       all(c("barcode","tissue","imagerow","imagecol"
             )%in%colnames(slide_info))){
        warning("Input slide info does not contain image.\n")
        slide_info <- list(slide=slide_info)
    }

    rownames(slide_info$slide) <- slide_info$slide$barcode
    if(!identical(sort(colnames(count_mat)),
                 sort(slide_info$slide$barcode)
                 )){
        stop("Barcodes in count matrix do not match ",
             "barcodes in slide information.")
    }

    # rearrange barcodes
    slide_info$slide <- slide_info$slide[colnames(count_mat),]
    count_ts_mat <- count_mat[,slide_info$slide$tissue==1]

    # filter genes
    gene_cutoff <- max(gene_cutoff,0)
    good_gene <- rowMeans(count_ts_mat)>gene_cutoff

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
