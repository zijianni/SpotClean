#' @title Calculate the ambient RNA contamination (ARC) score
#'
#' @description The ARC score is intended for measuring the level of
#' contamination caused by ambient RNAs for a given sptial transcriptomics
#' or droplet-based single-cell RNA-seq data. Intuitively, this score is a lower
#' bound of the average proportion of contaminated expressions in tissue
#' spots or cell droplets.
#'
#' @param count_mat (matrix of num) Input gene-by-barcode count matrix.
#' Ideally this should be the raw count matrix with all genes and barcodes
#' without any filtering.
#' Can be either standard matrix format or sparse matrix format.
#'
#' @param background_idx (vector of int) Indices of background barcodes in
#' \code{count_mat}. For spatial data, these are background spots not
#' covered by tissue. For single-cell data, these are empty droplets not
#' containing cells.
#'
#' @return (num) The ARC score of given data.
#'
#' @examples
#'
#' data(MbrainSmall)
#' mbrain_raw <- mbrain_raw[,mbrain_slide_info$slide$barcode]
#' background_index <- which(mbrain_slide_info$slide$tissue==0)
#' ARCScore(mbrain_raw, background_index)

#' @import Matrix

#' @export

ARCScore <- function(count_mat, background_idx){

    if(!all(background_idx%in%seq_len(ncol(count_mat)))){
        stop("Invalid background indices.")
    }

    total_counts <- colSums(count_mat)
    # underestimated contamination per barcode
    cont <- sum(total_counts[background_idx])/ncol(count_mat)
    # expression per tissue/cell barcode
    expr <- mean(total_counts[-background_idx])
    # proportion
    arc <- cont/expr

    return(arc)
}
