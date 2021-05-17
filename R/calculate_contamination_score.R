#' @title Calculate the ambient RNA contamination (ARC) score
#'
#' @description The ARC score is intended for subjectively measuring the level
#' of contamination caused by ambient RNAs for a given sptial transcriptomics
#' or droplet-based single-cell RNA-seq data. Intuitively, this score is a lower
#' bound of the average proportion of contaminated expressions in tissue
#' spots or cell droplets.
#'
#' This lower bound is very conservative as it only
#' accounts for contamination in background spots or empty droplets, neglecting
#' contamination in tissue spots or cell droplets. However, it is totally
#' subjective, not relying on any complicated assumptions and modelling.
#' ARC score serves as a relative contamination measurement for comparison
#' among different datasets and protocols.
#'
#' @param object A gene-by-barcode count matrix or a slide object created or
#' inherited from \code{CreateSlide()}.
#' Ideally this should give the raw count matrix with all genes and barcodes
#' without any filtering.
#'
#' @param ... Arguments passed to other methods
#'
#' @return (num) The ARC score of given data.
#'
#' @examples
#'
#' data(mbrain_raw)
#' data(mbrain_slide_info)
#' background_bcs <- dplyr::filter(mbrain_slide_info$slide, tissue==0)$barcode
#' ARCScore(mbrain_raw, background_bcs)
#'
#' mbrain_obj <- CreateSlide(mbrain_raw, mbrain_slide_info)
#' ARCScore(mbrain_obj)
#'
#' @rdname ARCScore
#'
#' @export

ARCScore <- function(object, ...) {
    UseMethod(generic = "ARCScore", object = object)
}


#' @param background_bcs (vector of chr) Background barcodes in
#' \code{count_mat}. For spatial data, these are background spots not
#' covered by tissue. For single-cell data, these are empty droplets not
#' containing cells.
#'
#' @importFrom Matrix colSums
#' @importFrom methods as
#'
#' @method ARCScore default
#' @rdname ARCScore
#' @export
#'
ARCScore.default <- function(object, background_bcs, ...){

    if (!inherits(x = object, 'Matrix')) {
        object <- as(object = as.matrix(x = object), Class = 'Matrix')
    }
    if (!inherits(x = object, what = 'dgCMatrix')) {
        object <- as(object = object, Class = 'dgCMatrix')
    }

    if(!all(background_bcs%in%colnames(object))){
        stop("Invalid background barcodes.")
    }

    total_counts <- colSums(object)
    # underestimated contamination per barcode
    cont <- sum(total_counts[background_bcs])/ncol(object)
    # expression per tissue/cell barcode
    ts_bcs <- setdiff(colnames(object),background_bcs)
    expr <- mean(total_counts[ts_bcs])
    # proportion
    arc <- cont/expr

    return(arc)
}


#' @method ARCScore SummarizedExperiment
#' @rdname ARCScore
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @importFrom dplyr filter
#'
#' @export
#'
ARCScore.SummarizedExperiment <- function(object, ...){

    # junk code... get rid of R CMD check notes
    tissue <- NULL

    background_bcs <- filter(metadata(object)$slide,tissue==0)$barcode
    count_mat <- assay(object)

    total_counts <- colSums(count_mat)
    # underestimated contamination per barcode
    cont <- sum(total_counts[background_bcs])/ncol(count_mat)
    # expression per tissue/cell barcode
    ts_bcs <- setdiff(colnames(count_mat),background_bcs)
    expr <- mean(total_counts[ts_bcs])
    # proportion
    arc <- cont/expr

    return(arc)
}
