#' @title Calculate the ambient RNA contamination (ARC) score
#'
#' @description The ARC score is intended for subjectively measuring the level
#' of contamination caused by ambient RNAs for a given spatial transcriptomics
#' or droplet-based single-cell RNA-seq data. Intuitively, this score is a lower
#' bound of the average proportion of contaminated expressions in tissue
#' spots or cell droplets in the observed contaminated data.
#'
#' ARC score is calculated as follows: (1) estimate the average amount of
#' contaminated UMI counts received per spot (droplet) using the total UMI
#' counts in background spots (empty droplets) divided by total number of
#' spots (droplets). This is an underestimation as it does not account for
#' contamination inside tissue spots or cell droplets. (2) estimate the average
#' amount of observed UMI counts per tissue spot (cell droplet) using the
#' total UMI counts in tissue spots (cell droplets) divided by number of
#' tissue spots (cell droplets). (3) Divide the value in (1) by the value in
#' (2) to get the ARC score.
#'
#' This lower bound is very conservative as it only
#' accounts for observable contamination in background spots or empty droplets,
#' neglecting contamination in tissue spots or cell droplets. However, it
#' is totally subjective, not relying on any complicated assumptions and
#' modelling. ARC score serves as a relative contamination measurement for
#' comparison among different datasets and protocols.
#'
#' @param object A gene-by-barcode count matrix or a slide object created or
#' inherited from \code{createSlide()}.
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
#' spatial_dir <- system.file(file.path("extdata",
#'                                      "V1_Adult_Mouse_Brain_spatial"),
#'                            package = "SpotClean")
#' mbrain_slide_info <- read10xSlide(tissue_csv_file=file.path(spatial_dir,
#'                                        "tissue_positions_list.csv"),
#'              tissue_img_file = file.path(spatial_dir,
#'                                        "tissue_lowres_image.png"),
#'              scale_factor_file = file.path(spatial_dir,
#'                                        "scalefactors_json.json"))
#'
#' background_bcs <- dplyr::filter(mbrain_slide_info$slide, tissue==0)$barcode
#' arcScore(mbrain_raw, background_bcs)
#'
#' mbrain_obj <- createSlide(mbrain_raw, mbrain_slide_info)
#' arcScore(mbrain_obj)
#'
#' @rdname arcScore
#'
#' @export

arcScore <- function(object, ...) {
    UseMethod(generic = "arcScore", object = object)
}


#' @param background_bcs (vector of chr) Background barcodes in
#' \code{count_mat}. For spatial data, these are background spots not
#' covered by tissue. For single-cell data, these are empty droplets not
#' containing cells.
#'
#' @importFrom Matrix colSums
#' @importFrom methods as
#'
#' @method arcScore default
#' @rdname arcScore
#' @export
#'
arcScore.default <- function(object, background_bcs, ...){

    if (!inherits(x = object, 'Matrix')) {
        object <- as(object = as.matrix(x = object), Class = 'Matrix')
    }
    if (!inherits(x = object, what = 'dgCMatrix')) {
        object <- as(object = object, Class = 'CsparseMatrix')
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


#' @method arcScore SummarizedExperiment
#' @rdname arcScore
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @export
#'
arcScore.SummarizedExperiment <- function(object, ...){

    background_bcs <- filter(metadata(object)$slide, .data$tissue==0)$barcode
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
