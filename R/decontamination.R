#' @title Decontaminate spot swapping effect in spatial transcriptomics data
#'
#' @description This is the main function implementing the TBD method
#' for decontaminating spot swapping effect in spatial transcriptomics data.
#'
#' @param slide_obj A slide object created or inherited from
#' \code{CreateSlide()}.
#'
#' @param verbose (logical) Whether print progress information
#' Default: \code{TRUE}
#'
#' @return A slide object where the decontaminated expression matrix is in the
#' "decont" assay slot and the contamination statistics are in the "cont_stat"
#' metadata slot.
#'
#' @details Briefly, the contamination level for the slide is estimated based on
#' the total counts of all spots. UMI counts travelling around the slide are
#' assumed to follow Poisson distributions and modeled by a mixture of
#' Gaussian and uniform kernels and adjusted for the mRNA density in each spot.
#' The underlying uncontaminated gene expressions are estimated by EM algorithm
#' to maximize the data likelihood.
#'
#'
#' @examples
#'
#' data(MbrainSmall)
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#'
#' @import Matrix
#'
#' @export
