#' SpotClean: a computational method to adjust for spot swapping in
#' spatial transcriptomics data
#'
#' SpotClean is a computational method to adjust for spot swapping in
#' spatial transcriptomics data. Recent spatial transcriptomics experiments
#' utilize slides containing thousands of spots with spot-specific barcodes
#' that bind mRNA. Ideally, unique molecular identifiers at a spot measure
#' spot-specific expression, but this is often not the case due to bleed from
#' nearby spots, an artifact we refer to as spot swapping. SpotClean is able
#' to estimate the contamination rate in observed data and decontaminate the
#' spot swapping effect, thus increase the sensitivity and precision of
#' downstream analyses.
#'
#' To learn more about \code{SpotClean}, see the vignette
#' using \code{browseVignettes(package = "SpotClean")}.
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
