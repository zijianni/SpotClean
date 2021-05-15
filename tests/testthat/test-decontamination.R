data(mbrain_raw)
data(mbrain_slide_info)

mbrain_obj <- CreateSlide(mbrain_raw,
                          mbrain_slide_info)

# test output

test_that("Decontamination",{
    expect_silent(mbrain_decont_obj <- SpotClean(mbrain_obj, candidate_radius=20,
                                           maxit = 3, verbose = FALSE))
    expect_s4_class(mbrain_decont_obj,"SummarizedExperiment")
    expect_identical(names(mbrain_decont_obj@assays),"decont")
    expect_equal(metadata(mbrain_decont_obj)$ARC_score, 0.05160659)
    expect_equal(metadata(mbrain_decont_obj)$bleeding_rate, 0.3301657,
                 tolerance = 4e-8)
    expect_equal(metadata(mbrain_decont_obj)$distal_rate, 0.2154095,
                 tolerance = 4e-8)
    expect_equal(metadata(mbrain_decont_obj)$weight_matrix[1,1], 4.916e-05,
                 tolerance = 4e-8)
    expect_equal(metadata(mbrain_decont_obj)$loglh,
                 c(167030816, 167314060, 167379802))
    expect_equal(unname(metadata(mbrain_decont_obj)$contamination_rate[1:3]),
                 c(0.06273360, 0.08079689, 0.09040354))

    expect_equal(rowSums(assay(mbrain_decont_obj)),rowSums(assay(mbrain_obj)))
    expect_equal(assay(mbrain_decont_obj)[1,1],2645.022778)

})


# test internal computations

test_that("Internal computations", {
    small_mat <- matrix(1:4,2,2)
    small_edist <- .calculate_euclidean_weight(small_mat)
    small_gdist <- .gaussian_kernel(small_edist,10)

    expect_equal(small_edist[1,2],sqrt(2))
    expect_equal(small_edist,t(small_edist))
    expect_equal(small_gdist[1,2],0.9900498,tolerance = 4e-8)
    expect_equal(.points_to_sdv(10,5), 25)

})

# test input parameters
test_that("Data loading", {
    expect_error(SpotClean(mbrain_obj, gene_keep = "foo"),
                 "Specified genes not found")
})

names(mbrain_obj@assays) <- "decont"

test_that("Wrong assay", {
    expect_error(SpotClean(mbrain_obj, gene_keep = "foo"),
                 "Cannot find raw data")
})
