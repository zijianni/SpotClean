data(mbrain_raw)
spatial_dir <- system.file(file.path("extdata",
                                     "V1_Adult_Mouse_Brain_spatial"),
                           package = "SpotClean")
mbrain_slide_info <- read10xSlide(
    tissue_csv_file=file.path(spatial_dir,
                              "tissue_positions_list.csv"),
    tissue_img_file = file.path(spatial_dir,
                                "tissue_lowres_image.png"),
    scale_factor_file = file.path(spatial_dir,
                                  "scalefactors_json.json"))

mbrain_obj <- createSlide(mbrain_raw,
                          mbrain_slide_info)

# test output

test_that("Decontamination",{
    expect_no_error(mbrain_decont_obj <- spotclean(
        mbrain_obj, candidate_radius=20, maxit = 3, verbose = FALSE))
    expect_s4_class(mbrain_decont_obj,"SummarizedExperiment")
    expect_identical(names(mbrain_decont_obj@assays),"decont")
    expect_equal(metadata(mbrain_decont_obj)$ARC_score, 0.05160659)
    expect_equal(metadata(mbrain_decont_obj)$bleeding_rate, 0.3301657,
                 tolerance = 4e-8)
    expect_equal(metadata(mbrain_decont_obj)$distal_rate, 0.2154095,
                 tolerance = 4e-8)
    expect_equal(metadata(mbrain_decont_obj)$weight_matrix[
        "AAACAAGTATCTCCCA-1","AAACAACGAATAGTTC-1"
    ], 0.00004315095, tolerance = 4e-8)
    expect_equal(metadata(mbrain_decont_obj)$loglh,
                 c(167030816, 167314060, 167379802))
    expect_equal(unname(metadata(mbrain_decont_obj)$contamination_rate[
        c("AAACAAGTATCTCCCA-1", "AAACAATCTACTAGCA-1", "AAACACCAATAACTGC-1")
        ]), c(0.2794837, 0.3102042, 0.2128133), tolerance = 4e-8)

    expect_equal(rowSums(assay(mbrain_decont_obj)),rowSums(assay(mbrain_obj)))
    expect_equal(assay(mbrain_decont_obj)["Bc1","AAACAAGTATCTCCCA-1"],
                 1457.840109)

})


# test internal computations

test_that("Internal computations", {
    small_mat <- matrix(1:4,2,2)
    small_edist <- .calculate_euclidean_weight(small_mat)
    small_gdist <- .gaussian_kernel(small_edist,10)
    small_ldist <- .linear_kernel(small_edist,20)
    small_lapdist <- .laplace_kernel(small_edist,30)
    small_cdist <- .cauchy_kernel(small_edist,40)

    expect_equal(small_edist[1,2],sqrt(2))
    expect_equal(small_edist,t(small_edist))
    expect_equal(small_gdist[1,2],0.99004984)
    expect_equal(small_ldist[1,2],18.5857865)
    expect_equal(small_lapdist[1,2],0.9539534)
    expect_equal(small_cdist[1,2],0.99875156)

    expect_equal(.points_to_sdv(10,5), 25)

})

# test input parameters
test_that("Data loading", {
    expect_error(spotclean(mbrain_obj, gene_keep = "foo"),
                 "Specified genes not found")
})

names(mbrain_obj@assays) <- "decont"

test_that("Wrong assay", {
    expect_error(spotclean(mbrain_obj, gene_keep = "foo"),
                 "Cannot find raw data")
})
