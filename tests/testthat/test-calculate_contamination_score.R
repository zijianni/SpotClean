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

mbrain_obj <- createSlide(mbrain_raw, mbrain_slide_info)

background_bcs <- dplyr::filter(mbrain_slide_info$slide,tissue==0)$barcode

test_that("Invalid background indices",{
    expect_error(arcScore(mbrain_raw, "foo"),
                 "Invalid background barcodes")
})

test_that("Calculation",{
    expect_equal(arcScore(mbrain_raw, background_bcs),
                 0.05160659)
    expect_equal(arcScore(mbrain_raw, background_bcs[1:100]),
                 0.00414822)
    expect_equal(arcScore(mbrain_obj),
                 0.05160659)
})
