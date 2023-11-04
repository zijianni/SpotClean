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


test_that("Barcodes not matching", {
    expect_warning(createSlide(mbrain_raw[,1:10], mbrain_slide_info),
                 "Barcodes in count matrix do not match")
})

test_that("Rearrange barcode order", {
    slide <- mbrain_slide_info$slide
    mbrain_slide_info_2 <- mbrain_slide_info
    mbrain_slide_info_2$slide <- slide[sample(nrow(slide)),]

    expect_identical(createSlide(mbrain_raw, mbrain_slide_info),
                     createSlide(mbrain_raw, mbrain_slide_info_2))
})

slide_obj <- createSlide(mbrain_raw, mbrain_slide_info, gene_cutoff = 100)

test_that("Object class", {
    expect_s4_class(slide_obj,"SummarizedExperiment")
})

test_that("Gene filtering", {
    expect_true(nrow(slide_obj)==14)
})

test_that("input data frame",{
    expect_warning(createSlide(mbrain_raw, mbrain_slide_info$slide),
                   "Input slide info does not contain image")
})
