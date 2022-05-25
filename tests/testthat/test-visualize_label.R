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

mbrain_slide_info2 <- mbrain_slide_info
mbrain_slide_info2$slide <- head(mbrain_slide_info2$slide, 100)
mbrain_obj2 <- createSlide(mbrain_raw[,mbrain_slide_info2$slide$barcode],
                           mbrain_slide_info2)

test_that("Non-existing column",{
    expect_error(visualizeLabel(mbrain_obj, label="foo"),
                 "Label name does not exist")
    expect_error(visualizeLabel(mbrain_slide_info$slide, label="foo"),
                 "Label name does not exist")

})

test_that("Invalid input label vector",{
    expect_error(visualizeLabel(mbrain_obj, label=rnorm(10)),
                 "Invalid label input")
    expect_error(visualizeLabel(mbrain_slide_info$slide, label=rnorm(10)),
                 "Invalid label input")

})

gp1 <- visualizeLabel(mbrain_obj, subset_barcodes = colnames(mbrain_obj2))
gp2 <- visualizeLabel(mbrain_obj2)

test_that("Subsetting barcodes", {
    expect_identical(gp1$data, gp2$data[match(gp1$data$barcode,
                                              gp2$data$barcode),])

})

test_that("Object class", {
    expect_s3_class(gp1,c("gg","ggplot"))

})
