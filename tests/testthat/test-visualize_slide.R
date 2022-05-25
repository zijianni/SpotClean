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

gp <- visualizeSlide(mbrain_obj)

test_that("Object class", {
    expect_s3_class(gp,c("gg","ggplot"))

})

S4Vectors::metadata(mbrain_obj)$grob <- NULL

test_that("Missing image grob", {
    expect_error(visualizeSlide(mbrain_obj),
                 "No valid image information")

})
