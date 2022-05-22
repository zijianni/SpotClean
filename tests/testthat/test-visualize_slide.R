data(mbrain_raw)
data(mbrain_slide_info)

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
