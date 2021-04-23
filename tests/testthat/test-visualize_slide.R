data(MbrainSmall)
mbrain_obj <- CreateSlide(mbrain_raw, mbrain_slide_info)

gp <- VisualizeSlide(mbrain_obj)

test_that("Object class", {
    expect_identical(class(gp),c("gg","ggplot"))

})

metadata(mbrain_obj)$grob <- NULL

test_that("Missing image grob", {
    expect_error(VisualizeSlide(mbrain_obj),
                 "No valid image information")

})
