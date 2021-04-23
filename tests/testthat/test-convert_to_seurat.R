data(MbrainSmall)

mbrain_obj <- CreateSlide(mbrain_raw, mbrain_slide_info)

test_that("Invalid assay",{
    expect_error(ConvertToSeurat(mbrain_obj, "foo", "a"),
                 "Specified assay not found in the slide object")
})

test_that("Invalid image directory",{
    expect_error(ConvertToSeurat(mbrain_obj, "foo", "raw"),
                 "unable to open foo/tissue_lowres_image.png")
})
