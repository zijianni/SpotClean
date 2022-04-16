data(mbrain_raw)
data(mbrain_slide_info)

mbrain_obj <- CreateSlide(mbrain_raw, mbrain_slide_info)

mbrain_slide_info2 <- mbrain_slide_info
mbrain_slide_info2$slide <- head(mbrain_slide_info2$slide, 100)
mbrain_obj2 <- CreateSlide(mbrain_raw[,mbrain_slide_info2$slide$barcode], mbrain_slide_info2)

test_that("Non-existing column",{
    expect_error(VisualizeLabel(mbrain_obj, label="foo"),
                 "Label name does not exist")
    expect_error(VisualizeLabel(mbrain_slide_info$slide, label="foo"),
                 "Label name does not exist")

})

test_that("Invalid input label vector",{
    expect_error(VisualizeLabel(mbrain_obj, label=rnorm(10)),
                 "Invalid label input")
    expect_error(VisualizeLabel(mbrain_slide_info$slide, label=rnorm(10)),
                 "Invalid label input")

})

gp1 <- VisualizeLabel(mbrain_obj, subset_barcodes = colnames(mbrain_obj2))
gp2 <- VisualizeLabel(mbrain_obj2)

test_that("Subsetting barcodes", {
    expect_identical(gp1$data, gp2$data[match(gp1$data$barcode,
                                              gp2$data$barcode),])

})

test_that("Object class", {
    expect_s3_class(gp1,c("gg","ggplot"))

})
