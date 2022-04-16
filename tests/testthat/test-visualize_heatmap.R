data(mbrain_raw)
data(mbrain_slide_info)

mbrain_obj <- CreateSlide(mbrain_raw, mbrain_slide_info)

mbrain_slide_info2 <- mbrain_slide_info
mbrain_slide_info2$slide <- head(mbrain_slide_info2$slide, 100)
mbrain_obj2 <- CreateSlide(mbrain_raw[,mbrain_slide_info2$slide$barcode], mbrain_slide_info2)


test_that("Non-existing gene",{
    expect_error(VisualizeHeatmap(mbrain_obj, value="foo"),
                 "Specified gene does not exist")
    expect_error(VisualizeHeatmap(mbrain_slide_info$slide, value="foo",mbrain_raw),
                 "Specified gene does not exist")

})

test_that("Missing input matrix",{
    expect_error(VisualizeHeatmap(mbrain_slide_info$slide, value="foo"),
                 "You must provide an input expression")

})

test_that("Invalid input value vector",{
    expect_error(VisualizeHeatmap(mbrain_obj, value=rnorm(10)),
                 "Invalid value input")
    expect_error(VisualizeHeatmap(mbrain_obj, value=rnorm(10), mbrain_raw),
                 "Invalid value input")

})

test_that("Matrix not match slide barcodes",{
    gp1 <- VisualizeHeatmap(mbrain_slide_info$slide, value="Bc1",
                     mbrain_raw[,mbrain_slide_info2$slide$barcode])
    gp2 <- VisualizeHeatmap(mbrain_slide_info2$slide, value="Bc1",
                            mbrain_raw)
    expect_identical(gp1$data, gp2$data)

})

gp1 <- VisualizeHeatmap(mbrain_obj, value="Bc1",
                      subset_barcodes = colnames(mbrain_obj2))
gp2 <- VisualizeHeatmap(mbrain_obj2, value="Bc1")

gp3 <- VisualizeHeatmap(mbrain_slide_info$slide, value="Bc1", mbrain_raw,
                        subset_barcodes = colnames(mbrain_obj2))
gp4 <- VisualizeHeatmap(mbrain_slide_info2$slide, value="Bc1", mbrain_raw)

test_that("Subsetting barcodes", {
    expect_identical(gp1$data, gp2$data[match(gp1$data$barcode,
                                              gp2$data$barcode),])
    expect_identical(gp3$data, gp4$data[match(gp3$data$barcode,
                                              gp4$data$barcode),])
})

test_that("Object class", {
    expect_s3_class(gp1,c("gg","ggplot"))
    expect_s3_class(gp3,c("gg","ggplot"))
})
