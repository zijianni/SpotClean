data(MbrainSmall)
mbrain_obj <- CreateSlide(mbrain_raw, mbrain_slide_info)

mbrain_slide_info2 <- mbrain_slide_info
mbrain_slide_info2$slide <- head(mbrain_slide_info2$slide, 100)
mbrain_obj2 <- CreateSlide(mbrain_raw[,mbrain_slide_info2$slide$barcode], mbrain_slide_info2)


test_that("Non-existing gene",{
    expect_error(VisualizeHeatmap(mbrain_obj, value="foo"),
                 "Specified gene does not exist")

})

test_that("Invalid input value vector",{
    expect_error(VisualizeHeatmap(mbrain_obj, value=rnorm(10)),
                 "Invalid value input")

})

# test_that("Invalid assay",{
#     expect_error(VisualizeHeatmap(mbrain_obj, value="Bc1", assay_name = "foo"),
#                  "Specified assay name does not exist")
#
# })

gp1 <- VisualizeHeatmap(mbrain_obj, value="Bc1",
                      subset_barcodes = colnames(mbrain_obj2))
gp2 <- VisualizeHeatmap(mbrain_obj2, value="Bc1", )

test_that("Subsetting barcodes", {
    expect_identical(gp1$data, gp2$data)

})

test_that("Object class", {
    expect_identical(class(gp1),c("gg","ggplot"))

})
