data(mbrain_raw)
data(mbrain_slide_info)


test_that("Barcodes not matching", {
    expect_error(createSlide(mbrain_raw[,1:10], mbrain_slide_info),
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
