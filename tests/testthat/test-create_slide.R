data(MbrainSmall)

test_that("Barcodes not matching", {
    expect_error(CreateSlide(mbrain_raw[,1:10], mbrain_slide_info),
                 "Barcodes in count matrix do not match")
})

test_that("Rearrange barcode order", {
    expect_identical(CreateSlide(mbrain_raw[,sample(colnames(mbrain_raw))], mbrain_slide_info),
                     CreateSlide(mbrain_raw, mbrain_slide_info))
})

slide_obj <- CreateSlide(mbrain_raw, mbrain_slide_info, gene_cutoff = 100)

test_that("Object class", {
    expect_true(class(slide_obj)=="SummarizedExperiment")
})

test_that("Gene filtering", {
    expect_true(nrow(slide_obj)==11)
})

