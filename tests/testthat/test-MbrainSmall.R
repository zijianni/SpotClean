data(MbrainSmall)

test_that("Data dimension",{
    expect_equal(dim(mbrain_raw),c(100,4992))
})

test_that("Slide info",{
    expect_true(all(names(mbrain_slide_info)==c("slide","grob")))
})

test_that("Gene and barcode",{
    expect_identical(rownames(mbrain_raw)[1],"Bc1")
    expect_identical(colnames(mbrain_raw)[999],"ATATCGTTCCTCGAAC-1")
    expect_identical(sort(colnames(mbrain_raw)),
                     sort(mbrain_slide_info$slide$barcode))
})

test_that("Values",{
    expect_equal(mbrain_raw[10,10],14)
    expect_equal(mbrain_slide_info$slide$total_counts[5], 695)
})
