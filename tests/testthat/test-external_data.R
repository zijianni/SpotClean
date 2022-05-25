data(mbrain_raw)

test_that("Data dimension",{
    expect_equal(dim(mbrain_raw),c(100,4992))
})

test_that("Gene and barcode",{
    expect_identical(rownames(mbrain_raw)[1],"Bc1")
    expect_identical(colnames(mbrain_raw)[999],"ATATCGTTCCTCGAAC-1")

})

test_that("Values",{
    expect_equal(mbrain_raw[10,10],14)
})
