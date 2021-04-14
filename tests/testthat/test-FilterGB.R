data(mbrainSub)

test_that("Input data format", {
    expect_error(FilterGB(1), "at least two dimensions")
    expect_error(FilterGB("1"), "at least two dimensions")
})

test_that("Gene filtering threshold", {
    for(k in sample(0:100,5)){
        dat_temp <- FilterGB(mbrainSub,g_threshold = k)
        expect_true(all(Matrix::rowSums(dat_temp) > k))
    }
})

test_that("Barcode filtering threshold", {
    for(k in sample(0:100,5)){
        dat_temp <- FilterGB(mbrainSub,b_threshold = k)
        expect_true(all(Matrix::colSums(dat_temp) > k))
    }
})
