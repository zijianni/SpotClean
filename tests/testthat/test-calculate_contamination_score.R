data(mbrain_raw)
data(mbrain_slide_info)

mbrain_raw <- mbrain_raw[,mbrain_slide_info$slide$barcode]
background_index <- which(mbrain_slide_info$slide$tissue==0)

test_that("Invalid background indices",{
    expect_error(ARCScore(mbrain_raw, 10000),
                 "Invalid background indices")
})

test_that("Calculation",{
    expect_equal(ARCScore(mbrain_raw, background_index),
                 0.05160659)
    expect_equal(ARCScore(mbrain_raw, 1:100),
                 0.00543154)
})
