hgmmh5 <- "/path/to/h5/file.h5"
hg19 <- "/path/to/folder/"


test_that("Both directory and HDF5 are used", {
    expect_error(QuickCB2(dir=hg19, h5file=hgmmh5),"dir and h5file can't be specified together.")
})

