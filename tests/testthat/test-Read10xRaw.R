data("mbrainSub")

data.dir <- file.path(tempdir(), "CB2_testthat")
H5.dir <- file.path(tempdir(), "CB2_testthat_H5")

if (dir.exists(data.dir)) {
    unlink(data.dir, recursive = TRUE)
}
if (file.exists(H5.dir)) {
    file.remove(H5.dir)
}

DropletUtils::write10xCounts(data.dir, mbrainSub, type = "sparse",
                             version = "2")
DropletUtils::write10xCounts(H5.dir, mbrainSub, type = "HDF5",
                             version = "3")


test_that("Directory/files existence", {
    expect_error(Read10xRaw(paste0(data.dir,"/foo")), 
                 "Directory does not exist")
    expect_error(Read10xRaw(tempdir()), "No 10x output file detected")
})

test_that("Data loading", {
    expect_identical(Read10xRaw(data.dir), Read10xRawH5(H5.dir))
})


test_that("Metadata loading", {
    expect_true(is.list(Read10xRaw(data.dir, meta = TRUE)))
})
