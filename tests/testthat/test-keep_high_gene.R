data(mbrain_raw)
data(mbrain_slide_info)

mbrain_filter <- KeepHighGene(mbrain_raw)

test_that("Default parameters", {
    expect_true(nrow(mbrain_filter)==100)

})

test_that("Object class", {
    expect_s4_class(mbrain_filter, "dgCMatrix")
})

test_that("Top high gene cutoff", {
    expect_true(nrow(KeepHighGene(mbrain_raw, top_high = 50, mean_cutoff = 1))==50)

})

test_that("Mean expression cutoff", {
    expect_true(nrow(KeepHighGene(mbrain_raw, mean_cutoff = Inf))==7)
    expect_true(nrow(KeepHighGene(mbrain_raw, mean_cutoff = 500))==10)
})
