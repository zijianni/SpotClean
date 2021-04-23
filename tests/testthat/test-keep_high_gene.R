data(MbrainSmall)

test_that("Default parameters", {
    expect_true(nrow(KeepHighGene(mbrain_raw))==100)

})

test_that("Top high gene cutoff", {
    expect_true(nrow(KeepHighGene(mbrain_raw, top_high = 50, mean_cutoff = 1))==50)

})

test_that("Mean expression cutoff", {
    expect_true(nrow(KeepHighGene(mbrain_raw, mean_cutoff = Inf))==7)
    expect_true(nrow(KeepHighGene(mbrain_raw, mean_cutoff = 500))==10)
})
