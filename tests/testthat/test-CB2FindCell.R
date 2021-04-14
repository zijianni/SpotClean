# Note: this test was deleted due to Bioconductor's time limit.
# The example of CB2FindCell.R runs almost the same code.
# So we are good if everything is OK when checking examples.

# data(mbrainSub)

# Only tested one upper and lower threshold due to 
# Bioconductor time limits

# for(lower in c(100)){
#     for(upper in c(1000)){
#         CBOut_temp <- CB2FindCell(mbrainSub,
#                                 lower = lower, 
#                                 upper = upper,
#                                 Ncores = 2, verbose = FALSE)
#         
#         test_that("Output format", {
#             expect_true(all(c("cluster_matrix","cell_matrix")%in%names(CBOut_temp)))
#         })
#         
#         cell_b <- colnames(cbind(CBOut_temp$cluster_matrix,
#                                  CBOut_temp$cell_matrix))
#     
#         test_that("Lower threshold", {
#             expect_true(all(Matrix::colSums(mbrainSub[,cell_b])>lower))
#         })
#     
#         upper_b <- colnames(mbrainSub)[Matrix::colSums(mbrainSub)>=upper]
#         test_that("upper threshold", {
#             expect_true(all(upper_b%in%cell_b))
#         })
#         
#     }
# }

