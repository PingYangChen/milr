context("test - softmax")

set.seed(5)
beta <- c(-0.5, 0.7, -0.9, 1.1)
miData <- DGP(70, 3, beta)
softmax_result <- softmax(miData$Z, miData$X, miData$ID, alpha = 0)
test_that("softmax", {
  expect_is(softmax_result, "softmax")
  expect_equal(softmax_result$alpha, 0)
  expect_equal(length(coef(softmax_result)), 4L)
  expect_equal(length(fitted(softmax_result)), 70L)
})

