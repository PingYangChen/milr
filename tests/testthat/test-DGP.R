context("test - DGP")

beta <- c(-0.5, 0.7, -0.9, 1.1)
miData <- DGP(70, 3, beta)
test_that("DGP", {
  expect_equal(length(miData$Z), 210L)
  expect_equal(dim(miData$X), c(210L, 3L))
  expect_equal(length(miData$ID), 210L)
  expect_equal(miData$ID, rep(1L:70L, each = 3))
  expect_true(all(abs(colMeans(miData$X)) < 1e-6))
  expect_true(all(tapply(miData$Z, miData$ID, function(x) all(diff(x) == 0))))
})
