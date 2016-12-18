context("test - milr")

set.seed(5)
beta <- c(-0.5, 0.7, -0.9, 1.1)
miData <- DGP(70, 3, beta)
milr_result_1 <- milr(miData$Z, miData$X, miData$ID)
milr_result_2 <- milr(miData$Z, miData$X, miData$ID, lambda = -1)
milr_result_3 <- milr(miData$Z, miData$X, miData$ID, lambda = c(0.5, 0.7, 0.9, 1.1))
milr_result_4 <- milr(miData$Z, miData$X, miData$ID, lambda = -1, numLambda = 30L)
test_that("milr", {
  expect_message(milr(miData$Z, miData$X, miData$ID, lambda = 0), "Lasso-penalty is not used.")
  expect_message(milr(miData$Z, miData$X, miData$ID, lambda = -1), "The penalty term is selected automatically with 20 candidates.")
  expect_message(milr(miData$Z, miData$X, miData$ID, lambda = c(0.5, 0.7, 0.9, 1.1)), "Use the user-defined lambdas")
  expect_message(milr(miData$Z, miData$X, miData$ID, lambda = -1, numLambda = 30L), 
                 "The penalty term is selected automatically with 30 candidates.")
  
  expect_is(milr_result_1, "milr")
  expect_equal(length(milr_result_1$beta), 4L)
  expect_equal(length(coef(milr_result_1)), 4L)
  expect_equal(length(fitted(milr_result_1)$bag), 70L)
  expect_equal(length(fitted(milr_result_1)$instance), 210L)
  
  expect_is(milr_result_2, "milr")
  expect_equal(dim(milr_result_2$beta), c(4L, 20L))
  expect_equal(length(coef(milr_result_2)), 4L)
  expect_equal(length(fitted(milr_result_2)$bag), 70L)
  expect_equal(length(fitted(milr_result_2)$instance), 210L)
  
  expect_is(milr_result_3, "milr")
  expect_equal(dim(milr_result_3$beta), c(4L, 4L))
  expect_equal(length(coef(milr_result_3)), 4L)
  expect_equal(length(fitted(milr_result_3)$bag), 70L)
  expect_equal(length(fitted(milr_result_3)$instance), 210L)
  
  expect_is(milr_result_4, "milr")
  expect_equal(dim(milr_result_4$beta), c(4L, 30L))
  expect_equal(length(coef(milr_result_4)), 4L)
  expect_equal(length(fitted(milr_result_4)$bag), 70L)
  expect_equal(length(fitted(milr_result_4)$instance), 210L)
})

