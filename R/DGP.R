#' DGP: data generation
#'
#' Simulate multiple-instance regression data.
#' 
#' @param n An integer. The number of subjects (bags).
#' @param m An integer or vector. When \code{m} is integer, subjects have the same number of instance. 
#'   When \code{m} is vector, subjects have corresponding number of instance. 
#' @param beta A vector. The true regression coefficients.
#' @return An list includes status of subjects (bags) \code{Z}, design matrix \code{X}, subject id for all components (instances) \code{ID}.
#' @examples
#' data1 <- DGP(50, 3, runif(10, -5, 5))
#' data2 <- DGP(50, sample(3:5, 50, TRUE), runif(10, -5, 5))
#' @importFrom purrr map
#' @export
DGP <- function(n, m, beta){
  assert_that(length(n) == 1, is.numeric(n), is.finite(n), n > 0, abs(n - floor(n)) < 1e-6)
  assert_that(all(is.numeric(m)), all(is.finite(m)), all(m > 0), all(abs(m - floor(m)) < 1e-6))
  assert_that(all(is.numeric(beta)), all(is.finite(beta)), n > 0)
  
  p <- length(beta)
  if (length(m) < n)
    m <- rep(m, length = n)
  X <- scale(matrix(rnorm(sum(m)*p),sum(m),p)) %>% 
    matrix(nrow(.), ncol(.))  # remove attributes
  X[ ,1] <- rep(1, nrow(X))
  Y <- logit(X, beta) %>% rbinom(sum(m), 1, .)
  ID <- rep(1:n, m)
  Z <- split(Y, ID) %>% map(~rep(any(. == 1), length(.))) %>% unlist %>% as.integer
  if(all(Z == 1))
    Z[ID == sample(1:n, 1)] <- 0
  if(all(Z == 0))
    Z[ID == sample(1:n, 1)] <- 1
  return(list(Z = Z, X = X, ID = ID))
}
