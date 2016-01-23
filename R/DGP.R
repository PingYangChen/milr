#' milr function:
#'
#' description here
#' 
#' @param n 123
#' @param m 123
#' @param beta 123
#' @return An list involves coefficients and fitted values.
#' @examples
#' data <- DGP(50, 3, runif(sample(5:21, 1), -5, 5))
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
  X[ ,1] <- 1
  pr <- logit(X, beta)
  Y <- rbinom(sum(m), 1, pr)
  ID <- rep(1:n, m)
  Z <- split(Y, ID) %>% map(~rep(any(.==1), length(.))) %>% 
    unlist %>% as.integer
  if(all(Z==1))
    Z[ID==sample(1:n,1)] <- 0
  if(all(Z==0))
    Z[ID==sample(1:n,1)] <- 1
  return(list(Z = Z, X = X, ID = ID))
}
