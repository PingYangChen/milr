#' DGP: data generation
#'
#' Generating the multiple-instance data set.
#' 
#' @param n an integer. The number of bags.
#' @param m an integer or vector of length \code{n}. If \code{m} is an integer, each bag has the identical number of instances, \code{m}. 
#'   If \code{m} is a vector, the \code{i}th bag has \code{m[i]} instances. 
#' @param beta a vector. The true regression coefficients.
#' @return a list including (1) bag-level labels, \code{Z}, (2) the design matrix, \code{X}, and (3) bag ID of each instance, \code{ID}.
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
  X <- scale(matrix(rnorm(sum(m)*p),sum(m),p)) %>>%
    `attr<-`("scaled:center", NULL) %>>% `attr<-`("scaled:scale", NULL) # remove attributes
  X[ ,1] <- rep(1, nrow(X))
  Y <- logit(X, beta) %>>% (rbinom(sum(m), 1, .))
  ID <- rep(1:n, m)
  Z <- split(Y, ID) %>>% map(~rep(any(. == 1), length(.))) %>>% unlist %>>% as.integer
  if (all(Z == 1))
    Z[ID == sample(1:n, 1)] <- 0
  if (all(Z == 0))
    Z[ID == sample(1:n, 1)] <- 1
  return(list(Z = Z, X = X[ , 2:ncol(X)], ID = ID))
}