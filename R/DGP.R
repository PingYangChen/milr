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
#' @importFrom stats rnorm rbinom
#' @export
DGP <- function(n, m, beta){
  assert_that(length(n) == 1, is.numeric(n), is.finite(n), n > 0, abs(n - floor(n)) < 1e-6)
  assert_that(all(is.numeric(m)), all(is.finite(m)), all(m > 0), all(abs(m - floor(m)) < 1e-6))
  assert_that(all(is.numeric(beta)), all(is.finite(beta)))
  
  p <- length(beta)
  if (length(m) < n)
    m <- rep(m, length = n)
  
  # get IDs
  ID <- rep(1L:n, m)
  # generate X
  X <- cbind(1, scale(matrix(rnorm(sum(m) * (p - 1L)), sum(m), p - 1L)))
  # generate output of each instances
  Y <- rbinom(sum(m), 1L, logit(X, beta))
  # find the bag response
  Z <- tapply(Y, ID, function(x){
    if (any(x == 1L))
      return(rep(1L, length(x)))
    return(rep(0L, length(x)))
  }) %>>% unlist %>>% `names<-`(NULL)
  
  # prevent only one class of Z
  if (all(Z == 1L))
    Z[ID == sample(1L:n, 1)] <- 0
  if (all(Z == 0L))
    Z[ID == sample(1L:n, 1)] <- 1
  
  return(list(Z = Z, X = X[ , 2:ncol(X)], ID = rep(1L:n, m)))
}