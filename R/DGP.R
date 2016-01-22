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
#' @importFrom dplyr data_frame group_by mutate
#' @export
DGP <- function(n, m, beta){
  assert_that(length(n) == 1, is.numeric(n), is.finite(n), n > 0, abs(n - floor(n)) < 1e-6)
  assert_that(length(m) == 1, is.numeric(m), is.finite(m), m > 0, abs(m - floor(m)) < 1e-6)
  assert_that(all(is.numeric(beta)), all(is.finite(beta)), n > 0)
  
  p <- length(beta)
  if(length(m) < n){
    m <- rep(m,length=n)
  }
  X <- scale(matrix(rnorm(sum(m)*p),sum(m),p))
  X[,1] <- 1
  pr <- logit(X, beta)
  # mu <- X %*% beta
  # pr <- 1/(1+exp(-mu))
  Y <- rbinom(sum(m), 1, pr)
  ID <- rep(1:n,m)
  Z <- data_frame(Y=Y,ID=ID) %>%
    group_by(ID) %>%
    mutate(Z = ifelse(all(Y==0),0,1)) %>% .$Z
  # ungroup() %>%
  # select(Z) %>%
  # as.matrix() %>%
  # as.numeric()
  if(all(Z==1)) Z[ID==sample(1:n,1)] <- 0
  if(all(Z==0)) Z[ID==sample(1:n,1)] <- 1
  return(list(Z=Z,X=X,ID=ID))
}
