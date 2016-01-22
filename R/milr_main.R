#' @export
#' @method coef milr
coef.milr <- function(object, ...){
  return(object$coeffiecents)
}

#' @export
#' @method fitted milr
fitted.milr <- function(object, ...){
  return(object$fitted)
}

#' milr function:
#'
#' description here
#' 
#' @param y 123
#' @param x 123
#' @param bag 123
#' @param lambda 123
#' @param alpha 123
#' @param maxit An integer, the maximum iteration for EM algorithm.
#' @return An list involves coefficients and fitted values.
#' @examples
#' testData <- DGP(50, 3, runif(sample(5:21, 1), -5, 5))
#' a <- milr(testData$Z, testData$X, testData$ID)
#' coef(a)
#' @importFrom purrr map_int
#' @export
milr <- function(y, x, bag = NULL, lambda = 0, alpha = 1, maxit = 500) {
  init_beta <- glm(y~x-1)$coefficients
  beta <- CLR_lasso(y, x, bag, init_beta, lambda, alpha, maxit)
  fit_y <- beta %>% {split(logit(x, .) > 0.5, bag)} %>%
    map_int(~ifelse(sum(.) > 1, 1L, 0L))
  out <- list(coeffiecents = as.vector(beta), fitted = fit_y)
  class(out) <- 'milr'
  return(out)
}