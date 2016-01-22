#' @export
#' @method coef milr
coef.milr <- function(object){
  return(object$coeffiecents)
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
#' @importFrom dplyr summarize
#' @export
milr <- function(y, x, bag = NULL, lambda = 0, alpha = 1, maxit = 500) {
  init_beta <- glm(y~x-1)$coefficients
  beta <- CLR_lasso(y, x, bag, init_beta, lambda)
  
  p_instance <- ifelse(logit(x, beta) > 0.5, 1, 0)
  fit_y <- data.frame(p_instance, bag) %>%
    group_by(bag) %>%
    summarize(p_bag = ifelse(sum(p_instance) > 1, 1, 0)) %>%
    select(p_bag) %>%
    as.matrix()
  out <- list(coeffiecents = as.vector(beta), fitted = as.vector(fit_y))
  class(out) <- 'milr'
  return(out)
}