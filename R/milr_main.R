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

#' @export
#' @method predict milr
predict.milr <- function(object, newdata, bag_newdata, ...){
  return(coef(object) %>% {split(logit(cbind(1, newdata), .) > 0.5, bag_newdata)} %>%
         map_int(~ifelse(sum(.) > 1, 1L, 0L)))
}

#' milr function:
#'
#' description here
#' 
#' @param y A vector, its elements contain two different values.
#' @param x The design matrix. The number of rows of x must be equal to the length of y.
#' @param bag bag variable.
#' @param lambda The panelty for lasso. Default is exp(seq(log(0.01), log(50),length = 30)).
#'   The panelty is chosen by BIC.
#' @param alpha 123
#' @param maxit An integer, the maximum iteration for EM algorithm.
#' @return An list includes BIC, chosen lambda, coefficients and fitted values.
#' @examples
#' set.seed(100)
#' beta <- runif(10, -5, 5)
#' trainData <- DGP(70, 3, beta)
#' testData <- DGP(30, 3, beta)
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID, 
#'   exp(seq(log(0.01), log(50),length = 30)))
#' coef(milr_result)      # coefficients
#' fitted(milr_result)    # fitted values
#' predict(milr_result, testData$X, testData$ID) # predicted label
#' @importFrom magrittr set_names
#' @importFrom purrr map map_int map2_dbl
#' @importFrom logistf logistf
#' @export
milr <- function(y, x, bag, lambda = 0, alpha = 1, maxit = 500) {
  # if x is vector, transform it to matrix
  if (is.vector(x))
    x <- matrix(x, ncol = 1)
  # input check
  assert_that(length(unique(y)) == 2, length(y) == nrow(x),
              all(is.finite(y)), is.numeric(y), all(is.finite(x)), is.numeric(x), 
              length(lambda) >= 1, all(is.finite(lambda)), is.numeric(lambda), 
              is.finite(alpha), is.numeric(alpha), is.finite(maxit), is.numeric(maxit), 
              abs(maxit - floor(maxit)) < 1e-4)
  
  # initial value for coefficients
  init_beta <- coef(logistf(y~x))
  if (length(lambda) > 1)
  {
    cat("Using BIC to choose optimized panelty.\n")
    n_bag <- length(unique(bag))
    loglik <- function(b, Z, X, ID){
      pii <- split(as.data.frame(X), ID) %>% 
        purrr::map(~1-prod(1-logit(as.matrix(.), b)))
      return(purrr::map2_dbl(split(Z, ID) %>% purrr::map(unique), pii, 
        ~.x*log(.y)+(1-.x)*log(1-.y)) %>% sum(na.rm = TRUE))
    }
    BIC <- vector('numeric', length(lambda))
    for (i in seq_along(lambda))
    {
      beta_select <- CLR_lasso(y, cbind(1, x), bag, init_beta, lambda[i])
      BIC[i] <- -2 * loglik(beta_select, y, cbind(1, x), bag) + 
        sum(beta_select != 0) * log(n_bag)
    }
    lambda_out <- lambda[which.min(BIC)]
    beta <- CLR_lasso(y, cbind(1, x), bag, init_beta, lambda_out, alpha, maxit)  
  } else if (lambda == -1)
  {
    cat("Using auto-tuning to choose panelty.\n")
    cat("This part is undone.\n")
    beta <- init_beta
    BIC <- NULL
    lambda_out <- NULL
  } else
  {
    cat("Using auto-tuning to choose panelty.\n")
    beta <- CLR_lasso(y, cbind(1, x), bag, init_beta, lambda, alpha, maxit)  
    BIC <- -2 * loglik(beta, y, cbind(1, x), bag) + sum(beta != 0) * log(n_bag)
    lambda_out <- lambda
  }
  beta %<>% as.vector %>% set_names(c("intercept", colnames(x)))
  fit_y <- beta %>% {split(logit(cbind(1, x), .) > 0.5, bag)} %>%
    purrr::map_int(~ifelse(sum(.) > 1, 1L, 0L))
  out <- list(BIC = BIC, lambda_chosen = lambda_out, 
              coeffiecents = beta, fitted = fit_y)
  class(out) <- 'milr'
  return(out)
}
