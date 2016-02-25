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

#' @export
#' @method summary milr
summary.milr <- function(object, ...){
  summary <- list(loglik = object$loglik, beta = object$coeffiecents, se = sqrt(object$var),
                  z = object$coeffiecents / sqrt(object$var))
  summary$pvalue <- pnorm(abs(summary$z), 0, 1, FALSE) * 2
  class(summary) <- "summary.milr"
  return(summary)
}

#' @export
#' @method print summary.milr
#' @importFrom magrittr set_colnames
print.summary.milr <- function(x, digits = max(3L, getOption("digits") - 2L), ...){
  cat(sprintf("Log-Likelihood: %.3f\n", x$loglik))
  outMat <- cbind(x$beta, x$se, x$z, x$pvalue) %>% 
    magrittr::set_colnames(c("Estimate", "Std.Err", "Z value", "Pr(>z)"))
  cat("Estimates:\n")
  printCoefmat(outMat, digits = digits)
}

#' Maximum likelihood estimation of multiple-instance logistic regression with LASSO penalty
#' 
#' Please refer to \link{milr-package}.
#' 
#' @param y A vector. Binay response.
#' @param x The design matrix. The number of rows of x must be equal to the length of y.
#' @param bag A vector, bag id.
#' @param lambda The penalty for LASSO. Default is 0 (not use LASSO). If \code{lambda} is vector, the penalty will be chosen by BIC.
#'   If \code{lambda} = 0, then the penalty will be chosen automatically.
#' @param maxit An integer, the maximum iteration for EM algorithm.
#' @return An list includes BIC, chosen lambda, coefficients, fitted values, log-likelihood and variances of coefficients.
#' @examples
#' set.seed(100)
#' beta <- runif(5, -5, 5)
#' trainData <- DGP(70, 3, beta)
#' testData <- DGP(30, 3, beta)
#' # default (not use LASSO)
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID)
#' coef(milr_result)      # coefficients
#' fitted(milr_result)    # fitted values
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID) # predicted label
#' 
#' # use BIC to choose penalty
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID,
#'   exp(seq(log(0.01), log(50),length = 30)))
#' coef(milr_result)      # coefficients
#' fitted(milr_result)    # fitted values
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID) # predicted label
#' 
#' # use auto-tuning
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID, lambda = -1)
#' coef(milr_result)      # coefficients
#' fitted(milr_result)    # fitted values
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID) # predicted label
#' @importFrom magrittr set_names
#' @importFrom purrr map map_int map2_dbl
#' @importFrom logistf logistf
#' @importFrom numDeriv hessian
#' @name milr
#' @rdname milr
#' @export
milr <- function(y, x, bag, lambda = 0, maxit = 500) {
  # if x is vector, transform it to matrix
  if (is.vector(x))
    x <- matrix(x, ncol = 1)
  if (!is.matrix(x))
    x %<>% as.matrix
  bag %<>% factor %>% as.integer
  # input check
  alpha <- 1
  assert_that(length(unique(y)) == 2, length(y) == nrow(x),
              all(is.finite(y)), is.numeric(y), all(is.finite(x)), is.numeric(x), 
              length(lambda) >= 1, all(is.finite(lambda)), is.numeric(lambda), 
              is.finite(alpha), is.numeric(alpha), is.finite(maxit), is.numeric(maxit), 
              abs(maxit - floor(maxit)) < 1e-4)
  loglik <- function(b, Z, X, ID){
    pii <- split(as.data.frame(X), ID) %>% 
      purrr::map(~1-prod(1-logit(as.matrix(.), b)))
    return(purrr::map2_dbl(split(Z, ID) %>% purrr::map(unique), pii, 
                           ~.x*log(.y)+(1-.x)*log(1-.y)) %>% sum(na.rm = TRUE))
  }

  if (length(lambda) == 1 && all(lambda == -1))
  {
    cat("Lambda is selected automatically.\n")
    m <- table(bag)
    zi <- tapply(y, bag, function(x) sum(x) > 0) %>% as.numeric
    lambdaMax <- sqrt(sum(m-1)) * sqrt(sum(m**(1-2*zi)))
    lambda <- exp(seq(log(lambdaMax/1000), log(lambdaMax), length = 20))
  }
  
  # initial value for coefficients
  # init_beta <- coef(logistf(y~x))
  init_beta <- coef(glm(y~x, family = binomial(link="logit")))
  n_bag <- length(unique(bag))
  if (length(lambda) > 1)
  {
    cat("Using BIC to choose optimized penalty.\n")
    BIC <- vector('numeric', length(lambda))
    for (i in seq_along(lambda))
    {
      beta_select <- CLR_lasso(y, cbind(1, x), bag, init_beta, lambda[i])
      BIC[i] <- -2 * loglik(beta_select, y, cbind(1, x), bag) + 
        sum(beta_select != 0) * log(n_bag)
    }
    lambda_out <- lambda[which.min(BIC)]
    beta <- CLR_lasso(y, cbind(1, x), bag, init_beta, lambda_out, alpha, maxit)  
  } else
  {
    beta <- CLR_lasso(y, cbind(1, x), bag, init_beta, lambda, alpha, maxit)  
    BIC <- -2 * loglik(beta, y, cbind(1, x), bag) + sum(beta != 0) * log(n_bag)
    lambda_out <- lambda
  }
  beta %<>% as.vector %>% set_names(c("intercept", colnames(x)))
  fit_y <- beta %>% {split(logit(cbind(1, x), .) > 0.5, bag)} %>%
    purrr::map_int(~ifelse(sum(.) > 1, 1L, 0L))
  out <- list(BIC = BIC, lambda_chosen = lambda_out, 
              coeffiecents = beta, fitted = fit_y, loglik = loglik(beta, y, cbind(1, x), bag),
              var = -diag(solve(numDeriv::hessian(function(b) loglik(b, y, cbind(1, x), bag), beta))))
  class(out) <- 'milr'
  return(out)
}
