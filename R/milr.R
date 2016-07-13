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
  return(coef(object) %>>% (split(logit(cbind(1, newdata), .), bag_newdata)) %>>%
           purrr::map_int(~1-prod(1-.) > 0.5))
}

#' @export
#' @method summary milr
summary.milr <- function(object, ...){
  if (object$lambda_chosen == 0)
  {
    summary <- list(loglik = object$loglik, beta = object$coeffiecents, se = sqrt(object$var),
                    z = object$coeffiecents / sqrt(object$var), lambda = object$lambda_chosen)
    summary$pvalue <- pnorm(abs(summary$z), 0, 1, FALSE) * 2
  } else
  {
    summary <- list(loglik = object$loglik, beta = object$coeffiecents, 
                    lambda = object$lambda_chosen)
  }
  class(summary) <- "summary.milr"
  return(summary)
}

#' @export
#' @method print summary.milr
print.summary.milr <- function(x, digits = max(3L, getOption("digits") - 2L), ...){
  message(sprintf("Log-Likelihood: %.3f.", x$loglik))
  if (x$lambda == 0)
  {
    outMat <- cbind(x$beta, x$se, x$z, x$pvalue) %>>% 
      `colnames<-`(c("Estimate", "Std.Err", "Z value", "Pr(>z)"))
    message("Estimates:")
    printCoefmat(outMat, digits = digits)
  } else
  {
    message(sprintf("Chosen Penalty: %.3f.", x$lambda))
    outMat <- cbind(x$beta) %>>% `colnames<-`("Estimate")
    message("Estimates:")
    printCoefmat(outMat, digits = digits)
  }
}

# function of generation of cv index
cvIndex_f <- function(n, fold){
  fold_n <- n %/% fold
  rem <- n - fold_n * fold
  index <- rep(1:fold, fold_n)
  if (rem > 0)
    index <- c(index, 1:rem)
  return(sample(index, length(index)))
}


#' Maximum likelihood estimation of multiple-instance logistic regression with LASSO penalty
#' 
#' Please refer to \link{milr-package}.
#' 
#' @param y a vector. Bag-level binary labels.
#' @param x the design matrix. The number of rows of \code{x} must be equal to the length of \code{y}.
#' @param bag a vector, bag id.
#' @param lambda the tuning parameter for LASSO-penalty.  If \code{lambda} is a real value number, then the \code{milr} 
#'  fits the model based on this lambda value.  Second, if \code{lambda} is vector, then the optimal lambda value would be
#'  be chosen based on the optimality criterion, \code{lambdaCriterion}.  
#'  Finally, if \code{lambda = -1}, then the optimal lambda value would be chosen automatically.
#'  The default is 0. 
#' @param lambdaCriterion a string, the used optimality criterion for tuning the \code{lambda} value.
#'  It can be specified with \code{lambdaCriterion = "BIC"} or \code{lambdaCriterion = "deviance"}.
#' @param nfold an integer, the number of fold for cross-validation to choose the optimal \code{lambda} when
#'  \code{lambdaCriterion = "deviance"}.
#' @param maxit an integer, the maximum iteration for the EM algorithm.
#' @return a list including deviance (not cv deviance), BIC, chosen lambda, coefficients, 
#'  fitted values, log-likelihood and variances of coefficients.
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
#' 
#' # use cv in auto-tuning
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID, 
#'                     lambda = -1, lambdaCriterion = "deviance")
#' coef(milr_result)      # coefficients
#' fitted(milr_result)    # fitted values
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID) # predicted label
#' @importFrom purrr map map_int map2_dbl map_dbl
#' @importFrom numDeriv hessian
#' @name milr
#' @rdname milr
#' @export
milr <- function(y, x, bag, lambda = 0, lambdaCriterion = "BIC", nfold = 10, maxit = 500) {
  # if x is vector, transform it to matrix
  if (is.vector(x))
    x <- matrix(x, ncol = 1)
  if (!is.matrix(x))
    x <- as.matrix(x)
  if (!all(y %in% c(0, 1)))
    stop('y must be 0 and 1.')
  bag <- factor(bag) %>>% as.integer
  # input check
  alpha <- 1
  assert_that(length(unique(y)) == 2, length(y) == nrow(x),
              all(is.finite(y)), is.numeric(y), all(is.finite(x)), is.numeric(x), 
              length(lambda) >= 1, all(is.finite(lambda)), is.numeric(lambda), 
              is.finite(alpha), is.numeric(alpha), is.finite(maxit), is.numeric(maxit), 
              abs(maxit - floor(maxit)) < 1e-4, nfold < nrow(x),
              lambdaCriterion %in% c("deviance", "BIC"))
  loglik <- function(b, Z, X, ID){
    pii <- split(as.data.frame(X), ID) %>>% purrr::map(~1-prod(1-logit(as.matrix(.), b)))
    return(purrr::map2_dbl(split(Z, ID) %>>% purrr::map(unique), pii, 
                           ~.x*log(.y)+(1-.x)*log(1-.y)) %>>% sum(na.rm = TRUE))
  }
  
  if (length(lambda) == 1 && all(lambda == -1))
  {
    message("The penalty term is selected automatically.")
    m <- table(bag)
    zi <- tapply(y, bag, function(x) sum(x) > 0) %>>% as.numeric
    lambdaMax <- sqrt(sum(m-1)) * sqrt(sum(m**(1-2*zi)))
    lambda <- exp(seq(log(lambdaMax/1000), log(lambdaMax), length = 20))
  }
  
  # initial value for coefficients
  init_beta <- coef(glm(y~x))
  unique_bag <- unique(bag)
  n_bag <- length(unique_bag)
  if (length(lambda) > 1)
  {
    # save the history of betas along lambda
    message(sprintf("Using the criterion %s to choose optimized penalty.", lambdaCriterion))
    beta_history <- matrix(NA, ncol(x) + 1, length(lambda) + 1)
    beta_history[ , 1] <- init_beta
    
    if (lambdaCriterion == "BIC")
    {
      dev <- vector('numeric', length(lambda))
      BIC <- vector('numeric', length(lambda))
      for (i in seq_along(lambda))
      {
        beta_history[, i + 1] <- CLR_lasso(y, cbind(1, x), bag, beta_history[, i], lambda[i])
        dev[i] <- -2 * loglik(beta_history[, i + 1], y, cbind(1, x), bag)
        BIC[i] <- dev[i] + sum(beta_history[, i + 1] != 0) * log(n_bag)
      }
      locMin <- which.min(BIC)
    } else if (lambdaCriterion == "deviance")
    {
      # calculate the bag label
      bagLvl <- tapply(y, bag, sum) %>>% (as.integer(. > 0))
      # avoid the sample size of one class is less than nfold to retrun error in cvIndex_f.
      if (nfold > min(table(bagLvl)))
      {
        message("nfold is bigger than the number of size of one class.")
        message(sprintf("Therefore, nfold is changed to %i.", min(table(bagLvl))))
        nfold <- min(table(bagLvl))
      }
      # stratified sampling for cv on bag labels
      cvIndex_bag_p <- cvIndex_f(sum(bagLvl == 1), nfold)
      cvIndex_bag_n <- cvIndex_f(sum(bagLvl == 0), nfold)
      # convert the cv index on bag labels to cv index on samples
      cvIndex_sample <- purrr::map(1:nfold, ~c(which(bagLvl == 1)[which(cvIndex_bag_p == .)],
                                               which(bagLvl == 0)[which(cvIndex_bag_n == .)])) %>>%
        purrr::map(~which(bag %in% unique_bag[.]))
      dev_cv <- matrix(NA, length(lambda), nfold)
      dev <- vector('numeric', length(lambda))
      BIC <- vector('numeric', length(lambda))
      cv_betas <- purrr::map(1:nfold, ~beta_history)
      for (i in seq_along(lambda))
      {
        # CV to calculate the deviances
        for (j in 1:nfold)
        {
          trainSetIndex <- do.call(c, cvIndex_sample[setdiff(1:nfold, j)])
          cv_betas[[j]][ , i + 1] <- CLR_lasso(y[trainSetIndex], cbind(1, x[trainSetIndex, ]), 
                                               bag[trainSetIndex], cv_betas[[j]][ , i], lambda[i])
          dev_cv[i, j] <- -2 * loglik(cv_betas[[j]][ , i + 1], y[cvIndex_sample[[j]]], 
                                      cbind(1, x[cvIndex_sample[[j]], ]), bag[cvIndex_sample[[j]]])
        }
        # calculate the beta under lambda
        beta_history[, i + 1] <- CLR_lasso(y, cbind(1, x), bag, beta_history[, i], lambda[i])
        dev[i] <- -2 * loglik(beta_history[, i + 1], y, cbind(1, x), bag)
        BIC[i] <- dev[i] + sum(beta_history[, i + 1] != 0) * log(n_bag)
      }
      locMin <- which.min(rowMeans(dev_cv))
    }
    lambda_out <- lambda[locMin]
    message(sprintf("The chosen penalty throuth %s is %.4f.", lambdaCriterion, lambda[locMin]))
    beta <- beta_history[ , locMin + 1]
  } else
  {
    beta <- CLR_lasso(y, cbind(1, x), bag, init_beta, lambda, alpha, maxit)  
    dev <- -2 * loglik(beta, y, cbind(1, x), bag)
    BIC <- dev + sum(beta != 0) * log(n_bag)
    lambda_out <- lambda
  }
  # return beta with names
  beta <- as.vector(beta) %>>% `names<-`(c("intercept", colnames(x)))
  # fitted response
  fit_y <- beta %>>% (split(logit(cbind(1, x), .), bag)) %>>%
    purrr::map_int(~1-prod(1-.) > 0.5)
  # calculate the hessian matrix when without lasso
  if (lambda_out == 0)
    bateVar = -diag(solve(numDeriv::hessian(function(b) loglik(b, y, cbind(1, x), bag), beta)))
  else
    bateVar = NULL
  out <- list(deviance = dev, BIC = BIC, lambda_chosen = lambda_out, 
              coeffiecents = beta, fitted = fit_y, loglik = loglik(beta, y, cbind(1, x), bag),
              var = bateVar)
  class(out) <- 'milr'
  return(out)
}