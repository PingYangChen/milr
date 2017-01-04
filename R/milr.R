#' @export
#' @method coef milr
coef.milr <- function(object, ...) {
  return(object$best_model$coeffiecents)
}

#' Fitted Response of milr Fits
#' 
#' @param object A fitted obejct of class inheriting from \code{"milr"}.
#' @param type The type of fitted response required. Default is \code{"bag"}, the fitted labels of bags.
#'   The \code{"instance"} option returns the fitted labels of instances.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method fitted milr
fitted.milr <- function(object, type = "bag", ...) {
  stopifnot(length(type) == 1)
  if (type == "bag") {
    return(object$best_model$fitted$bag)
  } else if (type == "instance") {
    return(object$best_model$fitted$instance)
  }
}

#' Predict Method for milr Fits
#' 
#' @param object A fitted obejct of class inheriting from \code{"milr"}.
#' @param newdata A matrix with variables to predict.
#' @param bag_newdata A vector. The labels of instances to bags.
#' @param type The type of prediction required. Default is \code{"bag"}, the predicted labels of bags.
#'   The \code{"instance"} option returns the predicted labels of instances.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method predict milr
predict.milr <- function(object, newdata, bag_newdata, type = "bag", ...) {
  stopifnot(length(type) == 1)
  if (type == "bag") {
    return(coef(object) %>>% getMilrProb(cbind(1, newdata), bag_newdata) %>>%
             `>`(0.5) %>>% as.numeric)
  } else if (type == "instance") {
    return(logit(cbind(1, newdata), coef(object)) %>>% `>`(0.5) %>>% as.numeric)
  }
}

#' @export
#' @method print milr
print.milr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (length(coef(x)) > 0) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  }
  cat("Residual Deviance:   ", format(signif(x$best_model$deviance, digits)), "\n")
  cat("              BIC:   ", format(signif(x$best_model$BIC, digits)), "\n")
  invisible(x)
}

#' @export
#' @method summary milr
summary.milr <- function(object, ...) {
  if (object$best_model$lambda_chosen == 0) {
    summary <- list(loglik = object$best_model$loglik, 
										beta = object$best_model$coeffiecents, 
										se = sqrt(object$best_model$var),
                    z = object$best_model$coeffiecents / sqrt(object$best_model$var), 
										lambda = object$best_model$lambda_chosen)
    summary$pvalue <- pnorm(abs(summary$z), 0, 1, FALSE) * 2
  } else {
    summary <- list(loglik = object$best_model$loglik, 
                    beta = object$best_model$coeffiecents, 
                    lambda = object$best_model$lambda_chosen)
  }
  class(summary) <- "summary.milr"
  return(summary)
}

#' @export
#' @method print summary.milr
print.summary.milr <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  message(sprintf("Log-Likelihood: %.3f.", x$loglik))
  if (x$lambda == 0) {
    outMat <- cbind(x$beta, x$se, x$z, x$pvalue) %>>% 
      `colnames<-`(c("Estimate", "Std.Err", "Z value", "Pr(>z)"))
    message("Estimates:")
    printCoefmat(outMat, digits = digits)
  } else {
    message(sprintf("Chosen Penalty: %.3f.", x$lambda))
    outMat <- cbind(x$beta) %>>% `colnames<-`("Estimate")
    message("Estimates:")
    printCoefmat(outMat, digits = digits)
  }
}

# function of generation of cv index
cvIndex_f <- function(n, fold) {
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
#' @param numLambda An integer, the maximum length of LASSO-penalty. in atuo-tunning mode 
#'  (\code{lambda = -1}). The default is 20.
#' @param lambdaCriterion a string, the used optimality criterion for tuning the \code{lambda} value.
#'  It can be specified with \code{lambdaCriterion = "BIC"} or \code{lambdaCriterion = "deviance"}.
#' @param nfold an integer, the number of fold for cross-validation to choose the optimal \code{lambda} when
#'  \code{lambdaCriterion = "deviance"}.
#' @param maxit an integer, the maximum iteration for the EM algorithm. The default is 1000.
#' @param tolerance Convergence threshold for coordinate descent. The default is 1e-5.
#' @return An object with S3 class "milr".
#' \itemize{
#' \item{lambda}{a vector of candidate lambda values.}
#' \item{cv}{a vector of predictive deviance via \code{nfold}-fold cross validation
#'  when \code{lambdaCriterion = "deviance"}.}
#' \item{deviance}{a vector of deviance of candidate model for each candidate lambda value.}
#' \item{BIC}{a vector of BIC of candidate model for each candidate lambda value.}
#' \item{best_index}{an integer, indicates the index of the best model among candidate lambda values.}
#' \item{best_model}{a list of the information for the best model including deviance (not cv deviance), 
#'  BIC, chosen lambda, coefficients, fitted values, log-likelihood and variances of coefficients.}
#' }
#' @examples
#' set.seed(100)
#' beta <- runif(5, -5, 5)
#' trainData <- DGP(70, 3, beta)
#' testData <- DGP(30, 3, beta)
#' # default (not use LASSO)
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID)
#' coef(milr_result)      # coefficients
#' fitted(milr_result)                    # fitted bag labels
#' fitted(milr_result, type = "instance") # fitted instance labels
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID)                    # predicted bag labels
#' predict(milr_result, testData$X, testData$ID, type = "instance") # predicted instance labels
#' 
#' # use BIC to choose penalty
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID,
#'   exp(seq(log(0.01), log(50), length = 30)))
#' coef(milr_result)      # coefficients
#' fitted(milr_result)                    # fitted bag labels
#' fitted(milr_result, type = "instance") # fitted instance labels
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID)                    # predicted bag labels
#' predict(milr_result, testData$X, testData$ID, type = "instance") # predicted instance labels
#' 
#' # use auto-tuning
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID, lambda = -1, numLambda = 20)
#' coef(milr_result)      # coefficients
#' fitted(milr_result)                    # fitted bag labels
#' fitted(milr_result, type = "instance") # fitted instance labels
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID)                    # predicted bag labels
#' predict(milr_result, testData$X, testData$ID, type = "instance") # predicted instance labels
#' 
#' # use cv in auto-tuning
#' milr_result <- milr(trainData$Z, trainData$X, trainData$ID, 
#'                     lambda = -1, numLambda = 20, lambdaCriterion = "deviance")
#' coef(milr_result)      # coefficients
#' fitted(milr_result)                    # fitted bag labels
#' fitted(milr_result, type = "instance") # fitted instance labels
#' summary(milr_result)   # summary milr
#' predict(milr_result, testData$X, testData$ID)                    # predicted bag labels
#' predict(milr_result, testData$X, testData$ID, type = "instance") # predicted instance labels
#' @importFrom numDeriv hessian
#' @importFrom glmnet glmnet
#' @name milr
#' @rdname milr
#' @export
milr <- function(y, x, bag, lambda = 0, numLambda = 20L, lambdaCriterion = "BIC", nfold = 10L, maxit = 1000L, tolerance = 1e-5) {
  # if x is vector, transform it to matrix
  if (is.vector(x))
    x <- matrix(x, ncol = 1)
  if (!is.matrix(x))
    x <- as.matrix(x)
  # if column names of x is missing, assign xi
  if (is.null(colnames(x)))
    colnames(x) <- paste0("x", 1L:ncol(x))
  if (!all(y %in% c(0, 1)))
    stop("y must be 0 and 1.")
  bag <- as.integer(factor(bag))
  # input check
  alpha <- 1
  assert_that(length(unique(y)) == 2, length(y) == nrow(x),
              all(is.finite(y)), is.numeric(y), all(is.finite(x)), is.numeric(x), 
              length(lambda) >= 1, all(is.finite(lambda)), is.numeric(lambda), 
              is.finite(alpha), is.numeric(alpha), is.finite(maxit), is.numeric(maxit), 
              abs(maxit - floor(maxit)) < 1e-4, nfold < nrow(x),
              lambdaCriterion %in% c("deviance", "BIC"))
  
  if (length(lambda) == 1 && lambda == -1) {
    message("The penalty term is selected automatically with ", numLambda, " candidates.")
    m <- table(bag)
    zi <- tapply(y, bag, function(x) sum(x) > 0) %>>% as.numeric
    lambdaMax <- sqrt(sum(m-1)) * sqrt(sum(m**(1-2*zi)))
    lambda <- exp(seq(log(lambdaMax/1000), log(lambdaMax), length = numLambda))
  } else if (length(lambda) == 1 && lambda == 0) {
    message("Lasso-penalty is not used.")
  } else {
    message("Use the user-defined lambda vector.")
  }
  
  # initial value for coefficients
  #init_beta <- coef(glm(y~x))
  init_beta <- glmnet(x, y, standardize = T, alpha = 0, lambda = lambda[1]) %>>% coef %>>% as.vector
  beta_history <- matrix(NA, ncol(x) + 1, length(lambda) + 1)
  beta_history[ , 1] <- init_beta
  unique_bag <- unique(bag)
  n_bag <- length(unique_bag)
	cv <- NULL
  if (length(lambda) > 1) {
    # save the history of betas along lambda
    message(sprintf("Using the criterion %s to choose optimized penalty.", lambdaCriterion))

    if (lambdaCriterion == "BIC") {
      dev <- vector("numeric", length(lambda))
      BIC <- vector("numeric", length(lambda))
      for (i in seq_along(lambda)) {
        beta_history[, i + 1] <- milr_cpp(y, cbind(1, x), bag, beta_history[, i], lambda[i], tolerance, alpha, maxit)
        dev[i] <- -2 * getLogLikMilr(beta_history[, i + 1], y, cbind(1, x), bag)
        BIC[i] <- dev[i] + sum(beta_history[, i + 1] != 0) * log(n_bag)
      }
      locMin <- which.min(BIC)
    } else if (lambdaCriterion == "deviance") {
      # calculate the bag label
      bagLvl <- tapply(y, bag, function(x) any(x > 0))
      # avoid the sample size of one class is less than nfold to retrun error in cvIndex_f.
      if (nfold > min(table(bagLvl))) {
        message("nfold is bigger than the number of size of one class.")
        message(sprintf("Therefore, nfold is changed to %i.", min(table(bagLvl))))
        nfold <- min(table(bagLvl))
      }
      # stratified sampling for cv on bag labels
      cvIndex_bag_p <- cvIndex_f(sum(bagLvl == 1), nfold)
      cvIndex_bag_n <- cvIndex_f(sum(bagLvl == 0), nfold)
      # convert the cv index on bag labels to cv index on samples
      cvIndex_sample <- lapply(1:nfold, function(x){
        tmp <- c(which(bagLvl == 1)[which(cvIndex_bag_p == x)], which(bagLvl == 0)[which(cvIndex_bag_n == x)])
        which(bag %in% unique_bag[tmp])
      })
      
      dev_cv <- matrix(NA, length(lambda), nfold)
      dev <- vector("numeric", length(lambda))
      BIC <- vector("numeric", length(lambda))
      cv_betas <- replicate(nfold, beta_history, simplify = FALSE)
      for (i in seq_along(lambda)) {
        # CV to calculate the deviances
        for (j in 1:nfold) {
          trainSetIndex <- do.call(c, cvIndex_sample[setdiff(1:nfold, j)])
          cv_betas[[j]][ , i + 1] <- milr_cpp(y[trainSetIndex], cbind(1, x[trainSetIndex, ]), 
                                              bag[trainSetIndex], cv_betas[[j]][ , i], lambda[i], 
                                              tolerance, alpha, maxit)
          dev_cv[i, j] <- -2 * getLogLikMilr(cv_betas[[j]][ , i + 1], y[cvIndex_sample[[j]]], 
                                      cbind(1, x[cvIndex_sample[[j]], ]), bag[cvIndex_sample[[j]]])
        }
        # calculate the beta under lambda
        beta_history[, i + 1] <- milr_cpp(y, cbind(1, x), bag, beta_history[, i], lambda[i], tolerance, alpha, maxit)
        dev[i] <- -2 * getLogLikMilr(beta_history[, i + 1], y, cbind(1, x), bag)
        BIC[i] <- dev[i] + sum(beta_history[, i + 1] != 0) * log(n_bag)
      }
			cv <- rowMeans(dev_cv)
      locMin <- which.min(cv)
    }
    lambda_out <- lambda[locMin]
    message(sprintf("The chosen penalty throuth %s is %.4f.", lambdaCriterion, lambda[locMin]))
    beta <- beta_history[ , locMin + 1]
  } else {
    beta_history[,2] <- milr_cpp(y, cbind(1, x), bag, init_beta, lambda, tolerance, alpha, maxit)  
		beta <- beta_history[,2]
    dev <- -2 * getLogLikMilr(beta, y, cbind(1, x), bag)
    BIC <- dev + sum(beta != 0) * log(n_bag)
    lambda_out <- lambda
		locMin <- 1
  }
  # return beta with names
  beta <- as.vector(beta) %>>% `names<-`(c("intercept", colnames(x)))
  rownames(beta_history) <- c("intercept", colnames(x))
	
  # fitted response
	fit_yij <- (beta %>>% (logit(cbind(1, x), .)) > 0.5) %>>% as.numeric(.)
  fit_y <- getMilrProb(beta, cbind(1, x), bag) %>>% `>`(0.5) %>>% as.numeric
  
  # calculate the hessian matrix when without using lasso
  if (lambda_out == 0)
    bateVar <- -diag(solve(hessian(function(b) getLogLikMilr(b, y, cbind(1, x), bag), beta)))
  else
    bateVar <- NULL
  out <- list(lambda = lambda, cv = cv, deviance = dev, BIC = BIC, 
							beta = beta_history[, 2L:ncol(beta_history)],
							best_index = locMin,
							best_model = list(lambda_chosen = lambda_out, 
																deviance = dev[locMin], 
																BIC = BIC[locMin],
																coeffiecents = beta, 
																fitted = list(bag = fit_y, instance = fit_yij), 
																loglik = getLogLikMilr(beta, y, cbind(1, x), bag),
																var = bateVar)
							)
  class(out) <- "milr"
  return(out)
}