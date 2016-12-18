#' @export
#' @method coef softmax
coef.softmax <- function(object, ...) {
  return(object$coeffiecents)
}

#' @export
#' @method fitted softmax
fitted.softmax <- function(object, ...) {
  return(object$fitted)
}

#' Predict Method for softmax Fits
#' 
#' @param object A fitted obejct of class inheriting from \code{"softmax"}.
#' @param newdata A matrix with variables to predict.
#' @param bag_newdata A vector. The labels of instances to bags.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method predict softmax
predict.softmax <- function(object, newdata, bag_newdata, ...) {
  return(getSoftmaxBag(cbind(1, newdata), coef(object), bag_newdata, object$alpha))
}

#' @export
#' @method print softmax
print.softmax <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (length(coef(x)) > 0) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  }
  cat("Log-Likelihood:   ", format(signif(x$loglik, digits)), "\n")
  invisible(x)
}

#' Multiple-instance logistic regression via softmax function
#'
#' This function calculates the alternative maximum likelihood estimation for 
#' multiple-instance logistic regression
#' through a softmax function (Xu and Frank, 2004; Ray and Craven, 2005).
#'
#' @param y a vector. Bag-level binary labels.
#' @param x the design matrix. The number of rows of \code{x} must be equal to the length of \code{y}.
#' @param bag a vector, bag id.
#' @param alpha A non-negative realnumber, the softmax parameter. 
#' @param ... arguments to be passed to the \code{optim} function.
#' @return a list including coefficients and fitted values.
#' @examples
#' set.seed(100)
#' beta <- runif(10, -5, 5)
#' trainData <- DGP(70, 3, beta)
#' testData <- DGP(30, 3, beta)
#' # Fit softmax-MILR model S(0)
#' softmax_result <- softmax(trainData$Z, trainData$X, trainData$ID, alpha = 0)
#' coef(softmax_result)      # coefficients
#' fitted(softmax_result)    # fitted values
#' predict(softmax_result, testData$X, testData$ID) # predicted label
#' # Fit softmax-MILR model S(3)
#' softmax_result <- softmax(trainData$Z, trainData$X, trainData$ID, alpha = 3)
#' @references
#' \enumerate{
#'	 \item S. Ray, and M. Craven. (2005) Supervised versus multiple instance learning: 
#'	 An empirical comparsion. in Proceedings of the 22nd International Conference on 
#'	 Machine Learnings, ACM, 697--704.
#'	 \item X. Xu, and E. Frank. (2004) Logistic regression and boosting for labeled bags 
#'	 of instances. in Advances in Knowledge Discovery and Data Mining, Springer, 272--281.
#' }
#' @export
softmax <- function(y, x, bag, alpha = 0, ...) {
  # if x is vector, transform it to matrix
  if (is.vector(x))
    x <- matrix(x, ncol = 1)
  if (!is.matrix(x))
    x <- as.matrix(x)
  # if column names of x is missing, assign xi
  if (is.null(colnames(x)))
    colnames(x) <- paste0("x", 1L:ncol(x))
  if (!all(y %in% c(0, 1)))
    stop('y must be 0 and 1.')
  # input check
  assert_that(length(unique(y)) == 2L, length(y) == nrow(x),
              all(is.finite(y)), is.numeric(y), all(is.finite(x)), is.numeric(x),  
              alpha >= 0, is.finite(alpha), is.numeric(alpha))
  
  # initial value for coefficients
  init_beta <- coef(glm(y ~ x))
  # find the bag response
  y_bag <- tapply(y, bag, function(z) any(z > 0)) %>>% as.integer
  bagTmp <- as.integer(as.factor(bag))
  # optimize coefficients
  beta <- optim(par = init_beta, fn = function(b){
    softmaxlogL(bagTmp, cbind(1, x), y_bag, b, alpha)
  }, ...)$par %>>% `names<-`(c("intercept", colnames(x)))
  
  # get fitted bag response
  fit_y <- getSoftmaxBag(cbind(1, x), beta, bag, alpha)
  out <- structure(list(alpha = alpha, coeffiecents = beta, fitted = fit_y, 
                        loglik = -softmaxlogL(bagTmp, cbind(1, x), y_bag, beta, alpha)), 
                   class = "softmax")
  return(out)
}
