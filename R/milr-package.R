#' The milr package: multiple-instance logistic regression with lasso penalty
#'
#' The multiple instance data set consists of many independent subjects (called bags) 
#' and each subject is composed of several components (called instances). The outcomes 
#' of such data set are binary or multinomial, and, we can only observe the subject-level 
#' outcomes. For example, in manufactory processes, a subject is labeled as "defective" 
#' if at least one of its own components is defective, and otherwise, is labeled as 
#' "non-defective".  The milr package focuses on the predictive model for the multiple 
#' instance data set with binary outcomes and performs the maximum likelihood estimation 
#' with the Expectation-Maximization algorithm under the framework of logistic regression.  
#' Moreover, the LASSO penalty is attached to the likelihood function for simultaneous parameter 
#' estimation and variable selection.
#'
#' @references 
#' \enumerate{
#'   \item Chen, R.-B., Cheng, K.-H., Chang, S.-M., Jeng, S.-L., Chen, P.-Y., Yang, C.-H., 
#'  and Hsia, C.-C. (2016). Multiple-Instance Logistic Regression with LASSO Penalty. arXiv:1607.03615 [stat.ML].
#' }
#'
#' @docType package
#' @name milr-package
#' @useDynLib milr
#' @import assertthat
#' @importFrom Rcpp cppFunction sourceCpp
#' @importFrom pipeR %>>%
#' @importFrom utils globalVariables 
#' @importFrom stats coef glm optim pnorm printCoefmat rbinom rnorm
utils::globalVariables(c(".", "%>>%"))
