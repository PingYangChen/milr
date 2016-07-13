#' The milr package: an implementation for multiple-instance logistic regression with lasso penalty
#'
#' This package performs maximum likelihood estimation for multiple-instance logistic regression  
#' utilizing the EM algorithm. Also, the LASSO penalty is implemented for variable selection. 
#'
#' @references 
#' \enumerate{
#'   \item Ray-Bing Chen, Kuang-Hung Cheng, Sheng-Mao Chang, Shuen-Lin Jeng, 
#'    Ping-Yang Chen, Chun-Hao Yang, and Chi-Chun Hsia. (2016+) Multiple-Instance Logistic Regression with LASSO Penalty. (Submitted)
#' }
#'
#' @docType package
#' @name milr-package
#' @useDynLib milr
#' @import assertthat
#' @importFrom Rcpp cppFunction sourceCpp
#' @importFrom pipeR %>>%
#' @importFrom utils globalVariables
utils::globalVariables(c(".", "%>>%"))
