#' The milr package: a implementation for multiple-instance logistic regression with lasso penalty
#'
#' Description here
#'
#' @section Reference:
#' \enumerate{
#'   \item to be filled in.
#' }
#'
#' @docType package
#' @name milr
#' @useDynLib milr
#' @import assertthat
#' @importFrom Rcpp cppFunction sourceCpp
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom utils globalVariables
utils::globalVariables(c(".", "%>%", "%<>%"))
