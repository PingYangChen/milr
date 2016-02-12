#' Benchmark dataset for multiple instance experiments
#' 
#' The MUSK data sets are the most widely used benchmark data sets when comparing 
#' different MI learning methods. The MUSK datasets consist of conformations. Each 
#' conformation is represented by 166 features. In MUSK1, the average conformation
#' in one bag is 5. For more detailed descriptions, please refer to the paper in 'Source'.
#' 
#' @format A data.frame with 476 observations on 168 variables.
#' \itemize{
#'   \item \code{bag}: bag id
#'   \item \code{label}: label
#'   \item \code{x1 ~ x166}: coveriates
#' }
#'
#' @source T. G. Dietterich, R. H. Lathrop, and T. Lozano-P`erez. (1997) Solving the multiple instance problem with axis-parallel rectangles. Artificial Intelligence, 89, 31 - 71.
#' @docType data
#' @name MUSK1
#' @usage MUSK1
NULL
