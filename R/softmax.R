#' softmax function:
#'
#' description here
#' 
#' @param y A vector, its elements contain two different values.
#' @param x The design matrix. The number of rows of x must be equal to the length of y.
#' @param bag bag variable.
#' @param alpha 123
#' @param maxit An integer, the maximum iteration for EM algorithm.
#' @return An list includes coefficients and fitted values.
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
#' @importFrom magrittr set_names
#' @importFrom purrr map map_int map2_dbl
#' @importFrom logistf logistf
#' @rdname milr_main
#' @export
softmax <- function(y, x, bag, alpha = 0, maxit = 500) {
  # if x is vector, transform it to matrix
  if (is.vector(x))
    x <- matrix(x, ncol = 1)
  # input check
  assert_that(length(unique(y)) == 2, length(y) == nrow(x),
              all(is.finite(y)), is.numeric(y), all(is.finite(x)), is.numeric(x),  
              is.finite(alpha), is.numeric(alpha), is.finite(maxit), is.numeric(maxit), 
              abs(maxit - floor(maxit)) < 1e-4)
  
  # initial value for coefficients
  init_beta <- coef(logistf(y~x))
	
	n_bag <- length(unique(bag))
  nloglik <- function(b, Z, X, ID){
    pii <- (split(as.data.frame(X), ID) %>% 
			purrr::map(~logit(as.matrix(.), b)*exp(alpha*logit(as.matrix(.), b))) %>% sum(na.rm = TRUE) )/
      (split(as.data.frame(X), ID) %>% 
				purrr::map(~exp(alpha*logit(as.matrix(.), b))) %>% sum(na.rm = TRUE) )
    return(-1*(purrr::map2_dbl(split(Z, ID) %>% purrr::map(unique), pii, 
      ~.x*log(.y)+(1-.x)*log(1-.y)) %>% sum(na.rm = TRUE)))
  }
	cat("Using auto-tuning to choose panelty.\n")
	beta <- optim(init_beta, nloglik(b, y, cbind(1, x), bag), control = list(maxit = maxit))
	
  beta %<>% as.vector %>% set_names(c("intercept", colnames(x)))
  fit_y <- beta %>% {split(logit(cbind(1, x), .) > 0.5, bag)} %>%
    purrr::map_int(~ifelse(sum(.) > 1, 1L, 0L))
  out <- list(coeffiecents = beta, fitted = fit_y)
  class(out) <- 'milr'
  return(out)
}
