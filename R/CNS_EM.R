#' CBR_FS function
#'
#' description here
#' 
#' @param Z 123
#' @param X 123
#' @param ID 123
#' @return An list involves coefficients and fitted values.
#' @examples
#' # testData <- DGP(400, 10, runif(12, -5, 5))
#' # a <- CBR_FS(testData$Z, testData$X, testData$ID)
#' @importFrom numDeriv hessian
#' @export
CBR_FS <- function(Z, X, ID){
  n <- nrow(X)
  p <- ncol(X)
  d <- names(table(ID))
  n.d <- length(d)
  
  out <- summary(glm(Z~X-1, family=binomial(link="logit")))$coefficients
  temp.Lp <- out[,4]
  temp.L <- out[,1]
  #b <- temp.L
  b <- c(-2,1,-1,0 )
  ## i: iteration for Fisher Scoring
  for(i in 1:20){
    # i2: iteration for experiment unit
    P <- exp(X%*%b); P <- P/(1+P)
    Xbar <- matrix(0,n,p) 
    Pi <- rep(1, n)
    z <- rep(0, n)
    for(i2 in 1:n.d)
    {
      sel <- which(ID == d[i2])
      for(j in 1:length(sel))
      {
        Xbar[i2,] <- Xbar[i2,] + X[sel[j],]*P[sel[j]] 
        Pi[i2] <- Pi[i2]*(1 - 1/(1+exp(-sum(X[sel[j],]*b)))) 	
      }
      Pi[i2] <- 1 - Pi[i2]
      z[i2] <- Z[sel[1]]
    }
    O <- diag((1-Pi)/Pi)
    D <- diag(1/Pi)
    V <- solve(t(Xbar)%*%O%*%Xbar)
    
    b.new <- b + V%*%t(Xbar)%*%D%*%(z-Pi)
    sse <- crossprod(b.new - b)
    b <- b.new
    if(sse < 1e-6)
      break
  }
  
  test <- -abs(b/sqrt(diag(V)))
  Lc.p <- 2*pnorm(test)
  return(list(L.b=temp.L, L.p=temp.Lp, Lc.b=b, Lc.p=Lc.p))
}

#' CBR_EM function:
#'
#' description here
#' 
#' @param Z 123
#' @param X 123
#' @param ID 123
#' @return An list involves coefficients and fitted values.
#' @examples
#' # testData <- DGP(400, 10, runif(12, -5, 5))
#' # a <- CBR_EM(testData$Z, testData$X, testData$ID)
#' @export
CBR_EM <- function(Z, X, ID){
  ## log-likelihood of z
  ## used for variance calculation
  
  loglik <- function(b,Z,X,ID){
    pii <- split(as.data.frame(X), ID) %>% 
      map(~1-prod(1-logit(as.matrix(.), b)))
    return(map2_dbl(split(Z, ID) %>% map(unique), pii, 
                    ~.x*log(.y)+(1-.x)*log(1-.y)) %>% sum)
  }
  
  ## solving logistic-like regression with MM algorithm
  MMLogit <- function(y, X, b)
  {
    for(i in 1:200)
    {
      pv <- exp(X%*%b) %>% {./(1+.)}
      b.new <- b + t(X) %*% (y - pv) / XXi
      sse <- crossprod(b.new-b)
      b <- b.new
      if(sse < 10^(-6))
        break
    }
    if(i > 199) cat("Do not converge in the inner loop after 200 iteration!\n")
    return(b)
  }
  
  ## prepare data
  d <- names(table(ID))
  n.d <- length(d)
  
  out <- summary(glm(Z~X-1, family=binomial(link="logit")))$coefficients
  temp.Lp <- out[,4]
  temp.L <- out[,1]
  b <- temp.L
  XXi <- .25*(nrow(X)-1)
  # EM iterations
  for(i in 1:500)
  {
    pv <- exp(X%*%b) %>% {./(1+.)}
    qv <- pv
    for(i2 in 1:n.d)
    {
      sel <- which(ID == d[i2])
      temp <- 1- prod(1-pv[sel])
      qv[sel] <- qv[sel] / temp
    }
    y <- Z*qv
    b.new <- MMLogit(y, X, b)
    
    sse <- crossprod(b.new-b)
    b <- b.new
    if (sse < 10^(-6)) break
  }
  if (i > 499)
    cat("Do not converge after 500 iterations!\n")
  
  b.var <- -diag(solve(hessian(function(b) loglik(b,Z,X,ID), b)))
  test <- -abs(b/sqrt(b.var))
  b.p <- 2*pnorm(test)
  return(list(NMLE=temp.L, NTest=temp.Lp, CMLE=b, CTest=b.p))
}
