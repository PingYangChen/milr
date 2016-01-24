#' CBR_FS function
#'
#' description here
#' 
#' @param Z 123
#' @param X 123
#' @param ID 123
#' @return An list involves coefficients and fitted values.
#' @examples
#' testData <- DGP(50, 3, runif(sample(5:21, 1), -5, 5))
#' a <- milr(testData$Z, testData$X, testData$ID)
#' coef(a)
#' @importFrom numDeriv hessian
#' @export
CBR_FS <- function(Z, X, ID){
  n <- nrow(X)
  p <- ncol(X)
  d <- names(table(ID))
  n.d <- length(d)
  
  out <- summary(glm(A$Z~A$X-1, family=binomial(link="logit")))$coefficients
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
      n.sel <- length(sel)
      for(j in 1:n.sel)
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
    test <- b.new - b; test <- sum(test*test)
    b <- b.new
    if(test < 1e-6) break
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
#' testData <- DGP(50, 3, runif(sample(5:21, 1), -5, 5))
#' a <- milr(testData$Z, testData$X, testData$ID)
#' coef(a)
#' @export
CBR_EM <- function(Z, X, ID){
  ## log-likelihood of z
  ## used for variance calculation
  loglik <- function(b)
  {
    out <- 0
    for(i in 1:n.d)
    {
      sel <- which(ID == d[i])
      n.sel <- length(sel)
      temp <- 1
      for(j in 1:n.sel)
        temp <- temp*(1 - 1/(1+exp(-sum(X[sel[j],]*b)))) 	
      pii <- 1 - temp
      out <- out + Z[sel[1]]*log(pii) + (1-Z[sel[1]])*log(1-pii) 
      #cat(i, out,"\n")
    }
    return(out)
  }
  
  ## solving logistic-like regression with MM algorithm
  MMLogit <- function(y, X, b)
  {
    for(i in 1:200)
    {
      pv <- exp(X%*%b); pv <- pv/(1+pv)
      b.new <- b + XXi%*%t(X)%*%(y - pv)
      test <- b.new - b; test <- sum(test*test)
      b <- b.new
      if(test < 10^(-6)) break
    }
    if(i > 199) cat("Do not converge in the inner loop after 200 iteration!\n")
    return(b)
  }
  
  ## prepare data
  d <- names(table(ID))
  n.d <- length(d)
  
  out <- summary(glm(A$Z~A$X-1, family=binomial(link="logit")))$coefficients
  temp.Lp <- out[,4]
  temp.L <- out[,1]
  b <- temp.L
  
  XXi <- solve(t(X)%*%X/4)
  
  # EM iterations
  for(i in 1:500)
  {
    pv <- exp(X%*%b); pv <- pv/(1+pv)		
    qv <- pv
    for(i2 in 1:n.d){
      sel <- which(ID == d[i2])
      n.sel <- length(sel)
      temp <- 1- prod(1-pv[sel])
      qv[sel] <- qv[sel]/temp
    }
    y <- Z*qv
    b.new <- MMLogit(y, X, b)
    
    test <- b.new-b; test <- sum(test*test)
    b <- b.new
    if (test < 10^(-6)) break
  }
  if (i > 499)
    cat("Do not converge after 500 iterations!\n")
  
  b.var <- -diag(solve(hessian(loglik, b)))
  test <- -abs(b/sqrt(b.var))
  b.p <- 2*pnorm(test)
  return(list(NMLE=temp.L, NTest=temp.Lp, CMLE=b, CTest=b.p))
}
