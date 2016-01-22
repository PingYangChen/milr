
library(numDeriv)
library(plyr)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(snowfall); sfInit(parallel=TRUE, cpus=6, type='SOCK')

sfLibrary(numDeriv);sfLibrary(plyr);sfLibrary(dplyr);
sfLibrary(Rcpp);sfLibrary(RcppArmadillo);

setwd("D:\\Ping_Yang\\Google Drive\\PYChen_Statistics_NCKU\\Researches\\2014 Project ITRI 2014Nov\\R code\\package-milr\\v1")

sourceCpp("CLR_lasso.cpp")
source("milr.r");source("CBR_EM.r");source("DGP.r")

sfSource("milr.r");sfSource("CBR_EM.r");sfSource("DGP.r")

nSim <- 100
#diffAB <- numeric(nSim)
sfExport('nSim')

diffAB <- sapply(1:nSim, function(k) {
#diffAB <- sfSapply(1:nSim, function(k) {
	pL <- round(runif(1, 5, 21))
	testData <- DGP(50,3,runif(pL, -5,5))
	setwd("D:\\Ping_Yang\\Google Drive\\PYChen_Statistics_NCKU\\Researches\\2014 Project ITRI 2014Nov\\R code\\package-milr\\v1")
	sourceCpp("CLR_lasso.cpp")
	a <- milr(testData$Z, testData$X, testData$ID)$coeffiecents
	#a <- 0
	bfit <- CBR_EM(testData$Z, testData$X, testData$ID)
	if (length(bfit) < 4) {
		return(NA)
	} else {
		b <- bfit$CMLE
		return(mean(a-b)^2)
	}
})

sfStop()
diffAB

for (i in 1:nSim) {
	pL <- round(runif(1, 5, 21))
	testData <- DGP(50,3,runif(pL, -5,5))
	#sourceCpp("CLR_lasso.cpp")
	a <- milr(testData$Z, testData$X, testData$ID)$coeffiecents
	#a <- 0
	bfit <- CBR_EM(testData$Z, testData$X, testData$ID)
	if (length(bfit) < 4) {
		diffAB[i] <- NA
	} else {
		b <- bfit$CMLE
		diffAB[i] <- mean(a-b)^2
	}
}
