###############################Function source code##########################
#File with functions that are used in the code for the analyses of the datasets
#############################################################################




##Transformation - Smithson and Verkuilen (2006)
#Transform response variable with the linear transformation
#data - response variable
tsv <- function(data) {
  n <- length(data)
  tsv_data <- ((n-1)*data + 0.5)/n
  return(tsv_data)
}

##Replacing by smallest/largest observations in (0,1)
#data - response variable
rep_mm <- function(data) {
  rep_data <- data  
  for(i in 1:length(data)){
    if (data[i] == 0) rep_data[i] <- min(data[data>0])
    if (data[i] == 1) rep_data[i] <- max(data[data<1])
  }
  return(rep_data)
}

##Adding/subtracting epsilon>0
#data - response variable 
#e - epsilon>0 to be added or subtracted from boundary observations
rep_e <- function(data, e) {
  rep_data <- data  
  for(i in 1:length(data)){
    if (data[i] == 0) rep_data[i] <- data[i]+e
    if (data[i] == 1) rep_data[i] <- data[i]-e
  }
  return(rep_data)
}


##Removing boundary observations
#data - dataset with variables
#response - column number of the response variable

rem_bd <- function(data, response) {
  rem <- (data[,response]==0)|(data[,response]==1)
  new_data <- data[!rem,]
  return(new_data)
}


##Calculates mu_i for a matrix X of covariate values (first column a vector of 1's)
mu_mle <- function(x, Beta) {
  l <- dim(as.matrix(x))[1]
  x1 <- cbind(rep(1, l),x)
  Eta <- x1%*%Beta
  mu_estim <- exp(Eta)/(1+exp(Eta))
  return(mu_estim)
}


mu_mle2 <- function(x, Beta) {
  l <- dim(as.matrix(x))[1]
  x1 <- cbind(rep(1, l),x)
  Eta <- x1%*%Beta
  mu_estim <- exp(Eta)/(1+exp(Eta))
  return(as.vector(mu_estim))
}


#Randomized quantile residuals for (robust) beta regression models

rqbeta <- function(y, X, Z, beta, gamma) {
  eta1 <- X%*%beta
  eta2 <- Z%*%gamma
  mu   <- exp(eta1)/(1+exp(eta1))
  phi  <- exp(eta2)
  res <- numeric(length(y))
  for (i in 1:length(y)) {
    res[i] <- qnorm(pbeta(y[i], mu[i]*phi[i], (1-mu[i])*phi[i]))
  }
  return(res)
}


#*************************************************************************Description*************************************************************************************************#
#Create normal probability plots of residuals with simulated envelope to assess the goodness inflated beta regression fit
#***ARGUMENTS***#
# y - response variable.
# X -  regressor matrix for the mean submodel.
# Z -  regressor matrix for the precision submodel.
# type - "BEZI" if zero-inflated beta regression; "BEOI" for  one-inflated beta regression
# main.title - main title for the plot. Default is "Envelope".
# faixa.fixed - range of residuals values (optional). Default is NULL.
# labels.fixed - labels of the observations used to create the plot (optional). Default is NULL meaning that all observations are used.
#**************************************************************************************************************************************************************************************#

envelope_BEI <- function(object, type = NULL, main.title = "Envelope", faixa.fixed = NULL, labels.fixed = NULL) { 
  y <- object$y
  X <- object$mu.x
  Z <- object$sigma.x
  
  
  B <- 100; #number of replicates
  n <- nrow(X)
  beta_p  <- coefficients(object)
  gama_p  <- object$sigma.coefficients
  
  
  mu    <- exp(X%*%beta_p)/(1+exp(X%*%beta_p))
  phi   <- exp(Z%*%gama_p)
  alpha <- object$nu.coefficients
  
  RP2 <- residuals(object, what = "z-scores")
  #kk1 <- ncol(X); kk2 <- ncol(Z); n <- nrow(X)
  #***parameters for parametric bootstrap***#
  #beta_p <- theta[1:kk1]
  #gama_p <- theta[(kk1+1.0):(kk1+kk2)]    	   
  #etahat <- X%*%beta_p
  #deltahat <- Z%*%gama_p 	
  
  #mu <- exp(X%*%beta)/(1+exp(X%*%beta))
  #phi <- exp(Z%*%gama)
  #***************link functions for mean submodel**********************#
  #if(linkmu == "logit") muhat <- exp(etahat)/(1.0+exp(etahat))
  #if(linkmu == "probit") muhat <- pnorm(etahat) 
  #if(linkmu == "cloglog") muhat <- 1.0 - exp(-exp(etahat)) 
  #if(linkmu == "log") muhat <- exp(etahat) 
  #if(linkmu == "loglog") muhat <- exp(-exp(etahat)) 
  #if(linkmu == "cauchit") muhat <- (pi^(-1))*atan(etahat) + 0.5 
  #************************************************************#
  #***************link functions for precision submodel**********************#
  #if(linkphi == "log") phihat <- exp(deltahat) 
  #if(linkphi == "identify") phihat <- deltahat 
  #if(linkphi == "sqrt") phihat <- deltahat^2
  #************************************************************#
  
  Menvelope_rp2 <- matrix(numeric(0),nrow=n,ncol=B)
  
  #------------> residuals for the observed sample<--------------#
  
  if(type == "BEZI") {
    
    
    set.seed(c(1994,1991), kind="Marsaglia-Multicarry")
  
   for(j in 1:B){		
    ygen <- rBEZI(n, mu, phi, alpha)
    inf_fit_b <- gamlss(ygen ~ X[,-1], 
                        sigma.formula = ~Z[,-1],
                        family = BEZI(mu.link = "logit", sigma.link = "log",
                                      nu.link = "identity"), 
                        control=gamlss.control(trace=FALSE))
   RP2_b <- residuals(inf_fit_b, what = "z-scores")
    Menvelope_rp2[,j] = RP2_b
    }
  }
  
  else {
    
    set.seed(c(1994,1991), kind="Marsaglia-Multicarry")
    
    for(j in 1:B){		
      ygen <- rBEOI(n, mu, phi, alpha)
      inf_fit_b <- gamlss(ygen ~ X[,-1], 
                          sigma.formula = ~Z[,-1],
                          family = BEOI(mu.link = "logit", sigma.link = "log",
                                        nu.link = "identity"), 
                          control=gamlss.control(trace=FALSE))
      RP2_b <- residuals(inf_fit_b, what = "z-scores")
      Menvelope_rp2[,j] = RP2_b
    }
  }
  
  Menvelope_rp2 <- apply(Menvelope_rp2,2,sort);          
  res_rp2 <-    RP2;    
  res_min_rp2  <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.05)));         
  res_mean_rp2 <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.5)));                              
  res_max_rp2  <-    as.numeric(t(apply(t(Menvelope_rp2), 2,quantile, probs =0.95)));           
  faixa <- range(res_rp2,res_min_rp2,res_max_rp2)
  if(is.vector(faixa.fixed)) faixa <- faixa.fixed
  if(is.vector(labels.fixed)) labels <- labels.fixed
  par(mar=c(5.0,5.0,4.0,2.0))
  v <- qqnorm(res_rp2, main=main.title, xlab="Normal quantiles", ylab=expression(r[q]), ylim=faixa, pch=16, cex=1, cex.lab=2.0, cex.axis=2, cex.main=2.0)
  identify(v$x,v$y,labels,cex =1.3) #identify points in the plot
  #identify(v$x[c(15,16,72)],v$y[c(15,16,72)],cex=1.3,labels=c("15","16","72"), cex=1.3) #Only for the the firm cost data
  par(new=T)
  #
  qqnorm(res_min_rp2,axes=F,main = "",xlab="",ylab="",type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_max_rp2,axes=F,main = "",xlab="",ylab="", type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_mean_rp2,axes=F,xlab="",main = "", ylab="", type="l",ylim=faixa,lty=2,lwd=2.0)
}#ends function




