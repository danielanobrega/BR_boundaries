#*************************************************************************Description*************************************************************************************************#
#Create normal probability plots of residuals with simulated envelope to assess the goodness of beta regression fit under MLE and SMLE
#***ARGUMENTS***#
# y - response variable.
# X -  regressor matrix for the mean submodel.
# Z -  regressor matrix for the precision submodel.
# theta - vector of parameter estimates (via SMLE or MLE) of the observed sample. 
# linkmu - character specification of the link function in the mean submodel. Currently, "logit", "probit", "cloglog", "cauchit", "log", "loglog" are supported. Default is "logit".
# linkphi - character specification of the link function in the precision submodel. Currently, "identity", "log", "sqrt" are supported. Default is "log".
# SMLE - logical; if TRUE, the SMLE is used to create the plot; if FALSE, the MLE is used to create the plot. Default is TRUE.
# main.title - main title for the plot. Default is "Envelope".
# faixa.fixed - range of residuals values (optional). Default is NULL.
# labels.fixed - labels of the observations used to create the plot (optional). Default is NULL meaning that all observations are used.
#**************************************************************************************************************************************************************************************#

envelope_SMLE <- function(y, X, Z, theta, linkmu="logit", linkphi="log", SMLE=T, main.title = "Envelope", faixa.fixed = NULL, labels.fixed = NULL) { 
source("Resfunction.r") 
source("SMLE.r") 
B <- 100; #number of replicates
kk1 <- ncol(X); kk2 <- ncol(Z); n <- nrow(X)
#***parameters for parametric bootstrap***#
beta_p <- theta[1:kk1]
gama_p <- theta[(kk1+1.0):(kk1+kk2)]    	   
etahat <- X%*%beta_p
deltahat <- Z%*%gama_p 	

#***************link functions for mean submodel**********************#
if(linkmu == "logit") muhat <- exp(etahat)/(1.0+exp(etahat))
if(linkmu == "probit") muhat <- pnorm(etahat) 
if(linkmu == "cloglog") muhat <- 1.0 - exp(-exp(etahat)) 
if(linkmu == "log") muhat <- exp(etahat) 
if(linkmu == "loglog") muhat <- exp(-exp(etahat)) 
if(linkmu == "cauchit") muhat <- (pi^(-1))*atan(etahat) + 0.5 
#************************************************************#
#***************link functions for precision submodel**********************#
if(linkphi == "log") phihat <- exp(deltahat) 
if(linkphi == "identify") phihat <- deltahat 
if(linkphi == "sqrt") phihat <- deltahat^2
#************************************************************#

Menvelope_rp2 <- matrix(numeric(0),nrow=n,ncol=B)

#------------> residuals for the observed sample<--------------#
RP2 <- residuals_beta(y, X, Z, c(beta_p,gama_p), linkmu= linkmu, linkphi= linkphi)
set.seed(c(1994,1991), kind="Marsaglia-Multicarry")

    for(j in 1:B){		
        ygen <- rbeta(n, muhat*phihat, (1.0-muhat)*phihat)
        while(any(round(ygen,5)==0|round(ygen,5)==1)){
        ygen <- rbeta(n, muhat*phihat, (1.0-muhat)*phihat)
    }
  if(SMLE==T){
     fit <- SMLE_BETA(y=ygen, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac =0.02, method="BFGS", startV="CP", linkmu=linkmu, linkphi=linkphi)
  }
  else{
     fit <- SMLE_BETA(y=ygen, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu=linkmu, linkphi=linkphi)
  }

RP2_b <- residuals_beta(ygen,X,Z,c(fit$beta,fit$gama), linkmu=linkmu, linkphi=linkphi)
Menvelope_rp2[,j] = RP2_b
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
v <- qqnorm(res_rp2, main=main.title, xlab="Normal quantiles", ylab="Residuals", ylim=faixa, pch=16, cex=1, cex.lab=2.0, cex.axis=2, cex.main=2.0)
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
