#*************************************************************************Description*************************************************************************************************#
#Residuals of beta regression fit under MLE and SMLE
#***ARGUMENTS***#
# y - response variable.
# X -  regressor matrix for the mean submodel.
# Z -  regressor matrix for the precision submodel.
# theta - vector of parameter estimates (via SMLE or MLE) of the observed sample. 
# linkmu - character specification of the link function in the mean submodel. Currently, "logit", "probit", "cloglog", "cauchit", "log", "loglog" are supported. Default is "logit".
# linkphi - character specification of the link function in the precision submodel. Currently, "identity", "log", "sqrt" are supported. Default is "log".
#**************************************************************************************************************************************************************************************#

residuals_beta <- function(y, X, Z, theta, linkmu="logit", linkphi="log") { 
kk1 <- ncol(X)
kk2 <- ncol(Z)
n <- nrow(X)
beta_p <- theta[1:kk1]
gama_p <- theta[(kk1+1.0):(kk1+kk2)]    	   
etahat <- as.vector(X%*%beta_p)	
deltahat <- as.vector(Z%*%gama_p) 		
#***mean link functions***#
if(linkmu == "logit"){
muhat <- exp(etahat)/(1.0+exp(etahat))
T_1hat <- diag((muhat*(1.0-muhat)) )
}
if(linkmu == "probit"){
muhat <- pnorm(etahat) 
T_1hat <- diag(dnorm(qnorm(muhat)) )
}
if(linkmu == "cloglog"){
muhat <- 1.0 - exp(-exp(etahat)) 
T_1hat <- diag((muhat - 1)*log(1 - muhat))
}
if(linkmu == "log"){
muhat <- exp(etahat) #funcao ligacao log
T_1hat <- diag(muhat);
}
if(linkmu == "loglog"){
muhat <- exp(-exp(etahat)) 
T_1hat <- diag(muhat*log(muhat))
}
if(linkmu == "cauchit"){
muhat <- (pi^(-1))*atan(etahat) + 0.5 
T_1hat <- diag((1/pi)*((cospi(muhat-0.5))^2))                                 
}
#***precision link functions***#
if(linkphi == "log"){
phihat <- exp(deltahat) 
T_2hat <- diag(phihat)
}
if(linkphi == "identify"){
phihat <- deltahat 
T_2hat <- diag(rep(1,n))
}
if(linkphi == "sqrt"){
phihat <- deltahat^2 
T_2hat <- diag(2*sqrt(phihat))
}
PhiM <- diag(phihat)
mustarhat <- psigamma(muhat*phihat, 0) - psigamma((1.0-muhat)*phihat, 0)	
psi1hat <- psigamma(muhat*phihat, 1.0) 
psi2hat <- psigamma((1.0-muhat)*phihat, 1.0)  
What <- diag(phihat*(psi1hat+psi2hat))%*%(T_1hat^2.0)
tempinvhat <- solve(t(X)%*%PhiM%*%What%*%X)
Hhat <- diag(((What%*%PhiM)^(0.5))%*%X%*%tempinvhat%*%t(X)*((PhiM%*%What)^(0.5)))
vhat <- psi1hat + psi2hat		 
r_p1 <- (log(y/(1-y)) - mustarhat)/sqrt(vhat) 	
r_p2_aux <- 1.0/sqrt(1.0- Hhat)
r_p2 <- r_p1*r_p2_aux 
return(r_p2)
}#ends function

