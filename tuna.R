#############################Tuna application##################################
#Application based on Ribeiro and Ferrari(2020)
#For more information check https://github.com/terezinharibeiro/RobustBetaRegression
###################################################################################

#Loading libraries
require(betareg)
require(gamlss)
require(mcglm)
require(ggplot2)

#Loading functions
source("functions.R")
source("SMLE.r") 
source("Resfunction.r")
source("envelope_function.R")

#Data
data <- read.table("dados_2000_Indian.txt",h=T)
head(data)
attach(data)


y <-  trop/100
ysul <- y[lat<0];
SSTsul <- SST[lat<0]

#Data w/o observation 46
data1 <- data[lat<0,]
data1$ysul <- data1$trop/100
data1_46 <- rem_bd(data1, 7)
names(data1_46)[6]<- "SSTsul"

##Transforming response variables
ysul_tsv <- tsv(ysul)  #Transforming using method proposed by Smithson and Verkuilen
ysul_0.01 <- rep_e(ysul, 0.01)  #Subtracting 0.01 from observation 46
ysul_0.001 <- rep_e(ysul, 0.001) #Subtracting 0.001 from observation 46
ysul_mm <- rep_mm(ysul) #Replacing by largest observation in (0,1)


####Beta regression

#Model estimation
mle_fit_tsv <- betareg(ysul_tsv ~ SSTsul|1) 
summary(mle_fit_tsv)


mle_fit_0.01 <- betareg(ysul_0.01 ~ SSTsul|1)
summary(mle_fit_0.01)


mle_fit_0.001 <- betareg(ysul_0.001 ~ SSTsul|1)
summary(mle_fit_0.001)


mle_fit_mm <- betareg(ysul_mm ~ SSTsul|1)
summary(mle_fit_mm)


mle_fit_bd <- betareg(data1_46$ysul ~ data1_46$SST|1)
summary(mle_fit_bd)


#Diagnostics 


par(mar = c(5.6, 4.6, 4.6, 2.6))
plot(mle_fit_tsv, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)
identify(f_values,pr_ql,lab,cex =1.3)

plot(mle_fit_tsv, which = 2, type = "sweighted2", cex.lab=1.5, cex.axis=1.5,
     caption = "", sub.caption = "")

par(mfrow=c(2,2))
plot(mle_fit_0.01, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)
plot(mle_fit_0.01, which = 2, type = "sweighted2", cex.lab=1.5, cex.axis=1.5,
     caption = "", sub.caption = "")

par(mfrow=c(2,2))
plot(mle_fit_0.001, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)
plot(mle_fit_0.001, which = 2, type = "sweighted2", cex.lab=1.5, cex.axis=1.5,
     caption = "", sub.caption = "")

par(mfrow=c(2,2))
plot(mle_fit_mm, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)

plot(mle_fit_mm, which = 2, type = "sweighted2", cex.lab=1.5, cex.axis=1.5,
     caption = "", sub.caption = "")

par(mfrow=c(2,2))
plot(mle_fit_bd, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)
plot(mle_fit_bd, which = 2, type = "sweighted2", cex.lab=1.5, cex.axis=1.5,
     caption = "", sub.caption = "")

####Robust beta regression

#Model matrices
X <- matrix(c(rep(1,77), SSTsul), ncol=2, byrow=F); 
Z <- matrix(c(rep(1,77)), ncol=1, byrow=F);
y46 <- ysul[-c(46)]; X46<- X[-c(46),];Z46 <- as.matrix(Z[-c(46),])

##Model estimation
fitMLE_tsv<- SMLE_BETA(y=ysul_tsv, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                       spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")
fitMLE_tsv


fitMLE_0.01<- SMLE_BETA(y=ysul_0.01, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                        spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")
fitMLE_0.01


fitMLE_0.001<- SMLE_BETA(y=ysul_0.001, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                         spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_0.001


fitMLE_mm<- SMLE_BETA(y=ysul_mm, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                      spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_mm

fitMLE_bd<- SMLE_BETA(y=y46, X=X46, Z=Z46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                      spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_bd


#Not robust (equivalent to betareg)

fitMLE1_tsv   <- SMLE_BETA(y=ysul_tsv, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                       spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")

fitMLE1_0.01  <- SMLE_BETA(y=ysul_0.01, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                        spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")
fitMLE1_0.001 <- SMLE_BETA(y=ysul_0.001, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                         spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE1_mm    <- SMLE_BETA(y=ysul_mm, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                      spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE1_bd    <- SMLE_BETA(y=y46, X=X46, Z=Z46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                      spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")



####Inflated beta regression

#Model estimation
beinf_fit <- gamlss(ysul ~ SSTsul, 
                    family = BEOI(mu.link = "logit", sigma.link = "log",
                                    nu.link = "identity"))
summary(beinf_fit)
beta_inf <- coefficients(beinf_fit)

#Diagnostics

wp(beinf_fit, line = FALSE, cex.lab = 2, cex.axis = 2)


####Quasi-likelihood (Bonat et al.(2019))

#Model estimation

#QL while estimating p
form1 <- ysul ~ SSTsul

Z0_1 <- mc_id(data1)
fit_mc1 <- mcglm(c(form1), list(Z0_1), link = "logit", 
                 variance = "binomialP", power_fixed = FALSE,
                 control_algorithm = list(tunning = 0.75, correct = FALSE, 
                                          verbose = FALSE),
                 data = data1)
summary(fit_mc1)
beta_ql <- coefficients(fit_mc1) [1:2,1]

#QL with fixed p
fit_mc1p <- mcglm(c(form1), list(Z0_1), link = "logit", 
                 variance = "binomialP", power_fixed = TRUE,
                 control_algorithm = list(tunning = 0.75, correct = FALSE, 
                                          verbose = FALSE),
                 data = data1)
summary(fit_mc1p)
form2 <- ysul ~ SSTsul

#QL without obs. 46 in the sample
Z0_2 <- mc_id(data1_46)
fit_mc2 <- mcglm(c(form2), list(Z0_2), link = "logit", 
                 variance = "binomialP", power_fixed = FALSE,
                 control_algorithm = list(tunning = 0.75, correct = FALSE),
                 data = data1_46)
summary(fit_mc2)
beta_ql2 <- coefficients(fit_mc2) [1:2,1]


#Diagnostics QL- Estimating p
 

pr_ql <- residuals(fit_mc1, type = "pearson")
pr_ql <- as.numeric(pr_ql)

f_values <- as.numeric(fitted(fit_mc1))

op <- par(mfrow = c(2,2))
par(op)
plot( f_values, pr_ql, ylab = "Pearson Residuals", xlab = "Fitted values", 
     pch=16, cex=1, cex.lab=1.5, cex.axis=1.5, ylim = c(-4, 7))
identify(f_values,pr_ql,lab,cex =1.3)
abline(a = 0, b = 0, lty = 3, col = "grey")


#Diagnostics Ql - p = 1
pr_qlp <- residuals(fit_mc1p, type = "pearson")
pr_qlp <- as.numeric(pr_qlp)

f_valuesp <- as.numeric(fitted(fit_mc1p))



lab <- 1:77
plot( f_valuesp, pr_qlp, ylab = "Pearson Residuals", xlab = "Fitted values", 
      pch=16, cex=1, cex.lab=1.5, cex.axis=1.5, ylim = c(-4, 8))
identify(f_valuesp,pr_qlp,lab,cex =1.3)
abline(a = 0, b = 0, lty = 3, col = "grey")




###Comparison of models

estim <- matrix(c(beta_tsv, beta_0.01, beta_0.001, beta_mm, beta_bd, beta_inf, beta_ql),
                ncol = 7, byrow = F)

values1 <- seq(15,30,length.out=1000)


mu_mtx <- mu_mle(values1, estim)

#Graph
par(mar=c(5.0,5.0,4.0,2.0))

par(mfrow=c(1,1))
plot(SSTsul,ysul,pch=16,xlab="SST",ylab="TTP",ylim=c(0,1),xlim=c(15,30),
     cex=1.5,cex.lab=1.0,cex.axis=1.5,cex.main=2.0)
lines(values1, mu_mtx[,1],       lwd=2,col=2,lty=1)
lines(values1, mu_mtx[,2],       lwd=2,col=3,lty=2)
lines(values1, mu_mtx[,3],       lwd=2,col=4,lty=3)
lines(values1, mu_mtx[,4],       lwd=2,col=5,lty=4)
lines(values1, mu_mtx[,5],       lwd=2,col=6,lty=5)
lines(values1, mu_mtx[,6],       lwd=2,col=7,lty=6)
lines(values1, mu_mtx[,7],       lwd=2,col=8,lty=7)
identify(SSTsul,ysul, cex=1.3)

legend(16,1,c("Transformation (Smithson & Verkuilen)","Subtracting 0.01 from obs. 46", 
                "Subtracting 0.001 from obs. 46", "Replacing by max.", 
              "Removing obs. 46", "One-inflated beta", "Quasi-likelihood"),
       col=2:8,lty=1:7,cex=0.75,lwd= rep(2,7))

##ggplot

d <- cbind(SSTsul, ysul)
d <- as.data.frame(d)
d1 <- cbind(values1, t(mu_mtx))
d1 <- as.data.frame(d1)
nm <- c("tsv", "epsilon0.01","epsilon0.001", "mm", "bd", "Beta_inf", "QL")
names(d1)[2:8] <- nm

data_gg <- data.frame(x = rep(values1, 7), y = as.vector(mu_mtx), 
                      method = rep(nm, each = 1000))
g <- ggplot(d, aes(SSTsul, ysul)) + geom_point()

    
g+geom_line(data = data_gg, aes(x = x,y = y, colour = method, linetype = method))


##Comparing quasi-likelihood with and without obs. 46

mu_ql_46 <- mu_mle(values1, beta_ql2)

par(mfrow = c(1,1))
plot(SSTsul,ysul,pch=16,xlab="SST",ylab="TTP",ylim=c(0,1),xlim=c(15,30),
     cex=1.5,cex.lab=1.0,cex.axis=1.5,cex.main=2.0)
lines(values1, mu_mtx[,7],       lwd=2,col=2,lty=1)
lines(values1, mu_ql_46,       lwd=2,col=3,lty=2)
legend(16,1,c("Quasi-likelihood","Quasi-likelihood w/o obs. 46"),
       col=2:3,lty=1:2,cex=0.75,lwd= rep(2,2))




##############################Diagnostics#####################################

#Residuals - Standard weighted 2
r_tsv <- residuals(mle_fit_tsv, type = "sweighted2")
r_0.01 <- residuals(mle_fit_0.01, type = "sweighted2")
r_0.001 <- residuals(mle_fit_0.001, type = "sweighted2")
r_mm <- residuals(mle_fit_mm, type = "sweighted2")
r_bd <- residuals(mle_fit_bd, type = "sweighted2")


#Randomized quantile residuals
rq_tsv   <- rqbeta(ysul_tsv, X=X, Z=Z, fitMLE_tsv$beta, fitMLE_tsv$gama)
rq_0.01  <- rqbeta(ysul_0.01, X=X, Z=Z, fitMLE_0.01$beta, fitMLE_0.01$gama)
rq_0.001 <- rqbeta(ysul_0.001, X=X, Z=Z, fitMLE_0.001$beta, fitMLE_0.001$gama)
rq_mm    <- rqbeta(ysul_mm, X=X, Z=Z, fitMLE_mm$beta, fitMLE_mm$gama)
rq_bd   <- rqbeta(y46, X=X46, Z=Z46, fitMLE_bd$beta, fitMLE_bd$gama)




#Standardized weighted residual 2 - SMLE
RP_tsv   <- residuals_beta(ysul_tsv, X, Z,c(fitMLE_tsv$beta,fitMLE_tsv$gama), 
                         linkmu="logit", linkphi="log")
RP_0.01  <- residuals_beta(ysul_0.01, X, Z,c(fitMLE_0.01$beta,fitMLE_0.01$gama), 
                         linkmu="logit", linkphi="log")
RP_0.001 <- residuals_beta(ysul_0.001, X, Z,c(fitMLE_0.001$beta,fitMLE_0.001$gama), 
                         linkmu="logit", linkphi="log")
RP_mm    <- residuals_beta(ysul_mm, X, Z,c(fitMLE_mm$beta,fitMLE_mm$gama), 
                         linkmu="logit", linkphi="log")
RP_bd    <- residuals_beta(y46, X46, Z46,c(fitMLE_bd$beta,fitMLE_bd$gama), 
                         linkmu="logit", linkphi="log")
###Worm Plots

fit_tsv2   <- gamlss(ysul_tsv ~ SSTsul, 
                    family = BE(mu.link = "logit", sigma.link = "log"))
fit_0.01_2  <- gamlss(ysul_0.01 ~ SSTsul, 
                   family = BE(mu.link = "logit", sigma.link = "log"))
fit_0.001_2 <- gamlss(ysul_0.001 ~ SSTsul, 
                   family = BE(mu.link = "logit", sigma.link = "log"))
fit_mm2     <- gamlss(ysul_mm ~ SSTsul, 
                   family = BE(mu.link = "logit", sigma.link = "log"))
fit_bd2     <- gamlss(ysul ~ SSTsul, 
                   family = BE(mu.link = "logit", sigma.link = "log"), data = data1_46)

wp(fit_tsv2, line = FALSE, cex.lab = 2, cex.axis = 2)  #transformation
wp(fit_0.01_2, line = FALSE, cex.lab = 2, cex.axis = 2) #adding 0.01
wp(fit_0.001_2, line = FALSE, cex.lab = 2, cex.axis = 2) #adding 0.001
wp(fit_mm2, line = FALSE, cex.lab = 2, cex.axis = 2) #Replacing by max. in (0,1)
wp(fit_bd2, line = FALSE, cex.lab = 2, cex.axis = 2) #Removing obs. 46
wp(resid = rq_tsv, line = FALSE, cex.lab = 2, cex.axis = 2) #transformation (RobBR)
wp(resid = rq_0.01, line = FALSE, cex.lab = 2, cex.axis = 2) #adding 0.01 (RobBR)
wp(resid = rq_0.001, line = FALSE, cex.lab = 2, cex.axis = 2)  #adding 0.001 (RobBR)
wp(resid = rq_mm, line = FALSE, cex.lab = 2, cex.axis = 2) #Replacing by max. in (0,1) (RobBR)
wp(resid = rq_bd, line = FALSE, cex.lab = 2, cex.axis = 2) #Removing obs. 46 (RobBR)
   


##BR plots

#Cook's distance

cd_tsv   <- betareg:::cooks.distance.betareg(mle_fit_tsv)
cd_0.01  <- betareg:::cooks.distance.betareg(mle_fit_0.01)
cd_0.001 <- betareg:::cooks.distance.betareg(mle_fit_0.001)
cd_mm    <- betareg:::cooks.distance.betareg(mle_fit_mm)
cd_bd    <- betareg:::cooks.distance.betareg(mle_fit_bd)

#Standard weighted residual type 2

r_tsv   <- residuals(mle_fit_tsv, type = "sweighted2")
r_0.01  <- residuals(mle_fit_0.01, type = "sweighted2")
r_0.001 <- residuals(mle_fit_0.001, type = "sweighted2")
r_mm    <- residuals(mle_fit_mm, type = "sweighted2")
r_bd    <- residuals(mle_fit_bd, type = "sweighted2")

#Generalized leverage

gl_tsv   <- gleverage(mle_fit_tsv)
gl_0.01  <- gleverage(mle_fit_0.01)
gl_0.001 <- gleverage(mle_fit_0.001)
gl_mm    <- gleverage(mle_fit_mm)
gl_bd    <- gleverage(mle_fit_bd)


#Fitted values 
fv_tsv   <- mle_fit_tsv$fitted.values
fv_0.01  <- mle_fit_0.01$fitted.values
fv_0.001 <- mle_fit_0.001$fitted.values
fv_mm    <- mle_fit_mm$fitted.values
fv_bd    <- mle_fit_bd$fitted.values

#Linear predictor values
lp_tsv   <- log(fv_tsv/(1-fv_tsv))
lp_0.01  <- log(fv_0.01/(1-fv_0.01))
lp_0.001 <- log(fv_0.001/(1-fv_0.001))
lp_mm    <- log(fv_mm/(1-fv_mm))
lp_bd    <- log(fv_bd/(1-fv_bd))


par(mar = c(5.6, 4.6, 4.6, 2.6))
lab <- 1:77
#BR after linear tranformation

plot(1:length(ysul_tsv), r_tsv, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(ysul_0.01), r_tsv, lab, cex = 1.2)

plot(1:length(ysul_tsv), cd_tsv, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_tsv, gl_tsv, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_tsv, gl_tsv, lab, cex = 1.2)

plot(lp_tsv, r_tsv, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(lp_tsv, r_tsv, lab, cex = 1.2)

#BR after adding 0.01

plot(1:length(ysul_0.01), r_0.01, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(ysul_0.01), r_0.01, lab, cex = 1.2)

plot(1:length(ysul_0.01), cd_0.01, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_0.01, gl_0.01, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_0.01, gl_0.01, lab, cex = 1.2)

plot(lp_0.01, r_0.01, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(lp_0.01, r_0.01, lab, cex = 1.2)

#BR after adding 0.001
plot(1:length(ysul_0.001), r_0.001, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(ysul_0.001), r_0.001, lab, cex = 1.2)

plot(1:length(ysul_0.001), cd_0.001, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_0.001, gl_0.001, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_0.001, gl_0.001, lab, cex = 1.2)

plot(lp_0.001, r_0.001, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(lp_0.001, r_0.001, lab, cex = 1.2)

#BR after replacing bound. obs. with min in (0,1)

plot(1:length(ysul_mm), r_mm, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(ysul_mm), r_mm, lab, cex = 1.2)

plot(1:length(ysul_mm), cd_mm, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_mm, gl_mm, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_mm, gl_mm, lab, cex = 1.2)

plot(lp_mm, r_mm, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")

#BR after excluding boundary obs. 

plot(1:length(data1_46$ysul), r_bd, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(data1_46$ysul), r_bd, lab, cex = 1.2)

plot(1:length(data1_46$ysul), cd_bd, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_bd, gl_bd, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_bd, gl_bd, lab, cex = 1.2)

plot(lp_bd, r_bd, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
##Envelopes

#MLE 
envelope_SMLE(y=ysul_tsv, X=X, Z=Z, theta=c(fitMLE1_tsv$beta,fitMLE1_tsv$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=ysul_0.01, X=X, Z=Z, theta=c(fitMLE1_0.01$beta,fitMLE1_0.01$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=ysul_0.001, X=X, Z=Z, theta=c(fitMLE1_0.001$beta,fitMLE1_0.001$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=ysul_mm, X=X, Z=Z, theta=c(fitMLE1_mm$beta,fitMLE1_mm$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=y46, X=X46, Z=Z46, theta=c(fitMLE1_bd$beta,fitMLE1_bd$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:76)

#SMLE
envelope_SMLE(y=ysul_tsv, X=X, Z=Z, theta=c(fitMLE_tsv$beta,fitMLE_tsv$gama), 
              linkmu="logit", linkphi="log", SMLE=T,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=ysul_0.01, X=X, Z=Z, theta=c(fitMLE_0.01$beta,fitMLE_0.01$gama), 
              linkmu="logit", linkphi="log", SMLE=T,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=ysul_0.001, X=X, Z=Z, theta=c(fitMLE_0.001$beta,fitMLE_0.001$gama), 
              linkmu="logit", linkphi="log", SMLE=T,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=ysul_mm, X=X, Z=Z, theta=c(fitMLE_mm$beta,fitMLE_mm$gama), 
              linkmu="logit", linkphi="log", SMLE=T,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:77)
envelope_SMLE(y=y46, X=X46, Z=Z46, theta=c(fitMLE_bd$beta,fitMLE_bd$gama), 
              linkmu="logit", linkphi="log", SMLE=T,
              main.title = "", faixa.fixed = c(-4,4), labels.fixed =1:76)

