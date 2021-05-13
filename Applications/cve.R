################################CVE application###############################
#Data from Morrison et al. (2020) [https://doi.org/10.1371/journal.pmed.1003049]
#Response variable is the median CVE percentage in each county in the state of
#Texas
############################################################################



library(dplyr)
library(broom)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(selectiveInference)
library(betareg)
library(glmnet)
library(lmtest)



source("functions.R")
source("make_functions.R")
source("SMLE.R") #File found in RobustBetaRegression repository
source("Resfunction.r")   #File found in RobustBetaRegression repository
source("envelope_function.R") #File found in RobustBetaRegression repository
load("District_Data_New2.RData")
load("County_Data_New.RData")

dictionary <-read.csv("S1_File.csv")


#Creating dataset
calculateVIF <- function(predictors) {
  vif.df = data.frame(matrix(nrow = ncol(predictors), ncol = 2))
  colnames(vif.df) = c("predictor", "VIF")
  for (i in 1:ncol(predictors)) {
    target = colnames(predictors)[i]
    covariates = paste(colnames(predictors)[-i], collapse = "+")
    formula = paste0(target, "~", covariates)
    r2 = summary(lm(formula, data = predictors))$r.squared
    vif = 1 / (1 - r2)
    vif.df[i,] = c(target, round(vif, digits = 2))
  }
  return(vif.df)
}
standard.remove = c(
  "County",
  "District",
  "district",
  "Metro_area",
  "Exemp_12.13",
  "Exemp_13.14",
  "Exemp_14.15",
  "Exemp_15.16",
  "Exemp_16.17",
  "Exemp_17.18",
  "Exemp_12.13_bc",
  "Exemp_13.14_bc",
  "Exemp_14.15_bc",
  "Exemp_15.16_bc",
  "Exemp_16.17_bc",
  "Exemp_17.18_bc",
  "Type"
)
theoretical.dist.variables = c(
  "MedIncome",
  "perc_white",
  "move_diff_state_est",
  "ESL",
  "DA00AR01S16R",
  "insured_kids_perc",
  "private_perc"
)
theoretical.county.variables = c(
  "MedIncome", 
  "PercWhite", 
  "AnnPopChng", 
  "ESL", 
  "PropPrivate",
  "Bachelors"
)




###County CVE
### Statewide  

county.data %>%
  mutate(
    metro.indicator = ifelse(is.na(Metro_area), "no", "yes"),
    Exemp_17.18.trans =  Exemp_17.18 / 100,
    target = log(Exemp_17.18.trans / (1 - Exemp_17.18.trans))
  ) -> dataset
remove.cols = which(
  colnames(dataset) %in%
    c(
      standard.remove,
      "Prop.risk3",
      "Prop.risk3_16.17",
      "Prop.risk2",
      "Prop.risk5",
      "Prop.risk4",
      "Prop.risk2",
      "Prop.risk1",
      "Urbanization",
      "target"
    )
)
dataset = dataset[, -remove.cols]
dataset = dataset %>% filter(!is.na(PropPrivate)) ## have to remove terrell --
dataset_wo <- rem_bd(dataset, 29)

attach(dataset)

#Response variable boxplot

par(mar = c(5.6, 4.6, 4.6, 2.6))

adjbox(Exemp_17.18.trans, ylab = "Median % of children with a CVE ", pch = 16,
       cex.lab = 1.4, cex.axis = 1.4)
identify(rep(1, 239), Exemp_17.18.trans, 1:239, cex = 1.2)
# Perform beta regression
y_tsv <- tsv(dataset$Exemp_17.18.trans)  #Transforming using method proposed by Smithson and Verkuilen
y_0.01 <- rep_e(dataset$Exemp_17.18.trans, 0.0001)  #Subtracting 0.01 from observation 46
y_0.001 <- rep_e(dataset$Exemp_17.18.trans, 0.00001) #Subtracting 0.001 from observation 46
y_mm <- rep_mm(dataset$Exemp_17.18.trans) #Replacing by smallest observation in (0,1)



###Beta regression


mle_fit_tsv <- betareg(y_tsv ~ ESL+Bachelors+PercWhite+MedIncome|PropUnder5Poverty+
                         metro.indicator) 
summary(mle_fit_tsv)


mle_fit_0.01 <- betareg(y_0.01 ~ ESL+Bachelors+PercWhite+MedIncome|PropUnder5Poverty+
                          metro.indicator)
summary(mle_fit_0.01)


mle_fit_0.001 <- betareg(y_0.001 ~ ESL+Bachelors+PercWhite+MedIncome|PropUnder5Poverty+
                           metro.indicator)
summary(mle_fit_0.001)


mle_fit_mm <- betareg(y_mm ~ ESL+Bachelors+PercWhite+MedIncome|PropUnder5Poverty+
                        metro.indicator)
summary(mle_fit_mm)


mle_fit_bd <- betareg(Exemp_17.18.trans ~ ESL+Bachelors+PercWhite+MedIncome|PropUnder5Poverty+
                        metro.indicator, data = dataset_wo)
summary(mle_fit_bd)


#Diagnostics
par(mar = c(5.6, 4.6, 4.6, 2.6))

plot(mle_fit_tsv, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = " ", sub.caption = "", pch = 16)


plot(mle_fit_0.01, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)


plot(mle_fit_0.001, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)


plot(mle_fit_mm, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)


plot(mle_fit_bd, which = 1:4, type = "sweighted2", cex.lab=2, cex.axis=2,
     caption = "", sub.caption = "", pch = 16)



#Robust estimation in Beta regression

#Defining model matrices
X <- matrix(c(rep(1,239), ESL, Bachelors, PercWhite, MedIncome), ncol=5, byrow=F);
Z <- matrix(c(rep(1,239), PropUnder5Poverty, as.numeric(factor(metro.indicator, labels = c(0,1)))-1), 
            ncol=3, byrow=F);


y_wo <- dataset_wo$Exemp_17.18.trans; 
Xwo  <- matrix(c(rep(1,227), dataset_wo$ESL, dataset_wo$Bachelors, dataset_wo$PercWhite, 
              dataset_wo$MedIncome), ncol=5, byrow=F)
Zwo  <- matrix(c(rep(1,227), dataset_wo$Bachelors, dataset_wo$PropUnder5Poverty), 
              ncol=3, byrow=F)
##Model estimation
fitMLE_tsv<- SMLE_BETA(y=y_tsv, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                       spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")
fitMLE_tsv

fitMLE_0.01<- SMLE_BETA(y=y_0.01, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                        spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")
fitMLE_0.01

fitMLE_0.001<- SMLE_BETA(y=y_0.001, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                         spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_0.001


fitMLE_mm<- SMLE_BETA(y=y_mm, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                      spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_mm


fitMLE_bd<- SMLE_BETA(y=y_wo, X=Xwo, Z=Zwo, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, 
                      spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_bd


#Not robust estimation(equivalent to betareg)

#Linear transformation
fitMLE1_tsv<- SMLE_BETA(y=y_tsv, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                        spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")
#Adding 0.0001
fitMLE1_0.01<- SMLE_BETA(y=y_0.01, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                         spac =0.02,method="BFGS",startV="CP", linkmu="logit", linkphi="log")
#Adding 0.00001
fitMLE1_0.001<- SMLE_BETA(y=y_0.001, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                          spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")
#Replacing boundaries with min. value in (0,1)
fitMLE1_mm<- SMLE_BETA(y=y_mm, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                       spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")

#Excluding boundary obs.
fitMLE1_bd<- SMLE_BETA(y=y_wo, X=Xwo, Z=Zwo, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, 
                       spac =0.02,method="BFGS", startV="CP", linkmu="logit", linkphi="log")

###Inflated beta regression

#Model estimation
beinf_fit <- gamlss(Exemp_17.18.trans ~ ESL+Bachelors+PercWhite+MedIncome, 
                    sigma.formula = ~PropUnder5Poverty+
                      metro.indicator,
                    family = BEZI(mu.link = "logit", sigma.link = "log",
                                  nu.link = "identity"), control=gamlss.control(trace=FALSE))
summary(beinf_fit)
wp(beinf_fit)

##Diagnostics

#Worm plot
wp(beinf_fit, line = FALSE, cex.lab = 2, cex.axis = 2)

plot(rq_bd, rqzibeta_wo, ylim = c(-3, 3), pch = 16, cex.lab = 2, cex.axis = 2, 
     xlab = "BR residuals", ylab = "Zero-inflated BR residuals")
abline(0,1, lty = 2, col = "red")

###Quasi-likelihood (Bonat et al.(2019))

#Model estimation
form1 <- Exemp_17.18.trans ~ ESL+Bachelors+PercWhite+MedIncome

Z0_1 <- mc_id(dataset)
fit_mc1 <- mcglm(c(form1), list(Z0_1), link = "logit", 
                 variance = "binomialP", power_fixed = FALSE, 
                 control_algorithm = list(tunning = 0.75, correct = FALSE),
                 data = dataset)
summary(fit_mc1)

fit_mc1p <- mcglm(c(form1), list(Z0_1), link = "logit", 
                 variance = "binomialP", power_fixed = TRUE, 
                 control_algorithm = list(tunning = 0.75, correct = FALSE),
                 data = dataset)
summary(fit_mc1p)



pr_ql <- residuals(fit_mc1, type = "pearson")
pr_ql <- as.numeric(pr_ql)


#Comparing residuals with fitted values
f_values <- as.numeric(fitted(fit_mc1))
plot( f_values, pr_ql, ylab = "Pearson Residuals", xlab = "Fitted values", 
      pch=16, cex=1, cex.lab=2, cex.axis=2, ylim = c(-4, 4))

abline(a = 0, b = 0, lty = 2, col = "lightgrey")

##Diagnostic plots
#Randomized quantile residuals
rq_tsv   <- rqbeta(y_tsv, X=X, Z=Z, fitMLE_tsv$beta, fitMLE_tsv$gama)
rq_0.01  <- rqbeta(y_0.01, X=X, Z=Z, fitMLE_0.01$beta, fitMLE_0.01$gama)
rq_0.001 <- rqbeta(y_0.001, X=X, Z=Z, fitMLE_0.001$beta, fitMLE_0.001$gama)
rq_mm    <- rqbeta(y_mm, X=X, Z=Z, fitMLE_mm$beta, fitMLE_mm$gama)
rq_bd   <- rqbeta(y_wo, X=Xwo, Z=Zwo, fitMLE_bd$beta, fitMLE_bd$gama)

rq1_tsv   <- rqbeta(y_tsv, X=X, Z=Z, fitMLE1_tsv$beta, fitMLE1_tsv$gama)
rq1_0.01  <- rqbeta(y_0.01, X=X, Z=Z, fitMLE1_0.01$beta, fitMLE1_0.01$gama)
rq1_0.001 <- rqbeta(y_0.001, X=X, Z=Z, fitMLE1_0.001$beta, fitMLE1_0.001$gama)
rq1_mm    <- rqbeta(y_mm, X=X, Z=Z, fitMLE1_mm$beta, fitMLE1_mm$gama)
rq1_bd   <- rqbeta(y_wo, X=Xwo, Z=Zwo, fitMLE1_bd$beta, fitMLE1_bd$gama)

rqzibeta <- residuals(beinf_fit, what = "z-scores")
rqzibeta_wo <- rqzibeta[-which(Exemp_17.18.trans == 0)]


#Worm plots
wp(resid = rq_tsv, line = FALSE, cex.lab = 2, cex.axis = 2, ylab = "Residuals")
wp(resid = rq_0.01, line = FALSE)
wp(resid = rq_0.001, line = FALSE)
wp(resid = rq_mm, line = FALSE)
wp(resid = rq_bd, line = FALSE)
wp(resid = scale(fit_mc1$residuals), line = FALSE) 

wp(resid = rq1_tsv, line = FALSE, cex.lab = 2, cex.axis = 2)
wp(resid = rq1_0.01, line = FALSE, cex.lab = 2, cex.axis = 2)
wp(resid = rq1_0.001, line = FALSE, cex.lab = 2, cex.axis = 2)
wp(resid = rq1_mm, line = FALSE, cex.lab = 2, cex.axis = 2)
wp(resid = rq1_bd, line = FALSE, cex.lab = 2, cex.axis = 2)


###Envelopes

##Beta regression

#Linear transformation
envelope_SMLE(y=y_tsv, X=X, Z=Z, theta=c(fitMLE1_tsv$beta,fitMLE1_tsv$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-5,5), labels.fixed =1:235)
#Adding 0.0001
envelope_SMLE(y=y_0.01, X=X, Z=Z, theta=c(fitMLE1_0.01$beta,fitMLE1_0.01$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-5,5), labels.fixed =1:235)

#Adding 0.00001
envelope_SMLE(y=y_0.001, X=X, Z=Z, theta=c(fitMLE1_0.001$beta,fitMLE1_0.001$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-5,5), labels.fixed =1:235)

#Replacing boundaries with min. value in (0,1)
envelope_SMLE(y=y_mm, X=X, Z=Z, theta=c(fitMLE1_mm$beta,fitMLE1_mm$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-5,5), labels.fixed =1:235)

#Excluding boundary obs.
envelope_SMLE(y=y_wo, X=Xwo, Z=Zwo, theta=c(fitMLE1_bd$beta,fitMLE1_bd$gama), 
              linkmu="logit", linkphi="log", SMLE=F,
              main.title = "", faixa.fixed = c(-5,5), labels.fixed =1:234)


#Inflated beta regression

envelope_BEI(beinf_fit, type = "BEZI", main.title = "", faixa.fixed = c(-4,2.5),
             labels.fixed = 1:239)


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
lab <- 1:239
#BR after linear tranformation

plot(1:length(y_tsv), r_tsv, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")

plot(1:length(y_tsv), cd_tsv, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_tsv, gl_tsv, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_tsv, gl_tsv, lab, cex = 1.2)

plot(lp_tsv, r_tsv, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")

#BR after adding 0.01

plot(1:length(y_0.01), r_0.01, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(y_0.01), r_0.01, lab, cex = 1.2)

plot(1:length(y_0.01), cd_0.01, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_0.01, gl_0.01, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_0.01, gl_0.01, lab, cex = 1.2)

plot(lp_0.01, r_0.01, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")

#BR after adding 0.001
plot(1:length(y_0.001), r_0.001, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(y_0.001), r_0.001, lab, cex = 1.2)

plot(1:length(y_0.001), cd_0.001, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_0.001, gl_0.001, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_0.001, gl_0.001, lab, cex = 1.2)

plot(lp_0.001, r_0.001, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")

#BR after replacing bound. obs. with min in (0,1)

plot(1:length(y_mm), r_mm, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(y_mm), r_mm, lab, cex = 1.2)

plot(1:length(y_mm), cd_mm, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_mm, gl_mm, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_mm, gl_mm, lab, cex = 1.2)

plot(lp_mm, r_mm, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")

#BR after excluding boundary obs. 

plot(1:length(y_wo), r_bd, xlab = "Obs. number", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")
identify(1:length(y_wo), r_bd, lab, cex = 1.2)

plot(1:length(y_wo), cd_bd, xlab = "Obs. number", ylab = "Cook's distance",
     cex.lab = 2.0, cex.axis = 2.0, type = "h")

plot(fv_bd, gl_bd, xlab = "Predicted values", ylab = "Generalized leverage",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
identify(fv_bd, gl_bd, lab, cex = 1.2)

plot(lp_bd, r_bd, xlab = "Linear predictor", ylab = "Residual",
     cex.lab = 2.0, cex.axis = 2.0, pch = 16)
abline(a = 0, b = 0, lty = 3, col = "grey")



###Comparison graphs

plot(fitted(mle_fit_bd), fitted(mle_fit_tsv)[-which(Exemp_17.18.trans == 0)], 
     xlab = "BR w/o boundary obs.", ylab = "BR after linear trans.", cex.axis = 2,
     cex.lab = 2, pch = 16)
abline(0, 1, lty = 2, col = "red")

plot(fitted(mle_fit_0.01), fitted(mle_fit_tsv), 
     xlab = "BR after adding 0.0001", ylab = "BR after linear trans.", cex.axis = 2,
     cex.lab = 2, pch = 16)
abline(0, 1, lty = 2, col = "red")

plot(fitted(mle_fit_mm), fitted(mle_fit_tsv), 
    xlab = "BR after replacing boundary obs.", ylab = "BR after linear trans.", 
    cex.axis = 2,
    cex.lab = 2, pch = 16)
abline(0, 1, lty = 2, col = "red")



###Simulation - Worm plot



G <- 9
beta <- coefficients(beinf_fit)[1:5]
gama <- beinf_fit$sigma.coefficients

mu <- exp(X%*%beta)/(1+exp(X%*%beta))
phi <- exp(Z%*%gama)
est <- matrix(0, nrow = 20, ncol =9)
par(mfrow = c(3, 3))
for (i in 1:G) {
  sim_data <- data.frame(numeric(length(mu)),ESL, Bachelors, PercWhite, MedIncome, 
                         PropUnder5Poverty, metro.indicator)
  names(sim_data)[1] <- "y"
  for (j in 1:length(sim_data$y)) {
    sim_data$y[j] <- rBEZI(1, mu[j], phi[j], beinf_fit$nu.coefficients)
  }
  bezi_fit <- gamlss(y ~ ESL+Bachelors+PercWhite+MedIncome, 
                      sigma.formula = ~PropUnder5Poverty+
                        metro.indicator,
                      family = BEZI(mu.link = "logit", sigma.link = "log",
                                    nu.link = "identity"), data = sim_data,
                     control=gamlss.control(trace=FALSE))
  est[i,] <- c(coefficients(bezi_fit), bezi_fit$sigma.coefficients, 
               bezi_fit$nu.coefficients)
  wp(bezi_fit, cex.lab = 2, cex.axis = 2)
  }

qqnorm(residuals(bezi_fit, what = "z-scores"), cex.lab = 2,
       cex.axis = 2, ylim = c(-3, 2.5), main = "")
abline(0,1, lty = 2)


#Simulation - worm plot - Beta regression
G <- 9
beta <- coefficients(mle_fit_0.01)[1:5]
gama <- coefficients(mle_fit_0.01)[6:8]

mu <- rep(exp(X%*%beta)/(1+exp(X%*%beta)),4)
phi <- rep(exp(Z%*%gama),4)
est <- matrix(0, nrow = 9, ncol =8)
par(mfrow = c(3, 3))
for (i in 1:G) {
  sim_data <- data.frame(numeric(length(mu)),rep(ESL,4), rep(Bachelors,4)
                         , rep(PercWhite,4), rep(MedIncome,4), 
                         rep(PropUnder5Poverty,4), rep(metro.indicator,4))
  names(sim_data)<- c("y", "ESL", "Bachelors", "PercWhite", "MedIncome", 
                      "PropUnder5Poverty", "metro.indicator")
  sim_data$y <- rbeta(length(mu), mu*phi, (1-mu)*phi)
  
  be_fit <- gamlss(y ~ ESL+Bachelors+PercWhite+MedIncome, 
                     sigma.formula = ~PropUnder5Poverty+
                       metro.indicator,
                     family = BE(mu.link = "logit", sigma.link = "log"), data = sim_data,
                     control=gamlss.control(trace=FALSE))
  est[i,] <- c(coefficients(be_fit), be_fit$sigma.coefficients)
  wp(be_fit, cex.lab = 2, cex.axis = 2)
}

qqnorm(residuals(bezi_fit, what = "z-scores"), cex.lab = 2,
       cex.axis = 2, ylim = c(-3, 2.5), main = "")
abline(0,1, lty = 2)




