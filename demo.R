library('nprobust')
library("np")
library(KernSmooth)
library(ATE.ncb)
library(foreach)
library(doParallel)
#######################
n = 100
variance = 1

################################################
################# generate data ################
################################################
source("generate_data_1.R")


###########################################
################# IPW  ####################
###########################################

##### use logistic regression to get propensity score #####
lg = glm(treat ~ X, family = binomial(link = "logit"))
prophat = fitted(lg)
ipw0 <- 1 / prophat * treat + 1 / (1 - prophat) * (1 - treat) # inverse propensity
Yadj <- ipw0 * treat * Y - ipw0 * (1 - treat) * Y

###### local constant regression #######
h = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3] * n ^ (1 / 5) * n ^ (-2 / 7)
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = h, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
####### L2 loss ######
l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("ipw:", l2loss, "\n"))



################################################
##### kernel-based covariate balancing ATE  ####
################################################
OUTstd <- transform.sob(X)
Xstd <- OUTstd$Xstd
Xlim <- OUTstd$Xlim
K <- getGram(Xstd) # get Gram matrix using Sobolev kernel
Vstd <- Xstd[, 1]
Vstd <- as.matrix(Vstd)
Vnewstd = (Vnew - Xlim[1, 1]) / diff(Xlim[, 1])
Vnewstd <- as.matrix(Vnewstd)

nlam <- 20
lams <- exp(seq(log(1e-6), log(1), len = nlam))
eta = 0.001

fit1_ate <- ATE.ncb.SN(treat, K, lam1s = lams, lam2s = eta * lams, traceit = FALSE)
fit0_ate <- ATE.ncb.SN(1 - treat, K, lam1s = lams, lam2s = eta * lams, traceit = FALSE)

Yadj <- fit1_ate$w * treat * Y - fit0_ate$w * (1 - treat) * Y

hate = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3]
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hate, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("ate:", l2loss, "\n"))



################################################
##### Proposed estimator PCATE balancing  ####
################################################
CATE_f = function(h) {
  lams <- exp(seq(log(1e-06 / h), log(1 / h), len = nlam))
  fit1 <- ATE.ncb.SN2(V = Vstd, Vg = Vnewstd, h, treat, K, lam1s = lams, lam2s = eta *
    h * lams, traceit = FALSE, method = 2)
  fit0 <- ATE.ncb.SN2(Vstd, Vg = Vnewstd, h, 1 - treat, K, lam1s = lams, lam2s = eta *
    h * lams, traceit = FALSE, method = 2)
  return(list(fit1 = fit1, fit0 = fit0))
}

Rcpp::sourceCpp("kernel_cpp.cpp")
source("core3_g_new.R")
hdpi = hate / diff(Xlim[, 1])
cate_noaug = CATE_f(hdpi)

Yadj <- cate_noaug$fit1$w * treat * Y - cate_noaug$fit0$w * (1 - treat) * Y
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal", bandwidth = hate,
  gridsize = 300, range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
dpi_loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("proposed:", dpi_loss, "\n"))



###########################################
############# LM augmentation #############
###########################################
X = as.data.frame(X)
names(X) = as.character(c(1:p))
data = as.data.frame(cbind(Y, X)[treat == 1,])
model1 = lm(Y ~ ., data = data)
Y1hat = predict(model1, newdata = as.data.frame(X))

data = as.data.frame(cbind(Y, X)[treat == 0,])
model0 = lm(Y ~ ., data = data)
Y0hat = predict(model0, newdata = as.data.frame(X))

phi1 = Y1hat + treat * (Y - Y1hat) / prophat
phi0 = Y0hat + (1 - treat) * (Y - Y0hat) / (1 - prophat)
Yadj = phi1 - phi0


hdpi_lc = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3] * n ^ (1 / 5) * n ^ (-2 / 7)
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hdpi_lc, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
l2loss = mean((Ytemp$y - real) ^ 2)

cat(paste0("aipw_LM:", l2loss, "\n"))


############ LM fit ############
Yadj = Y1hat - Y0hat
hdpi_lc = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3] * n ^ (1 / 5) * n ^ (-2 / 7)
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hdpi_lc, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("LMfit:", l2loss, "\n"))


phi1 = Y1hat + treat * (Y - Y1hat) * fit1_ate$w
phi0 = Y0hat + (1 - treat) * (Y - Y0hat) * fit0_ate$w
Yadj = phi1 - phi0
hate_LM = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3] * n ^ (1 / 5) * n ^ (-2 / 7)
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hate_LM, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("ate_LM:", l2loss, "\n"))


hdpi = hate_LM / diff(Xlim[, 1])
cate_LM = CATE_f(hdpi)
phi1 = Y1hat + treat * (Y - Y1hat) * cate_LM$fit1$w
phi0 = Y0hat + (1 - treat) * (Y - Y0hat) * cate_LM$fit0$w
Yadj = phi1 - phi0

Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hate_LM, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
dpi_l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("proposed_LM:", dpi_l2loss, "\n"))




#########################################
###########  KRR  #######################
#########################################

source(file = "KRR.R")
lamseq = exp(seq(log(1e-06), log(1e-02), length.out = 20))
######### KRR fit ############
data = as.data.frame(cbind(Y, Xstd)[treat == 1,])
KRR_out = KRR.CV(location = data[, 2:(p + 1)], X = data[, 1], lamseq = lamseq)

Knew = K[, treat == 1]
Y1hat = Knew %*% KRR_out$alpha


data = as.data.frame(cbind(Y, Xstd)[treat == 0,])
KRR_out = KRR.CV(location = data[, 2:(p + 1)], X = data[, 1], lamseq = lamseq)

Knew = K[, treat == 0]
Y0hat = Knew %*% KRR_out$alpha

Yadj = Y1hat - Y0hat

hKRR = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3] * n ^ (1 / 5) * n ^ (-2 / 7)
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hKRR, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
KRR_l2loss = mean((Ytemp$y - real) ^ 2)

cat(paste0("KRRfit:", KRR_l2loss, "\n"))


phi1 = Y1hat + treat * (Y - Y1hat) / prophat
phi0 = Y0hat + (1 - treat) * (Y - Y0hat) / (1 - prophat)
Yadj = phi1 - phi0

hdpi_lc = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3] * n ^ (1 / 5) * n ^ (-2 / 7)
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hdpi_lc, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("aipw_KRR:", l2loss, "\n"))




phi1 = Y1hat + treat * (Y - Y1hat) * fit1_ate$w
phi0 = Y0hat + (1 - treat) * (Y - Y0hat) * fit0_ate$w
Yadj = phi1 - phi0
hate_KRR = lpbwselect(y = Yadj, x = V, kernel = "gau", p = 0, bwselect = 'imse-dpi', bwcheck = NULL)$bws[1, 3] * n ^ (1 / 5) * n ^ (-2 / 7)
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hate_KRR, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("ate_KRR:", l2loss, "\n"))


hdpi = hate_KRR / diff(Xlim[, 1])
cate_KRR = CATE_f(hdpi)
phi1 = Y1hat + treat * (Y - Y1hat) * cate_KRR$fit1$w
phi0 = Y0hat + (1 - treat) * (Y - Y0hat) * cate_KRR$fit0$w
Yadj = phi1 - phi0
Ytemp = locpoly(x = V, y = Yadj, drv = 0, degree = 0, kernel = "normal",
        bandwidth = hate_KRR, gridsize = 300,
        range.x = c(min(Vnew), max(Vnew)), binned = FALSE, truncate = FALSE)
dpi_l2loss = mean((Ytemp$y - real) ^ 2)
cat(paste0("cate_KRR:", dpi_l2loss, "\n"))


