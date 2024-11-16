library(dplyr)
library(tidyr)
library(MASS)
library(ggplot2)
library(forcats)
library(reshape2)
library(Matrix)

# Golden method which uses whole data assuming no missing data
get_gold = function(X1, X2, D){
  z_alpha = qnorm(0.975)
  
  logit_gold = glm(D ~ X1 + X2, family = binomial)
  
  d = dim(X1)[2]
  d2 = dim(X2)[2]
  
  if(is.null(d2)){
    d2=1
  }
  
  est = c()
  sd = c()
  pval = c()
  or = c()
  ciL = c()
  ciU = c()
  
  for(i in 1:d){
    est = c(est,logit_gold$coefficients[1+i])
    sd = c(sd, sqrt(diag(vcov(logit_gold)[2:(1+d),2:(1+d)]))[i])
    pval = c(pval, 2*min(pnorm(logit_gold$coefficients[1+i], 0, sqrt(diag(vcov(logit_gold)[2:(1+d),2:(1+d)]))[i]),
                         1-pnorm(logit_gold$coefficients[1+i], 0, sqrt(diag(vcov(logit_gold)[2:(1+d),2:(1+d)]))[i])))
    or = c(or, exp(logit_gold$coefficients[1+i]))
    ciL = c(ciL, exp(logit_gold$coefficients[1+i]-z_alpha*sqrt(diag(vcov(logit_gold)[2:(1+d),2:(1+d)]))[i]))
    ciU = c(ciU, exp(logit_gold$coefficients[1+i]+z_alpha*sqrt(diag(vcov(logit_gold)[2:(1+d),2:(1+d)]))[i]))
  }
  #print("##########")
  
  for(j in 1:d2){
    est = c(est, logit_gold$coefficients[1+d+j])
    sd = c(sd, sqrt(vcov(logit_gold)[1+d+j,1+d+j]))
    pval = c(pval, 2*min(pnorm(logit_gold$coefficients[1+d+j], 0, sqrt(vcov(logit_gold)[1+d+j,1+d+j])),
                         1-pnorm(logit_gold$coefficients[1+d+j], 0, sqrt(vcov(logit_gold)[1+d+j,1+d+j]))))
    or = c(or, exp(logit_gold$coefficients[1+d+j]))
    ciL = c(ciL, exp(logit_gold$coefficients[1+d+j]-z_alpha*sqrt(vcov(logit_gold)[1+d+j,1+d+j])))
    ciU = c(ciU, exp(logit_gold$coefficients[1+d+j]+z_alpha*sqrt(vcov(logit_gold)[1+d+j,1+d+j])))
    
  }
  
  df = data.frame(est, sd, pval, or, ciL, ciU)
  
  return(df)
}

# Main method which uses only fully observed (main) data
get_main = function(X1_m, X2_m, D_m){
  z_alpha = qnorm(0.975)
  
  logit_main = glm(D_m ~ X1_m + X2_m, family = binomial)
  
  d = dim(X1_m)[2]
  d2 = dim(X2_m)[2]
  
  if(is.null(d2)){
    d2=1
  }
  
  est = c()
  sd = c()
  pval = c()
  or = c()
  ciL = c()
  ciU = c()
  
  for(i in 1:d){
    est = c(est,logit_main$coefficients[1+i])
    sd = c(sd, sqrt(diag(vcov(logit_main)[2:(1+d),2:(1+d)]))[i])
    pval = c(pval, 2*min(pnorm(logit_main$coefficients[1+i], 0, sqrt(diag(vcov(logit_main)[2:(1+d),2:(1+d)]))[i]),
                         1-pnorm(logit_main$coefficients[1+i], 0, sqrt(diag(vcov(logit_main)[2:(1+d),2:(1+d)]))[i])))
    or = c(or, exp(logit_main$coefficients[1+i]))
    ciL = c(ciL, exp(logit_main$coefficients[1+i]-z_alpha*sqrt(diag(vcov(logit_main)[2:(1+d),2:(1+d)]))[i]))
    ciU = c(ciU, exp(logit_main$coefficients[1+i]+z_alpha*sqrt(diag(vcov(logit_main)[2:(1+d),2:(1+d)]))[i]))
  }
  #print("##########")
  
  for(j in 1:d2){
    est = c(est, logit_main$coefficients[1+d+j])
    sd = c(sd, sqrt(vcov(logit_main)[1+d+j,1+d+j]))
    pval = c(pval, 2*min(pnorm(logit_main$coefficients[1+d+j], 0, sqrt(vcov(logit_main)[1+d+j,1+d+j])),
                         1-pnorm(logit_main$coefficients[1+d+j], 0, sqrt(vcov(logit_main)[1+d+j,1+d+j]))))
    or = c(or, exp(logit_main$coefficients[1+d+j]))
    ciL = c(ciL, exp(logit_main$coefficients[1+d+j]-z_alpha*sqrt(vcov(logit_main)[1+d+j,1+d+j])))
    ciU = c(ciU, exp(logit_main$coefficients[1+d+j]+z_alpha*sqrt(vcov(logit_main)[1+d+j,1+d+j])))
    
  }
  
  df = data.frame(est, sd, pval, or, ciL, ciU)
  
  return(df)
}

# Naive method which considers surrogate covariates as partially missing covariates
get_naive = function(X1_m, X2_m, D_m, X1_v, Z_v, D_v){
  
  z_alpha = qnorm(0.975)
  
  Rcpp::cppFunction('
  arma::mat inverse(const arma::mat& X) {
    return inv(X);
  }
', depends = "RcppArmadillo")
  
  d = dim(X1_m)[2]
  d2 = dim(X2_m)[2]
  
  if(is.null(d2)){
    d2=1
  }
  
  logit_m = glm(D_m ~ X1_m + X2_m, family = binomial)
  beta0_hat = logit_m$coefficients[1]
  beta1_hat = logit_m$coefficients[2:(1+d)]
  beta2_hat = logit_m$coefficients[(2+d):(1+d+d2)] 
  var_beta1_hat = vcov(logit_m)[2:(1+d),2:(1+d)]
  var_beta2_hat = as.matrix(vcov(logit_m)[(2+d):(1+d+d2),(2+d):(1+d+d2)])
  
  logit_v = glm(D_v ~ X1_v + Z_v, family = binomial)
  
  beta0_s_tilde = logit_v$coefficients[1]
  beta1_s_tilde = logit_v$coefficients[2:(1+d)]
  beta2_s_tilde = rep(logit_v$coefficients[2+d], d2)
  
  
  var_beta_s_tilde = vcov(logit_v)[2:(1+d+d2),2:(1+d+d2)]
  
  var_beta2_hat[!diag(1, nrow(var_beta2_hat), ncol(var_beta2_hat))] <- 0
  
  var_beta2_tilde = as.matrix(var_beta_s_tilde[(d+1):(d+d2),(d+1):(d+d2)])
  var_beta2_tilde[!diag(1, nrow(var_beta2_tilde), ncol(var_beta2_tilde))] <- 0
  
  beta1_naive = (t(as.matrix(beta1_hat, ncol=1))%*%inverse(diag(diag(var_beta1_hat)))+
                   t(as.matrix(beta1_s_tilde, ncol=1))%*%inverse(diag(diag(var_beta_s_tilde[1:d,1:d]))))%*%inverse((inverse(diag(diag(var_beta1_hat)))+inverse(var_beta_s_tilde[1:d,1:d])))
  beta2_naive = (t(as.matrix(beta2_hat, ncol=1))%*%inverse(var_beta2_hat)+
                   t(as.matrix(beta2_s_tilde, ncol=1))%*%inverse(var_beta2_tilde))%*%inverse((inverse(var_beta2_hat)+inverse(var_beta2_tilde)))
  var_beta1_naive = diag(diag(var_beta1_hat))%*%diag(diag(var_beta_s_tilde[1:d,1:d]))%*%inverse(diag(diag(var_beta1_hat))+diag(diag(var_beta_s_tilde[1:d,1:d])))
  
  var_beta2_naive = var_beta2_hat%*%var_beta2_tilde%*%inverse(var_beta2_hat+var_beta2_tilde)
  
  est = c()
  sd = c()
  pval = c()
  or = c()
  ciL = c()
  ciU = c()
  
  for(i in 1:d){
    est = c(est, beta1_naive[i])
    sd = c(sd, sqrt(diag(var_beta1_naive)[i]))
    pval = c(pval, 2*min(pnorm(beta1_naive[i], 0, sqrt(diag(var_beta1_naive)[i])),
                         1-pnorm(beta1_naive[i], 0, sqrt(diag(var_beta1_naive)[i]))))
    or = c(or, exp(beta1_naive[i]))
    ciL = c(ciL, exp(beta1_naive[i]-z_alpha*sqrt(diag(var_beta1_naive)[i])))
    ciU = c(ciU, exp(beta1_naive[i]+z_alpha*sqrt(diag(var_beta1_naive)[i])))
  }
  
  #print("##########")
  
  beta2_naive = as.vector(beta2_naive)
  var_beta2_naive = as.matrix(var_beta2_naive)
  
  for(j in 1:d2){
    est = c(est, beta2_naive[j])
    sd = c(sd, sqrt(var_beta2_naive[j,j]))
    pval = c(pval, 2*min(pnorm(beta2_naive[j], 0, sqrt(var_beta2_naive[j,j])),
                         1-pnorm(beta2_naive[j], 0, sqrt(var_beta2_naive[j,j]))))
    or = c(or, exp(beta2_naive[j]))
    ciL = c(ciL, exp(beta2_naive[j]-z_alpha*sqrt(var_beta2_naive[j,j])))
    ciU = c(ciU, exp(beta2_naive[j]+z_alpha*sqrt(var_beta2_naive[j,j])))
  }
  
  df = data.frame(est, sd, pval, or, ciL, ciU)
  
  return(df)
}

# Calibration method which adopts calibration model to handle missing covariates
get_calibration = function(X1_m, X2_m, Z_m, X1_v, Z_v, D_v){
  
  z_alpha = qnorm(0.975)
  
  Rcpp::cppFunction('
  arma::mat inverse(const arma::mat& X) {
    return inv(X);
  }
', depends = "RcppArmadillo")
  
  d = dim(X1_m)[2]
  d2 = dim(X2_m)[2]
  
  if(is.null(d2)){
    d2=1
  }
  
  logit_v = glm(D_v ~ X1_v + Z_v, family = binomial)
  
  beta0_s_tilde = logit_v$coefficients[1]
  beta1_s_tilde = logit_v$coefficients[2:(1+d)]
  beta2_s_tilde = logit_v$coefficients[(2+d):(1+d+d2)]
  
  var_beta_s_tilde = vcov(logit_v)[2:(1+d+d2),2:(1+d+d2)]
  
  
  linear_m = lm(X2_m ~ Z_m + X1_m)
  
  if (d2==1){
    gamma0_hat= linear_m$coefficients[1]
    gamma_z_hat = linear_m$coefficients[2]
    gamma_x_hat = linear_m$coefficients[3:(2+d)]
    
    linear_hat = diag(c(gamma_z_hat, rep(1, d)))
    linear_hat[1,2:(1+d)]=gamma_x_hat
  }
  else{
    gamma0_hat= linear_m$coefficients[1,]
    gamma_z_hat = linear_m$coefficients[2:(1+d2),]
    gamma_x_hat = linear_m$coefficients[(2+d2):(1+d+d2),]
    
    linear_hat = diag(rep(1,d+d2))
    linear_hat[1:d2,1:d2] = gamma_z_hat
    linear_hat[1:d2,(1+d2):(d+d2)] = t(gamma_x_hat)
  }
  
  
  #A is inverse of linear_hat
  A= inverse(linear_hat)
  
  var_linear_hat = vcov(linear_m)[-seq(1, 1+(d2-1)*(d+d2+1), by=(d+d2+1)),
                                  -seq(1, 1+(d2-1)*(d+d2+1), by=(d+d2+1))]
  
  var_linear_list = vector("list", d2)
  for (t in 1:d2) {
    var_linear_list[[t]] <- vector("list", d2)
  }
  
  for (t1 in 1:d2) {
    for (t2 in 1:d2) {
      var_linear_list[[t1]][[t2]] = var_linear_hat[(1+(t1-1)*(d+d2)):(t1*(d+d2)),(1+(t2-1)*(d+d2)):(t2*(d+d2))]
    }
    
  }
  
  
  covlist <- vector("list", d+d2)
  for (i in 1:(d+d2)) {
    covlist[[i]] <- vector("list", d+d2)
  }
  
  for (j1 in 1:(d+d2)) {
    for (j2 in 1:(d+d2)) {
      for (t1 in 1:d2) {
        for (t2 in 1:d2) {
          covlist[[j1]][[j2]] = diag(rep(0, d+d2))
          
          A_part_1 = A[, t1] %*% t(A[, j1]) # equivalent to A[i1, 1] * A[s, j1]
          A_part_2 = A[, t2] %*% t(A[, j2]) # equivalent to A[i2, 1] * A[u, j2]
          
          # Now, the outer products of A_part_1 and A_part_2 can be computed
          covlist[[j1]][[j2]] = covlist[[j1]][[j2]] + 
            A_part_1 %*% (var_linear_list[[t1]][[t2]] %*% A_part_2)
        }
      }
    }
  }
  
  
  
  var_logit_hat_calibrated = diag(rep(0, d+d2))
  
  for (j1 in 1:(d+d2)){
    for (j2 in 1:(d+d2)){
      var_logit_hat_calibrated[j1,j2] = (t(A)%*%var_beta_s_tilde%*%A)[j1, j2] +
        matrix(c(beta2_s_tilde, beta1_s_tilde),ncol=(d+d2))%*%covlist[[j1]][[j2]]%*%t(matrix(c(beta2_s_tilde, beta1_s_tilde),ncol=(d+d2)))
    }
  }
  
  beta_calibrated = matrix(c(beta2_s_tilde, beta1_s_tilde), nrow=1)%*%A
  
  # print(beta_calibrated[2:(1+d)])
  # print(var_logit_hat_calibrated[2:(1+d),2:(1+d)])
  # 
  # print(beta_calibrated[1])
  # print(var_logit_hat_calibrated[1,1])
  
  est = c()
  sd = c()
  pval = c()
  or = c()
  ciL = c()
  ciU = c()
  
  for(i in 1:d){
    est = c(est, beta_calibrated[d2+i])
    sd = c(sd, sqrt(diag(var_logit_hat_calibrated)[d2+i]))
    pval = c(pval, 2*min(pnorm(beta_calibrated[1+i], 0, sqrt(diag(var_logit_hat_calibrated)[d2+i])),
                         1-pnorm(beta_calibrated[1+i], 0, sqrt(diag(var_logit_hat_calibrated)[d2+i]))))
    or = c(or, exp(beta_calibrated[1+i]))
    ciL = c(ciL, exp(beta_calibrated[1+i]-z_alpha*sqrt(diag(var_logit_hat_calibrated)[d2+i])))
    ciU = c(ciU, exp(beta_calibrated[1+i]+z_alpha*sqrt(diag(var_logit_hat_calibrated)[d2+i])))
  }
  #print("##########")
  for(j in 1:d2) {
    est = c(est, beta_calibrated[j])
    sd = c(sd, sqrt(var_logit_hat_calibrated[j,j]))
    pval = c(pval, 2*min(pnorm(beta_calibrated[j], 0, sqrt(var_logit_hat_calibrated[j,j])),
                         1-pnorm(beta_calibrated[j], 0, sqrt(var_logit_hat_calibrated[j,j]))))
    or = c(or, exp(beta_calibrated[j]))
    ciL = c(ciL, exp(beta_calibrated[j]-z_alpha*sqrt(var_logit_hat_calibrated[j,j])))
    ciU = c(ciU, exp(beta_calibrated[j]+z_alpha*sqrt(var_logit_hat_calibrated[j,j])))
  }
  
  
  df = data.frame(est, sd, pval, or, ciL, ciU)
  
  return(df)
}

# Updating method which corrects omitted covariate bias to handle missing covariates
get_updating = function(X1_m, X2_m, D_m, X1_v, D_v){
  
  z_alpha = qnorm(0.975)
  
  Rcpp::cppFunction('
  arma::mat inverse(const arma::mat& X) {
    return inv(X);
  }
', depends = "RcppArmadillo")
  
  d = dim(X1_m)[2]
  d2 = dim(X2_m)[2]
  
  if(is.null(d2)){
    d2=1
  }
  
  
  logit_v1 = glm(D_v ~ X1_v, family = binomial)
  b1_tilde = logit_v1$coefficients[1:(1+d)]
  
  logit_m = glm(D_m ~ X1_m + X2_m, family = binomial)
  beta1_hat = logit_m$coefficients[1:(1+d)]
  beta2_hat = logit_m$coefficients[(2+d):(1+d+d2)] 
  var_beta1_hat = vcov(logit_m)[1:(1+d),1:(1+d)]
  var_beta2_hat = vcov(logit_m)[(2+d):(1+d+d2),(2+d):(1+d+d2)]
  
  
  logit_m1 = glm(D_m ~ X1_m, family = binomial)
  b1_hat = logit_m1$coefficients[1:(1+d)]
  
  beta1_tilde = b1_tilde-(b1_hat-beta1_hat)
  
  ###
  
  var_b1_tilde = vcov(logit_v1)
  
  var_delta_hat = 0;
  
  d1 = dim(X1_m)[1]
  
  mu_R = exp(as.matrix(cbind(rep(1, d1), X1_m))%*%as.matrix(b1_hat))/(1+exp(as.matrix(cbind(rep(1, d1), X1_m))%*%as.matrix(b1_hat)))
  mu_F = exp(as.matrix(cbind(rep(1, d1), X1_m))%*%as.matrix(beta1_hat)+X2_m%*%as.matrix(beta2_hat))/(1+exp(as.matrix(cbind(rep(1, d1), X1_m))%*%as.matrix(beta1_hat)+X2_m%*%as.matrix(beta2_hat)))
  
  const_R = diag(rep(0, 1+d))
  
  for (i in 1:length(mu_R)){
    const_R = const_R + mu_R[i]*(1-mu_R[i])*(matrix(c(1, X1_m[i,]), ncol=1)%*%t(matrix(c(1, X1_m[i,]), ncol=1)))
  }
  
  const_R = inverse(const_R/length(mu_R))
  
  const_F = diag(rep(0, 1+d+d2))
  
  for (i in 1:length(mu_F)){
    const_F = const_F + mu_F[i]*(1-mu_F[i])*(matrix(c(1, X1_m[i,], X2_m[i,]), ncol=1)%*%t(matrix(c(1, X1_m[i,], X2_m[i,]), ncol=1)))
  }
  
  const_F = inverse(const_F/length(mu_F))
  
  for (i in 1:length(mu_F)){
    ji = (const_R%*%c(1,X1_m[i,])*(D_m[i]-mu_R[i]))-
      (const_F%*%c(1,X1_m[i,], X2_m[i,])*(D_m[i]-mu_F[i]))[1:(1+d)]
    
    var_delta_hat= var_delta_hat+ji%*%t(ji)
  }
  
  var_delta_hat = var_delta_hat/length(mu_F)^2
  
  ###
  
  cov_beta1_hat_tilde = 0
  
  for (i in 1:length(mu_F)){
    cov_beta1_hat_tilde = cov_beta1_hat_tilde +
      (const_R%*%c(1,X1_m[i,])*(D_m[i]-mu_R[i]))%*%
      t((const_F%*%c(1,X1_m[i,], X2_m[i,])*(D_m[i]-mu_F[i]))[1:d1])
  }
  
  cov_beta1_hat_tilde = cov_beta1_hat_tilde/length(mu_F)^2
  
  ###
  
  var_beta1_tilde = var_delta_hat+var_b1_tilde
  
  W1 = inverse(var_beta1_tilde[2:(1+d),2:(1+d)])%*%inverse(inverse(var_beta1_tilde[2:(1+d),2:(1+d)])+inverse(vcov(logit_m)[2:(1+d),2:(1+d)]))
  W2 = inverse(vcov(logit_m)[2:(1+d),2:(1+d)])%*%inverse(inverse(var_beta1_tilde[2:(1+d),2:(1+d)])+inverse(vcov(logit_m)[2:(1+d),2:(1+d)]))
  
  
  beta1_updating = (matrix(beta1_tilde[2:(1+d)], nrow=1)%*%inverse(var_beta1_tilde[2:(1+d),2:(1+d)])+
                      matrix(beta1_hat[2:(1+d)], nrow=1)%*%inverse(var_beta1_hat[2:(1+d),2:(1+d)]))%*%inverse(inverse(var_beta1_tilde[2:(1+d),2:(1+d)])+inverse(vcov(logit_m)[2:(1+d),2:(1+d)]))
  var_beta1_updating = W1%*%var_beta1_tilde[2:(1+d),2:(1+d)]%*%t(W1)+W2%*%vcov(logit_m)[2:(1+d),2:(1+d)]%*%t(W2)+
    W1%*%(vcov(logit_m)[2:(1+d),2:(1+d)]-cov_beta1_hat_tilde[2:(1+d),2:(1+d)])%*%t(W2)+
    t(W1%*%(vcov(logit_m)[2:(1+d),2:(1+d)]-cov_beta1_hat_tilde[2:(1+d),2:(1+d)])%*%t(W2))
  
  # print(beta1_updating)
  # print(var_beta1_updating)
  # print(beta2_hat)
  # print(var_beta2_hat)
  
  est = c()
  sd = c()
  pval = c()
  or = c()
  ciL = c()
  ciU = c()
  
  for(i in 1:d){
    est = c(est, beta1_updating[i])
    sd = c(sd, sqrt(diag(var_beta1_updating)[i]))
    pval = c(pval, 2*min(pnorm(beta1_updating[i], 0, sqrt(diag(var_beta1_updating)[i])),
                         1-pnorm(beta1_updating[i], 0, sqrt(diag(var_beta1_updating)[i]))))
    or = c(or, exp(beta1_updating[i]))
    ciL = c(ciL, exp(beta1_updating[i]-z_alpha*sqrt(diag(var_beta1_updating)[i])))
    ciU = c(ciU, exp(beta1_updating[i]+z_alpha*sqrt(diag(var_beta1_updating)[i])))
  }
  #print("##########")
  
  var_beta2_hat = as.matrix(var_beta2_hat)
  
  for(j in 1:d2){
    est = c(est, beta2_hat[j])
    sd = c(sd, sqrt(var_beta2_hat[j,j]))
    pval = c(pval, 2*min(pnorm(beta2_hat[j], 0, sqrt(var_beta2_hat[j,j])),
                         1-pnorm(beta2_hat[j], 0, sqrt(var_beta2_hat[j,j]))))
    or = c(or, exp(beta2_hat[j]))
    ciL = c(ciL, exp(beta2_hat[j]-z_alpha*sqrt(var_beta2_hat[j,j])))
    ciU = c(ciU, exp(beta2_hat[j]+z_alpha*sqrt(var_beta2_hat[j,j])))
  }
  
  df = data.frame(est, sd, pval, or, ciL, ciU)
  
  return(df)
}

# New method which adopts calibration model on updating method
get_new = function(X1_m, X2_m, Z_m, D_m, X1_v, Z_v, D_v){
  
  z_alpha = qnorm(0.975)
  
  Rcpp::cppFunction('
  arma::mat inverse(const arma::mat& X) {
    return inv(X);
  }
', depends = "RcppArmadillo")
  
  d = dim(X1_m)[2]
  d1 = dim(X1_m)[1]
  d2 = dim(X2_m)[2]
  
  if(is.null(d2)){
    d2=1
  }
  
  linear_m = lm(X2_m ~ Z_m + X1_m)
  
  if (d2==1){
    gamma0_hat= linear_m$coefficients[1]
    gamma_z_hat = linear_m$coefficients[2]
    gamma_x_hat = linear_m$coefficients[3:(2+d)]
    
    linear_hat = diag(c(gamma_z_hat, rep(1, d)))
    linear_hat[1,2:(1+d)]=gamma_x_hat
  }
  else{
    gamma0_hat= linear_m$coefficients[1,]
    gamma_z_hat = linear_m$coefficients[2:(1+d2),]
    gamma_x_hat = linear_m$coefficients[(2+d2):(1+d+d2),]
    
    linear_hat = diag(rep(1,d+d2))
    linear_hat[1:d2,1:d2] = gamma_z_hat
    linear_hat[1:d2,(1+d2):(d+d2)] = t(gamma_x_hat)
  }
  
  
  #A is inverse of linear_hat
  A= inverse(linear_hat)
  
  var_linear_hat = vcov(linear_m)[-seq(1, 1+(d2-1)*(d+d2+1), by=(d+d2+1)),
                                  -seq(1, 1+(d2-1)*(d+d2+1), by=(d+d2+1))]
  
  var_linear_list = vector("list", d2)
  for (t in 1:d2) {
    var_linear_list[[t]] <- vector("list", d2)
  }
  
  for (t1 in 1:d2) {
    for (t2 in 1:d2) {
      var_linear_list[[t1]][[t2]] = var_linear_hat[(1+(t1-1)*(d+d2)):(t1*(d+d2)),(1+(t2-1)*(d+d2)):(t2*(d+d2))]
    }
    
  }
  
  
  covlist <- vector("list", d+d2)
  for (i in 1:(d+d2)) {
    covlist[[i]] <- vector("list", d+d2)
  }
  
  for (j1 in 1:(d+d2)) {
    for (j2 in 1:(d+d2)) {
      for (t1 in 1:d2) {
        for (t2 in 1:d2) {
          covlist[[j1]][[j2]] = diag(rep(0, d+d2))
          
          A_part_1 = A[, t1] %*% t(A[, j1]) # equivalent to A[i1, 1] * A[s, j1]
          A_part_2 = A[, t2] %*% t(A[, j2]) # equivalent to A[i2, 1] * A[u, j2]
          
          # Now, the outer products of A_part_1 and A_part_2 can be computed
          covlist[[j1]][[j2]] = covlist[[j1]][[j2]] + 
            A_part_1 %*% (var_linear_list[[t1]][[t2]] %*% A_part_2)
        }
      }
    }
  }
  
  
  resid_x2_z = lm(X2_m ~ X1_m+ Z_m)$residuals
  resid_x2_z = as.matrix(resid_x2_z)
  
  logit_m = glm(D_m ~ X1_m + X2_m, family = binomial)
  beta1_hat = logit_m$coefficients[1:(1+d)]
  var_beta1_hat = vcov(logit_m)[1:(1+d),1:(1+d)]
  
  beta2_hat = logit_m$coefficients[(2+d):(1+d+d2)]
  var_beta2_hat = vcov(logit_m)[(2+d):(1+d+d2),(2+d):(1+d+d2)]
  
  logit_m0 = glm(D_m ~ X1_m + Z_m + resid_x2_z, family=binomial)
  beta1_hat0 = logit_m0$coefficients[1:(1+d+d2)]
  
  logit_m1 = glm(D_m ~ X1_m + Z_m, family = binomial)
  b1_hat = logit_m1$coefficients
  
  
  logit_v1 = glm(D_v ~ X1_v + Z_v, family = "binomial")
  b1_tilde = logit_v1$coefficients
  
  
  beta1_tilde = b1_tilde-(b1_hat-beta1_hat0)
  
  ###
  var_b1_tilde = vcov(logit_v1)
  
  var_delta_hat = 0;
  
  
  mu_R = exp(as.matrix(cbind(rep(1, d1), X1_m, Z_m))%*%as.matrix(b1_hat))/(1+exp(as.matrix(cbind(rep(1, d1), X1_m, Z_m))%*%as.matrix(b1_hat)))
  mu_F = exp(as.matrix(cbind(rep(1, d1), X1_m, Z_m))%*%as.matrix(beta1_hat0)+resid_x2_z%*%as.matrix(beta2_hat))/(1+exp(as.matrix(cbind(rep(1, d1), X1_m, Z_m))%*%as.matrix(beta1_hat0)+resid_x2_z%*%as.matrix(beta2_hat)))
  
  const_R=matrix(rep(0,(1+d+d2)^2), ncol=1+d+d2)
  
  for (i in 1:length(mu_R)){
    const_R = const_R + mu_R[i]*(1-mu_R[i])*(matrix(c(1, X1_m[i,], Z_m[i,]), ncol=1)%*%t(matrix(c(1, X1_m[i,], Z_m[i,]), ncol=1)))
  }
  
  
  const_R = inverse(const_R/length(mu_R))
  
  const_F=matrix(rep(0,(1+d+2*d2)^2), ncol=1+d+2*d2)
  
  for (i in 1:length(mu_F)){
    const_F = const_F + mu_F[i]*(1-mu_F[i])*(matrix(c(1, X1_m[i,], Z_m[i,], resid_x2_z[i,]), ncol=1)%*%t(matrix(c(1, X1_m[i,], Z_m[i,], resid_x2_z[i,]), ncol=1)))
  }
  
  const_F = inverse(const_F/length(mu_F))
  
  for (i in 1:length(mu_F)){
    ji = (const_R%*%c(1,X1_m[i,],Z_m[i,])*(D_m[i]-mu_R[i]))-
      (const_F%*%c(1,X1_m[i,], Z_m[i,], resid_x2_z[i,])*(D_m[i]-mu_F[i]))[1:(1+d+d2)]
    
    var_delta_hat= var_delta_hat+ji%*%t(ji)
  }
  
  var_delta_hat = var_delta_hat/length(mu_F)^2
  
  var_beta1_tilde = var_delta_hat+var_b1_tilde
  
  beta_calibrated = matrix(c(beta1_tilde[(2+d):(1+d+d2)], beta1_tilde[2:(1+d)]), nrow=1)%*%A
  
  beta1_tilde1 = beta_calibrated[(1+d2):(d+d2)]
  beta2_tilde1 = beta_calibrated[1:d2]
  
  var_beta1_tilde1_calibrated = diag(rep(0,d+d2))
  
  for (j1 in 1:(d+d2)){
    for (j2 in 1:(d+d2)){
      var_beta1_tilde1_calibrated[j1,j2] = (t(A)%*%var_beta1_tilde[c((2+d):(1+d+d2),2:(1+d)),c((2+d):(1+d+d2),2:(1+d))]%*%A)[j1, j2] +
        matrix(c(beta1_tilde[(2+d):(1+d+d2)], beta1_tilde[2:(1+d)]), nrow=1)%*%covlist[[j1]][[j2]]%*%t(matrix(c(beta1_tilde[(2+d):(1+d+d2)], beta1_tilde[2:(1+d)]), nrow=1))
    }
  }
  
  W1 = var_beta1_tilde1_calibrated[(1+d2):(d+d2),(1+d2):(d+d2)]
  W2 = as.matrix(var_beta1_tilde1_calibrated[1:d2,1:d2])
  
  
  return(c((matrix(beta1_tilde1, nrow=1)%*%inverse(W1)+
              matrix(beta1_hat[2:(1+d)], nrow=1)%*%inverse(var_beta1_hat[2:(1+d),2:(1+d)]))
           %*%inverse(inverse(W1)+inverse(var_beta1_hat[2:(1+d),2:(1+d)])),
           (matrix(beta2_tilde1, nrow=1)%*%inverse(W2)+matrix(beta2_hat, nrow=1)%*%inverse(as.matrix(var_beta2_hat)))
           %*%inverse(inverse(W2)+inverse(as.matrix(var_beta2_hat)))
    )
  )
  
}

function_oh = function(X1, X2, Z, D, method){
  main = !is.na(X2[,1])
  
  covname = c(colnames(X1), colnames(X2))
  
  d = dim(X1)[2]
  d2 = dim(X2)[2]
  
  
  X1_m = as.matrix(X1[main,])
  X1_v = as.matrix(X1[!main,])
  
  X2_m = as.matrix(X2[main,])
  
  Z_m = as.matrix(Z[main,])
  Z_v = as.matrix(Z[!main,])
  
  D_m = D[main]
  D_v = D[!main]
  
  if (method=="gold"){
    X1 = as.matrix(X1)
    X2 = as.matrix(X2)
    df_gold = get_gold(X1, X2, D)
    rownames(df_gold) = covname
    
    return(df_gold)
  }
  else if (method=="main"){
    df_main = get_main(X1_m, X2_m, D_m)
    rownames(df_main) = covname
    
    return(df_main)
  }
  else if (method=="naive"){
    df_naive = get_naive(X1_m, X2_m, D_m, X1_v, Z_v, D_v)
    rownames(df_naive) = covname
    
    return(df_naive)
  }
  else if (method=="calibration"){
    df_calibration = get_calibration(X1_m, X2_m, Z_m, X1_v, Z_v, D_v)
    rownames(df_calibration) = covname
    
    return(df_calibration)
  }
  else if (method=="updating"){
    df_updating = get_updating(X1_m, X2_m, D_m, X1_v, D_v)
    rownames(df_updating) = covname
    
    return(df_updating)
  }
  else if (method=="new"){
    
    z_alpha = qnorm(0.975)
    
    result_new = get_new(X1_m, X2_m, Z_m, D_m, X1_v, Z_v, D_v)
    
    est_new = result_new
    
    beta1_new = result_new[1:d]
    beta2_new = result_new[(1+d):(d+d2)]
    
    jack_knife1 = matrix(rep(0, length(D)*d), ncol= d)
    jack_knife2 = matrix(rep(0, length(D)*d2), ncol= d2)
    
    pb = txtProgressBar(min = 0, max = dim(X1)[1], style = 3)
    
    for (k in 1:dim(X1_m)[1]){
      result_jack = get_new(X1_m[-k,], X2_m[-k,], as.matrix(Z_m[-k,]), D_m[-k], 
                            X1_v, Z_v, D_v)
      
      jack_knife1[k,] = result_jack[1:d]
      jack_knife2[k,] = result_jack[(1+d):(d+d2)]
      setTxtProgressBar(pb, k)
    }
    
    for (l in 1:dim(X1_v)[1]){
      result_jack = get_new(X1_m, X2_m, Z_m, D_m, 
                            X1_v[-l,], as.matrix(Z_v[-l,]), D_v[-l])
      jack_knife1[l+dim(X1_m)[1],] = result_jack[1:d]
      jack_knife2[l+dim(X1_m)[1],] = result_jack[(1+d):(d+d2)]
      setTxtProgressBar(pb, l+dim(X1_m)[1])
    }
    
    # trim_ratio <- 0.05
    # 
    # var_beta1_new <- apply(jack_knife1, 2, function(x) {
    #   # 상하위 trim_ratio 비율로 trimming
    #   trimmed_x <- sort(x)[(ceiling(length(x)*trim_ratio) + 1):(floor(length(x)*(1 - trim_ratio)))]
    #   # trimmed 데이터로 분산 계산
    #   return(var(trimmed_x) * (length(trimmed_x) - 1)^2 / length(trimmed_x))
    # })
    # 
    # # jackknife2에 대해 trimmed standard deviation을 계산
    # var_beta2_new <- apply(jack_knife2, 2, function(x) {
    #   # 상하위 trim_ratio 비율로 trimming
    #   trimmed_x <- sort(x)[(ceiling(length(x)*trim_ratio) + 1):(floor(length(x)*(1 - trim_ratio)))]
    #   # trimmed 데이터로 분산 계산
    #   return(var(trimmed_x) * (length(trimmed_x) - 1)^2 / length(trimmed_x))
    # })
    
    var_beta1_new = apply(jack_knife1, 2, var)*(dim(jack_knife1)[1]-1)^2/dim(jack_knife1)[1]
    var_beta2_new = apply(jack_knife2, 2, var)*(dim(jack_knife2)[1]-1)^2/dim(jack_knife2)[1]
    
    
    
    sd_new = c(sqrt(var_beta1_new), sqrt(var_beta2_new))
    
    pval_new = c()
    or_new = c()
    ciL_new = c()
    ciU_new = c()
    
    for(i in 1:d){
      pval_new = c(pval_new, 2*min(pnorm(beta1_new[i], 0, sqrt(var_beta1_new[i])),
                                   1-pnorm(beta1_new[i], 0, sqrt(var_beta1_new[i]))))
      or_new = c(or_new, exp(beta1_new[i]))
      ciL_new = c(ciL_new, exp(beta1_new[i]-z_alpha*sqrt(var_beta1_new[i])))
      ciU_new = c(ciU_new, exp(beta1_new[i]+z_alpha*sqrt(var_beta1_new[i])))
    }
    
    for(j in 1:d2){
      pval_new = c(pval_new, 2*min(pnorm(beta2_new[j], 0, sqrt(var_beta2_new[j])),
                                   1-pnorm(beta2_new[j], 0, sqrt(var_beta2_new[j]))))
      or_new = c(or_new, exp(beta2_new[j]))
      ciL_new = c(ciL_new, exp(beta2_new[j]-z_alpha*sqrt(var_beta2_new[j])))
      ciU_new = c(ciU_new, exp(beta2_new[j]+z_alpha*sqrt(var_beta2_new[j])))
    }
    
    
    df_new = data.frame(est = est_new, sd = sd_new, 
                        pval = pval_new, or = or_new, 
                        ciL = ciL_new, ciU = ciU_new)
    rownames(df_new) = covname
    close(pb)
    
    return(df_new)
  }
  else{
    print("Error : wrong method")
  }
}


# X1, X2 and Z should be form of data frame with column names. 
# D should be either form of data frame or vector.
# Dimension of X2 and Z should be same.

# Output of the function is (number of covariates) by 6 matrix,
# containing estimation, standard deviation, p-value, odds ratio,
# and confidence interval of odds ratio.
