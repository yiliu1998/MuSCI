#### Multi-Source Conformal Inference (MuSCI)
#### Data analysis R script 

#### load R packages
library(MASS)
library(numDeriv)
library(rootSolve)
library(osqp)
library(densratio)
library(dplyr)
library(glmnet)
library(quantreg)
library(randomForest)
library(SuperLearner)

#### Main function MuSIC()
MuSCI <- function(data.train, 
                  data.test, 
                  alpha=0.1, 
                  conf.score="ASR",
                  Y.name, 
                  T.name, 
                  tgt.site.name, 
                  ps.Library="SL.glm",  
                  m.Library="SL.glm", 
                  p1=0.5, 
                  seed=10012) {
  
  #### argument description
  # --  data.train: training data, a data.frame() objective with defined variable names
  # --  data.test: testing data, which should have the same structure as the training data
  # --  alpha: desired coverage level of the prediction intervals
  # --  Y.name: name of outcome in the training and testing data
  # --  T.name: name of the site variable in training data
  # --  tgt.site.name: name of the target site
  # --  ps.Library: SuperLearner library for fitting the propensity score
  # --              use SuperLearner::listWrappers() to see all possible choices
  # --  m.Library: SuperLearner library for fitting the putative cdf of conformal score
  # --              use SuperLearner::listWrappers() to see all possible choices
  # --  p1: proportion of training data used for nuisance function training, default value 0.5
  # --  seed: seed for random sampling when splitting data
  
  #### Before you input the data, we assume
  # -- the is.na() function in base R can help you to find missing outcome values
  
  ###############################################################################
  ##################            Data pre-processing           ###################
  ###############################################################################
  colnames(data.train)[colnames(data.train)==T.name] <- "Site"
  site.names <- unique(data.train$Site) # extract names of sites
  s <- length(site.names) # number of sites
  non.tgt.site.names <- site.names[site.names%in%tgt.site.name==F] # names of non-target sites
  
  data.train$Site[data.train$Site==tgt.site.name] <- 999
  for(k in 2:s) {
    data.train$Site[data.train$Site==non.tgt.site.names[k-1]] <- k*999 
  }
  data.train$Site <- data.train$Site/999
  covar.names <- names(data.train)[names(data.train)%in%c(Y.name, "Site")==F]
  
  D1 <- D2 <- data.frame()
  set.seed(seed)
  
  for(k in 1:s) {
    dat <- data.train[data.train$Site==k,]
    n <- nrow(dat)
    n.d1 <- floor(n*p1)
    n.d2 <- n-n.d1
    I1 <- sample(1:n, n.d1, replace=F)
    I2 <- (1:n)[-I1] 
    D1 <- rbind(D1, dat[I1, ])
    D2 <- rbind(D2, dat[I2, ])
  }
    
  ###############################################################################
  ##################            Pooled sample method           ##################
  ###############################################################################
  D1$R <- as.numeric(is.na(D1[,Y.name])==F) # define the missingness indicator
  D2$R <- as.numeric(is.na(D2[,Y.name])==F)
  D11 <- D1[D1$R==1,]
  D21 <- D2[D2$R==1,]
  
  y = D1$R
  x = D1[,covar.names]
  ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
  pred = as.numeric(predict(ps.fit, D21[,covar.names])$library.predict)
  prob = pmax(pmin(pred, 0.9), 0.1) 
  eta = (1-prob)/prob 
  
  if(conf.score=="ASR") {
    CFS <- ASR(W=data.matrix(D11[,covar.names]), 
               Y=D11[,Y.name],
               W.pred=data.matrix(D21[,covar.names]),
               Y.pred=D21[,Y.name])
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }
  if(conf.score=="localASR") {
    CFS <- localASR(W=data.matrix(D11[,covar.names]), 
                    Y=D11[,Y.name],
                    W.pred=data.matrix(D21[,covar.names]),
                    Y.pred=D21[,Y.name])
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    CFS <- CQR(data=D11, Y=D11[,Y.name],
               data.pred=D21, Y.pred=D21[,Y.name],
               rq.form=rq.form, alpha=alpha)
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }
  
  m <- function(theta) {
    D11 <- D11 %>% mutate(score=as.numeric(CFS1<=theta)) 
    y = D11$score
    x = D11[,covar.names]
    ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=m.Library)
    p_CFS = as.numeric(predict(ps.fit, D2[,covar.names])$library.predict)
    return(p_CFS)
  }
  f <- function(theta){ 
    sum(m(theta)[D2$R==0]-(1-alpha))/nrow(D2[D2$R==0,]) + 
       sum(eta*(I(CFS2<=theta)-m(theta)[D2$R==1]))/nrow(D2[D2$R==1,])
  }
  r_0 <- c(CFS1, CFS2)
  rhat_pool = try(uniroot(f,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
  if(class(rhat_pool)!= "try-error"){
    rhat_pool = rhat_pool
  } else if (f(max(r_0)/2)>0) {
    rhat_pool = try(uniroot(f,c(0,max(r_0)/2),tol=1e-5)$root)
  } else {
    rhat_pool = try(uniroot(f,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
  }
  
  ###############################################################################
  ##################         Use target site data only         ##################
  ###############################################################################
  D1_s1 <- D1[D1$Site==1, ]
  D2_s1 <- D2[D2$Site==1, ]
  D11_s1 <- D1_s1[D1_s1$R==1, ]
  D21_s1 <- D2_s1[D2_s1$R==1, ]
  
  y = D1_s1$R
  x = D1_s1[,covar.names]
  ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
  pred = as.numeric(predict(ps.fit, D21_s1[,covar.names])$library.predict)
  prob = pmax(pmin(pred, 0.9), 0.1) 
  eta_s1 = (1-prob)/prob 
  
  if(conf.score=="ASR") {
    CFS <- ASR(W=data.matrix(D11_s1[,covar.names]), 
               Y=D11_s1[,Y.name],
               W.pred=data.matrix(D21_s1[,covar.names]),
               Y.pred=D21_s1[,Y.name])
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }
  if(conf.score=="localASR") {
    CFS <- localASR(W=data.matrix(D11_s1[,covar.names]), 
                    Y=D11_s1[,Y.name],
                    W.pred=data.matrix(D21_s1[,covar.names]),
                    Y.pred=D21_s1[,Y.name])
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    CFS <- CQR(data=D11_s1, Y=D11_s1[,Y.name],
               data.pred=D21_s1, Y.pred=D21_s1[,Y.name],
               rq.form=rq.form, alpha=alpha)
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }
  
  m_s1 <- function(theta) {
    D11_s1 <- D11_s1 %>% mutate(score=as.numeric(CFS1_s1<=theta)) 
    y = D11_s1$score
    x = D11_s1[,covar.names]
    ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=m.Library)
    p_CFS = as.numeric(predict(ps.fit, D2_s1[,covar.names])$library.predict)
    return(p_CFS)
  }
  
  f_s1 <- function(theta){ 
    sum(m_s1(theta)[D2_s1$R==0]-(1-alpha))/nrow(D2_s1[D2_s1$R==0,])  + 
       sum(eta_s1*(I(CFS2_s1<=theta)-m_s1(theta)[D2_s1$R==1]))/nrow(D2_s1[D2_s1$R==1,])  
  }
  r_0 <- c(CFS1_s1, CFS2_s1)
  rhat_tgSt = try(uniroot(f_s1,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
  if(class(rhat_tgSt)!= "try-error"){
    rhat_tgSt = rhat_tgSt
  } else if (f_s1(max(r_0)/2)>0) {
    rhat_tgSt = try(uniroot(f_s1,c(0,max(r_0)/2),tol=1e-5)$root)
  } else {
    rhat_tgSt = try(uniroot(f_s1,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
  }
  
  ###############################################################################
  ##############       Federated weights and equal weights        ###############
  ###############################################################################
  chi <- rhat_ntgt <- c()
  omega <- eta <- phi <- CFS1 <- CFS2 <- list()
  D10_s1 <- D1_s1[D1_s1$R==0, ]
  D20_s1 <- D2_s1[D2_s1$R==0, ]
  for(k in 2:s) {
    D1_sk <- D1[D1$Site==k, ]
    D2_sk <- D2[D2$Site==k, ]
    D11_sk <- D11[D11$Site==k, ]
    D21_sk <- D21[D21$Site==k, ]
    
    y = D1_sk$R
    x = D1_sk[,covar.names]
    ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
    pred = as.numeric(predict(ps.fit, D21_sk[,covar.names])$library.predict)
    prob = pmax(pmin(pred, 0.9), 0.1) 
    eta_sk = (1-prob)/prob
    eta[[k-1]] <- eta_sk
    
    if(conf.score=="ASR") {
      CFS <- ASR(W=data.matrix(D11_sk[,covar.names]), 
                 Y=D11_sk[,Y.name],
                 W.pred=data.matrix(D21_sk[,covar.names]),
                 Y.pred=D21_sk[,Y.name])
      CFS1_sk <- CFS[[1]]; CFS2_sk <- CFS[[2]]
    }
    if(conf.score=="localASR") {
      CFS <- localASR(W=data.matrix(D11_sk[,covar.names]), 
                      Y=D11_sk[,Y.name],
                      W.pred=data.matrix(D21_sk[,covar.names]),
                      Y.pred=D21_sk[,Y.name])
      CFS1_sk <- CFS[[1]]; CFS2_sk <- CFS[[2]]
    }
    if(conf.score=="CQR") {
      rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
      CFS <- CQR(data=D11_sk, Y=D11_sk[,Y.name],
                 data.pred=D21_sk, Y.pred=D21_sk[,Y.name],
                 rq.form=rq.form, alpha=alpha)
      CFS1_sk <- CFS[[1]]; CFS2_sk <- CFS[[2]]
    }
    CFS1[[k-1]] <- CFS1_sk
    CFS2[[k-1]] <- CFS2_sk
    
    D10_s1M <- as.matrix(D10_s1)
    D11_skM <- as.matrix(D11_sk)
    D20_s1M <- as.matrix(D20_s1)
    D21_skM <- as.matrix(D21_sk)
    ratio <- Estimate_omega_np(x=D11_skM[,covar.names], 
                               x_target=D10_s1M[,covar.names],
                               x.pred=D21_skM[,covar.names], 
                               x_target.pred=D20_s1M[,covar.names]) %>%
      as.numeric()
    omega[[k-1]] <- ratio
    
    m_sk <- function(theta) {
      D11_sk <- D11_sk %>% mutate(score = as.numeric(CFS1_sk<=theta)) 
      y = D11_sk$score
      x = D11_sk[,covar.names]
      ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=m.Library)
      p_CFS_sk = as.numeric(predict(ps.fit, D2_sk[,covar.names])$library.predict)
      return(p_CFS_sk)
    }
    
    f_sk = function(theta){ 
      sum((m_s1(theta)[D2_s1$R==0]-(1-alpha)))/nrow(D2_s1[D2_s1$R==0,]) + 
        sum(ratio*(I(CFS2_sk<=theta)-m_sk(theta)[D2_sk$R==1]))/nrow(D2_sk[D2_sk$R==1,])
    }
    r_0 <- c(CFS1_sk, CFS2_sk) 
    rhat_sk = try(uniroot(f_sk,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
    if(class(rhat_sk)!= "try-error"){
      rhat_sk=rhat_sk
    } else if (f_sk(max(r_0)/2)>0) {
      rhat_sk = try(uniroot(f_sk,c(0,max(r_0)/2),tol=1e-5)$root)
    } else {
      rhat_sk = try(uniroot(f_sk,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
    }
    rhat_ntgt[k-1] <- rhat_sk
    chi[k-1] <- abs(rhat_sk - rhat_tgSt)
  }
  
  ### Calculating federated weights via optimization
  phi_s1 = function(theta) { # phi function of the target site
    ind1 <- which(D2$R==0 & D2$Site==1)
    ind2 <- which(D2$R==1 & D2$Site==1)
    phi <- rep(NA, length(D2$R))
    phi[ind1] <- m_s1(theta)[D2_s1$R==0]-(1-alpha)/mean(D2$Site==1 & D2$R==0)
    phi[ind2] <- eta_s1*(I(CFS2_s1<=theta)-m_s1(theta)[D2_s1$R==1])/mean(D2$Site==1 & D2$R==1)
    phi[-c(ind1,ind2)] <- 0
    return(matrix(phi, ncol=1))
  }
  phi_k <- function(theta, k){
    D1_sk <- D1 %>% filter(Site==k) 
    D2_sk <- D2 %>% filter(Site==k) 
    D11_sk <- D11 %>% filter(Site==k) 
    D21_sk <- D21 %>% filter(Site==k) 
    
    m_sk <- function(theta) {
      D11_sk <- D11_sk %>% mutate(score = as.numeric(CFS1[[k-1]]<=theta))
      y = D11_sk$score
      x = D11_sk[,covar.names]
      ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
      p_CFS_sk = as.numeric(predict(ps.fit, D2_sk[,covar.names])$library.predict)
      return(p_CFS_sk)
    }
    
    ind1 <- which(D2$R==0 & D2$Site==1)
    ind2 <- which(D2$R==1 & D2$Site==k)
    phi <- rep(NA, length(D2$R))
    phi[ind1] <- (m_s1(theta)[D2_s1$R==0]-(1-alpha))/mean(D2$Site==1 & D2$R==0)
    phi[ind2] <- (omega[[k-1]]*(I(CFS2[[k-1]]<=theta)-m_sk(theta)[D2_sk$R==1]))/mean(D2$Site==k & D2$R==1)
    phi[-c(ind1,ind2)] <- 0
    
    inds1 <- which(D2$Site==1)
    phi[inds1] <- as.numeric(phi_s1(theta))[inds1]
    return(matrix(phi, ncol=1))
  }
  
  IF_s1 <- phi_s1(rhat_tgSt)
  IF_diff <- matrix(0, length(IF_s1), s-1)
  for(k in 2:s) {
    IF_diff[,(k-1)] <- phi_k(rhat_tgSt, k)
  }
  
  cvfit = cv.glmnet(x = IF_diff, y = IF_s1) 
  lambda = cvfit$lambda.min
  fit = glmnet(x = IF_diff, y = IF_s1, 
               penalty.factor = chi^2,
               intercept = FALSE,
               alpha = 1,
               lambda = lambda,
               lower.limits = 0,
               upper.limits = 1)
  w1 = coef(fit, s = lambda)[2:s]; w1
  w2 = w1*(s-1)/s; w2 # to ensure at least 1/s weight assigned to target site
  
  fit3 = glmnet(x = IF_diff, y = IF_s1, 
                penalty.factor = chi^2,
                intercept = FALSE,
                alpha = 1,
                lambda = lambda,
                lower.limits = 0,
                upper.limits = 1/s)
  w3 = coef(fit3, s = lambda)[2:s]; w3
  
  ### Use federated weights
  rhat_fed1 = (1-sum(w1))*rhat_tgSt + sum(w1*rhat_ntgt)
  rhat_fed2 = (1-sum(w2))*rhat_tgSt + sum(w2*rhat_ntgt)
  rhat_fed3 = (1-sum(w3))*rhat_tgSt + sum(w3*rhat_ntgt)
  
  ### Use equal weights
  rhat_eqwt = mean(c(rhat_tgSt, rhat_ntgt))
  
  ###############################################################################
  ##############       Prediction on testing data point           ###############
  ###############################################################################
  data.test$R <- as.numeric(is.na(data.test[,Y.name])==F)
  data.test1 <- data.test %>% filter(R==1)
  
  if(conf.score=="ASR") {
    W=data.matrix(data.test1[,covar.names])
    W.pred=data.matrix(data.test[data.test$R==0,covar.names])
    Y=data.test1[,Y.name]
    lambda_seq <- 10^seq(2, -2, by = -.1)
    mu <- glmnet(W, Y, alpha=0, lambda=lambda_seq)
    ridge_cv <- cv.glmnet(W, Y, alpha=0, lambda=lambda_seq)
    pred_y0 <- predict(mu, s=ridge_cv$lambda.min, W.pred)
    
    PI.fed1 <- data.frame(lwr= -rhat_fed1+pred_y0, upr=rhat_fed1+pred_y0)
    PI.fed2 <- data.frame(lwr= -rhat_fed2+pred_y0, upr=rhat_fed2+pred_y0)
    PI.fed3 <- data.frame(lwr= -rhat_fed3+pred_y0, upr=rhat_fed3+pred_y0)
    PI.pool <- data.frame(lwr= -rhat_pool+pred_y0, upr=rhat_pool+pred_y0)
    PI.tgSt <- data.frame(lwr= -rhat_tgSt+pred_y0, upr=rhat_tgSt+pred_y0)
    PI.eqwt <- data.frame(lwr= -rhat_eqwt+pred_y0, upr=rhat_eqwt+pred_y0)
  }
  if(conf.score=="localASR") {
    W=data.matrix(data.test1[,covar.names])
    W.pred=data.matrix(data.test[data.test$R==0,covar.names])
    Y=data.test1[,Y.name]
    lambda_seq <- 10^seq(2, -2, by = -.1)
    mu <- glmnet(W, Y, alpha=0, lambda=lambda_seq)
    ridge_cv <- cv.glmnet(W, Y, alpha=0, lambda=lambda_seq)
    pred_y0 <- predict(mu, s=ridge_cv$lambda.min, W.pred)
    pred_y1 <- predict(mu, s=ridge_cv$lambda.min, W)
    
    ASR0 <- abs(pred_y1-Y)
    Rho <- glmnet(W, ASR0, alpha=0, lambda=lambda_seq)
    Rho_cv <- cv.glmnet(W, ASR0, alpha=0, lambda=lambda_seq)
    rho0 <- predict(Rho, s=Rho_cv$lambda.min, W.pred)
    
    PI.fed1 <- data.frame(lwr= -rhat_fed1*rho0+pred_y0, upr=rhat_fed1*rho0+pred_y0)
    PI.fed2 <- data.frame(lwr= -rhat_fed2*rho0+pred_y0, upr=rhat_fed2*rho0+pred_y0)
    PI.fed3 <- data.frame(lwr= -rhat_fed3*rho0+pred_y0, upr=rhat_fed3*rho0+pred_y0)
    PI.pool <- data.frame(lwr= -rhat_pool*rho0+pred_y0, upr=rhat_pool*rho0+pred_y0)
    PI.tgSt <- data.frame(lwr= -rhat_tgSt*rho0+pred_y0, upr=rhat_tgSt*rho0+pred_y0)
    PI.eqwt <- data.frame(lwr= -rhat_eqwt*rho0+pred_y0, upr=rhat_eqwt*rho0+pred_y0)
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    fit.lw <- rq(rq.form, data=data.test1, tau=alpha/2)
    fit.up <- rq(rq.form, data=data.test1, tau=1-alpha/2)
    
    W.pred=data.test[data.test$R==0,]
    pred_y0.lw <- predict(fit.lw, W.pred)
    pred_y0.up <- predict(fit.up, W.pred)
    
    PI.fed1 <- data.frame(lwr= pred_y0.lw-rhat_fed1, upr= pred_y0.up+rhat_fed1)
    PI.fed2 <- data.frame(lwr= pred_y0.lw-rhat_fed2, upr= pred_y0.up+rhat_fed2)
    PI.fed3 <- data.frame(lwr= pred_y0.lw-rhat_fed3, upr= pred_y0.up+rhat_fed3)
    PI.pool <- data.frame(lwr= pred_y0.lw-rhat_pool, upr= pred_y0.up+rhat_pool)
    PI.tgSt <- data.frame(lwr= pred_y0.lw-rhat_tgSt, upr= pred_y0.up+rhat_tgSt)
    PI.eqwt <- data.frame(lwr= pred_y0.lw-rhat_eqwt, upr= pred_y0.up+rhat_eqwt)
  }
  colnames(PI.fed1) <- colnames(PI.fed2) <- colnames(PI.fed3) <- 
    colnames(PI.pool) <- colnames(PI.tgSt) <- colnames(PI.eqwt) <- c("lower", "upper")
  return(list(PI.fed1=PI.fed1, 
              PI.fed2=PI.fed2, 
              PI.fed3=PI.fed3, 
              PI.pool=PI.pool, 
              PI.tgSt=PI.tgSt, 
              PI.eqwt=PI.eqwt,
              weights1=w1, 
              weights2=w2, 
              weights3=w3, 
              Chi=chi))
}

Estimate_omega_np = function(x, x_target, x.pred, x_target.pred) {
  x.all = rbind(x, x_target)
  z.all = c(rep(1,nrow(x)),rep(0,nrow(x_target)))
  
  fit <- cv.glmnet(x=x.all, y=z.all, nfolds=5, family='binomial')
  src.predict = predict(fit, x.pred, type="response", lambda=fit$lambda.min)
  
  omega = (1-src.predict)/src.predict*nrow(x.pred)/nrow(x_target.pred)
  omega = pmax(pmin(omega, 20), 0.05)
  return(omega)
}

#### Functions for conformal scores: ASR, local ASR and CQR
ASR <- function(W, Y, W.pred, Y.pred) {
  lambda_seq <- 10^seq(2, -2, by = -.1)
  mu <- glmnet(W, Y, alpha=0, lambda=lambda_seq)
  ridge_cv <- cv.glmnet(W, Y, alpha=0, lambda=lambda_seq)
  pred_y <- predict(mu, s=ridge_cv$lambda.min, W.pred)
  pred_y1 <- predict(mu, s=ridge_cv$lambda.min, W)
  
  CFS1 <- as.numeric(abs(pred_y1-Y))
  CFS2 <- as.numeric(abs(pred_y-Y.pred))
  return(CFS=list(CFS1=CFS1, CFS2=CFS2))
}

localASR <- function(W, Y, W.pred, Y.pred) {
  lambda_seq <- 10^seq(2, -2, by = -.1)
  mu <- glmnet(W, Y, alpha=0, lambda=lambda_seq)
  ridge_cv <- cv.glmnet(W, Y, alpha=0, lambda=lambda_seq)
  
  pred_y1 <- predict(mu, s=ridge_cv$lambda.min, W)
  pred_y <- predict(mu, s=ridge_cv$lambda.min, W.pred)
  
  ASR1 <- abs(pred_y1-Y)
  ASR <- abs(pred_y-Y.pred)
  
  Rho <- glmnet(W, ASR1, alpha=0, lambda=lambda_seq)
  Rho_cv <- cv.glmnet(W, ASR1, alpha=0, lambda=lambda_seq)
  rho1 <- predict(Rho, s=Rho_cv$lambda.min, W)
  rho <- predict(Rho, s=Rho_cv$lambda.min, W.pred)
  
  CFS1 <- as.numeric(ASR1/rho1)
  CFS2 <- as.numeric(ASR/rho)
  return(CFS=list(CFS1=CFS1, CFS2=CFS2))
}

CQR <- function(data, Y, data.pred, Y.pred, rq.form, alpha=0.1) {
  fit.lw <- rq(rq.form, data=data, tau=alpha/2)
  fit.up <- rq(rq.form, data=data, tau=1-alpha/2)
  
  pred_y1.lw <- predict(fit.lw, data)
  pred_y1.up <- predict(fit.up, data)
  
  pred_y.lw <- predict(fit.lw, data.pred)
  pred_y.up <- predict(fit.up, data.pred)
  
  CFS1 <- as.numeric(apply(cbind(pred_y1.lw-Y, Y-pred_y1.up), 1, max))
  CFS2 <- as.numeric(apply(cbind(pred_y.lw-Y.pred, Y.pred-pred_y.up), 1, max))
  return(CFS=list(CFS1=CFS1, CFS2=CFS2))
}
