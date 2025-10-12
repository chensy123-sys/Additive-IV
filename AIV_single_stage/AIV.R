rm(list = ls())
library(dplyr);library(tibble);library(randomForest)
library(splines)
library(SuperLearner)
# library(npiv)
# library(ranger)
library(mgcv)
n <- 5000
Looptime <- 1000



# fit_gam <- function(formula, data, family, max_try = 5) {
#   for (i in 1:max_try) {
#     result <- tryCatch({
#       m <- gam(formula, data = data, family = family, method="REML")
#       m  # 返回模型
#     }, warning = function(w) {
#       if (grepl("Fitting terminated with step failure", w$message)) {
#         message("Step failure detected, retrying... (attempt ", i, ")")
#         return(NULL)  # 返回 NULL，下一轮重试
#       } else {
#         warning(w)  # 其他警告直接抛出
#       }
#     }, error = function(e) {
#       message("Error encountered: ", e$message)
#       return(NULL)
#     })
#     if (!is.null(result)) return(result)
#   }
#   stop("Failed to fit GAM after ", max_try, " attempts")
# }




sim_AIV <- function(n,par,print_state=T){
  mydat <- tibble(
    L = runif(n,-1,1),
    U = runif(n,-1,1)
  ) %>%
    mutate(Z = par$f(L) + rnorm(n,sd=2),
           pZ=pnorm(par$pZ(Z,L)),
           pU=pnorm(par$pU(U,L)),
           A = rbinom(n,1,pZ*par$merge_prob+pU*(1-par$merge_prob)),
           Y0 = par$g0(U,L),
           Y1 = par$g1(U,L),
           Y = Y0 + A*(Y1-Y0) + rnorm(n,sd=1)
    )
  if(print_state){
    suppressMessages({
      attach(mydat)
      cat(cor(U,A),'\n',cor(Z,A))
      # par(mfrow=c(2,2))
      # plot(L[A==1],Y1[A==1],col='blue')
      # points(L[A==0],Y0[A==0],col='red')
      # plot(L,Y1,col='blue')
      # points(L,Y0,col='red')
      # plot(L,Z)
      par(mfrow=c(2,2))
      cat('\n True ATE:',mean(Y1-Y0),',',mean(Y1),',',mean(Y0))
      cat('\n Treat as RCT:',mean(Y[A==1])-mean(Y[A==0]))
      cat('\n Adjust for L:',mean(predict(lm(Y~A*ns(L,df = 10),data=mydat),
                                          newdata = data.frame(A=rep(1,n),L=L)))-
            mean(predict(lm(Y~A*ns(L,df = 10),data=mydat),
                         newdata = data.frame(A=rep(0,n),L=L))),'\n'
      )
      detach(mydat)
    })
    
  }
  return(mydat%>% select(L,Z,A,Y))
}


sim_AIV2 <- function(n,par,print_state=T){
  mydat <- tibble(
    L = runif(n,-1,1),
    U = runif(n,-1,1)
  ) %>%
    mutate(Z = par$f(L) + rnorm(n,sd=2),
          #  pZ=pnorm(par$pZ(Z,L)),
          #  pU=pnorm(par$pU(U,L)),
           A = rbinom(n,1,plogis(Z-L+U)),
           Y0 = par$g0(U,L),
           Y1 = par$g1(U,L),
           Y = Y0 + A*(Y1-Y0) + rnorm(n,sd=1)
    )
  if(print_state){
    suppressMessages({
      attach(mydat)
      cat(cor(U,A),'\n',cor(Z,A))
      # par(mfrow=c(2,2))
      # plot(L[A==1],Y1[A==1],col='blue')
      # points(L[A==0],Y0[A==0],col='red')
      # plot(L,Y1,col='blue')
      # points(L,Y0,col='red')
      # plot(L,Z)
      par(mfrow=c(2,2))
      cat('\n True ATE:',mean(Y1-Y0),',',mean(Y1),',',mean(Y0))
      cat('\n Treat as RCT:',mean(Y[A==1])-mean(Y[A==0]))
      cat('\n Adjust for L:',mean(predict(lm(Y~A*ns(L,df = 10),data=mydat),
                                          newdata = data.frame(A=rep(1,n),L=L)))-
            mean(predict(lm(Y~A*ns(L,df = 10),data=mydat),
                         newdata = data.frame(A=rep(0,n),L=L))),'\n'
      )
      detach(mydat)
    })
    
  }
  return(mydat%>% select(L,Z,A,Y))
}







AIV_unknown_spline <- function(data, K=2,type = 'ATE'){
  n <- nrow(data)
  cross_fit_split_index <- function(n, K = 5) {
    folds <- sample(rep(1:K, length.out = n))
    splits <- vector("list", K)
    
    for (k in 1:K) {
      test_idx <- which(folds == k)
      train_idx <- setdiff(1:n, test_idx)
      splits[[k]] <- list(train_idx = train_idx, test_idx = test_idx)
    }
    
    return(splits)
  }
  
  mydat <- data
  if(type == 'treat'){
    mydat <- mydat %>% mutate(Y=Y*A)
  }else if (type == 'control') {
    mydat <- mydat %>% mutate(Y=Y*(A-1))
  }

  EST_AUG <- matrix(0,nrow(mydat),1)
  EST_OR <- matrix(0,nrow(mydat),1)
  EST_IPW <- matrix(0,nrow(mydat),1)
  split_indices <- cross_fit_split_index(n = nrow(mydat), K)
  for (k in 1:K) {
    train_index <- split_indices[[k]]$train
    test_index <- split_indices[[k]]$test
    
    # superLearners <- c(
    #   "SL.gam",       # 广义线性模型
    #   # "SL.randomForest" # 随机森林
    #   "SL.ranger" # 随机森林
    # )
    # pi_model   <- SuperLearner(mydat[train_index,]$A,
    #                             mydat[train_index,] %>% select(Z,L), 
    #                             family = binomial(), SL.library = superLearners)
    # delta_model <- SuperLearner(mydat[train_index,]$A,
    #                             mydat[train_index,] %>% select(L), 
    #                             family = binomial(), SL.library = superLearners)
    # eta_model   <- SuperLearner(mydat[train_index,]$Y,
    #                             mydat[train_index,] %>% select(L), 
    #                             family = gaussian(), SL.library = superLearners)
    # xi_model   <- SuperLearner(mydat[train_index,]$Y,
    #                             mydat[train_index,] %>% select(Z,L), 
    #                             family = gaussian(), SL.library = superLearners)
    # 
    # pi_pred    <- predict(pi_model, newdata = mydat[test_index,] %>% select(-Y,-A),
    #                       type = 'response')$pred
    # delta_pred <- predict(delta_model, newdata = mydat[test_index,] %>% select(L),
    #                       type = 'response')$pred
    # xi_pred    <- predict(xi_model, newdata = mydat[test_index,] %>% select(-Y,-A))$pred
    # eta_pred   <- predict(eta_model, newdata = mydat[test_index,] %>% select(L))$pred
    # 
    # pi_pred_train <- predict(pi_model, newdata = mydat[train_index,] %>% select(-Y,-A),
    #                          type = 'response')$pred
    # xi_pred_train <- predict(xi_model, newdata = mydat[train_index,] %>% select(-Y,-A))$pred
    # 
    # 
    # kappa_model <- SuperLearner(mydat[train_index,]$A * pi_pred_train * 2 - pi_pred_train ^ 2,
    #                             mydat[train_index,] %>% select(L),
    #                             family = gaussian(), SL.library = superLearners)
    # zeta_model  <- SuperLearner(mydat[train_index,]$Y * pi_pred_train 
    #                             + xi_pred_train * (mydat[train_index,]$A - pi_pred_train),
    #                             mydat[train_index,] %>% select(L),
    #                             family = gaussian(), SL.library = superLearners)
    # 
    # kappa_pred <- predict(kappa_model, newdata = mydat[test_index,] %>% select(L))$pred
    # zeta_pred <- predict(zeta_model, newdata = mydat[test_index,] %>% select(L))$pred
    
    
    
    pi_model <- gam(A ~ s(Z, L, bs = "tp"),data = mydat[train_index,], family = binomial())
    delta_model <- gam(A ~ s(L, bs = "tp"),data = mydat[train_index,], family = binomial())
    xi_model <- gam(Y ~ s(Z, L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    eta_model <- gam(Y ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())

    pi_pred    <- as.vector(predict(pi_model, newdata = mydat[test_index,] %>% select(-Y,-A),
                                    type = 'response'))
    delta_pred <- as.vector(predict(delta_model, newdata = mydat[test_index,] %>% select(L),
                                    type = 'response'))
    xi_pred    <- as.vector(predict(xi_model, newdata = mydat[test_index,] %>% select(-Y,-A),type = 'response'))
    eta_pred   <- as.vector(predict(eta_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))

    pi_pred_train    <- as.vector(predict(pi_model, newdata = mydat[train_index,] %>% select(-Y,-A),
                                          type = 'response'))
    xi_pred_train    <- as.vector(predict(xi_model, newdata = mydat[train_index,] %>% select(-Y,-A),
                                          type = 'response'))



    kappa_model <- gam((A * pi_pred_train * 2 - pi_pred_train ^ 2) ~ s(L, bs = "tp"),
                       data = mydat[train_index,], family = gaussian())
    zeta_model <- gam((Y * pi_pred_train + xi_pred_train * (A - pi_pred_train)) ~ s(L, bs = "tp"),
                      data = mydat[train_index,], family = gaussian())

    kappa_pred <- as.vector(predict(kappa_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))
    zeta_pred <- as.vector(predict(zeta_model, newdata = mydat[test_index,] %>% select(L),type = 'response'))
    
    
    
    
    
    den <- kappa_pred - delta_pred ^2
    den <- pmax(den,quantile(abs(den),0.1))
    # hist(den, breaks=30)
    
    
    gamma_pred <- (zeta_pred - eta_pred * delta_pred)/den
    # est_aug <- gamma_pred + 1/den * (
    #   mydat[test_index,]$Y * pi_pred - zeta_pred
    #   +xi_pred * (mydat[test_index,]$A - pi_pred)
    #   -delta_pred * (mydat[test_index,]$Y - eta_pred)
    #   -eta_pred*(mydat[test_index,]$A - delta_pred)
    # ) - gamma_pred/(den)* (
    #   mydat[test_index,]$A * pi_pred - kappa_pred
    #   + pi_pred * (mydat[test_index,]$A - pi_pred)
    #   - 2 * delta_pred * (mydat[test_index,]$A - delta_pred)
    # )
    
    est_aug <- (pi_pred - delta_pred)/den * mydat[test_index,]$Y + gamma_pred/den * (
      den + (mydat[test_index,]$A -pi_pred)^2 - (mydat[test_index,]$A - delta_pred)^2
    ) + 1/den * (
      xi_pred * (mydat[test_index,]$A - pi_pred) - eta_pred * (mydat[test_index,]$A - delta_pred)
    )
    
    
    EST_AUG[test_index,]=est_aug
    EST_OR[test_index,]=gamma_pred
    EST_IPW[test_index,]=mydat[test_index,]$Y * (pi_pred-delta_pred)/den
    
  }
  RES = list(est_aug = mean(EST_AUG), sd_aug = sd(EST_AUG)/sqrt(n),
             est_or = mean(EST_OR), est_ipw = mean(EST_IPW))
  return(RES)
}





# fixed a pi function
# e.g., pi=function(Z,L){return(Z)}
AIV_known_spline <- function(data, K=2,pi=function(Z,L){return(Z)},type = 'ATE'){
  n <- nrow(data)
  cross_fit_split_index <- function(n, K = 5) {
    folds <- sample(rep(1:K, length.out = n))
    splits <- vector("list", K)
    
    for (k in 1:K) {
      test_idx <- which(folds == k)
      train_idx <- setdiff(1:n, test_idx)
      splits[[k]] <- list(train_idx = train_idx, test_idx = test_idx)
    }
    
    return(splits)
  }
  
  mydat <- data
  if(type == 'treat'){
    mydat <- mydat %>% mutate(Y=Y*A)
  }else if (type == 'control') {
    mydat <- mydat %>% mutate(Y=Y*(A-1))
  }
    
  mydat <- mydat %>% mutate(Z=pi(Z,L))
  EST_AUG <- matrix(0,nrow(mydat),1)
  EST_OR <- matrix(0,nrow(mydat),1)
  EST_IPW <- matrix(0,nrow(mydat),1)
  split_indices <- cross_fit_split_index(n = nrow(mydat), K)
  for (k in 1:K) {
    train_index <- split_indices[[k]]$train
    test_index <- split_indices[[k]]$test
    
    # superLearners <- c(
    #   "SL.glm",       # 广义线性模型
    #   # "SL.randomForest" # 随机森林
    #   "SL.ranger" # 随机森林
    # )
    # rho_model   <- SuperLearner(mydat[train_index,]$Z,
    #                              mydat[train_index,] %>% select(L), 
    #                              family = gaussian(), SL.library = superLearners)
    # delta_model <- SuperLearner(mydat[train_index,]$A,
    #                             mydat[train_index,] %>% select(L), 
    #                             family = binomial(), SL.library = superLearners)
    # eta_model   <- SuperLearner(mydat[train_index,]$Y,
    #                             mydat[train_index,] %>% select(L), 
    #                             family = gaussian(), SL.library = superLearners)
    # kappa_model <- SuperLearner(mydat[train_index,]$A * mydat[train_index,]$Z,
    #                             mydat[train_index,] %>% select(L), 
    #                             family = gaussian(), SL.library = superLearners)
    # zeta_model  <- SuperLearner(mydat[train_index,]$Y * mydat[train_index,]$Z,
    #                            mydat[train_index,] %>% select(L), 
    #                            family = gaussian(), SL.library = superLearners)
    # 
    # rho_pred   <- predict(rho_model, newdata = mydat[test_index,] %>% select(L))$pred
    # delta_pred <- predict(delta_model, newdata = mydat[test_index,] %>% select(L))$pred
    # eta_pred   <- predict(eta_model, newdata = mydat[test_index,] %>% select(L))$pred
    # kappa_pred <- predict(kappa_model, newdata = mydat[test_index,] %>% select(L))$pred
    # zeta_pred  <- predict(zeta_model, newdata = mydat[test_index,] %>% select(L))$pred
    
    
    
    
    
    rho_model <- gam(Z ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    delta_model <- gam(A ~ s(L, bs = "tp"),data = mydat[train_index,], family = binomial())
    eta_model <- gam(Y ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    # kappa_model <- gam(Z ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    kappa_model <- gam(Z ~ s(A,L, bs = "tp"), data = mydat[train_index,], family = gaussian())
    zeta_model <- gam(Y * Z ~ s(L, bs = "tp"),data = mydat[train_index,], family = gaussian())
    
    rho_pred    <- as.vector(predict(rho_model, newdata = mydat[test_index,] %>% select(L)))
    delta_pred <- as.vector(predict(delta_model, newdata = mydat[test_index,] %>% select(L),
                                    type = 'response'))
    eta_pred   <- as.vector(predict(eta_model, newdata = mydat[test_index,] %>% select(L)))
    kappa_pred <- as.vector(predict(kappa_model, newdata = mydat[test_index,] %>% select(L) %>% mutate(A=1))) * delta_pred
    zeta_pred <- as.vector(predict(zeta_model, newdata = mydat[test_index,] %>% select(L)))
    
    
    
    
    den <- kappa_pred - delta_pred * rho_pred
    den <- pmin(pmax(den,quantile(den,0.1)), quantile(den,0.9))
    # hist(den,breaks = 30)
    gamma_pred <- (zeta_pred - eta_pred * rho_pred)/den
    
    # est_aug <- gamma_pred + 1/den * (
    #   mydat[test_index,]$Y * mydat[test_index,]$Z - zeta_pred
    #   -rho_pred * (mydat[test_index,]$Y - eta_pred)
    #   -eta_pred*(mydat[test_index,]$Z - rho_pred)
    # ) - gamma_pred/den* (
    #   mydat[test_index,]$A * mydat[test_index,]$Z - kappa_pred
    #   - rho_pred * (mydat[test_index,]$A - delta_pred)
    #   - delta_pred * (mydat[test_index,]$Z - rho_pred)
    # )
    
    est_aug <- 1/den * (mydat[test_index,]$Z - rho_pred) * (mydat[test_index,]$Y - eta_pred) + (
      1 - (mydat[test_index,]$Z - rho_pred) * (mydat[test_index,]$A - delta_pred) / den
    ) * gamma_pred
    
    
    EST_AUG[test_index,]=est_aug
    EST_OR[test_index,]=gamma_pred
    EST_IPW[test_index,]=mydat[test_index,]$Y * (mydat[test_index,]$Z-rho_pred)/den
    
  }
  RES = list(est_aug = mean(EST_AUG), sd_aug = sd(EST_AUG)/sqrt(n),
             est_or = mean(EST_OR), est_ipw = mean(EST_IPW))
  return(RES)
}