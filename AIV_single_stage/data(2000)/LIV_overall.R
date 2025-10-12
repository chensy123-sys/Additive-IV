library(simcausal)
library(tibble)
library(dplyr)
library(SuperLearner)
library(purrr)
library(mgcv)
# library(tidyverse)
library(readr)
library(tidyr)
# install.packages("remotes")
# remotes::install_github("tedwestling/survSuperLearner")
# library(survSuperLearner)

# 模拟真实数据
simLIV_overall <- function(n,D){
  mydat <- tibble(sim(D,n=n,wide = F)) %>%
    select(-contains("OPA"), 
           -contains("deltaA"), 
           -contains("pA0"),
           -contains("U"))
  return(mydat)  
}








# K=2; H=1; mydat=data_LIV; superLearners=c("SL.glm", "SL.randomForest"); action_list=c(0,0); time=TM;k=2
# Debiased machine learning with longituidinal instrumental variable
DML_LIV_overall2 <- function(mydat,
                    action_list,
                    K=2, H=2, 
                    superLearners=c("SL.gam"),
                    print_state=F){
  TM <- max(mydat$t)-1
  n <- length(unique(mydat$ID))
  
  
  # 随机切分数据
  make_id_folds <- function(mydat, id_var = "ID", K = 2, seed = NULL) {
    
    if(K==1){
      ids <- unique(mydat[[id_var]])
      folds_id <- list(list(train_ids = 1:length(ids), test_ids = 1:length(ids)))
    }else{
      if (!is.null(seed)) set.seed(seed)
      
      ids <- unique(mydat[[id_var]])
      folds <- sample(rep(1:K, length.out = length(ids)))  # 随机打乱并分组
      id_folds <- split(ids, folds)  # list: fold 1, fold 2, ..., fold K
      
      folds_id <- lapply(seq_along(id_folds), function(k) {
        test_ids <- id_folds[[k]]
        train_ids <- setdiff(ids, test_ids)
        list(train_ids = train_ids, test_ids = test_ids)
      })
    }
    # a=list(list(a=11))
    # a
    return(folds_id)
  }
  
  build_fY_history <- function(mydat, time, H){
    fY <- mydat %>% filter(t == time) %>% select(L, Z, A)
    
    h<-1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fY <- bind_cols(fY, temp)
      } else {
        break
      }
      h <- h+1
    }  
    return(fY)
  }
  
  build_fA_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L,Z)
    
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    return(fA)
  }
  
  build_fZ_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L)
    
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    return(fA)
  }
  
  # smooth_all <- function(response, data, bs_type = "tp") {
  #   # 获取自变量名称（去掉响应变量）
  #   vars <- setdiff(names(data), response)
  #   if (length(vars) == 0) stop("No predictors found in data.")
  # 
  #   # 拼接平滑项
  #   # smooth_terms <- paste0("s(", vars, ", bs='", bs_type, "')", collapse = " + ")
  #   smooth_terms <- paste0("s(", vars, ", bs='", bs_type, "', k=5)", collapse = " + ")
  #   
  #   # 返回公式对象
  #   as.formula(paste(response, "~", smooth_terms))
  # }
  
  smooth_all <- function(response, data, bs_type = "tp") {
    # 获取自变量名称（去掉响应变量）
    vars <- setdiff(names(data), response)
    if (length(vars) == 0) stop("No predictors found in data.")

    # 拼接多维平滑项
    # 例如: s(X1, X2, X3, bs='tp')
    smooth_term <- paste0("s(", paste(vars, collapse = ", "), ", bs='", bs_type, "')")

    # 返回公式对象
    as.formula(paste(response, "~", smooth_term))
  }
  
  
  time <- max(mydat$t)-1 # start from the last time  
  
  
  ETA <- array(0,dim = c(n,K,TM+2))
  ETA[,,TM+2] <- matrix(rep(mydat[which(mydat$t==time+1),]$Y,K),n,K)
  # if(K>1){
  folds <- make_id_folds(mydat, id_var = "ID", K = K)
  # }
  while (time >=0){
    # time 时间的action
    
    for (k in 1:K) {
      cat(time,k,"\n")
      # if(K==1){
      #   train <- seq(1,nrow(mydat))
      #   test <- seq(1,nrow(mydat))
      # }else{
      train <- folds[[k]]$train_ids
      test <- folds[[k]]$test_ids
      # }
      oY <- ETA[,k,time+2]
      # fY <- build_fY_history(mydat,time,H)
      



      # fA <- build_fA_history(mydat,time,H)
      oZ <- mydat[which(mydat$t==time),]$Z
      # fZ <- build_fZ_history(mydat,time,H)
      his <- build_fZ_history(mydat,time,H)
      

      trt <- tryCatch(
        {
          action_list[[time + 1]](his)  # 尝试调用函数
          # print(action_list)
          # print(f)
          # return(f)
        },
        error = function(e) {
          # message("Function call failed, falling back to raw object.")
          action_list[time + 1]     # 如果失败，返回原始对象
        }
      )
      # print('hello!')
      # print(trt)
      # oA <- as.numeric(mydat[which(mydat$t==time),]$A == trt)
      oA <- mydat[which(mydat$t==time),]$A
      # print(table(oA))
      
      his_A <- his %>% mutate(oA=oA)
      his_Z <- his_A %>% mutate(oZ=oZ)
      
      
      
      famil <- if (length(unique(oZ[train])) <= 2) {
        binomial()
      } else {
        gaussian()
      }
      
      
      
      
      # fit_delta <- SuperLearner(oA[train],his[train,],
      #                           family = binomial(), SL.library = superLearners)
      # fit_rho <- SuperLearner(oZ[train],his[train,],
      #                         family = famil, SL.library = superLearners)
      # fit_kappa <- SuperLearner(oZ[train], his_A[train,],
      #                           family = famil, SL.library = superLearners)
      # fit_eta <- SuperLearner(oY[train],his_A[train,],
      #                         family = gaussian(), SL.library = superLearners)
      # # fit_zeta <- SuperLearner(oY[train], his_Z[train,],
      # #                          family = gaussian(), SL.library = superLearners)
      # fit_zeta <- SuperLearner(oY[train] * oZ[train], his_A[train,],
      #                          family = gaussian(), SL.library = superLearners)
      
      
      # delta <- predict(fit_delta, newdata = his,type='response')$pred
      # rho <- predict(fit_rho, newdata = his,type='response')$pred
      # kappa <- predict(fit_kappa, newdata = his_A %>% mutate(oA=1),type='response')$pred
      # # zeta <- predict(fit_zeta, newdata = his_Z %>% mutate(oA=1,oZ=1),
      # #                 type='response')$pred * delta * kappa
      # zeta <- predict(fit_zeta, newdata = his_A %>% mutate(oA=1),
      #                 type='response')$pred * delta
      # kappa <- kappa * delta - delta * rho
      # # kappa <- pmax(abs(kappa),quantile(abs(kappa),0.02)) * ifelse(kappa>=0,1,-1)
      # kappa <- pmin(pmax(kappa,quantile(kappa,0.1)), quantile(kappa,0.9))
      # eta <- predict(fit_eta, newdata = his_A %>% mutate(oA=1),
      #                          type='response')$pred * delta
      # gamma <- (zeta - rho * eta) / kappa
      
      
      
      

      # if(time==0){
      #   fom <- as.formula("o ~ L")
      #   fom_A <- as.formula("o ~ L*oA")
      #   fom_Z <- smooth_all("o",his_Z[train,]) # o~L,oA,Z
      # }else{
      fom <- smooth_all("o",his[train,])
      fom_A <- smooth_all("o",his_A[train,])
      fom_Z <- smooth_all("o",his_Z[train,])
      # }
      
      # print(fom)
      # print(fom_A)
      # print(fom_Z)

      fit_delta <- gam(fom, data = his[train,] %>%
                             mutate(o = oA[train]),family = binomial()) # A~L
      fit_rho   <- gam(fom, data = his[train,] %>%
                             mutate(o = oZ[train]),family = famil) # Z~L


      fit_kappa <- gam(fom_A, data = his_A[train,] %>%
                             mutate(o = oZ[train]),family = famil) # Z~L,A
      fit_eta   <- gam(fom_A, data = his_A[train,] %>%
                             mutate(o = oY[train]),family = gaussian()) # Y~L,A
      fit_zeta  <- gam(fom_A, data = his_A[train,] %>%
                             mutate(o = oY[train] * oZ[train]),family = gaussian()) # Y*Z~L,A
      # fit_zeta  <- gam(fom_Z, data = his_Z[train,] %>%
      #                    mutate(o = oY[train]),family = gaussian())




      delta <- as.vector(predict(fit_delta, newdata = his,type='response'))
      rho <- as.vector(predict(fit_rho, newdata = his,type='response'))
      kappa <- as.vector(predict(fit_kappa, newdata = his_A %>% mutate(oA=1),type='response'))
      zeta <- as.vector(predict(fit_zeta, newdata = his_A %>% mutate(oA=trt),
                                type='response')) * (delta + trt -1)
      # zeta <- as.vector(predict(fit_zeta, newdata = his_Z %>% mutate(oA=1,oZ=1),
      #                           type='response')) * delta * kappa
      kappa <- kappa * delta - delta * rho
      # kappa <- pmax(abs(kappa),quantile(abs(kappa),0.02)) * ifelse(kappa>=0,1,-1)
      kappa <- pmin(pmax(kappa,quantile(kappa,0.1)), quantile(kappa,0.9))

      eta <- as.vector(predict(fit_eta, newdata = his_A %>% mutate(oA=trt),
                               type='response')) * (delta +trt -1)
      gamma <- (zeta - rho * eta) / kappa

      
      if(print_state){
        # cat(names(his_A[train,]),'\n',names(his_Z[train,]),'\n')
        hist(kappa,breaks=30)
      }
      
      est_aug <- 1/kappa * (oZ - rho) * ((oA+trt-1) * oY - eta) + (
        1 - (oZ - rho) * (oA - delta) / kappa
      ) * gamma
      
      ETA[,k,time+1] <- est_aug
    }
    if(print_state){
      print(time) 
    }
    time =time-1
  }
  return(list(
    ETA = ETA,
    folds = folds,
    input_data = mydat,
    Action = action_list, 
    H = H, 
    SupersuperLearners = superLearners
  ))
}

























































# K=2; H=1; mydat=data_LIV; superLearners=c("SL.glm", "SL.randomForest"); action_list=c(0,0); time=TM;k=2
# Debiased machine learning with longituidinal instrumental variable
DML_LIV_overall <- function(mydat,
                    action_list,
                    K=2, H=2, 
                    superLearners=c("SL.gam"),
                    print_state=F){
  TM <- max(mydat$t)-1
  n <- length(unique(mydat$ID))
  
  
  # 随机切分数据
  make_id_folds <- function(mydat, id_var = "ID", K = 2, seed = NULL) {
    
    if(K==1){
      ids <- unique(mydat[[id_var]])
      folds_id <- list(list(train_ids = 1:length(ids), test_ids = 1:length(ids)))
    }else{
      if (!is.null(seed)) set.seed(seed)
      
      ids <- unique(mydat[[id_var]])
      folds <- sample(rep(1:K, length.out = length(ids)))  # 随机打乱并分组
      id_folds <- split(ids, folds)  # list: fold 1, fold 2, ..., fold K
      
      folds_id <- lapply(seq_along(id_folds), function(k) {
        test_ids <- id_folds[[k]]
        train_ids <- setdiff(ids, test_ids)
        list(train_ids = train_ids, test_ids = test_ids)
      })
    }
    # a=list(list(a=11))
    # a
    return(folds_id)
  }
  
  build_fY_history <- function(mydat, time, H){
    fY <- mydat %>% filter(t == time) %>% select(L, Z, A)
    
    h<-1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fY <- bind_cols(fY, temp)
      } else {
        break
      }
      h <- h+1
    }  
    return(fY)
  }
  
  build_fA_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L,Z)
    
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L, Z) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    return(fA)
  }
  
  build_fZ_history <- function(mydat, time, H) {
    fA <- mydat[which(mydat$t==time),] %>% select(L)
    
    h <- 1
    while (h <=H-1) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(L) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    h <- 1
    while(h <=H) {
      past_time <- time - h
      if (past_time >= 0) {
        temp <- mydat %>%
          filter(t == past_time) %>%
          select(Z, A) %>%
          rename_with(~ paste0(., "_lag", h))
        fA <- bind_cols(fA, temp)
      } else {
        break
      }
      h=h+1
    }
    
    return(fA)
  }
  
  # smooth_all <- function(response, data, bs_type = "tp") {
  #   # 获取自变量名称（去掉响应变量）
  #   vars <- setdiff(names(data), response)
  #   if (length(vars) == 0) stop("No predictors found in data.")
  # 
  #   # 拼接平滑项
  #   # smooth_terms <- paste0("s(", vars, ", bs='", bs_type, "')", collapse = " + ")
  #   smooth_terms <- paste0("s(", vars, ", bs='", bs_type, "', k=5)", collapse = " + ")
  #   
  #   # 返回公式对象
  #   as.formula(paste(response, "~", smooth_terms))
  # }
  
  smooth_all <- function(response, data, bs_type = "tp") {
    # 获取自变量名称（去掉响应变量）
    vars <- setdiff(names(data), response)
    if (length(vars) == 0) stop("No predictors found in data.")

    # 拼接多维平滑项
    # 例如: s(X1, X2, X3, bs='tp')
    smooth_term <- paste0("s(", paste(vars, collapse = ", "), ", bs='", bs_type, "')")

    # 返回公式对象
    as.formula(paste(response, "~", smooth_term))
  }
  
  
  time <- max(mydat$t)-1 # start from the last time  
  
  
  ETA <- array(0,dim = c(n,K,TM+2))
  ETA[,,TM+2] <- matrix(rep(mydat[which(mydat$t==time+1),]$Y,K),n,K)
  if(K>1){
    folds <- make_id_folds(mydat, id_var = "ID", K = K)
  }
  while (time >=0){
    # time 时间的action
    
    for (k in 1:K) {
      cat(time,k,"\n")
      if(K==1){
        train <- seq(1,nrow(mydat))
        test <- seq(1,nrow(mydat))
      }else{
        train <- folds[[k]]$train_ids
        test <- folds[[k]]$test_ids
      }
      oY <- ETA[,k,time+2]
      # fY <- build_fY_history(mydat,time,H)
      



      # fA <- build_fA_history(mydat,time,H)
      oZ <- mydat[which(mydat$t==time),]$Z
      # fZ <- build_fZ_history(mydat,time,H)
      his <- build_fZ_history(mydat,time,H)
      

      trt <- tryCatch(
        {
          action_list[[time + 1]](his)  # 尝试调用函数
          # print(action_list)
          # print(f)
          # return(f)
        },
        error = function(e) {
          # message("Function call failed, falling back to raw object.")
          action_list[time + 1]     # 如果失败，返回原始对象
        }
      )
      # print('hello!')
      # print(trt)
      # oA <- as.numeric(mydat[which(mydat$t==time),]$A == trt)
      oA <- mydat[which(mydat$t==time),]$A
      # print(table(oA))
      
      his_A <- his %>% mutate(oA=oA)
      his_Z <- his %>% mutate(oZ=oZ)
      his_ZA <- his_A %>% mutate(oZ=oZ)
      
      
      
      famil <- if (length(unique(oZ[train])) <= 2) {
        binomial()
      } else {
        gaussian()
      }
      
      
      
      
      
      

      # if(time==0){
      #   fom <- as.formula("o ~ L")
      #   fom_A <- as.formula("o ~ L*oA")
      #   fom_Z <- smooth_all("o",his_Z[train,]) # o~L,oA,Z
      # }else{
      fom <- smooth_all("o",his[train,])
      fom_A <- smooth_all("o",his_A[train,])
      fom_Z <- smooth_all("o",his_Z[train,])
      fom_ZA <- smooth_all("o",his_ZA[train,])
      # }
  
      fit_delta <- gam(fom, data = his[train,] %>%
                             mutate(o = oA[train]),family = binomial()) # A ~ L
      delta <- as.vector(predict(fit_delta, newdata = his,type='response'))


      fit_pi <- gam(fom_Z, data = his_Z[train,] %>%
                             mutate(o = oA[train]),family = binomial()) # A ~ L,Z
      pi <- as.vector(predict(fit_pi, newdata = his_Z,type='response'))


      fit_kappa <- gam(fom, data = his[train,] %>%
                             mutate(o = (pi[train]-delta[train])*(2*oA[train] - pi[train]-delta[train])),
                             family = gaussian) # (E[A|Z,L] - E[A|L])^2 ~ L,A
      kappa_orig <- as.vector(predict(fit_kappa, newdata = his,type='response'))
      
      # fit_kappa <- gam(fom, data = his[train,] %>%
      #                        mutate(o = (pi[train]-delta[train])**2),
      #                        family = Gamma(link = "log")) # (E[A|Z,L] - E[A|L])^2 ~ L
      # kappa_orig <- as.vector(predict(fit_kappa, newdata = his,type='response'))
      # print(summary(kappa))
      kappa <- pmax(kappa_orig,quantile(kappa_orig,0.1))
      
      
      fit_eta   <- gam(fom_A, data = his_A[train,] %>%
                             mutate(o = oY[train]),family = gaussian()) # Y ~ L,A
      eta <- as.vector(predict(fit_eta, newdata = his_A %>% mutate(oA=trt),
                               type='response')) * (delta + trt -1)


      fit_zeta  <- gam(fom_A, data = his_A[train,] %>%
                             mutate(o = oY[train] * pi[train]),family = gaussian()) # Y*pi ~ L,A
      zeta <- as.vector(predict(fit_zeta, newdata = his_A %>% mutate(oA=trt),
                                type='response')) * (delta + trt -1)
      
      
      fit_xi  <- gam(fom_ZA, data = his_ZA[train,] %>%
                             mutate(o = oY[train]),family = gaussian()) # Y~L,A,Z
      xi <- as.vector(predict(fit_zeta, newdata = his_ZA %>% mutate(oA=trt),
                                type='response')) * (pi + trt -1)
      
      gamma <- (zeta - delta * eta) / kappa
      # cat(
      #   'OR:',
      #   mean(gamma),
      #   mean((zeta - delta * eta) / kappa_orig),
      #   ',  IPW:',
      #   mean((pi - delta) / kappa * (oA + trt -1) * oY),
      #   mean((pi - delta) / kappa_orig * (oA + trt -1) * oY),'\n'
      # )

      if(print_state){
        # cat(names(his_A[train,]),'\n',names(his_Z[train,]),'\n')
        hist(kappa,breaks=30)
      }
      
      est_aug <- (pi - delta) / kappa * (oA + trt -1) * oY +
            (oA - pi) / kappa * xi - (oA - delta) / kappa * eta + 
            (kappa_orig + (oA - pi)^2 - (oA - delta)^2) * gamma / kappa
      # print(mean)
      ETA[,k,time+1] <- est_aug
    }
    if(print_state){
      print(time) 
    }
    time =time-1
  }
  return(list(
    ETA = ETA,
    folds = folds,
    input_data = mydat,
    Action = action_list, 
    H = H, 
    SupersuperLearners = superLearners
  ))
}













































predict_LIV2_overall <- function(fit_LIV){
  
  predict_LIV <- function(fit_LIV, time){
    n <- length(unique(fit_LIV$input_data$ID))
    # print(n)
    K <- length(fit_LIV$folds)
    # print(K)
    ETA_merge <- rep(0,n)
    # print(length(ETA_merge))
    for (k in 1:K) {
      index <- fit_LIV$folds[[k]]$test_ids
      # print(length(index))
      # print(k)
      # print(time+1)
      # print(dim(fit_LIV$ETA))
      # print(fit_LIV$ETA[index,,])
      ETA_merge[index] <- fit_LIV$ETA[index,k,time+1]
    }
    # print(length(ETA_merge))
    ATE <- mean(ETA_merge)
    SD <- sd(ETA_merge)/sqrt(n)
    return(list(time=time,
                est=ATE, 
                std=SD
    ))
  }
  
  res <- unlist(predict_LIV(fit_LIV,time=0))
  # print('nihao!')
  for (time in 1:(dim(fit_LIV$ETA)[3]-1)) {
    res <- rbind(res,unlist(predict_LIV(fit_LIV,time)))  
  }
  res <- as.data.frame(res)
  rownames(res) <- NULL
  res <- res %>%
    mutate(
      lower = est - 1.96 * std,
      upper = est + 1.96 * std
    )
  
  return(res)
}





ATE_overall <- function(D,TM,PROB,accuracy=1e5){
  generate_action_from_vector <- function(treatment_vec, varname = "A") {
    nodes=list()
    for (i in 1: length(treatment_vec)) {
      temp=treatment_vec[i]
      if(temp==1){
        nodes = c(nodes, node(varname, t = i-1, distr = "rconst", const = 1))  
      }else{
        nodes = c(nodes, node(varname, t = i-1, distr = "rconst", const = 0))  
      }
      
      # nodes = c(nodes, node(varname, t = i-1, distr = "rbern", prob = 0))
    }
    
    return(nodes)
  }
  
  PROB<-generate_action_from_vector(PROB)
  ate <- rep(0,TM+2) 
  for(s in 0:(TM)){
    prob = PROB[(s+1):length(PROB)]
    
    Dact <- D + action(name = "do", nodes = prob)
    dat_do <- sim(Dact, n = accuracy, actions = "do",wide = T)$do
    ate[s+1] <- mean(dat_do[[paste('Y',TM+1,sep = '_')]])
    # ate[s+1] <- eval.target(DAG = set.targetE(Dact, outcome=paste("Y", TM+1, sep = "_"), param="do"), 
    #                         n=accuracy,actions = "do")$res
  }
  ate[TM+2] <- mean(simLIV_overall(accuracy,D)$Y,na.rm = T)
  return(ate)
}


plot_LIV_overall <- function(res,D,PROB){
  TM<-nrow(res)-2
  res$ate=ATE_overall(D,TM,PROB)
  fig <- ggplot(res, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.2) +
    geom_line(aes(y = est, color = "Estimate"), linewidth = 1) +
    geom_point(aes(y = est, color = "Estimate")) +
    geom_line(aes(y = ate, color = "True ATE"), linewidth = 1) +
    geom_point(aes(y = ate, color = "True ATE")) +
    scale_color_manual(
      name = "Lines",
      values = c("Estimate" = "blue", "True ATE" = "red")
    ) +
    scale_fill_manual(
      name = "Confidence Band",
      values = c("95% CI" = "blue")
    ) +
    labs(
      x = "Time", y = "Estimate", 
      title = "Estimate with 95% Confidence Band"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  return(fig)
}
