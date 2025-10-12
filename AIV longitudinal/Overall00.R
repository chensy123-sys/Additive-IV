setwd("./AIV_longitudinal")
# library(tidyverse)
rm(list=ls())
source('par_overall1.R')
source('LIV_overall.R')

# par(mfrow = c(2, 2))
data_LIV <- simLIV_overall(n,D)
# print(data_LIV, n=100)
# hist((data_LIV %>% filter(t==0))$Z)
# hist((data_LIV %>% filter(t==1))$Z)
# data_LIV %>% filter(t==0) %>% summarise(corr=cor(A,Z), mean_A = mean(A))
# data_LIV %>% filter(t==1) %>% summarise(corr=cor(A,Z), mean_A = mean(A))


# fit_LIV <- DML_LIV_overall(data_LIV,action_list = c(0,0))
# fit_LIV_deg <- DML_LIV_overall_deg(data_LIV,action_list = c(0,0))
# predict_LIV2_overall(fit_LIV)
# predict_LIV2_overall(fit_LIV_deg)
# hist(fit_LIV$ETA[,1,1])


seed <- as.integer(2025)
set.seed(2025)

RES <- array(0,dim = c(TM+2,5,Looptime))
for(loop in 1:Looptime){
  print(loop)
  data_LIV <- simLIV_overall(n,D)
  fit_LIV <- DML_LIV_overall(data_LIV,action_list = c(0,0))
  RES[,,loop] <- as.matrix(predict_LIV2_overall(fit_LIV))
}
save(RES, file = "Overall_00.RData")