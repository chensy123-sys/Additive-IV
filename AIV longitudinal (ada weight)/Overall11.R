source('par_overall1.R')
source('LIV_overall.R')
seed <- as.integer(2025)
set.seed(2025)

RES <- array(0,dim = c(TM+2,5,Looptime))
for(loop in 1:Looptime){
  print(loop)
  data_LIV <- simLIV_overall(n,D)
  fit_LIV <- DML_LIV_overall(data_LIV,action_list = c(1,1))
  RES[,,loop] <- as.matrix(predict_LIV2_overall(fit_LIV))
}
save(RES, file = "Overall_11.RData")
