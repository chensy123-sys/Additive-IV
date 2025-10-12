source('AIV.R')
par <- list(
  f = function(L) {
    L+sin(3*L)
  },
  
  pZ = function(Z,L){
    -2*Z+2*L
  },
  
  pU = function(U,L){
    3*U-L
  },
  
  merge_prob=0.7,
  
  g0 = function(U,L){
    (U*2-2*L)
  },
  
  g1 = function(U,L){
    (U*2+2*L)
  }
)




# sim_AIV(n*100,par,print_state = T)
set.seed(2025)
REC_unknown <- matrix(0,Looptime,4)
REC_known <- matrix(0,Looptime,4)



i=1
while(i <= Looptime){
  tryCatch({
    sim_dat <- sim_AIV2(n,par,print_state = F)
    result_AIV_unknown <- AIV_unknown_spline(sim_dat,type = "treat")
    result_AIV_known <- AIV_known_spline(sim_dat,type = "treat")
    REC_unknown[i,] <- unlist(result_AIV_unknown)
    REC_known[i,] <- unlist(result_AIV_known)
    if(i==1){
      colnames(REC_unknown) <- names(unlist(result_AIV_unknown))
      colnames(REC_known) <- names(unlist(result_AIV_known))
    }
    print(i)
    i=i+1
  }, warning = function(w) {
    if (grepl("Fitting terminated with step failure", w$message)) {
      message("Step failure detected, retrying... (attempt ", i, ")")
    } else {
      warning(w)  # 其他警告直接抛出
    }
  })
  
}
save(REC_unknown, REC_known, file = "treat3.RData")
# load("REC_results.RData")







