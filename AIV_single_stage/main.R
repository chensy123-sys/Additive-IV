rm(list = ls());source('AIV.R')
n <- 5000
par(mfrow=c(2,2))
sim_dat <- sim_AIV(n,print_state = F,index = 5)
print(summary(sim_dat))
print(quantile(3*sim_dat$Z+sim_dat$L,c(1/3,2/3)))
# summary((sim_dat$Z-2)^2)
sim_dat %>% summarise(
  cor(Z,A==-1),cor(Z^2,A==0),cor(Z,A==1)
)

result_AIV_unknown <- AIV_unknown_spline(sim_dat,type = 1,K=5)
print(unlist(result_AIV_unknown))
result_AIV_known <- AIV_known_spline(sim_dat,type = 0,K=5,pi = function(Z,L){Z^2})
print(unlist(result_AIV_known))

