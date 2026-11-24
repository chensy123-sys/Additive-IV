TM <- 1 # maximal time
library(simcausal)
D <- DAG.empty() +
  node("L",t=0, distr = "runif", min = -3, max = 3) +
  node("U",t=0, distr = "runif", min = -3, max = 3) +
  node("Z",t=0, distr = "runif", min = 0.3*L[t]+sin(1.5*L[t]) - 4,
                                 max = 0.3*L[t]+sin(1.5*L[t]) + 4) +
  node("A",t=0, distr = "rbern", prob = pnorm(-3*Z[t]+L[t])*0.8 + pnorm(3*U[t]-L[t])*0.2)+


  node("L",t=1:TM, distr = "runif", min = (A[t-1]-0.5)+0.5*L[t-1]+0.3*U[t-1]-1,
                                    max = (A[t-1]-0.5)+0.5*L[t-1]+0.3*U[t-1]+1) +
  node("U",t=1:TM, distr = "runif", min = (A[t-1]-0.5)+0.5*U[t-1]+0.3*L[t]-1,
                                    max = (A[t-1]-0.5)+0.5*U[t-1]+0.3*L[t]+1) +
  node("Z",t=1:TM, distr = "runif", min = 0.5*L[t]-0.5*(A[t-1]-0.5)-0.3*Z[t-1]-4,
                                    max = 0.5*L[t]-0.5*(A[t-1]-0.5)-0.3*Z[t-1]+4) +
  node("A",t=1:TM, distr = "rbern", prob = pnorm(-3*Z[t]+L[t])*0.8 + pnorm(3*U[t]-L[t])*0.2)+
  node("Y",t=TM+1,distr = "runif",  min = (A[t-1]-0.5)+2*L[t-1]+U[t-1] + 1,
                                    max = (A[t-1]-0.5)+2*L[t-1]+U[t-1] + 1)

D <- set.DAG(D, vecfun=c("pnorm"))
# D <- set.DAG(D)
# plotDAG(D, xjitter = 1)
























# TM <- 1 # maximal time
# library(simcausal)
# library(tibble)
# library(dplyr)
# library(purrr)
# library(ggplot2)
# 
# D <- DAG.empty() +
#   node("L",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
#   node("U",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
#   node("Z",t=0, distr = "rnorm", mean = 0.3*L[t]+sin(1.5*L[t]), sd = 2) +
#   node("A",t=0, distr = "rbern", prob = pnorm(-2*Z[t]+0.6*L[t])*0.7 + pnorm(3*U[t]-L[t])*0.3)+
# 
#   node("L",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*L[t-1]+0.3*U[t-1], sd = 0.5) +
#   node("U",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*U[t-1]+0.3*L[t], sd = 0.5) +
#   node("Z",t=1:TM, distr = "rnorm", mean = 0.5*L[t]-0.5*(A[t-1]-0.5)-0.3*Z[t-1],sd=2) +
#   node("A",t=1:TM, distr = "rbern", prob = pnorm(-2*Z[t]+L[t])*0.7 + pnorm(3*U[t]-L[t])*0.3)+
#   node("Y",t=TM+1,distr = "rnorm", mean = (A[t-1]-0.5)+2*L[t-1]+U[t-1], sd = 0.5)
# 
# D <- set.DAG(D, vecfun=c("pnorm"))
# # D <- set.DAG(D)
# # plotDAG(D, xjitter = 1)


