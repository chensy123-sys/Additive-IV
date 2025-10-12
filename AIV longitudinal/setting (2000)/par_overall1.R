TM <- 1 # maximal time
n<-5000
Looptime<-500
library(simcausal)
library(tibble)
library(dplyr)
library(SuperLearner)
library(purrr)
library(ggplot2)
library(randomForest)
# install.packages("remotes")
# remotes::install_github("tedwestling/survSuperLearner")
# library(survSuperLearner)

# deltaA = function(L) {
#   # 生成 Pr(A=1|Z=1,X,U)-Pr(A=1|Z=0,X,U)
#   tanh(0.7 + 0.4 * sin(0.25 * L))
# }
# OPA = function(L, U) {
#   # 生成 Pr(A=1|Z=1,X,U)/(1-Pr(A=1|Z=1,X,U)) * Pr(A=1|Z=0,X,U)/(1-Pr(A=1|Z=0,X,U))
#   exp(-sin(U*0.25))
# }
# pi = function(OPA, deltaA) {
#   numerator1 <- OPA * (2 - deltaA) + deltaA
#   numerator2 <- (OPA * (deltaA - 2) - deltaA)^2 + 4 * OPA * (1 - OPA) * (1 - deltaA)
#   denumerator <- 2 * (OPA - 1)
#   return( (numerator1 - sqrt(numerator2)) / denumerator )
# }


D <- DAG.empty() +
  node("L",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
  node("U",t=0, distr = "rnorm", mean = 0, sd = 1.5) +
  node("Z",t=0, distr = "rnorm", mean = 0.3*L[t]+sin(1.5*L[t]), sd = 2) +
  node("A",t=0, distr = "rbern", prob = pnorm(-2*Z[t]+0.6*L[t])*0.7 + pnorm(3*U[t]-L[t])*0.3)+
  
  node("L",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*L[t-1]+0.3*U[t-1], sd = 0.5) +
  node("U",t=1:TM, distr = "rnorm", mean = (A[t-1]-0.5)+0.5*U[t-1]+0.3*L[t], sd = 0.5) +
  node("Z",t=1:TM, distr = "rnorm", mean = 0.5*L[t]-0.5*(A[t-1]-0.5)-0.3*Z[t-1],sd=2) +
  node("A",t=1:TM, distr = "rbern", prob = pnorm(-2*Z[t]+L[t])*0.7 + pnorm(3*U[t]-L[t])*0.3)+
  node("Y",t=TM+1,distr = "rnorm", mean = (A[t-1]-0.5)+2*L[t-1]+U[t-1], sd = 0.5)

D <- set.DAG(D, vecfun=c("pnorm"))
# D <- set.DAG(D)
# plotDAG(D, xjitter = 1)