# setwd('./data_K1_overall')
load('overall_00.RData')
print('treatment regimes: 00')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)

load('overall_01.RData')
print('treatment regimes: 01')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)

load('overall_02.RData')
print('treatment regimes: 02')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)

load('overall_03.RData')
print('treatment regimes: 03')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)




load('overall_10.RData')
print('treatment regimes: 10')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)

load('overall_11.RData')
print('treatment regimes: 11')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)

load('overall_12.RData')
print('treatment regimes: 12')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)

load('overall_13.RData')
print('treatment regimes: 13')
round(apply(RES[,,],c(1,2),mean),3)
round(apply(RES[,,],c(1,2),sd),3)

# apply(RES[,1,],2,mean)
# apply(RES[,2,],2,mean)
# apply(RES[,3,],2,mean)



