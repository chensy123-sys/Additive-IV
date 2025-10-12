rm(list=ls())
source('par_overall1.R')
source('LIV_overall.R')

alpha=0.000
trimmed_sd <- function(x, trim = 0.01) {
  q <- quantile(x, c(trim, 1 - trim))   # 得到上下分位数
  x_trim <- x[x >= q[1] & x <= q[2]]    # 去掉极端值
  sd(x_trim)
}


seed <- as.integer(2025)
set.seed(2025)

pdf('plot_overall.pdf', width = 6, height = 4)
par(mfrow=c(2,2))
load("Overall_00.RData")
Truth <- suppressMessages(ATE_overall(D,TM,c(0,0),accuracy = 5e5))
print(round(Truth,3))
print(round(apply(RES[,2,], 1, mean)-Truth,4))
print(round(apply(RES[,3,], 1, mean),4))
print(round(apply(RES[,2,], 1, trimmed_sd),4))
for(t in 1:(TM+2)){
  print(mean((RES[t,4,] <= Truth[t]) & (RES[t,5,] >= Truth[t]))  )
}
plot(0:(TM+1), apply(RES[,2,], 1, mean),'o',lwd=2)
lines(0:(TM+1),Truth,'o',lwd=2,col='red')



load("Overall_01.RData")
Truth <- ATE_overall(D,TM,c(0,1),accuracy = 5e5)
print(round(Truth,3))
print(round(apply(RES[,2,], 1, mean)-Truth,4))
print(round(apply(RES[,3,], 1, mean),4))
print(round(apply(RES[,2,], 1, trimmed_sd),4))
for(t in 1:(TM+2)){
  print(mean((RES[t,4,] <= Truth[t]) & (RES[t,5,] >= Truth[t]))  )
}
plot(0:(TM+1), apply(RES[,2,], 1, mean),'o',lwd=2)
lines(0:(TM+1),Truth,'o',lwd=2,col='red')


load("Overall_10.RData")
Truth <- ATE_overall(D,TM,c(1,0),accuracy = 5e5)
print(round(Truth,3))
print(round(apply(RES[,2,], 1, mean)-Truth,4))
print(round(apply(RES[,3,], 1, mean),4))
print(round(apply(RES[,2,], 1, trimmed_sd),4))
for(t in 1:(TM+2)){
  print(mean((RES[t,4,] <= Truth[t]) & (RES[t,5,] >= Truth[t]))  )
}
plot(0:(TM+1), apply(RES[,2,], 1, mean),'o',lwd=2)
lines(0:(TM+1),Truth,'o',lwd=2,col='red')


load("Overall_11.RData")
Truth <- ATE_overall(D,TM,c(1,1),accuracy = 5e5)
print(round(Truth,3))
print(round(apply(RES[,2,], 1, mean)-Truth,4))
print(round(apply(RES[,3,], 1, mean),4))
print(round(apply(RES[,2,], 1, trimmed_sd),4))
for(t in 1:(TM+2)){
  print(mean((RES[t,4,] <= Truth[t]) & (RES[t,5,] >= Truth[t]))  )
}
plot(0:(TM+1), apply(RES[,2,], 1, mean),'o',lwd=2)
lines(0:(TM+1),Truth,'o',lwd=2,col='red')
par(mfrow=c(1,1))
dev.off()
