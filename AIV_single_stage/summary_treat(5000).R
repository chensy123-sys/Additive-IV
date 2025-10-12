alpha=0.00
load("treat1.RData")
index <- (REC_unknown[,1] <= quantile(REC_unknown[,1], 1-alpha)) & 
         (REC_unknown[,1] >= quantile(REC_unknown[,1], alpha))
print("\n case 1")
print(round(apply(REC_unknown[index,], 2, mean),4))
print(round(apply(REC_unknown[index,], 2, sd),4))
print(round(mean(abs(REC_unknown[,1]/REC_unknown[,2])<qnorm(0.975)),4))
index <- (REC_known[,1] <= quantile(REC_known[,1], 1-alpha)) & 
         (REC_known[,1] >= quantile(REC_known[,1], alpha))
print(round(apply(REC_known[index,], 2, mean),4))
print(round(apply(REC_known[index,], 2, sd),4))
print(round(mean(abs(REC_known[,1]/REC_known[,2])<qnorm(0.975)),4))





load("treat2.RData")
index <- (REC_unknown[,1] <= quantile(REC_unknown[,1], 1-alpha)) & 
         (REC_unknown[,1] >= quantile(REC_unknown[,1], alpha))
print("\n case 2")
print(round(apply(REC_unknown[index,], 2, mean),4))
print(round(apply(REC_unknown[index,], 2, sd),4))
print(round(mean(abs(REC_unknown[,1]/REC_unknown[,2])<qnorm(0.975)),4))
index <- (REC_known[,1] <= quantile(REC_known[,1], 1-alpha)) & 
         (REC_known[,1] >= quantile(REC_known[,1], alpha))
print(round(apply(REC_known[index,], 2, mean),4))
print(round(apply(REC_known[index,], 2, sd),4))
print(round(mean(abs(REC_known[,1]/REC_known[,2])<qnorm(0.975)),4))





load("treat3.RData")
index <- (REC_unknown[,1] <= quantile(REC_unknown[,1], 1-alpha)) & 
         (REC_unknown[,1] >= quantile(REC_unknown[,1], alpha))
print("\n case 3")
print(round(apply(REC_unknown[index,], 2, mean),4))
print(round(apply(REC_unknown[index,], 2, sd),4))
print(round(mean(abs(REC_unknown[,1]/REC_unknown[,2])<qnorm(0.975)),4))
index <- (REC_known[,1] <= quantile(REC_known[,1], 1-alpha)) & 
         (REC_known[,1] >= quantile(REC_known[,1], alpha))
print(round(apply(REC_known[index,], 2, mean),4))
print(round(apply(REC_known[index,], 2, sd),4))
print(round(mean(abs(REC_known[,1]/REC_known[,2])<qnorm(0.975)),4))






load("treat4.RData")
index <- (REC_unknown[,1] <= quantile(REC_unknown[,1], 1-alpha)) & 
         (REC_unknown[,1] >= quantile(REC_unknown[,1], alpha))
print("\n case 4")
print(round(apply(REC_unknown[index,], 2, mean),4))
print(round(apply(REC_unknown[index,], 2, sd),4))
print(round(mean(abs(REC_unknown[,1]/REC_unknown[,2])<qnorm(0.975)),4))
index <- (REC_known[,1] <= quantile(REC_known[,1], 1-alpha)) & 
         (REC_known[,1] >= quantile(REC_known[,1], alpha))
print(round(apply(REC_known[index,], 2, mean),4))
print(round(apply(REC_known[index,], 2, sd),4))
print(round(mean(abs(REC_known[,1]/REC_known[,2])<qnorm(0.975)),4))