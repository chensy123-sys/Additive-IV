rm(list=ls())
# setwd('../')
source('LIV_overall.R')
set.seed(1234)
my_data <- read.csv("jpta_han.csv")[,2:8]
# tibble(my_data)
print(head(my_data))
Looptime <- 500
RES <- array(0,dim = c(3,5,Looptime))
for(loop in 1:Looptime){
  print(loop)
  my_data_boot <- my_data[sample(nrow(my_data),replace = TRUE),]
  my_data_boot_long <- my_data_boot %>%
    mutate(id = row_number()) %>% 
    pivot_longer(
      cols = c(A0, Z0, L0, A1, Z1, L1, Y2),
      names_to = c(".value", "time"),
      names_pattern = "([A-Z]+)([0-9]+|Y)?"
    ) %>%
    arrange(id, time) %>% rename(ID=id, t=time) %>%
    type_convert()
  # print(my_data_boot_long, n=100)
  # print(head(my_data_boot_long))
  # print(head(my_data_boot))
  fit_LIV <- DML_LIV_overall(my_data_boot_long,action_list = c(0,0),K=1)
  # fit_LIV2 <- DML_LIV_overall(my_data_boot_long,action_list = c(1,0),K=1)
  RES[,,loop] <- as.matrix(predict_LIV2_overall(fit_LIV))
}

save(RES, file = "overall_00.RData")

# predict_LIV2_overall(fit_LIV)

# fit_LIV <- DML_LIV_overall(my_data_boot_long,action_list = list(
#     function(fY) {return(0)},
#     function(fY) {
#       # print(fY)
#       # print(quantile(fY$L,0.8))
#       # print(table(fY$L>=0.8368623))
#       return(as.numeric(fY$L<=0.8368623)) # 0.8 quantile 
#     }
#   ),K=1)









