rm(list=ls())
source('LIV_overall.R')
set.seed(2025)
my_data <- read.csv("jpta_han.csv")[,2:8]

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

  fit_LIV <- DML_LIV_overall(my_data_boot_long,action_list = c(1,0),K=1)
  RES[,,loop] <- as.matrix(predict_LIV2_overall(fit_LIV))
}

save(RES, file = "overall_10.RData")

# predict_LIV2_overall(fit_LIV)
