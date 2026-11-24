library(dplyr);library(tibble);library(tidyr);library(readr)
my_data <- read.table("results/jtpa_han.tab", header = TRUE, sep = "\t") 
head(my_data)
attach(my_data)
my_data <- as_tibble(my_data %>% 
  select(-recid,-n_hs2_abv) %>% 
  mutate(
    edu = as.numeric(edu >=12), 
    prevearn = as.vector(scale(log(prevearn+1))),
    n_hs2 = as.vector(scale((1+pmin(n_hs2,40))**0.5)),
    earnings = as.numeric(earnings>median(earnings))
  ))
print(head(my_data))
colnames(my_data) <- c('A1', 'Z1', 'Y2', 'L1', 'A0', 'Z0', 'L0')
write.csv(my_data,'results/jpta_han.csv')