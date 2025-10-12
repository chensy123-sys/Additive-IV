library(dplyr);library(tibble);library(tidyr);library(readr)
# sep = "\t" 明确指定使用制表符作为分隔符
# setwd('./jtpa_han')
# getwd()
my_data <- read.table("jtpa_han.tab", header = TRUE, sep = "\t") 

head(my_data)
attach(my_data)
# cor(earnings,n_hs2)

# earnings - prevearn




my_data <- as_tibble(my_data %>% 
  select(-recid,-n_hs2_abv) %>% 
  mutate(
    edu = as.numeric(edu >=12), 
    prevearn = as.vector(scale(log(prevearn+1))),
    # prevearn = as.vector(scale(log(prevearn+1))),
    n_hs2 = pmin(n_hs2, 40),
    # n_hs2 = as.vector(scale(log(1+n_hs2))), # 0    1    2    5    9   13   15   16   17   19   28   37   39   40  346 
    n_hs2 = as.vector(scale((1+n_hs2)**0.5)),
    # n_hs2 = as.numeric(n_hs2>=13),
    # earnings = as.numeric(earnings - prevearn>median(earnings - prevearn))
    earnings = as.numeric(earnings>median(earnings))
  ))
print(head(my_data))
colnames(my_data) <- c('A1', 'Z1', 'Y2', 'L1', 'A0', 'Z0', 'L0')

# cor(my_data$A0, my_data$L1)
# cor(my_data$Z0, my_data$L1)

# cor(my_data$A0, my_data$Y2)
# cor(my_data$Z0, my_data$Y2)
# cor(my_data$Z0, my_data$A0)

write.csv(my_data,'jpta_han.csv')