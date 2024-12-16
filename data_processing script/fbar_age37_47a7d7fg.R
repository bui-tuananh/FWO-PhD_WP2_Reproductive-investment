library(tidyverse)
library(readxl)
dir_ices <- "./data/ices"

# f-at-age range vs stock assessment Fbar: 
## 4: 1-10 (Fbar 2-6)
## 7a: 2-8 (Fbar 4-7)
## 7d: 1-11 (Fbar 3-7)
## 7fg: 1-10 (Fbar 3-8)
# -> new Fbar 3-7

f4 <- read_xlsx(file.path(dir_ices, "sol_4_f-at-age.xlsx")) 
f7a <- read_xlsx(file.path(dir_ices, "sol_7a_f-at-age.xlsx")) 
f7d <- read_xlsx(file.path(dir_ices, "sol_7d_f-at-age.xlsx")) 
f7fg <- read_xlsx(file.path(dir_ices, "sol_7fg_f-at-age.xlsx")) 

f4 <- f4 %>% 
  pivot_longer(cols = age1:age10, names_to = "age", values_to = "f") %>%
  mutate(age = as.integer(substr(age, 4, 5))) %>%
  mutate(pop = "4bc")

f7a <- f7a %>% 
  pivot_longer(cols = age2:age8, names_to = "age", values_to = "f") %>%
  mutate(age = as.integer(substr(age, 4, 4))) %>%
  mutate(pop = "7a")

f7d <- f7d %>% 
  pivot_longer(cols = age1:age11, names_to = "age", values_to = "f") %>%
  mutate(age = as.integer(substr(age, 4, 5))) %>%
  mutate(pop = "7d")

f7fg <- f7fg %>% 
  pivot_longer(cols = age1:age10, names_to = "age", values_to = "f") %>%
  mutate(age = as.integer(substr(age, 4, 5))) %>%
  mutate(pop = "7fg")

# combine all f together
data_sol_f <- bind_rows(f4, f7a, f7d, f7fg)

# calculate Fbar for age 3-7
data_sol_fbar <- data_sol_f %>%
  filter(age %in% seq(3,7)) %>%
  group_by(pop, year) %>%
  summarize(f = mean(f))

write_rds(data_sol_fbar, file.path(dir_ices, "sol_47a7d7fg_fbar_age3-7_stock-assessment_2023.rds"))
