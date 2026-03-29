
# setup ----
#install.packages("icesSAG")
library(icesSAG)
library(tidyverse)
library(readxl)

# get stock assessment data ----
meta <- getListStocks(2023)
sol <- meta %>% filter(Purpose == "Advice",
                StockKeyLabel %in% c("sol.27.4", "sol.27.7a", "sol.27.8ab"))
# assessmentKey 
# North Sea               sol.27.4    18458
# Irish Sea               sol.27.7a   18059
# Bay of Biscay           sol.27.8ab  18024
# Eastern English Channel sol.27.7d   18120
# Celtic Sea              sol.27.7fg  18212


stock_4 <- getSummaryTable(assessmentKey = 18458)[[1]] %>%
  mutate(pop = "4bc") %>%
  rename(year = Year)

stock_7a <- getSummaryTable(assessmentKey = 18059)[[1]] %>%
  mutate(pop = "7a") %>%
  rename(year = Year)

stock_8ab <- getSummaryTable(assessmentKey = 18024)[[1]] %>%
  mutate(pop = "8ab") %>%
  rename(year = Year)

stock_7d <- getSummaryTable(assessmentKey = 18120)[[1]] %>%
  mutate(pop = "7d") %>%
  rename(year = Year)

stock_7fg <- getSummaryTable(assessmentKey = 18120)[[1]] %>%
  mutate(pop = "7fg") %>%
  rename(year = Year)

# merge all data
stock_all <- rbind(stock_4, stock_7a, stock_8ab, stock_7d, stock_7fg)

# save data
write_rds(stock_all, file.path("./data/ices", "stock-assessment_47adfg8ab_2023.rds"))
