# input: get_sol_gonads
# output: sol_gonads_2004_2023_sub.rds

# Note: 
# pre process data extracted from SmartFish to remove sensitive information and 
# make the after-processed data shareable

# SETUP
library(tidyverse)

sol <- read_rds("./data/sol_gonads_2004_2023.rds") 

sol_sub <- sol %>% 
  select(TRI.Year, 
         TRI.DepartureDate,
         TRI.ReturnDate,
         HAU.IcesArea, 
         SAM.SpeciesFaoCode,
         SpecimenID,
         SPE.Length,
         SPE.Weight,
         SPE.WeightGonads,
         SPE.Sex,
         SPE.MaturityStageDescription,
         SPA.Age)

write_rds(sol_sub, "./data/sol_gonads_2004_2023_sub.rds")


