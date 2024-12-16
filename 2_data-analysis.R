# BEFORE
# 0_data-exploration
# 1_data-processing

# AFTER
# output

# SETUP ----
library(tidyverse)  # data processing
library(glmmTMB)    # model 
library(lme4)       # model
library(lmerTest)   # model
library(AICcmodavg) # model comparison
library(MuMIn)      # model comparison
library(effects)    # effect visualisation
library(sf)         # spatial data processing
library(patchwork)  # visualisation

# set theme
theme_set(theme_bw())

# set directory
dir_output <- "./output"

# DATA ----
## gonad ----
dir_data <- "./data"
data_obs <- read_rds(file.path(dir_data, "sol_gonads_2004_2023_sub_ready-for-analysis.rds"))

## temperature ----
dir_temp <- "./data/temp"
data_temp <- read_rds(file.path(dir_temp, "oras5_datras.rds")) 

# winter temperature 2004-2022
data_temp <- data_temp %>%
  mutate(pop = if_else(IcesArea == "4abc", "4bc", IcesArea),
         month = month(Date),
         year = Year) %>%
  filter(month >= 1 | month <= 3) %>%
  filter(year >= 2004, year <= 2022) %>%
  group_by(pop, year) %>%
  summarize(temp = mean(oras_sbt))

# calculate temp anomaly for each population
data_temp <- data_temp %>%
  group_by(pop) %>%
  mutate(ave.temp = mean(temp),
         c.temp = temp - ave.temp)

## fishing ----
dir_ices <- "./data/ices"

# sole distribution area - survey datras
datras <- read_sf(file.path(dir_ices, "hl_loc_4abc7adefg8ab.gpkg"))
datras <- as_tibble(datras) %>%
  mutate(pop = if_else(Area_27 == "4abc", "4bc", Area_27)) %>%
  select(pop, area_km2)

# sole stock assessment
data_sol <- read_rds(file.path(dir_ices, "stock-assessment_47adfg8ab_2023.rds")) %>%
  rename(ssb = SSB,
         f = `F`) %>%
  select(pop, year, ssb, f, recruitment) %>%
  left_join(datras) %>%
  mutate(ssb.i = ssb/area_km2,
         recruitment.i = recruitment/area_km2) %>%
  select(-f) #remove fbar from stock assessment (different fbar age range across populations)

# add fbar from f-at-age 3-7 (consistent fbar age range across all populations)
data_sol_fbar <- read_rds(file.path(dir_ices, "sol_47a7d7fg_fbar_age3-7_stock-assessment_2023.rds"))
data_sol <- data_sol %>% left_join(data_sol_fbar)

# DATA ANALYSIS ----
# model weight -----  
data <- data_obs %>%
  left_join(data_sol) %>%
  left_join(data_temp)

data <- data %>%
  mutate(pop.week = paste0(pop, ":", week),
         pop.year = paste0(pop, ":", year))

## intrinsic -----
m2 <- lmer(log.gonad_weight ~ 1 + log.body_weight + condition + maturity_stage + 
                pop*log.body_weight + 
                (1 | pop.week) + (1 | pop.year), 
              data = data,
              REML = F,
              na.action = "na.fail")

model_compare <- dredge(m2, subset = `log.body_weight`)
write_rds(model_compare, file.path(dir_output, "model.weight_intrinsic_comparison.rds"))

m2_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + condition + maturity_stage + 
             pop*log.body_weight + 
             (1 | pop.week) + (1 | pop.year), 
           data = data,
           REML = T)

write_rds(m2_reml, file.path(dir_output, "model.weight_intrinsic.rds"))

## extrinsic ----
m3 <- lmer(log.gonad_weight ~ 1 + log.body_weight + condition + maturity_stage + 
                  pop*log.body_weight + c.temp*log.body_weight + c.temp*pop + f + ssb.i + 
                    (1 | pop.week) + (1 | pop.year), 
                  data = data,
                  REML = F,
                  na.action = "na.fail")

model_compare <- dredge(m3, subset = `log.body_weight`)
write_rds(model_compare, file.path(dir_output, "model.weight_extrinsic_comparison.rds"))

m3_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + condition + maturity_stage + 
             pop*log.body_weight + c.temp*log.body_weight + 
             (1 | pop.week) + (1 | pop.year), 
           data = data,
           REML = T)

write_rds(m3_reml, file.path(dir_output, "model.weight_extrinsic.rds"))

# model length -----
## intrinsic -----
m2 <- lmer(log.gonad_weight ~ 1 + log.length + condition + maturity_stage + 
             pop*log.length + 
             (1 | pop.week) + (1 | pop.year), 
           data = data,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m2, subset = `log.length`)
write_rds(model_compare, file.path(dir_output, "model.length_intrinsic_comparison.rds"))

m2_reml <- lmer(log.gonad_weight ~ 1 + log.length + condition + maturity_stage + 
                  pop*log.length + 
                  (1 | pop.week) + (1 | pop.year), 
                data = data,
                REML = T)

write_rds(m2_reml, file.path(dir_output, "model.length_intrinsic.rds"))

## extrinsic ----
m3 <- lmer(log.gonad_weight ~ 1 + log.length + condition + maturity_stage + 
             pop*log.length + c.temp*log.length + c.temp*pop + f + ssb.i + 
             (1 | pop.week) + (1 | pop.year), 
           data = data,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m3, subset = `log.length`)
write_rds(model_compare, file.path(dir_output, "model.length_extrinsic_comparison.rds"))

m3_reml <- lmer(log.gonad_weight ~ 1 + log.length + condition + maturity_stage + 
                  pop*log.length + c.temp*log.length + 
                  (1 | pop.week) + (1 | pop.year), 
                data = data,
                REML = T)

write_rds(m3_reml, file.path(dir_output, "model.length_extrinsic.rds"))




