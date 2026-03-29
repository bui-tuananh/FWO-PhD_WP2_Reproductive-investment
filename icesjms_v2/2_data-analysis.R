# BEFORE
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

# K_rel (instead of Fulton K)
data_obs <- data_obs %>%
  mutate(somatic_weight = body_weight - gonad_weight)

data_obs_new <- tibble()
for(p in c("4bc", "7a", "7fg", "7d")) {
  data_pop <-  data_obs %>% filter(pop == p)
  lm = lm(log(somatic_weight) ~ 1 + log(length), data = data_pop)
  data_pop <- data_pop %>%
    mutate(pred_weight = exp(predict(lm, data_pop, type = "response")),
           K_rel = somatic_weight/pred_weight)
  
  data_obs_new <- bind_rows(data_obs_new, data_pop)
}

## temp autumn winter ----
dir_temp <- "./data/temp"
data_temp <- read_rds(file.path(dir_temp, "oras5_datras.rds")) 

# temp_autumn
data_temp_autumn <- data_temp %>%
  mutate(pop = if_else(IcesArea == "4abc", "4bc", IcesArea),
         month = month(Date),
         year = Year)  %>%
  filter(month >= 10 | month <= 12) %>%
  #filter(year >= 2002) %>%
  group_by(pop, year) %>%
  summarize(temp_aut = mean(oras_sbt)) %>%
  mutate(year = year + 1) #autumn 2018 - corresponding to data 2019

# temp_winter
data_temp_winter <- data_temp %>%
  mutate(pop = if_else(IcesArea == "4abc", "4bc", IcesArea),
         month = month(Date),
         year = Year) %>%
  filter(month >= 1 | month <= 3) %>%
  #filter(year >= 2002) %>%
  group_by(pop, year) %>%
  summarize(temp_win = mean(oras_sbt)) 

# aut_win
data_temp_aut_win <- data_temp_autumn %>%
  left_join(data_temp_winter) %>%
  drop_na() %>% 
  mutate(temp = (temp_aut + temp_win)/2 )

data_temp_aut_win <- data_temp_aut_win %>%
  group_by(pop) %>%
  mutate(ave.temp = mean(temp),
         c.temp = temp - ave.temp)

## ssb lag1 ----
dir_ices <- "./data/ices"

# sole stock assessment
data_sol <- read_rds(file.path(dir_ices, "stock-assessment_47adfg8ab_2023.rds")) %>%
  rename(ssb = SSB,
         f = `F`) 

data_ssb_lag1 <- data_sol %>%
  mutate(year_lag1 = year + 1,
         ssb_lag1 = ssb) %>%
  select(pop, year_lag1, ssb_lag1) %>%
  mutate(ssb_lag1 = ssb_lag1/1000) # kilotonnes instead of tonnes

## data all ----
data <- data_obs_new %>%
  left_join(data_temp_aut_win, by = join_by(pop, year == year)) %>%
  left_join(data_ssb_lag1, by = join_by(pop, year == year_lag1)) 

# DATA ANALYSIS ----
# Individual-level variation and temporal trend ----
## 4bc ----
# data
data_sub <- data %>%
  filter(pop == "4bc")

# model comparison
m1 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
             (1 | week) + (1 | year), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m1, subset = `log.body_weight`)

# best model
m1_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m1_reml, file.path(dir_output, "m1_4bc.rds"))

## 7a ----
# data
data_sub <- data %>%
  filter(pop == "7a")

# model comparison
m1 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
             (1 | week) + (1 | year) + (1 | cohort), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m1, subset = `log.body_weight`)

# best model
m1_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m1_reml, file.path(dir_output, "m1_7a.rds"))

## 7fg ----
# data
data_sub <- data %>%
  filter(pop == "7fg")

# model comparison
m1 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
             (1 | week) + (1 | year) + (1 | cohort), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m1, subset = `log.body_weight`)

# best model
m1_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m1_reml, file.path(dir_output, "m1_7fg.rds"))

## 7d ----
# data
data_sub <- data %>%
  filter(pop == "7d")

# model comparison
m1 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
             (1 | week) + (1 | year) + (1 | cohort), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m1, subset = `log.body_weight`)

# best model
m1_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m1_reml, file.path(dir_output, "m1_7d.rds"))

# Size-specific effects of warming ----
## 4bc ----
# data
data_sub <- data %>%
  filter(pop == "4bc")

# model comparison
m2 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
             c.temp*log.body_weight + c.temp*ssb_lag1 + 
             (1 | week) + (1 | year), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m2, subset = `log.body_weight` & `maturity_stage`)

# best model
m2_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
                  c.temp*log.body_weight +
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m2_reml, file.path(dir_output, "m2_4bc.rds"))

## 7a ----
# data
data_sub <- data %>%
  filter(pop == "7a")

# model comparison
m2 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
             c.temp*log.body_weight + c.temp*ssb_lag1 + 
             (1 | week) + (1 | year), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m2, subset = `log.body_weight` & `maturity_stage`)

# best model
m2_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
                  c.temp*log.body_weight + c.temp*ssb_lag1 +
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m2_reml, file.path(dir_output, "m2_7a.rds"))

## 7fg ----
# data
data_sub <- data %>%
  filter(pop == "7fg")

# model comparison
m2 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
             c.temp*log.body_weight + c.temp*ssb_lag1 + 
             (1 | week) + (1 | year), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m2, subset = `log.body_weight` & `maturity_stage`)

# best model
m2_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + 
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m2_reml, file.path(dir_output, "m2_7fg.rds"))

## 7d ----
# data
data_sub <- data %>%
  filter(pop == "7d")

# model comparison
m2 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
             c.temp*log.body_weight + c.temp*ssb_lag1 + 
             (1 | week) + (1 | year), 
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m2, subset = `log.body_weight` & `maturity_stage`)

# best model
m2_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
                  c.temp*log.body_weight +
                  (1 | week) + (1 | year), 
                data = data_sub,
                REML = T)
write_rds(m2_reml, file.path(dir_output, "m2_7d.rds"))

# Population-specific effects of warming ---- 
# data
data_sub <- data 

# model comparison
m3 <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
             pop*log.body_weight +
             c.temp*log.body_weight + c.temp*pop + 
             (1 | pop.week) + (1 | pop.year), #note: pop.week and pop.year, not week and year
           data = data_sub,
           REML = F,
           na.action = "na.fail")

model_compare <- dredge(m3, subset = `log.body_weight` & `maturity_stage`)

# best model
m3_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
                  pop*log.body_weight +
                  c.temp*log.body_weight +  
                  (1 | pop.week) + (1 | pop.year),
                data = data_sub,
                REML = T)
write_rds(m3_reml, file.path(dir_output, "m3.rds"))
## model with temp ----
# model with temp
data_sub <- data
m3_temp <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
                  ave.temp*log.body_weight +
                  c.temp*log.body_weight + c.temp*ave.temp + 
                  (1 | pop.week) + (1 | pop.year), #note: pop.week and pop.year, not week and year
                data = data_sub,
                REML = F,
                na.action = "na.fail")

model_compare <- dredge(m3, subset = `log.body_weight` & `maturity_stage` & `K_rel`)

# best model
m3_temp_reml <- lmer(log.gonad_weight ~ 1 + log.body_weight + maturity_stage + K_rel +
                       ave.temp +
                       c.temp*log.body_weight +  
                       (1 | pop.week) + (1 | pop.year),
                     data = data_sub,
                     REML = T)
write_rds(m3_temp_reml, file.path(dir_output, "m3_temp.rds"))
