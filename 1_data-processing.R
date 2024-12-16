# BEFORE
# 0_data-exploration
# 1_data-processing

# SETUP ----
library(tidyverse)  # data processing
library(glmmTMB)    # model 

# DATA ----
## pre-processing obs ----
obs <- read_rds("./data/sol_gonads_2004_2023_sub.rds")

# select areas, remove NA, abnormal, skipped spawning stages, NA age, select female
obs_sub <- obs %>% 
  filter(HAU.IcesArea %in% c("4b", "4c", "7a", "7d", "7f", "7g"),
         is.na(SPE.MaturityStageDescription) == F,
         !SPE.MaturityStageDescription %in% c("Abnormal", "Skipped spawning"), #2 obs Abnormal, 4 obs Skipped spawning
         is.na(SPA.Age) == F,
         SPE.Sex == "F") %>%
  mutate(SPE.MaturityStageDescription = if_else(SPE.MaturityStageDescription %in% c("Maturing", "Developing"),
                                                "Developing",
                                                SPE.MaturityStageDescription))

# sampling week: majority of trips within 1 week, yet there are trips in 1-2 weeks
# -> get the middle value of week, month
# rename field and do some merging, factoring
obs_sub <- obs_sub %>%
  mutate(id = SpecimenID, 
         pop = if_else(HAU.IcesArea %in% c("4b", "4c"), "4bc",
                       if_else(HAU.IcesArea %in% c("7f", "7g"), "7fg", HAU.IcesArea)),
         age = SPA.Age,
         year = TRI.Year,
         cohort = year - age,
         month = round((month(TRI.DepartureDate) + month(TRI.ReturnDate))/2),
         week = round((week(TRI.DepartureDate) + week(TRI.ReturnDate))/2),
         week_diff = week(TRI.ReturnDate) - week(TRI.DepartureDate),
         length = SPE.Length,
         body_weight = SPE.Weight,
         gonad_weight = SPE.WeightGonads,
         gsi = gonad_weight/body_weight*100,
         maturity_stage = factor(SPE.MaturityStageDescription, level = c("Immature", "Developing", "Spawning", "Spent")),
         maturity = if_else(SPE.MaturityStageDescription %in% c("Immature"), 0, 1),
         maturity_desc = if_else(SPE.MaturityStageDescription %in% c("Immature"), "Immature", "Mature"),
  ) %>%
  select(id, pop, year, cohort, month, week, age, length, body_weight, gonad_weight, gsi, maturity_stage, maturity, maturity_desc)

## filter weeks within gsi >= 50% peak gsi ----
# remove gsi >= 100 (gutted weight include gonad weight -> cannot have gsi >= 100)
obs_sub <- obs_sub %>% 
  filter(gsi < 100, is.na(gsi) == F) %>%
  filter(maturity_stage %in% c("Developing", "Spawning")) %>%
  mutate(week2 = if_else(week >= 40, week - 53, week),
         year2 = if_else(week2 <0, year+1, year))

pred_gsi <- read_rds("./output/pred_gsi.rds")
pred_gsi_sub <- pred_gsi %>%
  group_by(pop) %>%
  mutate(gsi50 = max(gsi_pred)*0.5) %>%
  filter(gsi_pred >= gsi50) %>%
  summarize(week2_start = round(min(week2)),
            week2_end   = round(max(week2)))

pred_gsi <- pred_gsi %>% 
  left_join(pred_gsi_sub)

# week2_start 0 is the same as week2_start 1
obs_sub <- obs_sub %>%
  left_join(pred_gsi_sub) %>%
  filter(week2 >= week2_start,
         week2 <= week2_end)

## remove very large/small values ----
quantile(obs_sub$gsi)
# upper and lower whisker
11.99+1.5*(11.99-5.85) # 21.2
5.85-1.5*(11.99-5.85) # -3.36

ggplot(data = obs_sub, aes(x = gsi, y = id)) +
  geom_point(alpha = 0.3) +
  geom_vline(xintercept = 30) +
  geom_vline(xintercept = 21.2, linetype = "dashed") +
  labs(x = "Gonadosomatic index - GSI",
       y = "Observation ID") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# check distribution gsi
ggplot(data = obs_sub %>% filter(gsi <= 30), 
       aes(x = gsi, fill = maturity_stage)) +
  geom_histogram(bins = 200, position = "identity", alpha = 0.7) +
  facet_grid( ~ pop) +
  geom_vline(xintercept = 2)

# remove GSI > 30 and GSI < 2
nrow(obs_sub %>% filter(gsi > 30)) #27 observations

obs_sub <- obs_sub %>% filter(gsi >= 2, gsi <= 30)

### data for analysis ----
data_obs <- obs_sub %>%
  mutate(condition = 100*(body_weight)/(length*0.1)^3) %>%
  mutate(log.gonad_weight = log(gonad_weight),
         log.body_weight = log(body_weight),
         log.length = log(length),
         log.age = log(age),
         log.condition = log(condition),
         fweek = factor(week),
         fyear = factor(year),
         fcohort = factor(cohort)) 

# sample size vs year
# full data
View(data_obs %>%
       group_by(pop, year) %>%
       summarize(n = n()) %>%
       pivot_wider(names_from = pop, values_from = n) %>%
       arrange(year))

# sample size vs age
View(data_obs %>%
       group_by(pop, age) %>%
       summarize(n = n()) %>%
       arrange(age) %>%
       pivot_wider(names_from = age, values_from = n))
# only select age with at least 10 obs
# 4bc: 2-12
# 7a: 3-14
# 7d: 3-19
# 7fg: 3-16
data_obs <- data_obs %>%
  filter(age >= if_else(pop == "4bc", 2, 3), 
         age <= if_else(pop == "4bc", 12, 
                        if_else(pop == "7a", 14,
                                if_else(pop == "7d", 19, 16))))

# add pop:week, pop:year
data_obs <- data_obs %>%
  mutate(pop.week = paste0(pop, ":", week),
         pop.year = paste0(pop, ":", year))

write_rds(data_obs, './data/sol_gonads_2004_2023_sub_ready-for-analysis.rds')
