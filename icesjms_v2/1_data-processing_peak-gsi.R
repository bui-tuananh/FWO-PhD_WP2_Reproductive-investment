# setup ----
library(tidyverse)
library(glmmTMB)
library(patchwork)

theme_set(theme_bw())

# data ----
obs <- read_rds("./data/sol_gonads_2004_2023.rds")

obs_sub <- obs %>% 
  filter(HAU.IcesArea %in% c("4b", "4c", "7a", "7d", "7f", "7g"),
         is.na(SPE.MaturityStageDescription) == F,
         !SPE.MaturityStageDescription %in% c("Abnormal", "Skipped spawning"), #2 obs Abnormal, 4 obs Skipped spawning
         is.na(SPA.Age) == F,
         SPE.Sex == "F") %>%
  mutate(SPE.MaturityStageDescription = if_else(SPE.MaturityStageDescription %in% c("Maturing", "Developing"),
                                                "Maturing",
                                                SPE.MaturityStageDescription))

# sampling week: majority of trips within 1 week, yet there are trips in 1-2 weeks
# -> get the middle value of week; for trip 1 week, get the later week 
obs_sub <- obs_sub %>%
  mutate(id = SpecimenID, 
         pop = if_else(HAU.IcesArea %in% c("4b", "4c"), "4bc", 
                       if_else(HAU.IcesArea %in% c("7f", "7g"), "7fg", HAU.IcesArea)),
         age = SPA.Age,
         year = TRI.Year,
         cohort = year - age,
         week = round((week(TRI.DepartureDate) + week(TRI.ReturnDate))/2),
         week_diff = week(TRI.ReturnDate) - week(TRI.DepartureDate),
         length = SPE.Length,
         body_weight = SPE.Weight,
         gonad_weight = SPE.WeightGonads,
         gsi = gonad_weight/body_weight*100,
         maturity_stage = factor(SPE.MaturityStageDescription, level = c("Immature", "Maturing", "Spawning", "Spent")),
         maturity = if_else(SPE.MaturityStageDescription %in% c("Immature"), 0, 1),
         maturity_desc = if_else(SPE.MaturityStageDescription %in% c("Immature"), "Immature", "Mature"),
         ) %>%
  select(pop, year, cohort, week, age, length, body_weight, gonad_weight, gsi, maturity_stage, maturity, maturity_desc)

# data_strange
obs_strange <- obs_sub %>% 
  filter(gsi >= 100)

# data
obs_sub <- obs_sub %>% 
  filter(gsi < 30,          #only 40/25410 obs with gsi > 25 (among those with gsi < 100)
         is.na(gsi) == F) %>%  #3 obs with NA gsi
  mutate(log.gonad_weight = log(gonad_weight),
       log.body_weight = log(body_weight),
       log.length = log(length),
       week2 = if_else(week >= 40, week - 53, week),
       year2 = if_else(week2 <0, year+1, year))

# temporal 
ggplot(data = obs_sub, aes(x = week2, y = gsi, color = maturity_stage)) +
  geom_jitter() +
  facet_grid(pop ~ .)

# distribution gsi
# maturing but gsi < 2% - maturing: functionally mature
ggplot(data = obs_sub, aes(x = gsi, fill = maturity_stage)) +
  geom_histogram(bins = 200, position = "identity", alpha = 0.7) +
  facet_grid(pop ~ .) 

ggplot(data = obs_sub %>% filter(maturity_stage %in% c("Maturing", "Spawning")), 
       aes(x = gsi, fill = maturity_stage)) +
  geom_histogram(bins = 200, position = "identity", alpha = 0.7) +
  facet_grid(pop ~ .) +
  geom_vline(xintercept = 2)

# peak gsi, spawning ratio ----

## peak gsi ----
list_pop <- sort(unique(obs_sub$pop))

pred_gsi <- tibble()
for(pop_name in list_pop) {
  print(paste0("processing ", pop_name))
  
  # data
  obs_temp <- obs_sub %>% 
    filter(pop == pop_name) %>%
    mutate(fage = factor(age))
  
  # model
  m_gsi <- glmmTMB(gsi ~ 1 + s(week2) + (1 | year2), data = obs_temp)
  # wanted to use bs="cc" to be more elegant but is not supported yet in glmmTMB
  # bs="cc" specifies a cyclic cubic regression splines (see cyclic.cubic.spline). i.e. a penalized cubic regression splines whose ends match, up to second derivative.
  # the smooth is a cycle - the end of the smooth is the same as the starts 
  
  # pred
  range_week2 = range(obs_temp$week2)
  pred <- tibble(pop = pop_name,
                 #fage = 4, #before fage = 3 
                 week2 = seq(min(range_week2), max(range_week2), 0.1),
                 year2 = NA)
  
  pred <- pred %>%
    mutate(gsi_pred = predict(m_gsi, newdata = pred, type = "response", se = T)$fit,
           gsi_se = predict(m_gsi, newdata = pred, type = "response", se = T)$se.fit) 
  
  pred_gsi <- bind_rows(pred_gsi, pred)
}

# summary
pred_gsi_sub <- pred_gsi %>%
  group_by(pop) %>%
  mutate(gsi50 = max(gsi_pred)*0.5) %>%
  filter(gsi_pred >= gsi50) %>%
  summarize(week2_start = round(min(week2)),
            week2_end   = round(max(week2)))

# plot
ggplot() +
  geom_jitter(data = obs_sub, aes(x = week2, y = gsi, color = maturity_stage), alpha = 0.5) +
  geom_line(data = pred_gsi, aes(x = week2, y = gsi_pred)) +
  geom_ribbon(data = pred_gsi, aes(x = week2, 
                                       ymin = gsi_pred - 1.96*gsi_se, 
                                       ymax = gsi_pred + 1.96*gsi_se),
              alpha = 0.5) +
  facet_grid(pop ~ .)

ggplot() +
  geom_line(data = pred_gsi, aes(x = week2, y = gsi_pred, color = pop)) +
  geom_ribbon(data = pred_gsi, aes(x = week2, 
                                   ymin = gsi_pred - 1.96*gsi_se, 
                                   ymax = gsi_pred + 1.96*gsi_se,
                                   fill = pop),
              alpha = 0.5)


# save file
write_rds(pred_gsi, "./output/pred_gsi.rds")

## peak spawning ratio -----
# group by week then do beta (very limited data to group by year also);
# fincham 2013 - peak spawning is about week 20

list_pop <- sort(unique(obs_sub$pop))

obs_mature <- obs_sub %>%
  group_by(pop, week, week2, maturity_desc) %>%
  summarize(n_mature = n()) %>%
  filter(maturity_desc == "Mature") %>%
  select(pop, week, week2, n_mature)

obs_spawn <- obs_sub %>%
  group_by(pop, week, week2, maturity_stage) %>%
  summarize(n_spawn = n()) %>%
  filter(maturity_stage == "Spawning") %>%
  select(pop, week, week2, n_spawn)

obs_spawn <- obs_spawn %>% 
  left_join(obs_mature) %>%
  mutate(spawn_ratio = n_spawn/n_mature)

pred_spawn <- tibble()
for(pop_name in list_pop) {
  print(paste0("processing ", pop_name))
  
  # data
  obs_spawn_temp <- obs_spawn %>% filter(pop == pop_name)
  
  # model
  m_spawn <- glmmTMB(spawn_ratio ~ 1 + s(week2), 
                     family = beta_family(link = "logit"), 
                     data = obs_spawn_temp)
  
  # pred
  range_week2 = range(obs_spawn_temp$week2)
  pred <- tibble(pop = pop_name,
                 week2 = seq(min(range_week2), max(range_week2), 0.1))
  
  pred <- pred %>%
    mutate(spawn_ratio_pred = predict(m_spawn, newdata = pred, type = "response", se = T)$fit,
           spawn_ratio_se = predict(m_spawn, newdata = pred, type = "response", se = T)$se.fit) 
  
  pred_spawn <- bind_rows(pred_spawn, pred)
    
}

# plot
ggplot() +
  geom_jitter(data = obs_spawn, aes(x = week2, y = spawn_ratio, color = pop), alpha = 0.5) +
  geom_line(data = pred_spawn, aes(x = week2, y = spawn_ratio_pred, color = pop)) +
  geom_ribbon(data = pred_spawn, aes(x = week2,
                                   ymin = spawn_ratio_pred - 1.96*spawn_ratio_se,
                                   ymax = spawn_ratio_pred + 1.96*spawn_ratio_se,
                                   fill = pop),
              alpha = 0.5) +
  facet_grid(pop ~ .)

ggplot() +
  geom_line(data = pred_spawn, aes(x = week2, y = spawn_ratio_pred, color = pop)) +
  geom_ribbon(data = pred_spawn, aes(x = week2,
                                     ymin = spawn_ratio_pred - 1.96*spawn_ratio_se,
                                     ymax = spawn_ratio_pred + 1.96*spawn_ratio_se,
                                     fill = pop),
              alpha = 0.1)

# save file
write_rds(pred_spawn, "./output/pred_spawn.rds")
