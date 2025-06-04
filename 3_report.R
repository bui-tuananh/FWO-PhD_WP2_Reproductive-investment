# BEFORE
# 2_data-analysis

# AFTER
# report

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
library(car)        # calculate regression coefficients
library(sjPlot)     # table
library(GGally)     # pair plots

# set theme
theme_set(theme_bw())

# set directory
dir_output <- "./output"
dir_report <- "./report"

# population name
df_pop <- tibble(pop = c("4bc", "7a", "7fg", "7d"),
                 pop_name = factor(c("North Sea", "Irish Sea", "Bristol Channel, Celtic Sea North", "Eastern English Channel"),
                                   levels = c("North Sea", "Irish Sea", "Bristol Channel, Celtic Sea North", "Eastern English Channel")))

df_pop_2line <- tibble(pop = c("4bc", "7a", "7fg", "7d"),
                       pop_name = factor(c("North Sea", "Irish Sea", "Bristol Channel\nCeltic Sea North", "Eastern English Channel"),
                                         levels = c("North Sea", "Irish Sea", "Bristol Channel\nCeltic Sea North", "Eastern English Channel")))

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

## model data ----
# model 1 - Individual-level variation and temporal trend
m1_4bc <- read_rds("./output/m1_4bc.rds")
m1_7a <- read_rds("./output/m1_7a.rds")
m1_7fg <- read_rds("./output/m1_7fg.rds")
m1_7d <- read_rds("./output/m1_7d.rds")
list_m1 <- c("4bc" = m1_4bc,
             "7a" = m1_7a,
             "7fg" = m1_7fg,
             "7d" = m1_7d)

# model 2 - Size-specific effects of warming
m2_4bc <- read_rds("./output/m2_4bc.rds")
m2_7a <- read_rds("./output/m2_7a.rds")
m2_7fg <- read_rds("./output/m2_7fg.rds")
m2_7d <- read_rds("./output/m2_7d.rds")
list_m2 <- c("4bc" = m2_4bc,
             "7a" = m2_7a,
             "7fg" = m2_7fg,
             "7d" = m2_7d)

# model 3 - Population-specific effects of warming
m3 <- read_rds("./output/m3.rds")
m3_temp <- read_rds("./output/m3_temp.rds")

# FIGURE ----
## fig1 - sampling site + temp -----
source("./3_report_sampling-site.R")

## fig2 - individual effects ----
### body weight ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)

  if(pop_name != "7d") {
    pred_temp <- as.data.frame(Effect(c("log.body_weight"),
                                      model,
                                      xlevels = list("log.body_weight" = unique(data_sub$log.body_weight)))) %>%
      mutate(pop = pop_name)
  } else {
    pred_temp <- as.data.frame(Effect(c("log.body_weight", "K_rel"),
                                      model,
                                      xlevels = list("log.body_weight" = unique(data_sub$log.body_weight),
                                                     "K_rel" = c(0,1)))) %>%
      filter(K_rel == max(K_rel)) %>%
      mutate(pop = pop_name)
  }

  pred <- bind_rows(pred, pred_temp)
  
}

# est_ci95
est_ci95 <- rbind(deltaMethod(m1_4bc, "log.body_weight"),
                  deltaMethod(m1_7a, "log.body_weight"),
                  deltaMethod(m1_7fg, "log.body_weight"),
                  deltaMethod(m1_7d, "log.body_weight")
) %>%
  as.data.frame() %>%
  mutate(pop = c("4bc", "7a", "7fg", "7d")) %>%
  mutate(est_ci95 = sprintf("%.2f (%.2f-%.2f)", Estimate, `2.5 %`, `97.5 %`)) %>%
  select(pop, est_ci95)

# label pred
pred <- pred %>%
  left_join(df_pop) %>%
  left_join(est_ci95) %>%
  mutate(label = paste0("\u03B2 = ", est_ci95, " - ", pop_name)) 

unique(pred$label)
pred <- pred %>%
  mutate(label = factor(label, 
                        levels = c("β = 1.37 (1.33-1.42) - North Sea", 
                                   "β = 1.31 (1.27-1.34) - Irish Sea",
                                   "β = 1.41 (1.38-1.45) - Bristol Channel, Celtic Sea North",
                                   "β = 1.32 (1.29-1.34) - Eastern English Channel"
                        )))

# plot
ggplot() +
  geom_line(data = pred,
            aes(x = log.body_weight,
                y = fit,
                color = label)) +
  geom_ribbon(data = pred,
              aes(x = log.body_weight,
                  ymin = lower,
                  ymax = upper,
                  fill = label),
              alpha = 0.3) +
  labs(x = "ln(Gutted body weight) (g)",
       y = "ln(Gonad weight) (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  theme(legend.position = c(0.2, 0.8), #0.28, 0.88 
        legend.background = element_rect(fill = "transparent"))
p1 <- last_plot()

# histogram plot
ggplot(data = data, aes(x = log.body_weight)) +
  # inverse level + color values
  geom_histogram(bins = 200, aes(fill = fct_relevel(pop, "7d", "7fg", "7a", "4bc")), position = "identity", alpha = 0.7) +
  scale_color_manual(values = c("#d7191c","#fdae61","#abd9e9","#2c7bb6")) +
  scale_fill_manual(values = c("#d7191c","#fdae61","#abd9e9","#2c7bb6")) +
  theme_void() +
  theme(legend.position = "none") 
p2 <- last_plot()

### maturity_stage ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  if(pop_name != "7d") {
    pred_temp <- as.data.frame(Effect(c("maturity_stage", "log.body_weight"),
                                      model,
                                      xlevels = list("log.body_weight" = c(log(150), log(300)) ))) %>%
      filter(log.body_weight == max(log.body_weight)) %>%
      mutate(pop = pop_name)
  } else {
    pred_temp <- as.data.frame(Effect(c("maturity_stage", "log.body_weight", "K_rel"),
                                      model,
                                      xlevels = list("log.body_weight" = c(log(150), log(300)),
                                                     "K_rel" = c(0,1)))) %>%
      filter(log.body_weight == max(log.body_weight), 
             K_rel == max(K_rel)) %>%
      mutate(pop = pop_name)
  }
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop_2line)

# plot
ggplot() +
  geom_linerange(data = pred, 
                 aes(x = pop_name, 
                     y = fit,
                     ymin = lower, 
                     ymax = upper,
                     group = maturity_stage),
                 position = position_dodge(width = 0.5)
  ) +
  geom_point(data = pred, 
             aes(x = pop_name, 
                 y = fit,
                 fill = maturity_stage),
             shape = 21,
             position = position_dodge(width = 0.5)) +
  
  labs(x = "Population",
       y = "ln(Gonad weight) (g)",
       fill = "Maturity stage") +
  scale_fill_manual(values = c("white", "black")) +
  # theme(legend.position = "bottom",
  #       legend.title.position = "top",
  #       legend.title.align = 0.5) +
  theme(legend.position = c(0.8, 0.8), #0.28, 0.88 
        legend.background = element_rect(fill = "transparent"))

p3 <- last_plot()

### condition 7d ----
model <- m1_7d
data_sub <- data %>% filter(pop == "7d")

pred <- as.data.frame(Effect(c("log.body_weight", "K_rel"),
                                  model,
                                  xlevels = list("log.body_weight" = c(log(150), log(300)),
                                                 "K_rel" = unique(data_sub$K_rel)))) %>%
  filter(log.body_weight == max(log.body_weight)) %>%
  mutate(pop = "7d")

# plot
ggplot() +
  geom_line(data = pred, 
            aes(x = K_rel, 
                y = fit),
            color = "#d7191c") +
  geom_ribbon(data = pred, 
              aes(x = K_rel, 
                  ymin = lower, 
                  ymax = upper),
              fill = "#d7191c",
              alpha = 0.3) +
  labs(x = "Le Cren’s relative condition factor",
       y = "ln(Gonad weight) (g)") 
p4 <- last_plot()

### merge plot ----

layout <- "
AAAAAA
BBBBBB
BBBBBB
BBBBBB
CCCDDD
CCCDDD
"

# add tag
p2_ <- p2 + labs(tag = "A")
p3_ <- p3 + labs(tag = "B")
p4_ <- p4 + labs(tag = "C")

(p2_ + p1) + p3_ + p4_ + plot_layout(design = layout) 
# save plot
ggsave(last_plot(), file = file.path(dir_report, "fig2_individual-effects.png"),
       width =  17, height = 17, 
       units = "cm",
       dpi = 1200,
       scale = 1.7)
ggsave(last_plot(), file = file.path(dir_report, "fig2_individual-effects.pdf"),
       device = cairo_pdf,
       width =  17*1.7, height = 17*1.7,
       units = "cm")

## fig3 - random effect ----
### week ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- ranef(model)$`week` 
  pred_temp_se <- as.data.frame(arm::se.ranef(model)$`week`)
  
  pred_temp <- pred_temp %>%
    mutate(pop = pop_name,
           week = as.numeric(rownames(pred_temp)),
           intercept = `(Intercept)`,
           intercept_se = pred_temp_se$`(Intercept)`)
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop_2line)

ggplot(data = pred, aes(x = week, y = intercept)) +
  geom_linerange(data = pred, 
                 aes(ymin = intercept - 1.96*intercept_se, 
                     ymax = intercept + 1.96*intercept_se)) +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(pop_name ~ .) +
  labs(x = "Week",
       y = "Random effect") +
  scale_x_continuous(breaks = seq(2,18,2))

p_week <- last_plot()

### year ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- ranef(model)$`year` 
  pred_temp_se <- as.data.frame(arm::se.ranef(model)$`year`)
  
  pred_temp <- pred_temp %>%
    mutate(pop = pop_name,
           year = as.numeric(rownames(pred_temp)),
           intercept = `(Intercept)`,
           intercept_se = pred_temp_se$`(Intercept)`)
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop_2line)

ggplot(data = pred, aes(x = year, y = intercept)) +
  geom_linerange(data = pred, 
                 aes(ymin = intercept - 1.96*intercept_se, 
                     ymax = intercept + 1.96*intercept_se)) +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(pop_name ~ .) +
  labs(x = "Year",
       y = "Random effect") +
  scale_x_continuous(breaks = seq(2004,2022,4))

p_year <- last_plot()

### merge plot ----
(p_week | p_year) +
  plot_annotation(tag_levels = 'A') 

ggsave(last_plot(), file = file.path(dir_report, "fig3_random-effect.png"),
       width =  17, height = 17,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "fig3_random-effect.pdf"),
       device = cairo_pdf,
       width =  17, height = 17,
       units = "cm")

## fig 4 - temp*weight (0.5 deg) ----
pred <- tibble()
for(pop_name in c("4bc", "7a", "7d")) {
  
  print(paste("processing", pop_name))
  
  model <- list_m2[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  if(pop_name != "7d") {
    pred_temp <- as.data.frame(Effect(c("c.temp", "log.body_weight"),
                                      model,
                                      xlevels = list("log.body_weight" = log(c(150, 300, 600, 900)),
                                                     "c.temp" = c(0, 0.5)))) %>%
      mutate(pop = pop_name)
  } else {
    pred_temp <- as.data.frame(Effect(c("c.temp", "log.body_weight", "K_rel"),
                                      model,
                                      xlevels = list("log.body_weight" = log(c(150, 300, 600, 900)),
                                                     "c.temp" = c(0, 0.5),
                                                     "K_rel" = c(0,1)))) %>%
      filter(K_rel == max(K_rel)) %>%
      mutate(pop = pop_name)
  }
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop) %>%
  mutate(body_weight = round(exp(log.body_weight)))

ggplot(data = pred, aes(x = c.temp, y = fit)) +
  # geom_point(data = pred %>% filter(c.temp == 0), aes(x = c.temp, y = fit), color = "blue") +
  # geom_point(data = pred %>% filter(c.temp == 0.5), aes(x = c.temp, y = fit), color = "red") +
  geom_line(aes(linetype = factor(body_weight) )) +
  geom_ribbon(aes(ymin = lower, ymax = upper, linetype = factor(body_weight)), alpha = 0.1) +
  facet_grid(~ pop_name) +
  labs(x = "Temperature anomaly (°C)",
       y = "ln(Gonad weight) (g)",
       linetype = "Gutted body weight (g)") +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title.align = 0.5)

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig4_temp.png"),
       width =  17, height = 11,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "fig4_temp.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")

### summary ----
pred_sum <- pred %>%
  select(pop, body_weight, c.temp, fit) %>%
  pivot_wider(names_from = c.temp, values_from = fit) %>%
  mutate(diff = (`0.5` - `0`)/`0`*100) %>%
  mutate(diff_name = sprintf("%.1f", diff)) %>%
  mutate(diff_exp = (exp(`0.5`) - exp(`0`) )/exp(`0`)*100) %>%
  mutate(diff_exp_name = sprintf("%.1f", diff_exp))


## fig 5 - temp*weight + temp*ssb_lag1 7a ----
pred <- tibble()
for(pop_name in c("7a")) {
  
  print(paste("processing", pop_name))
  
  model <- list_m2[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  ssb_grp = tibble(ssb_name = factor(c("P25", "Mean", "P75"),
                                     levels = c("P25", "Mean", "P75")),
                   ssb_lag1 = c(quantile(data_sub$ssb_lag1, probs = c(0.25)), 
                                mean(data_sub$ssb_lag1), 
                                quantile(data_sub$ssb_lag1, probs = c(0.75)))
                   )
  
  pred_temp <- as.data.frame(Effect(c("c.temp", "log.body_weight", "ssb_lag1"),
                                    model,
                                    xlevels = list("log.body_weight" = log(c(150, 300, 600, 900)),
                                                   "c.temp" = c(0, 0.5, 1),
                                                   "ssb_lag1" = ssb_grp$ssb_lag1))) %>%
    mutate(pop = pop_name)
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop) %>%
  left_join(ssb_grp) %>%
  mutate(body_weight = round(exp(log.body_weight))) 

ggplot(data = pred, aes(x = c.temp, y = fit)) +
  geom_line(aes(linetype = factor(body_weight) )) +
  geom_ribbon(aes(ymin = lower, ymax = upper, linetype = factor(body_weight)), alpha = 0.1) +
  facet_grid(~ ssb_name) +
  labs(x = "Temperature anomaly (°C)",
       y = "ln(Gonad weight) (g)",
       linetype = "Gutted body weight (g)") +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title.align = 0.5)

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig5_temp_ssb_7a.png"),
       width =  17, height = 11,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "fig5_temp_ssb_7a.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")

# ssb 
ggplot(data = pred, aes(x = c.temp, y = fit)) +
  geom_line(aes(linetype = factor(body_weight), color = ssb_name )) +
  geom_ribbon(aes(ymin = lower, ymax = upper, linetype = factor(body_weight), fill = ssb_name), alpha = 0.1) +
  labs(x = "Temperature anomaly (°C)",
       y = "ln(Gonad weight) (g)",
       linetype = "Gutted body weight (g)") +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title.align = 0.5)

## fig6 - pop ----
model <- m3
data_sub <- data 
pred <- as.data.frame(Effect(c("pop", "log.body_weight", "K_rel"),
                                  model,
                                  xlevels = list("log.body_weight" = log(c(150, 300, 600, 900)),
                                                 "K_rel" = c(0,1)))) %>%
  filter(K_rel == max(K_rel))

# add average temp data
data_temp_aut_win_ave <- data_temp_aut_win %>%
  select(pop, ave.temp) %>%
  unique()

pred <- pred %>%
  left_join(data_temp_aut_win_ave) %>%
  left_join(df_pop)

# model with temperature
pred_temp <- as.data.frame(Effect(c("ave.temp", "log.body_weight", "K_rel"),
                                  m3_temp,
                                  xlevels = list("log.body_weight" = log(c(150, 300, 600, 900)),
                                                 "K_rel" = c(0,1)))) %>%
  filter(K_rel == max(K_rel))

# plot
# create annotation
df_text <- tibble(x = 11.75,
                  y = c(2.4, 3.4, 4.3, 4.9),
                  label = c("150 g", "300 g", "600 g", "900 g"))
ggplot() +
  geom_line(data = pred_temp, 
            aes(x = ave.temp, 
                y = fit, 
                group = log.body_weight),
            linetype = "dashed") +
  geom_linerange(data = pred, 
                 aes(x = ave.temp, 
                     y = fit, 
                     ymin = lower, 
                     ymax = upper, 
                     color = pop_name)) +
  geom_point(data = pred, 
             aes(x = ave.temp, 
                 y = fit, 
                 color = pop_name), shape = 21, fill = "white") + 
  geom_text(data = df_text, aes(x = x, y = y, label = label)) +
  labs(x = "Average temperature (°C)",
       y = "ln(Gonad weight) (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c") )+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title.align = 0.5) +
  guides(color = guide_legend(nrow = 2))

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig6_pop.png"),
       width =  8.5, height = 11,
       units = "cm",
       dpi = 1200,
       scale = 1.2)
ggsave(last_plot(), file = file.path(dir_report, "fig6_pop.pdf"),
       device = cairo_pdf,
       width =  8.5*1.2, height = 5.5*1.2,
       units = "cm")

# TABLE ----
## table 1 - in text ----
# SUPPLEMENTARY INFORMATION -----
## table S1-S3 - in text ----
## figS - peak gsi ----
# pred gsi 
pred_gsi <- read_rds("./output/pred_gsi.rds")
pred_gsi_sub <- pred_gsi %>%
  group_by(pop) %>%
  mutate(gsi50 = max(gsi_pred)*0.5) %>%
  filter(gsi_pred >= gsi50) %>%
  summarize(week2_start = round(min(week2)),
            week2_end   = round(max(week2)))

pred_gsi <- pred_gsi %>%
  left_join(df_pop_2line) %>%
  left_join(pred_gsi_sub)

# obs data (the same as in analysis)
obs_sub <- read_rds("./data/sol_gonads_2004_2023_sub.rds")

obs_sub <- obs_sub %>% 
  filter(HAU.IcesArea %in% c("4b", "4c", "7a", "7d", "7f", "7g"),
         is.na(SPE.MaturityStageDescription) == F,
         !SPE.MaturityStageDescription %in% c("Abnormal", "Skipped spawning"), #2 obs Abnormal, 4 obs Skipped spawning
         is.na(SPA.Age) == F,
         SPE.Sex == "F") %>%
  mutate(SPE.MaturityStageDescription = if_else(SPE.MaturityStageDescription %in% c("Maturing", "Developing"),
                                                "Developing",
                                                SPE.MaturityStageDescription))
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

obs_sub <- obs_sub %>%
  filter(gsi < 100,
         is.na(gsi) == F) %>%  #3 obs with NA gsi
  mutate(log.gonad_weight = log(gonad_weight),
         log.body_weight = log(body_weight),
         log.length = log(length),
         week2 = if_else(week >= 40, week - 53, week),
         year2 = if_else(week2 <0, year+1, year))

obs_sub_lt30 <- obs_sub %>%
  filter(gsi < 30)

obs_sub_lt30 <- obs_sub_lt30 %>%
  left_join(df_pop_2line)

# plot
df_month <- obs_sub_lt30 %>% 
  select(month, week2) %>%
  arrange(month, week2) %>%
  unique()

one_month <- 53/12
start_month <- -13
end_month <- 40
break_month <- seq(start_month, end_month,  one_month)
name_month <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "")

ggplot() +
  geom_jitter(data = obs_sub_lt30, aes(x = week2, y = gsi, color = maturity_stage), alpha = 0.3, size = rel(0.8)) +
  geom_line(data = pred_gsi, aes(x = week2, y = gsi_pred)) +
  geom_ribbon(data = pred_gsi, aes(x = week2, 
                                   ymin = gsi_pred - 1.96*gsi_se, 
                                   ymax = gsi_pred + 1.96*gsi_se),
              alpha = 0.5) +
  geom_vline(data = pred_gsi, aes(xintercept = week2_start), linetype = "dashed") +
  geom_vline(data = pred_gsi, aes(xintercept = week2_end), linetype = "dashed") +
  facet_grid(pop_name ~ .) +
  labs(x = NULL,
       y = "Gonadosomatic index - GSI",
       color = "Maturity stage") +
  scale_x_continuous(breaks = break_month[1:13], labels = name_month) +
  theme(axis.text.x = element_text(vjust = 0, hjust = -0.5),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title.align = 0.5)  +
  scale_color_manual(values = c("black", "#66c2a5", "#fc8d62", "#8da0cb"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
p1 <- last_plot()

# save plot 
(p1) &
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.95, 0.8),
        title = element_text(size = 9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text = element_text(size = 8),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.direction = "horizontal",
        legend.box="vertical",
        legend.title.position = "top",
        legend.title.align = 0.5
  ) &
  guides(color = guide_legend(override.aes = list(alpha = 1)))

ggsave(last_plot(), file = file.path(dir_report, "figS_peak-gsi.pdf"),
       device = cairo_pdf,
       width =  17, height = 17,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_peak-gsi.png"),
       width =  17, height = 17,
       units = "cm",
       dpi = 1200)

## figS - distribution gsi ----
obs_sub2 <- obs_sub2 %>%
  left_join(df_pop_2line)

ggplot(data = obs_sub2 %>% filter(gsi < 30), aes(x = gsi, fill = maturity_stage)) +
  geom_histogram(bins = 100, position = "identity") +
  geom_vline(xintercept = 2, linetype = "dashed") +
  facet_grid(. ~ pop_name) +
  labs(x = "Gonadosomatic index - GSI",
       y = "Count",
       fill = "Maturity stage") +
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "bottom", 
        legend.title.position = "top",
        legend.title.align = 0.5)

ggsave(last_plot(), file = file.path(dir_report, "figS_distribution-gsi.png"),
       width = 17, height = 11,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_distribution-gsi.pdf"),
       device = cairo_pdf,
       width = 17, height = 11,
       units = "cm")
## figS - outlier ----
obs_sub2 <- obs_sub %>%
  filter(maturity_stage %in% c("Developing", "Spawning")) %>%
  left_join(pred_gsi_sub) %>%
  filter(week2 >= week2_start,
         week2 <= week2_end)

quantile(obs_sub2$gsi)
# upper and lower whisker
11.99+1.5*(11.99-5.85) # 21.2
5.85-1.5*(11.99-5.85) # -3.36
ggplot(data = obs_sub2, aes(x = gsi, y = id)) +
  geom_point(alpha = 0.3) +
  geom_vline(xintercept = 30) +
  geom_vline(xintercept = 21.2, linetype = "dashed") +
  labs(x = "Gonadosomatic index - GSI",
       y = "Observation ID") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS_outlier.png"),
       width =  17, height = 11,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_outlier.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")

## figS - temp ----
data_temp_aut_win <- data_temp_aut_win %>% 
  left_join(df_pop) %>%
  filter(pop %in% c("4bc", "7a", "7fg", "7d"))

# c.temp
ggplot() +
  geom_rect(aes(xmin = 2004, xmax = 2022, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.2) +
  geom_line(data = data_temp_aut_win, aes(x = year, y = c.temp, color = pop_name)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year",
       y = "Autumn-winter temperature anomaly (°C)",
       color = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.title.align = 0.5)  +
  guides(color = guide_legend(nrow = 2))

ggsave(last_plot(), file = file.path(dir_report, "figS_temp.png"),
       width =  17, height = 12,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_temp.pdf"),
       device = cairo_pdf,
       width =  17, height = 12,
       units = "cm")

## figS - survey data ----
#### setup
dir_datras <- "./data/ices"

hh <- read_rds(file.path(dir_datras, "sf_HHflats_1985till2022_in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS"))
hl <- read_rds(file.path(dir_datras, "sf_Solea_soleaHL_withAbs_BTS+BTS-VIII+DYFS+SNS_in1985till2021_in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS"))

hh_loc <- hh %>% 
  select(Survey, Area_27) %>%
  filter(Area_27 %in% c("4.a", "4.b", "4.c", "7.a", "7.f", "7.g", "7.d")) %>%
  unique()

hl_loc <- hl %>% 
  filter(!is.na(HLNoAtLngt),
         #Quarter %in% c(1,3,4),
         Area_27 %in% c("4.a", "4.b", "4.c", "7.a", "7.f", "7.g", "7.d")) %>%
  select(Survey, Area_27) %>%
  unique()

# gis
dir_gis <- "./data/admin"
# continents
continents <- read_sf(file.path(dir_gis, "esri_continent.gpkg"))
# countries
countries <- read_sf(file.path(dir_gis, "esri_countries.gpkg"))
# ices area 
ices_area <- read_sf(file.path(dir_gis, "ices_areas_sub_group_4abc_4326_new.gpkg")) 
ices_area <- ices_area %>% 
  filter(Area_27 %in% c("4abc", "7a", "7fg", "7d")) %>% 
  st_simplify(dTolerance = 2000)

#### plot
p1 <- ggplot() +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.5) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 1) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = hh_loc, aes(color = Survey)) +
  coord_sf(xlim = c(-9.5, 9.5), ylim = c(48.5, 62.5), expand = FALSE) +
  theme_bw()

p2 <- ggplot() +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.5) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 1) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = hl_loc, aes(color = Survey)) +
  coord_sf(xlim = c(-9.5, 9.5), ylim = c(48.5, 62.5), expand = FALSE) +
  theme_bw()

(p1 + p2) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")  &
  theme(legend.position = "bottom")

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS_datras survey.png"),
       width = 17, height = 11,
       units = "cm",
       scale = 1.5,
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_datras survey.pdf"),
       device = cairo_pdf,
       width = 17, height = 11,
       units = "cm")

## appendixS - collinearity ----
# functions from Zuur et al., 2010
source("./ref/HighstatLib.r") 

df_vif <- tibble()
for(pop_name in c("4bc", "7a", "7fg", "7d")) {
  print(pop_name)
  
  data_sub <- data %>%
    filter(pop == pop_name) 
  MyVar <- c("log.body_weight", "K_rel", "c.temp", "ssb_lag1")
  df_vif_temp <- corvif(data_sub[,MyVar]) 
  df_vif_temp <- df_vif_temp %>% 
    mutate(var = row.names(df_vif_temp),
           VIF = sprintf("%.2f", GVIF),
           pop = pop_name) %>%
    select(var, VIF, pop)
  
  df_vif <- bind_rows(df_vif, df_vif_temp)
}

df_vif <- df_vif %>%
  pivot_wider(names_from = pop, values_from = VIF)

tab_df(df_vif,
       file = file.path("./report", "tableS_collinearity.html"))

# ggpairs 7a
data_sub <- data %>%
  filter(pop == "7a") 
ggpairs(data_sub, 
        columns = c("log.body_weight", "K_rel", "c.temp", "ssb_lag1"),
        columnLabels = c("ln(W)", "K", "T", "SSB"))


ggsave(last_plot(), file = file.path(dir_report, "figS_correlation_7a.png"),
       width =  17, height = 15,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_temp.pdf"),
       width =  17, height = 12,
       units = "cm",
       dpi = 1200)

## figS - gonad weight ~ body weight -----
data <- data_obs %>%
  left_join(df_pop_2line)

ggplot(data = data, aes(x = body_weight, y = gonad_weight)) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  facet_grid(. ~ pop_name) +
  labs(x = "Gutted body weight (g)",
       y = "Gonad weight (g)") +
  xlim(0, 2200)
p1 <- last_plot()

ggplot(data = data, aes(x = log(body_weight), y = log(gonad_weight) )) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  facet_grid(. ~ pop_name) +
  labs(x = "ln(Gutted body weight) (g)",
       y = "ln(Gonad weight) (g)") 
p2 <- last_plot()

(p1 / p2) & 
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

# save file
ggsave(last_plot(), file = file.path(dir_report, "figS_gonad-weight-vs-body-weight.pdf"),
       device = cairo_pdf,
       width =  17, height = 13,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_gonad-weight-vs-body-weight.png"),
       width =  17, height = 13,
       units = "cm",
       dpi = 1200)
## figS - model diagnostic ----
### m1 ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- data_sub %>%
    mutate(res = resid(model),
           pred = predict(model, type = "response")) 
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop)

list_plot <- list()
for(pop_name2 in unique(data$pop)) {
  
  print(pop_name2)
  
  pred_sub <- pred %>%
    filter(pop == pop_name2)
  
  # res vs pred
  ggplot(data = pred_sub, aes(x = pred, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "Fitted values",
         y = "Residuals") +
    facet_grid(. ~ pop_name)
  p1 <- last_plot()

  # res vs vars
  ggplot(data = pred_sub, aes(x = log.body_weight, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "ln(Gutted body weight) (g)",
         y = "Residuals") 
  p2 <- last_plot()
  
  ggplot(data = pred_sub, aes(x = maturity_stage, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "Maturity stage",
         y = "Residuals")
  p3 <- last_plot()
  
  ggplot(data = pred_sub, aes(x = K_rel, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "Le Cren's relative condition factor",
         y = "Residuals") 
  p4 <- last_plot()
  
  # qqplot
  ggplot(data = pred_sub, aes(sample = res)) +
    stat_qq(alpha = 0.3) +
    stat_qq_line() +
    labs(x  = "Theoretical quantiles",
         y = "Sample quantiles") 
  p5 <- last_plot()
  
  list_plot[[pop_name2]] <- (p1/p2/p3/p4/p5)
}
list_plot[["4bc"]] | list_plot[["7a"]] | list_plot[["7fg"]] | list_plot[["7d"]]

ggsave(last_plot(), file = file.path(dir_report, "figS_m1_diagnostic.png"),
       width =  17, height = 17,
       units = "cm",
       scale = 2,
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_m1_diagnostics.pdf"),
       device = cairo_pdf,
       width =  17*2, height = 17*2,
       units = "cm")

### m2 ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m2[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- data_sub %>%
    mutate(res = resid(model),
           pred = predict(model, type = "response")) 
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop)

list_plot <- list()
for(pop_name2 in unique(data$pop)) {
  
  print(pop_name2)
  
  pred_sub <- pred %>%
    filter(pop == pop_name2)
  
  # res vs pred
  ggplot(data = pred_sub, aes(x = pred, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "Fitted values",
         y = "Residuals") +
    facet_grid(. ~ pop_name)
  p1 <- last_plot()
  
  # res vs vars
  ggplot(data = pred_sub, aes(x = ssb_lag1, y = res)) +
    geom_point(alpha = 0.3, size = rel(0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "Spawning stock biomass (kilotonnes)",
         y = "Residuals") 
  p2 <- last_plot()
  
  ggplot(data = pred_sub, aes(x = c.temp, y = res)) +
    geom_point(alpha = 0.3, size = rel(0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "Temperature anomaly (°C)",
         y = "Residuals") 
  p3 <- last_plot()
  
  # qqplot
  ggplot(data = pred_sub, aes(sample = res)) +
    stat_qq(alpha = 0.3) +
    stat_qq_line() +
    labs(x  = "Theoretical quantiles",
         y = "Sample quantiles") 
  p4 <- last_plot()
  
  list_plot[[pop_name2]] <- (p1/p2/p3/p4)
}

list_plot[["4bc"]] | list_plot[["7a"]] | list_plot[["7fg"]] | list_plot[["7d"]]

ggsave(last_plot(), file = file.path(dir_report, "figS_m2_diagnostic.png"),
       width =  17, height = 17,
       units = "cm",
       scale = 2,
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_m2_diagnostics.pdf"),
       device = cairo_pdf,
       width =  17*2, height = 17*2,
       units = "cm")


### m3 ----
model <- m3
data_sub <- data 
pred <- data_sub %>%
  mutate(res = resid(model),
         pred = predict(model, type = "response")) 

pred <- pred %>%
  left_join(df_pop_2line)

pred_sub <- pred
# res vs pred
ggplot(data = pred_sub, aes(x = pred, y = res)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x  = "Fitted values",
       y = "Residuals") 
p1 <- last_plot()

ggplot(data = pred_sub, aes(x = pop_name, y = res)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x  = "Population",
       y = "Residuals")
p2 <- last_plot()

# qqplot
ggplot(data = pred_sub, aes(sample = res)) +
  stat_qq(alpha = 0.3) +
  stat_qq_line() +
  labs(x  = "Theoretical quantiles",
       y = "Sample quantiles") 
p3 <- last_plot()

p1/p2/p3

ggsave(last_plot(), file = file.path(dir_report, "figS_m3_diagnostic.png"),
       width =  7, height = 12.75,
       units = "cm",
       scale = 2,
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_m3_diagnostics.pdf"),
       device = cairo_pdf,
       width =  7*2, height = 12.75*2,
       units = "cm")


## appendixS - model summary ----
css_list = list(
  css.firsttablerow =  "font-weight:bold; font-style:normal; border:1px solid;",
  css.depvarhead = "font-style:normal; font-weight:bold; text-align:center; padding-top:.8em; border:1px solid;",
  css.firsttablecol = "text-align:left; border:1px solid;",
  css.centeralign = "text-align:center; border:1px solid;",
  css.randomparts = "font-weight:bold; text-align:left; padding-top:.8em; border:1px solid;",
  css.summary = "border:1px solid;",
  css.summarydata = "text-align:left;"
)

### tableS - m1 ----
pred_labels <- c(
  "Intercept",
  "ln(Gutted body weight)",
  "Maturity stage (Spawning)",
  "Relative condition factor"
)

# summary table
tab_model(m1_4bc, m1_7a, m1_7fg, m1_7d, 
          show.se = NULL, 
          #show.ci = TRUE, 
          collapse.ci = TRUE,
          show.p = TRUE,
          p.style = "stars", 
          #dv.labels = NULL,
          dv.labels = c("North Sea", "Irish Sea", "Bristol Channel and Celtic Sea North", "Eastern English Channel"),
          pred.labels = pred_labels,
          string.pred = "Fixed Effects",
          string.est = "Estimate (95CI)",
          string.ci = "95% CI",
          #collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "tableS_m1.html")
)

### tableS - m2 ----
# labels should be in the default order of vars in tab_model
pred_labels <- c(
  "Intercept",
  "ln(Gutted body weight)",
  "Maturity stage (Spawning)",
  "Temperature anomaly",
  "Temperature anomaly x ln(Gutted body weight)",
  "Spawning stock biomass",
  "Temperature anomaly x Spawning stock biomass",
  "Relative condition factor"
)

# summary table
tab_model(m2_4bc, m2_7a, m2_7fg, m2_7d, 
          show.se = NULL, 
          #show.ci = TRUE, 
          collapse.ci = TRUE,
          show.p = TRUE,
          p.style = "stars", 
          #dv.labels = NULL,
          dv.labels = c("North Sea", "Irish Sea", "Bristol Channel and Celtic Sea North", "Eastern English Channel"),
          pred.labels = pred_labels,
          order.terms = c(1,2,3,8,6,4,5,7),
          string.pred = "Fixed Effects",
          string.est = "Estimate (95CI)",
          string.ci = "95% CI",
          #collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "tableS_m2.html")
)

### tableS - m3 ----
pred_labels <- c(
  "Intercept",
  "ln(Gutted body weight)",
  "Maturity stage (Spawning)",
  "Relative condition factor",
  "Population (Irish Sea)",
  "Population (Eastern English Channel)",
  "Population (Bristol Channel and Celtic Sea North)",
  "Temperature anomaly",
  "ln(Gutted body weight) x Population (Irish Sea)",
  "ln(Gutted body weight) x Population (Eastern English Channel)",
  "ln(Gutted body weight) x Population (Bristol Channel and Celtic Sea North)",
  "Temperature anomaly x ln(Gutted body weight)"
)

# summary table
tab_model(m3, 
          show.se = NULL, 
          #show.ci = TRUE, 
          collapse.ci = TRUE,
          show.p = TRUE,
          p.style = "stars", 
          #dv.labels = NULL,
          dv.labels = c("Model 3"),
          pred.labels = pred_labels,
          order.terms = c(1,2,3,4,5,7,6,9,11,10,8,12),
          string.pred = "Fixed Effects",
          string.est = "Estimate (95CI)",
          string.ci = "95% CI",
          #collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "tableS_m3.html")
)
## appendixS - m3_temp ----

css_list = list(
  css.firsttablerow =  "font-weight:bold; font-style:normal; border:1px solid;",
  css.depvarhead = "font-style:normal; font-weight:bold; text-align:center; padding-top:.8em; border:1px solid;",
  css.firsttablecol = "text-align:left; border:1px solid;",
  css.centeralign = "text-align:center; border:1px solid;",
  css.randomparts = "font-weight:bold; text-align:left; padding-top:.8em; border:1px solid;",
  css.summary = "border:1px solid;",
  css.summarydata = "text-align:left;"
)

pred_labels <- c(
  "Intercept",
  "ln(Gutted body weight)",
  "Maturity stage (Spawning)",
  "Relative condition factor",
  "Population's average temperature",
  "Temperature anomaly",
  "Temperature anomaly x ln(Gutted body weight)"
)

# summary table
tab_model(m3_temp, 
          show.se = NULL, 
          #show.ci = TRUE, 
          collapse.ci = TRUE,
          show.p = TRUE,
          p.style = "stars", 
          #dv.labels = NULL,
          dv.labels = c("Model 3 extended"),
          pred.labels = pred_labels,
          string.pred = "Fixed Effects",
          string.est = "Estimate (95CI)",
          string.ci = "95% CI",
          #collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "tableS_m3_temp.html")
)
