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

# get only Developing and log K_rel and ssb_lag1
data <- data %>%
  filter(maturity_stage == "Developing") %>%
  mutate(log.K_rel = log(K_rel),
         log.ssb_lag1 = log(ssb_lag1)) %>%
  # temp
  mutate(log.temp = log(temp),
         log.temp_rel = log(temp/ave.temp),
         log.ave.temp = log(ave.temp))

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


# SUMMARY ----
# n and year self-sampling ----
sol_full <- read_rds("./data/sol_gonads_2004_2023.rds")
sol_full %>% group_by(TRI.ProjectTypeDescription) %>%
  summarise(range(TRI.Year))
table(sol_full$TRI.ProjectTypeDescription)
703/32629*100 # 2% seflsampling

# GSI range within data peak GSI (see 1_data_processing)

# FIGURE ----
## fig1 - sampling site + temp -----
source("./3_report_sampling-site.R")

## fig2 - peak gsi ----
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

ggsave(last_plot(), file = file.path(dir_report, "fig2_peak-gsi.pdf"),
       device = cairo_pdf,
       width =  17, height = 17,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig2_peak-gsi.png"),
       width =  17, height = 17,
       units = "cm",
       dpi = 1200)

## fig3 - individual effects ----
### body length ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- as.data.frame(Effect(c("log.length", "log.age", "log.K_rel"),
                                    model,
                                    xlevels = list("log.length" = unique(data_sub$log.length),
                                                   "log.age" = c(log(1), log(5)),
                                                   "log.K_rel" = c(log(0.5), log(1))))) %>%
    filter(log.age == max(log.age),
           log.K_rel == max(log.K_rel)) %>%
    mutate(pop = pop_name) 
  
  pred <- bind_rows(pred, pred_temp)
  
}

# est_ci95
est_ci95 <- rbind(deltaMethod(m1_4bc, "log.length"),
                  deltaMethod(m1_7a, "log.length"),
                  deltaMethod(m1_7fg, "log.length"),
                  deltaMethod(m1_7d, "log.length")
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
                        levels = c("β = 4.40 (3.94-4.85) - North Sea", 
                                   "β = 3.85 (3.60-4.10) - Irish Sea",
                                   "β = 4.29 (4.02-4.55) - Bristol Channel, Celtic Sea North",
                                   "β = 4.27 (4.10-4.44) - Eastern English Channel"
                        )))

# pred within range age 5
# get range length age 5
length_age5 <- data %>%
  filter(age == 5) %>%
  group_by(pop) %>%
  summarize(min_length = min(length),
            max_length = max(length))

pred_age5 <- pred %>%
  left_join(length_age5) %>%
  filter(exp(log.length) >= min_length, exp(log.length) <= max_length)


# plot
ggplot() +
  geom_line(data = pred_age5,
            aes(x = exp(log.length),
                y = exp(fit),
                color = label)) +
  # geom_line(data = pred,
  #           aes(x = exp(log.length),
  #               y = exp(fit),
  #               color = label,
  #               ),
  #           linetype = "dashed") +
  geom_ribbon(data = pred_age5,
              aes(x = exp(log.length),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = label),
              alpha = 0.3) +
  labs(x = "Total body length (mm)",
       y = "Ovary weight (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  theme(legend.position = c(0.2, 0.8), #0.28, 0.88 
        legend.background = element_rect(fill = "transparent"))
p1 <- last_plot()

# histogram plot
# ggplot(data = data, aes(x = length)) +
#   # inverse level + color values
#   geom_histogram(bins = 30, aes(fill = fct_relevel(pop, "7d", "7fg", "7a", "4bc")), position = "identity", alpha = 0.7) +
#   scale_color_manual(values = c("#d7191c","#fdae61","#abd9e9","#2c7bb6")) +
#   scale_fill_manual(values = c("#d7191c","#fdae61","#abd9e9","#2c7bb6")) +
#   #theme_void() +
#   theme(legend.position = "none") 
# p2 <- last_plot()

### age ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- as.data.frame(Effect(c("log.length", "log.age", "log.K_rel"),
                                    model,
                                    xlevels = list("log.length" = c(log(100), log(350)), 
                                                   "log.age" = unique(data_sub$log.age),
                                                   "log.K_rel" = c(log(0.5), log(1))))) %>%
    filter(log.length == max(log.length),
           log.K_rel == max(log.K_rel)) %>%
    mutate(pop = pop_name) 
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop)

# pred within range length 350
# get range length 350
age_length350 <- data %>%
  filter(length == 350) %>%
  group_by(pop) %>%
  summarize(min_age = min(age),
            max_age = max(age))

pred_length350 <- pred %>%
  left_join(age_length350) %>%
  filter(exp(log.age) >= min_age, exp(log.age) <= max_age)

# plot
ggplot() +
  geom_line(data = pred_length350,
            aes(x = exp(log.age),
                y = exp(fit),
                color = pop_name)) +
  # geom_line(data = pred,
  #           aes(x = exp(log.age),
  #               y = exp(fit),
  #               color = pop_name),
  #           linetype = "dashed") +
  geom_ribbon(data = pred_length350,
              aes(x = exp(log.age),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = pop_name),
              alpha = 0.3) +
  labs(x = "Age (year)",
       y = "Ovary weight (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  theme(legend.position = c(0.2, 0.8), #0.28, 0.88 
        legend.background = element_rect(fill = "transparent")) +
  theme(legend.position = "none")
p2 <- last_plot()

### condition ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- list_m1[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- as.data.frame(Effect(c("log.length", "log.age", "log.K_rel"),
                               model,
                               xlevels = list("log.body_weight" = c(log(150), log(350)),
                                              "log.age" = c(log(1), log(5)),
                                              "log.K_rel" = unique(data_sub$log.K_rel)))) %>%
    filter(log.length == max(log.length)) %>%
    filter(log.age == max(log.age)) %>%
    mutate(pop = pop_name) 

  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop)

# plot
ggplot() +
  geom_line(data = pred,
            aes(x = exp(log.K_rel),
                y = exp(fit),
                color = pop_name)) +
  geom_ribbon(data = pred,
              aes(x = exp(log.K_rel),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = pop_name),
              alpha = 0.3) +
  labs(x = "Le Cren’s relative condition factor",
       y = "Ovary weight (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  theme(legend.position = c(0.2, 0.8), #0.28, 0.88 
        legend.background = element_rect(fill = "transparent")) +
  theme(legend.position = "none")
p3 <- last_plot()

### merge plot ----

layout <- "
AAAAAA
AAAAAA
AAAAAA
AAAAAA
CCCDDD
CCCDDD
"

p1 + p2 + p3 + plot_layout(design = layout) 
# save plot
ggsave(last_plot(), file = file.path(dir_report, "fig3_individual-effects.png"),
       width =  17, height = 17, 
       units = "cm",
       dpi = 1200,
       scale = 1.7)
ggsave(last_plot(), file = file.path(dir_report, "fig3_individual-effects.pdf"),
       device = cairo_pdf,
       width =  17*1.7, height = 17*1.7,
       units = "cm")

## fig4 - random effect ----
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
           intercept_se = pred_temp_se$`(Intercept)`) %>%
    mutate(y = (exp(intercept)),
           ymin = (exp(intercept - intercept_se*1.96)),
           ymax = (exp(intercept + intercept_se*1.96)))
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop_2line)

ggplot(data = pred, aes(x = week, y = y)) +
  geom_linerange(data = pred, 
                 aes(ymin = ymin, 
                     ymax = ymax)) +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(pop_name ~ .) +
  labs(x = "Week",
       y = "Relative ovary weight") +
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
           intercept_se = pred_temp_se$`(Intercept)`) %>%
    mutate(y = (exp(intercept)),
           ymin = (exp(intercept - intercept_se*1.96)),
           ymax = (exp(intercept + intercept_se*1.96)))
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop_2line)

ggplot(data = pred, aes(x = year, y = y)) +
  geom_linerange(data = pred, 
                 aes(ymin = ymin, 
                     ymax = ymax)) +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(pop_name ~ .) +
  labs(x = "Year",
       y = "Relative ovary weight") +
  scale_x_continuous(breaks = seq(2004,2022,4))

p_year <- last_plot()

### merge plot ----
(p_week | p_year) +
  plot_annotation(tag_levels = 'A') 

ggsave(last_plot(), file = file.path(dir_report, "fig4_random-effect.png"),
       width =  17, height = 17,
       units = "cm",
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "fig4_random-effect.pdf"),
       device = cairo_pdf,
       width =  17, height = 17,
       units = "cm")

## fig5 - temp*length (min max) ----
pred <- tibble()
for(pop_name in c("4bc", "7a")) {
  
  print(paste("processing", pop_name))
  
  model <- list_m2[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  
  pred_temp <- as.data.frame(Effect(c("log.temp", "log.length", "log.age", "log.K_rel"),
                                    model,
                                    xlevels = list("log.temp" = seq(min(data_sub$log.temp), max(data_sub$log.temp), 0.001),
                                                   "log.length" = log(c(250, 350, 450)),
                                                   "log.age" = c(log(1), log(5)),
                                                   "log.K_rel" = c(log(0.5), log(1))
                                    ))) %>%
    filter(log.age == max(log.age),
           log.K_rel == max(log.K_rel)) %>%
    mutate(pop = pop_name) 
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop) %>%
  mutate(temp = exp(log.temp),
         length = round(exp(log.length)),
         fit_exp = exp(fit),
         lower_exp = exp(lower),
         upper_exp = exp(upper))

ggplot(data = pred, aes(x = temp, y = fit_exp )) +
  # geom_point(data = pred %>% filter(c.temp == 0), aes(x = c.temp, y = fit), color = "blue") +
  # geom_point(data = pred %>% filter(c.temp == 0.5), aes(x = c.temp, y = fit), color = "red") +
  geom_line(aes(linetype = factor(length), color = pop_name)) +
  geom_ribbon(aes(ymin = lower_exp, ymax = upper_exp, 
                  linetype = factor(length), fill = pop_name), alpha = 0.2) +
  #facet_grid(~ pop_name) +
  labs(x = "Autumn-Winter Temperature (°C)",
       y = "Ovary weight (g)",
       linetype = "Total body length (mm)",
       color = "Population",
       fill = "Population") +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title.align = 0.5) +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9"))

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig5_temp.png"),
       width =  8.5, height = 11,
       units = "cm",
       dpi = 1200,
       scale = 1.5)
ggsave(last_plot(), file = file.path(dir_report, "fig5_temp.pdf"),
       device = cairo_pdf,
       width =  8.5*1.5, height = 11*1.5,
       units = "cm")

### summary ----
pred <- tibble()
for(pop_name in c("4bc", "7a")) {
  
  print(paste("processing", pop_name))
  
  model <- list_m2[[pop_name]]
  
  data_sub <- data %>% filter(pop == pop_name)
  ave.temp = unique(data_sub$ave.temp)
  
  pred_temp <- as.data.frame(Effect(c("log.temp", "log.length", "log.age", "log.K_rel"),
                                    model,
                                    xlevels = list("log.temp" = c(log(ave.temp), log(ave.temp + 0.5)),
                                                   "log.length" = log(c(250, 350, 450)),
                                                   "log.age" = c(log(1), log(5)),
                                                   "log.K_rel" = c(log(0.5), log(1))
                                    ))) %>%
    filter(log.age == max(log.age),
           log.K_rel == max(log.K_rel)) %>%
    mutate(pop = pop_name)
  
  pred <- bind_rows(pred, pred_temp)
  
}

pred <- pred %>%
  left_join(df_pop) %>%
  group_by(pop) %>%
  mutate(temp = exp(log.temp),
         ave.temp = min(temp),
         length = round(exp(log.length)),
         fit_exp = exp(fit),
         lower_exp = exp(lower),
         upper_exp = exp(upper))

pred_sum <- pred %>%
  mutate(temp_scale = round((temp - ave.temp), 1)) %>%
  select(pop, length, temp_scale, fit_exp) %>%
  pivot_wider(names_from = temp_scale, values_from = fit_exp) %>%
  mutate(diff = (`0.5` - `0`)/`0`*100) %>%
  mutate(diff_name = sprintf("%.1f", diff)) 

# TABLE ----
## table 1 - in text ----
# SUPPLEMENTARY INFORMATION -----
## table S1-S3 - in text ----
## figS - distribution gsi ----
obs_sub2 <- obs_sub %>%
  filter(maturity_stage %in% c("Developing")) %>%
  left_join(pred_gsi_sub) %>%
  filter(week2 >= week2_start,
         week2 <= week2_end)

obs_sub2 <- obs_sub2 %>%
  left_join(df_pop_2line)

ggplot(data = obs_sub2 %>% filter(gsi < 30), aes(x = gsi)) +
  geom_histogram(bins = 100, position = "identity") +
  geom_vline(xintercept = 2, linetype = "dashed") +
  facet_grid(. ~ pop_name) +
  labs(x = "Gonadosomatic index - GSI",
       y = "Count") +
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
# mean gsi (before removing small and large)
obs_sub2 %>% 
  group_by(pop) %>%
  summarize(mean = mean(gsi),
            sd = sd(gsi)) %>%
  arrange(mean)

quantile(obs_sub2$gsi)
# upper and lower whisker
11.1+1.5*(11.1-4.9) # 20.4
4.9-1.5*(11.1-4.9) # -4.4
ggplot(data = obs_sub2, aes(x = gsi, y = id)) +
  geom_point(alpha = 0.3) +
  geom_vline(xintercept = 30) +
  geom_vline(xintercept = 20.4, linetype = "dashed") +
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
  MyVar <- c("log.length", "log.K_rel", "log.age", "log.temp", "log.ssb_lag1")
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

## figS - gonad weight ~ total body length -----
data <- data_obs %>%
  left_join(df_pop_2line)

ggplot(data = data, aes(x = length, y = gonad_weight)) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  facet_grid(. ~ pop_name) +
  labs(x = "Total body length (mm)",
       y = "Gonad weight (g)") #+
  #xlim(0, 2200)
p1 <- last_plot()

ggplot(data = data, aes(x = log(length), y = log(gonad_weight) )) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  facet_grid(. ~ pop_name) +
  labs(x = "ln(Total body length) (mm)",
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
ggsave(last_plot(), file = file.path(dir_report, "figS_gonad-weight-vs-body-length.pdf"),
       device = cairo_pdf,
       width =  17, height = 13,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_gonad-weight-vs-body-length.png"),
       width =  17, height = 13,
       units = "cm",
       dpi = 1200)
## figS - model diagnostic ----
### m1 ----
pred <- tibble()
for(pop_name2 in unique(data$pop)) {
  
  print(paste("processing", pop_name2))
  
  model <- list_m1[[pop_name2]]
  
  data_sub <- data %>% filter(pop == pop_name2)
  
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
  ggplot(data = pred_sub, aes(x = log.length, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "ln(Total body length) (g)",
         y = "Residuals") 
  p2 <- last_plot()

  ggplot(data = pred_sub, aes(x = log.age, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "ln(Age) (year)",
         y = "Residuals") 
  p3 <- last_plot()
  
  ggplot(data = pred_sub, aes(x = log.K_rel, y = res)) +
    geom_point(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "ln(Le Cren's relative condition factor)",
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
for(pop_name2 in unique(data$pop)) {
  
  print(paste("processing", pop_name2))
  
  model <- list_m2[[pop_name2]]
  
  data_sub <- data %>% filter(pop == pop_name2)
  
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
  ggplot(data = pred_sub, aes(x = temp, y = res)) +
    geom_point(alpha = 0.3, size = rel(0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x  = "ln(Autumn-winter temperature) (°C)",
         y = "Residuals") 
  p2 <- last_plot()
  
  # qqplot
  ggplot(data = pred_sub, aes(sample = res)) +
    stat_qq(alpha = 0.3) +
    stat_qq_line() +
    labs(x  = "Theoretical quantiles",
         y = "Sample quantiles") 
  p3 <- last_plot()
  
  list_plot[[pop_name2]] <- (p1/p2/p3)
}

list_plot[["4bc"]] | list_plot[["7a"]] | list_plot[["7fg"]] | list_plot[["7d"]]

ggsave(last_plot(), file = file.path(dir_report, "figS_m2_diagnostic.png"),
       width =  17, height = 12,
       units = "cm",
       scale = 2,
       dpi = 1200)
ggsave(last_plot(), file = file.path(dir_report, "figS_m2_diagnostics.pdf"),
       device = cairo_pdf,
       width =  17*2, height = 12*2,
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
  "ln(Total body length)",
  "ln(Age)",
  "ln(Relative condition factor)"
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
  "ln(Total body length)",
  "ln(Age)",
  "ln(Relative condition factor)",
  "ln(Autumn-winter temperature)",
  "ln(Autumn-winter temperature) x ln(Total body length)"
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
          #order.terms = c(1,2,3,8,6,4,5,7),
          string.pred = "Fixed Effects",
          string.est = "Estimate (95CI)",
          string.ci = "95% CI",
          #collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "tableS_m2.html")
)


## table sample size ----

df_summary <- data %>%
  count(year, pop) %>%              # count observations per year & pop
  pivot_wider(
    names_from = pop,
    values_from = n,
    values_fill = 0                 # fill missing combos with 0
  ) %>%
  select(year, `4bc`, `7a`, `7fg`, `7d`)

df_summary

table(data$pop)
tab_df(df_summary, 
          col.header = c("Year", 
                         "North Sea (n = 361)", 
                         "Irish Sea (n = 1078)", 
                         "Bristol Channel and Celtic Sea North (n = 1930)", 
                         "Eastern English Channel (n = 1268)"),
          file = file.path(dir_report, "tableS_analysed-data-summary.html")
)
