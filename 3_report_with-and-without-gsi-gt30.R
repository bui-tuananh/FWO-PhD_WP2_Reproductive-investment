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

## model data ----
dir_output_without <- "./output"
dir_output_with <- "./output/with-gsi-gt30"

# model 2 - Size-specific effects of warming
m2_4bc_without <- read_rds(file.path(dir_output_without, "m2_4bc.rds"))
m2_7a_without <- read_rds(file.path(dir_output_without, "m2_7a.rds"))
m2_7fg_without <- read_rds(file.path(dir_output_without, "m2_7fg.rds"))
m2_7d_without <- read_rds(file.path(dir_output_without, "m2_7d.rds"))

m2_4bc_with <- read_rds(file.path(dir_output_with, "m2_4bc.rds"))
m2_7a_with <- read_rds(file.path(dir_output_with, "m2_7a.rds"))
m2_7fg_with <- read_rds(file.path(dir_output_with, "m2_7fg.rds"))
m2_7d_with <- read_rds(file.path(dir_output_with, "m2_7d.rds"))


# check difference ----
## m2 ----
# without
df_m2_4bc_without <- as.data.frame(summary(m2_4bc_without)$coefficients) 
df_m2_4bc_without <- df_m2_4bc_without %>%
  select(1:2) %>%
  mutate(group = "Without GSI > 30%",
         pop_name = "North Sea") %>%
  mutate(pars = row.names(df_m2_4bc_without))

df_m2_7a_without <- as.data.frame(summary(m2_7a_without)$coefficients) 
df_m2_7a_without <- df_m2_7a_without %>%
  select(1:2) %>%
  mutate(group = "Without GSI > 30%",
         pop_name = "Irish Sea") %>%
  mutate(pars = row.names(df_m2_7a_without))

df_m2_7fg_without <- as.data.frame(summary(m2_7fg_without)$coefficients) 
df_m2_7fg_without <- df_m2_7fg_without %>%
  select(1:2) %>%
  mutate(group = "Without GSI > 30%",
         pop_name = "Bristol Channel \n Celtic Sea North") %>%
  mutate(pars = row.names(df_m2_7fg_without))

df_m2_7d_without <- as.data.frame(summary(m2_7d_without)$coefficients) 
df_m2_7d_without <- df_m2_7d_without %>%
  select(1:2) %>%
  mutate(group = "Without GSI > 30%",
         pop_name = "Eastern English Channel") %>%
  mutate(pars = row.names(df_m2_7d_without))

# with
df_m2_4bc_with <- as.data.frame(summary(m2_4bc_with)$coefficients) 
df_m2_4bc_with <- df_m2_4bc_with %>%
  select(1:2) %>%
  mutate(group = "With GSI > 30%",
         pop_name = "North Sea") %>%
  mutate(pars = row.names(df_m2_4bc_with))

df_m2_7a_with <- as.data.frame(summary(m2_7a_with)$coefficients) 
df_m2_7a_with <- df_m2_7a_with %>%
  select(1:2) %>%
  mutate(group = "With GSI > 30%",
         pop_name = "Irish Sea") %>%
  mutate(pars = row.names(df_m2_7a_with))

df_m2_7fg_with <- as.data.frame(summary(m2_7fg_with)$coefficients) 
df_m2_7fg_with <- df_m2_7fg_with %>%
  select(1:2) %>%
  mutate(group = "With GSI > 30%",
         pop_name = "Bristol Channel \n Celtic Sea North") %>%
  mutate(pars = row.names(df_m2_7fg_with))

df_m2_7d_with <- as.data.frame(summary(m2_7d_with)$coefficients) 
df_m2_7d_with <- df_m2_7d_with %>%
  select(1:2) %>%
  mutate(group = "With GSI > 30%",
         pop_name = "Eastern English Channel") %>%
  mutate(pars = row.names(df_m2_7d_with))

# merge
df_m2 <- bind_rows(df_m2_4bc_without, df_m2_7a_without, df_m2_7fg_without, df_m2_7d_without,
                   df_m2_4bc_with, df_m2_7a_with, df_m2_7fg_with, df_m2_7d_with)

df_pars <- tibble(pars = c("log.length",
                           "log.age",
                           "log.K_rel",
                           "log.temp",
                           "log.length:log.temp"),
                  pars_name = c("ln(Total body length)",
                                "ln(Age)",
                                "ln(Relative condition factor)",
                                "ln(Autumn-winter temperature)",
                                "ln(Autumn-winter temperature) x ln(Total body length)"))

df_m2 <- df_m2 %>%
  left_join(df_pars) %>%
  filter(is.na(pars_name) == F)

df_m2 <- df_m2 %>%
  mutate(pop_name = factor(pop_name, 
                           levels = c("North Sea", "Irish Sea", "Bristol Channel \n Celtic Sea North", "Eastern English Channel")))

ggplot() +
  geom_point(data = df_m2 %>% filter(pars %in% c("log.length", "log.temp", "log.length:log.temp") ), 
             aes(x = `Estimate`, y = pars_name, color = group),
             position = position_dodge(width = 1)) +
  geom_linerange(data = df_m2 %>% filter(pars %in% c("log.length", "log.temp", "log.length:log.temp")), 
                 aes(x = `Estimate`, y = pars_name,
                                   xmin = `Estimate` - `Std. Error`,
                                   xmax = `Estimate` + `Std. Error`,
                                   color = group),
                 position = position_dodge(width = 1)) +
  facet_grid(~ pop_name) +
  labs(x = "Parameter estimate",
       y = NULL,
       color = "Dataset") +
  theme(legend.position = "bottom")

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS_gsi30_m2.png"),
       width =  17, height = 6, 
       units = "cm",
       dpi = 1200,
       scale = 1.5)

ggsave(last_plot(), file = file.path(dir_report, "figS_gsi30_m2.pdf"),
       device = cairo_pdf,
       width =  17, height = 12, 
       units = "cm")

