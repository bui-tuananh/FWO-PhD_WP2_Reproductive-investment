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

## model data ----
data <- data_obs %>%
  left_join(data_sol) %>%
  left_join(data_temp)

# weigth model
m2_w <- read_rds("./output/model.weight_intrinsic.rds")
m3_w <- read_rds("./output/model.weight_extrinsic.rds")

# length model
m2_l <- read_rds("./output/model.length_intrinsic.rds")
m3_l <- read_rds("./output/model.length_extrinsic.rds")

# FIGURE ----
## fig1 - sampling site -----
source("./3_report_sampling-site.R")

## fig2 - intrinsic effect | body weight ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- m2_w
  
  data_sub <- data %>% filter(pop == pop_name)
  var_value <- sort(unique(data_sub$log.body_weight))
  
  pred_temp <- as.data.frame(Effect(c("log.body_weight", "pop"), 
                                    model, 
                                    xlevels = list("log.body_weight" = var_value))) %>%
    filter(pop == pop_name)
  
  pred <- bind_rows(pred, pred_temp)
  
}

# est_ci95
est_ci95 <- rbind(deltaMethod(m2_w, "log.body_weight"),
                  deltaMethod(m2_w, "log.body_weight + `log.body_weight:pop7a`"),
                  deltaMethod(m2_w, "log.body_weight + `log.body_weight:pop7fg`"),
                  deltaMethod(m2_w, "log.body_weight + `log.body_weight:pop7d`")
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

pred <- pred %>%
  mutate(label = factor(label, 
                        levels = c("β = 1.18 (1.13-1.22) - North Sea", 
                                   "β = 1.06 (1.02-1.10) - Irish Sea",
                                   "β = 1.16 (1.13-1.20) - Bristol Channel, Celtic Sea North",
                                   "β = 1.09 (1.06-1.11) - Eastern English Channel"
                        )))

### 4 isometric lines ----
df_label <- pred %>% 
  select(pop,
         label) %>%
  unique()

df_iso <- tibble(pop = c("4bc", "7a", "7fg", "7d"),
                 a = c(fixef(m2_w)["(Intercept)"],
                       fixef(m2_w)["(Intercept)"] + fixef(m2_w)["pop7a"],
                       fixef(m2_w)["(Intercept)"] + fixef(m2_w)["pop7fg"],
                       fixef(m2_w)["(Intercept)"] + fixef(m2_w)["pop7d"]),
                 b = 1,
                 b_c = fixef(m2_w)["condition"]
)

pred_iso <- data %>%
  left_join(df_iso) %>%
  left_join(df_label) %>%
  mutate(fit_iso = a + b*log.body_weight + b_c*mean(condition)) 

# plot 
ggplot() +
  geom_line(data = pred,
            aes(x = round(exp(log.body_weight)),
                y = exp(fit),
                color = label)) +
  geom_ribbon(data = pred,
              aes(x = round(exp(log.body_weight)),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = label),
              alpha = 0.3) +
  geom_line(data = pred_iso, 
            aes(x = round(exp(log.body_weight)),
                y = exp(fit_iso),
                color = label),
            linetype = "dashed") +
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p1 <- last_plot()

(p1) & 
  theme(legend.position = c(0.28, 0.88), 
        legend.background = element_rect(fill = "transparent"),
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

# save plot
ggsave(last_plot(), file = file.path(dir_report, "fig2_intrinsic-effect_body-weight.pdf"),
       device = cairo_pdf,
       width = 17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig2_intrinsic-effect_body-weight.png"),
       width = 17, height = 11,
       units = "cm",
       dpi = 1200)

### 1 isometric line ----
# isometric line - average intercept of all populations
a <- (fixef(m2_w)["(Intercept)"]*4 + fixef(m2_w)["pop7a"] + fixef(m2_w)["pop7d"] + fixef(m2_w)["pop7fg"])/4
b <- 1
b_c <- fixef(m2_w)["condition"]

pred_iso <- data %>%
  mutate(fit_iso = a + b*log.body_weight + b_c*mean(condition)) 

ggplot() +
  geom_line(data = pred,
            aes(x = round(exp(log.body_weight)),
                y = exp(fit),
                color = label)) +
  geom_ribbon(data = pred,
              aes(x = round(exp(log.body_weight)),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = label),
              alpha = 0.3) +
  geom_line(data = pred_iso, 
            aes(x = round(exp(log.body_weight)),
                y = exp(fit_iso)),
            color = "grey") +
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p1 <- last_plot()

(p1) & 
  theme(legend.position = c(0.28, 0.88), 
        legend.background = element_rect(fill = "transparent"),
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

# save plot
ggsave(last_plot(), file = file.path(dir_report, "fig2_intrinsic-effect_body-weight_1iso.pdf"),
       device = cairo_pdf,
       width = 17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig2_intrinsic-effect_body-weight_1iso.png"),
       width = 17, height = 11,
       units = "cm",
       dpi = 1200)

## fig3 - intrinsic effect | others -----
### condition ----
pred <- as.data.frame(Effect(c("condition"), 
                             model, 
                             xlevels = 100))

ggplot() +
  geom_line(data = pred, 
            aes(x = condition, 
                y = exp(fit))) +
  geom_ribbon(data = pred, 
              aes(x = condition, 
                  ymin = exp(lower), 
                  ymax = exp(upper)), 
              alpha = 0.3) +
  labs(x = "Fulton's condition factor",
       y = "Gonad weight (g)") 

p_cond <- last_plot()

### maturity_stage ----
pred <- as.data.frame(Effect(c("maturity_stage"), 
                             model, 
                             xlevels = 100))

ggplot() +
  geom_linerange(data = pred, 
                 aes(x = maturity_stage, 
                     y = exp(fit),
                     ymin = exp(lower), 
                     ymax = exp(upper))) +
  geom_point(data = pred, 
             aes(x = maturity_stage, 
                 y = exp(fit)),
             fill = "white") +
  labs(x = "Maturity stage",
       y = "Gonad weight (g)") 

p_mat <- last_plot()

### merge plot -----
(p_cond/p_mat) +
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 6.5),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.direction = "horizontal",
        legend.box="vertical",
        legend.title.position = "top",
        legend.title.align = 0.5
  ) &
  guides(color = guide_legend(nrow = 2))

ggsave(last_plot(), file = file.path(dir_report, "fig3_intrinsic-effect_others.pdf"),
       device = cairo_pdf,
       width =  8.5, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig3_intrinsic-effect_others.png"),
       width =  8.5, height = 11, 
       units = "cm",
       dpi = 1200)

## fig4 - random effect ----
### week ----
pred_week <- ranef(m2_w)$`pop.week` 
pred_week_se <- as.data.frame(arm::se.ranef(m2_w)$`pop.week`)

pred_week <- pred_week %>%
  mutate(pop = sub(":.*", "", rownames(pred_week)),
         week = as.numeric(sub(".*:", "", rownames(pred_week))),
         intercept = `(Intercept)`,
         intercept_se = pred_week_se$`(Intercept)`) %>%
  left_join(df_pop_2line)

ggplot(data = pred_week, aes(x = week, y = intercept)) +
  geom_point() +
  geom_linerange(data = pred_week, 
                 aes(ymin = intercept - 1.96*intercept_se, 
                     ymax = intercept + 1.96*intercept_se)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(pop_name ~ .) +
  labs(x = "Week",
       y = "Random effect")

p_week <- last_plot()

### year ----
pred_year <- ranef(m2_w)$`pop.year` 
pred_year_se <- as.data.frame(arm::se.ranef(m2_w)$`pop.year`)

pred_year <- pred_year %>%
  mutate(pop = sub(":.*", "", rownames(pred_year)),
         year = as.numeric(sub(".*:", "", rownames(pred_year))),
         intercept = `(Intercept)`,
         intercept_se = pred_year_se$`(Intercept)`) %>%
  left_join(df_pop_2line)

ggplot(data = pred_year, aes(x = year, y = intercept)) +
  geom_point() +
  geom_linerange(data = pred_year, 
                 aes(ymin = intercept - 1.96*intercept_se, 
                     ymax = intercept + 1.96*intercept_se)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(pop_name ~ .) +
  labs(x = "Year",
       y = "Random effect")

p_year <- last_plot()

### merge plot ----
(p_week | p_year) +
  #plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "bottom", 
        #plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.direction = "horizontal",
        #legend.box="vertical",
        #legend.title.position = "top",
        #legend.title.align = 0.5
  ) 

ggsave(last_plot(), file = file.path(dir_report, "fig4_random-effect.pdf"),
       device = cairo_pdf,
       width =  17, height = 17,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig4_random-effect.png"),
       width =  17, height = 17,
       units = "cm",
       dpi = 1200)

# fig 5 - extrinsic - temp*weight (0.5 deg) ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- m3_w
  
  data_sub <- data %>% filter(pop == pop_name)
  var_value <- sort(unique(data_sub$log.body_weight))
  var_value2 <- c(0,0.5)
  
  pred_temp <- as.data.frame(Effect(c("log.body_weight", "c.temp", "pop"), 
                                    model, 
                                    xlevels = list("log.body_weight" = var_value,
                                                   "c.temp" = var_value2 ))) %>%
    filter(pop == pop_name)
  
  pred_temp <- pred_temp %>% 
    mutate(pop = pop_name)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>%
  left_join(df_pop_2line)

# gonad_weight ~ body_weight by temp
pred <- pred %>% 
  mutate(temp_desc = if_else(c.temp == 0, "Average temperature (T0)", "Average temperature + 0.5°C (T0.5)"))

# est_ci95
est_ci95 <- rbind(deltaMethod(m3_w, "log.body_weight + `log.body_weight:c.temp`*0"),
                  deltaMethod(m3_w, "log.body_weight + `log.body_weight:c.temp`*0.5"),
                  deltaMethod(m3_w, "log.body_weight + `log.body_weight:pop7a` + `log.body_weight:c.temp`*0"),
                  deltaMethod(m3_w, "log.body_weight + `log.body_weight:pop7a` + `log.body_weight:c.temp`*0.5"),
                  deltaMethod(m3_w, "log.body_weight + `log.body_weight:pop7fg` + `log.body_weight:c.temp`*0"),
                  deltaMethod(m3_w, "log.body_weight + `log.body_weight:pop7fg` + `log.body_weight:c.temp`*0.5"),
                  deltaMethod(m3_w, "log.body_weight + `log.body_weight:pop7d` + `log.body_weight:c.temp`*0"),
                  deltaMethod(m3_w, "log.body_weight + `log.body_weight:pop7d` + `log.body_weight:c.temp`*0.5")) %>%
  as.data.frame() %>%
  mutate(pop = rep(c("4bc", "7a", "7fg", "7d"), each = 2),
         temp_range = rep(c("T0", "T0.5"), 4)) %>% 
  mutate(est_ci95 = sprintf("%.2f (%.2f-%.2f)", Estimate, `2.5 %`, `97.5 %`)) %>%
  left_join(df_pop_2line) %>%
  select(pop, pop_name, temp_range, est_ci95) %>%
  mutate(label = paste0("β", temp_range, " = ", est_ci95)) %>%
  mutate(label2 = paste0(" = ", est_ci95))

ggplot() +
  geom_line(data = pred,
            aes(x = exp(log.body_weight), 
                y = exp(fit), 
                color = temp_desc)) +
  geom_ribbon(data = pred,
              aes(x = exp(log.body_weight), 
                  ymin = exp(lower), 
                  ymax = exp(upper),
                  fill = temp_desc),
              alpha = 0.3) +
  geom_text(data = est_ci95 %>% filter(temp_range == "T0"), 
            aes(x = 50, y = 240, label = "beta[T0]"), size = 2, hjust = "left", parse = T) +
  geom_text(data = est_ci95 %>% filter(temp_range == "T0"), 
            aes(x = 250, y = 240, label = label2), size = 2, hjust = "left") + 
  geom_text(data = est_ci95 %>% filter(temp_range == "T0.5"), 
            aes(x = 50, y = 225, label = "beta[T0.5]"), size = 2, hjust = "left", parse = T) + 
  geom_text(data = est_ci95 %>% filter(temp_range == "T0.5"), 
            aes(x = 250, y = 225, label = label2), size = 2, hjust = "left") + 
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)",
       color = NULL,
       fill = NULL) +
  scale_color_manual(values = c("#00bfc4", "#f8766d")) +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  facet_grid(. ~ pop_name)  +
  theme(legend.position = "bottom", 
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.direction = "horizontal",
  ) 
p_main <- last_plot()


p_hist <- ggplot() +
  geom_histogram(data = data %>% left_join(df_pop_2line), aes(x = body_weight), bins = 200) +
  facet_grid(~ pop_name) +
  theme_void() +
  theme(strip.text = element_blank())

p1 <- p_hist/p_main

# save file
(p1) + plot_layout(nrow = 2, heights = c(1, 3))

ggsave(last_plot(), file = file.path(dir_report, "fig5_extrinsic-effect_temp.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig5_extrinsic-effect_temp.png"),
       width =  17, height = 11,
       units = "cm",
       dpi = 1200)

## figS - zoom in -----
ggplot() +
  geom_line(data = pred,
            aes(x = exp(log.body_weight), 
                y = exp(fit), 
                color = temp_desc)) +
  geom_ribbon(data = pred,
              aes(x = exp(log.body_weight), 
                  ymin = exp(lower), 
                  ymax = exp(upper),
                  fill = temp_desc),
              alpha = 0.3) +
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)",
       color = NULL,
       fill = NULL) +
  scale_color_manual(values = c("#00bfc4", "#f8766d")) +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  facet_grid(. ~ pop_name) +
  coord_cartesian(xlim = c(0, 500), ylim = c(0, 50))
p1 <- last_plot()
(p1) +  theme(legend.position = "bottom", 
              axis.text = element_text(size = 7),
              axis.title = element_text(size = 9),
              legend.key.size = unit(10, "pt"),
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 9),
              legend.direction = "horizontal")

ggsave(last_plot(), file = file.path(dir_report, "figS_extrinsic-effect_temp_zoomin.pdf"),
       device = cairo_pdf,
       width =  17, height = 9,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_extrinsic-effect_temp_zoomin.png"),
       width =  17, height = 9,
       units = "cm",
       dpi = 1200)

## summary size specific -----
pred_sum <- tibble(body_weight = seq(100, 2000, 50)) %>%
  mutate(diff = exp(fixef(m3)["c.temp"]*0.5 + fixef(m3)["log.body_weight:c.temp"]*0.5*log(body_weight))*100 - 100) %>%
  mutate(diff_name = sprintf("%.2f", diff))

# TABLE ----
## table1 - sample size ----
data_sum <- data %>% 
  group_by(pop, year) %>% 
  summarize(n = n()) %>% 
  arrange(year) %>%
  pivot_wider(names_from = pop, values_from = n) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(year, `4bc`, `7a`, `7fg`, `7d`)

data %>% group_by(pop) %>% summarize(n = n())

col.header <- c("Year", 
                "North Sea (age 2-12, n = 991)", 
                "Irish Sea (age 3-14, n = 1677)", 
                "Bristol Channel and Celtic Sea North (age 3-16, n = 2066)",
                "Eastern English Channel (age 3-19, n = 3343)")

tab_df(data_sum,
       col.header = col.header,
       file = file.path("./report", "table1_analysed-data-summary.html"))

## table2 - analysis summary ----
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
  "ln(Body weight)",
  "Condition",
  "Maturity stage (Spawning)",
  "Population (Irish Sea)",
  "Population (Eastern English Channel)",
  "Population (Bristol Channel and Celtic Sea North)",
  "ln(Body weight) x Population (Irish Sea)",
  "ln(Body weight) x Population (Eastern English Channel)",
  "ln(Body weight) x Population (Bristol Channel and Celtic Sea North)",
  "Temperature",
  "ln(Body weight) x Temperature"
)

# summary table
tab_model(m2_w, m3_w, 
          show.se = NULL, 
          #show.ci = TRUE, 
          collapse.ci = TRUE,
          show.p = TRUE,
          p.style = "stars", 
          #dv.labels = NULL,
          dv.labels = c("Intrinsic model", 
                        "Extrinsic model"),
          pred.labels = pred_labels,
          order.terms = c(1,2,3,4,5,7,6,8,10,9,11,12),
          string.pred = "Fixed Effects",
          string.est = "Estimate (95CI)",
          string.ci = "95% CI",
          #collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "table2_model-summary.html")
)

# SUPPLEMENTARY INFORMATION -----
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
ggsave(last_plot(), file = file.path(dir_report, "figS_outlier.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_outlier.png"),
       width =  17, height = 11,
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
  scale_fill_manual(values = c("grey", "black"))

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
  )

ggsave(last_plot(), file = file.path(dir_report, "figS_distribution-gsi.pdf"),
       device = cairo_pdf,
       width = 17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_distribution-gsi.png"),
       width = 17, height = 11,
       units = "cm",
       dpi = 1200)

#################################################### SI ###
## fig 
## figS - temp ----
data_temp <- data_temp %>% 
  left_join(df_pop) %>%
  filter(is.na(pop_name) == F)

# absolute temp
ggplot(data = data_temp, aes(x = year, y = temp, color = pop_name)) +
  geom_line() +
  labs(x = "Year",
       y = "Winter temperature (°C)",
       color = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p1 <- last_plot()

# c.temp
ggplot(data = data_temp, aes(x = year, y = c.temp, color = pop_name)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year",
       y = "Winter temperature anomaly (°C)",
       color = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p2 <- last_plot()

(p1 / p2) + 
  plot_layout(axis_title = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.title.align = 0.5) &
  guides(color = guide_legend(nrow = 2))

ggsave(last_plot(), file = file.path(dir_report, "figS_temp.pdf"),
       device = cairo_pdf,
       width =  10, height = 20,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_temp.png"),
       width =  10, height = 20,
       units = "cm",
       dpi = 1200)

## figS - f and ssb ----
data_sol <- data_sol %>%
  left_join(df_pop) %>%
  filter(is.na(pop_name) == F)

# f
ggplot(data = data_sol, aes(x = year, y = f, color = pop_name)) +
  geom_line() +
  labs(x = "Year",
       y = "Fishing mortality\nage 3-7",
       color = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p1 <- last_plot()

# ssb
ggplot(data = data_sol, aes(x = year, y = ssb/1000, color = pop_name)) +
  geom_line() +
  labs(x = "Year",
       y = "Spawning Stock Biomass\n(1000 tonne)",
       color = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p2 <- last_plot()

ggplot(data = data_sol, aes(x = year, y = ssb.i, color = pop_name)) +
  geom_line() +
  labs(x = "Year",
       y = "Spawning Stock Biomass\n(tonne/km²)",
       color = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p3 <- last_plot()


(p1 / p2 / p3) + 
  plot_layout(axis_title = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top",
        legend.title.align = 0.5) &
  guides(color = guide_legend(nrow = 2))

ggsave(last_plot(), file = file.path(dir_report, "figS_f_ssb.pdf"),
       device = cairo_pdf,
       width =  10, height = 20,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_f_ssb.png"),
       width =  10, height = 20,
       units = "cm",
       dpi = 1200)

## tableS - collinearity ----
source("./ref/HighstatLib.r") # support
MyVar <- c("log.body_weight", "log.age", "condition", "c.temp", "f", "ssb.i")
corvif(data[,MyVar]) 

corvif_full <- corvif(data[,MyVar]) 
corvif_full <- corvif_full %>% 
  mutate(var = row.names(corvif_full),
         VIF_full = sprintf("%.2f", GVIF)) %>%
  select(var, VIF_full)

MyVar <- c("log.body_weight", "condition", "c.temp", "f", "ssb.i")
corvif(data[,MyVar]) 
corvif_sub <- corvif(data[,MyVar]) 
corvif_sub <- corvif_sub %>% 
  mutate(var = row.names(corvif_sub),
         VIF_sub = sprintf("%.2f", GVIF)) %>%
  select(var, VIF_sub)

corvif <- corvif_full %>% 
  left_join(corvif_sub, by = join_by(var)) 

# save table
col.header <- c("Variable",
                "VIF (collinearity)",
                "VIF (collinearity removed)")

# table S4
tab_df(corvif,
       col.header = col.header,
       file = file.path("./report", "tableS_collinearity.html"))

# table S5 get manually from here
MyVar <- c("log.body_weight", "condition", "c.temp", "f", "ssb.i")
corvif(data[,MyVar]) 


## figS - gonad weight ~ body weight -----
data <- data_obs %>%
  left_join(df_pop_2line)

ggplot(data = data, aes(x = body_weight, y = gonad_weight)) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  facet_grid(. ~ pop_name) +
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)") +
  xlim(0, 2200)
p1 <- last_plot()

ggplot(data = data, aes(x = log(body_weight), y = log(gonad_weight) )) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  facet_grid(. ~ pop_name) +
  labs(x = "ln(Body weight) (g)",
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
### intrinsic model ----
list_model <- list("Intrinsic model (body weight)" = m2_w, 
                   "Intrinsic model (body length)" = m2_l)
pred <- tibble()
for(i in 1:2) {
  
  model <- list_model[[1]]
  data <- model@frame
  
  pred_temp <- data %>%
    mutate(res = resid(model),
           pred = predict(model, type = "response"))
  
  pred_temp <- pred_temp %>% 
    mutate(model = names(list_model[i]))
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>%
  mutate(model = factor(model, levels = c("Intrinsic model (body weight)", "Intrinsic model (body length)")))

ggplot(data = pred, aes(x = pred, y = res)) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x  = "Fitted values",
       y = "Residuals") +
  facet_grid(. ~ model)
p1 <- last_plot()

ggplot(data = pred, aes(sample = res)) +
  stat_qq(alpha = 0.3, size = rel(0.8)) +
  stat_qq_line() +
  labs(x  = "Theoretical quantiles",
       y = "Sample quantiles") +
  facet_grid(. ~ model)
p2 <- last_plot()

# save file
(p1/p2) & 
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

ggsave(last_plot(), file = file.path(dir_report, "figS_intrinsic-model-diagnostic.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_intrinsic-model-diagnostic.png"),
       width =  17, height = 11,
       units = "cm",
       dpi = 1200)

### extrinsic model ----
list_model <- list("Extrinsic model (body weight)" = m3_w, 
                   "Extrinsic model (body length)" = m3_l)
pred <- tibble()
for(i in 1:2) {
  
  model <- list_model[[1]]
  data <- model@frame
  
  pred_temp <- data %>%
    mutate(res = resid(model),
           pred = predict(model, type = "response"))
  
  pred_temp <- pred_temp %>% 
    mutate(model = names(list_model[i]))
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>%
  mutate(model = factor(model, levels = c("Extrinsic model (body weight)", "Extrinsic model (body length)")))

ggplot(data = pred, aes(x = pred, y = res)) +
  geom_point(alpha = 0.3, size = rel(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x  = "Fitted values",
       y = "Residuals") +
  facet_grid(. ~ model)
p1 <- last_plot()

ggplot(data = pred, aes(sample = res)) +
  stat_qq(alpha = 0.3, size = rel(0.8)) +
  stat_qq_line() +
  labs(x  = "Theoretical quantiles",
       y = "Sample quantiles") +
  facet_grid(. ~ model)
p2 <- last_plot()

# save file
(p1/p2) & 
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

ggsave(last_plot(), file = file.path(dir_report, "figS_extrinsic-model-diagnostic.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_extrinsic-model-diagnostic.png"),
       width =  17, height = 11,
       units = "cm",
       dpi = 1200)

## figS - intrinsic effect | length ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- m2_l
  
  data_sub <- data %>% filter(pop == pop_name)
  var_value <- sort(unique(data_sub$log.length))
  
  pred_temp <- as.data.frame(Effect(c("log.length", "pop"), 
                                    model, 
                                    xlevels = list("log.length" = var_value))) %>%
    filter(pop == pop_name)
  
  pred <- bind_rows(pred, pred_temp)
  
}

# est_ci95
est_ci95 <- rbind(deltaMethod(m2_l, "log.length"),
                  deltaMethod(m2_l, "log.length + `log.length:pop7a`"),
                  deltaMethod(m2_l, "log.length + `log.length:pop7fg`"),
                  deltaMethod(m2_l, "log.length + `log.length:pop7d`")
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

pred <- pred %>%
  mutate(label = factor(label, 
                        levels = c("β = 3.56 (3.42-3.71) - North Sea", 
                                   "β = 3.16 (3.04-3.28) - Irish Sea",
                                   "β = 3.51 (3.39-3.63) - Bristol Channel, Celtic Sea North",
                                   "β = 3.25 (3.16-3.33) - Eastern English Channel"
                        )))

### 4 isometric lines ----
df_label <- pred %>% 
  select(pop,
         label) %>%
  unique()

df_iso <- tibble(pop = c("4bc", "7a", "7fg", "7d"),
                 a = c(fixef(m2_l)["(Intercept)"],
                       fixef(m2_l)["(Intercept)"] + fixef(m2_l)["pop7a"],
                       fixef(m2_l)["(Intercept)"] + fixef(m2_l)["pop7fg"],
                       fixef(m2_l)["(Intercept)"] + fixef(m2_l)["pop7d"]),
                 b = 3,
                 b_c = fixef(m2_l)["condition"]
)

pred_iso <- data %>%
  left_join(df_iso) %>%
  left_join(df_label) %>%
  mutate(fit_iso = a + b*log.length + b_c*mean(condition)) 

# plot
ggplot() +
  geom_line(data = pred, 
            aes(x = round(exp(log.length)), 
                y = exp(fit), 
                color = label)) +
  geom_ribbon(data = pred, 
              aes(x = round(exp(log.length)), 
                  ymin = exp(lower), 
                  ymax = exp(upper),
                  fill = label), 
              alpha = 0.3) +
  geom_line(data = pred_iso, 
            aes(x = round(exp(log.length)),
                y = exp(fit_iso),
                color = label),
            linetype = "dashed") +
  labs(x = "Body length (mm)",
       y = "Gonad weight (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p1 <- last_plot()

(p1) & 
  theme(legend.position = c(0.28, 0.88), 
        legend.background = element_rect(fill = "transparent"),
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS_intrinsic-effect_body-length.pdf"),
       device = cairo_pdf,
       width = 17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_intrinsic-effect_body-length.png"),
       width = 17, height = 11,
       units = "cm",
       dpi = 1200)

### 1 isometric line ----
# average intercept 
a <- (fixef(m2_l)["(Intercept)"]*4 + fixef(m2_l)["pop7a"] + fixef(m2_l)["pop7d"] + fixef(m2_l)["pop7fg"])/4
b <- 3
b_c <- fixef(m2_l)["condition"]

pred_iso <- data %>%
  mutate(fit_iso = a + b*log.length + b_c*mean(condition))

ggplot() +
  geom_line(data = pred, 
            aes(x = round(exp(log.length)), 
                y = exp(fit), 
                color = label)) +
  geom_ribbon(data = pred, 
              aes(x = round(exp(log.length)), 
                  ymin = exp(lower), 
                  ymax = exp(upper),
                  fill = label), 
              alpha = 0.3) +
  geom_line(data = pred_iso, 
            aes(x = round(exp(log.length)),
                y = exp(fit_iso)),
            color = "grey") +
  labs(x = "Body length (mm)",
       y = "Gonad weight (g)",
       color = "Population", 
       fill = "Population") +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"))
p1 <- last_plot()

(p1) & 
  theme(legend.position = c(0.28, 0.88), 
        legend.background = element_rect(fill = "transparent"),
        plot.tag.position  = c(0.95, 0.9),
        #axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS_intrinsic-effect_body-length_1iso.pdf"),
       device = cairo_pdf,
       width = 17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_intrinsic-effect_body-length_1iso.png"),
       width = 17, height = 11,
       units = "cm",
       dpi = 1200)

## figS - extrinsic - temp*length (0.5 deg) ----

pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- m3_l
  
  data_sub <- data %>% filter(pop == pop_name)
  var_value <- sort(unique(data_sub$log.length))
  var_value2 <- c(0,0.5)
  
  pred_temp <- as.data.frame(Effect(c("log.length", "c.temp", "pop"), 
                                    model, 
                                    xlevels = list("log.length" = var_value,
                                                   "c.temp" = var_value2 ))) %>%
    filter(pop == pop_name)
  
  pred_temp <- pred_temp %>% 
    mutate(pop = pop_name)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>%
  left_join(df_pop_2line)

# gonad_weight ~ length by temp
pred <- pred %>% 
  mutate(temp_desc = if_else(c.temp == 0, "Average temperature (T0)", "Average temperature + 0.5°C (T0.5)"))

# est_ci95
est_ci95 <- rbind(deltaMethod(m3_l, "log.length + `log.length:c.temp`*0"),
                  deltaMethod(m3_l, "log.length + `log.length:c.temp`*0.5"),
                  deltaMethod(m3_l, "log.length + `log.length:pop7a` + `log.length:c.temp`*0"),
                  deltaMethod(m3_l, "log.length + `log.length:pop7a` + `log.length:c.temp`*0.5"),
                  deltaMethod(m3_l, "log.length + `log.length:pop7fg` + `log.length:c.temp`*0"),
                  deltaMethod(m3_l, "log.length + `log.length:pop7fg` + `log.length:c.temp`*0.5"),
                  deltaMethod(m3_l, "log.length + `log.length:pop7d` + `log.length:c.temp`*0"),
                  deltaMethod(m3_l, "log.length + `log.length:pop7d` + `log.length:c.temp`*0.5")) %>%
  as.data.frame() %>%
  mutate(pop = rep(c("4bc", "7a", "7fg", "7d"), each = 2),
         temp_range = rep(c("T0", "T0.5"), 4)) %>% 
  mutate(est_ci95 = sprintf("%.2f (%.2f-%.2f)", Estimate, `2.5 %`, `97.5 %`)) %>%
  left_join(df_pop_2line) %>%
  select(pop, pop_name, temp_range, est_ci95) %>%
  mutate(label = paste0("β", temp_range, " = ", est_ci95)) %>%
  mutate(label2 = paste0(" = ", est_ci95))

ggplot() +
  geom_line(data = pred,
            aes(x = exp(log.length), 
                y = exp(fit), 
                color = temp_desc)) +
  geom_ribbon(data = pred,
              aes(x = exp(log.length), 
                  ymin = exp(lower), 
                  ymax = exp(upper),
                  fill = temp_desc),
              alpha = 0.5) +
  geom_text(data = est_ci95 %>% filter(temp_range == "T0"), 
            aes(x = 50, y = 150, label = "beta[T0]"), size = 2, hjust = "left", parse = T) +
  geom_text(data = est_ci95 %>% filter(temp_range == "T0"), 
            aes(x = 100, y = 150, label = label2), size = 2, hjust = "left") + 
  geom_text(data = est_ci95 %>% filter(temp_range == "T0.5"), 
            aes(x = 50, y = 135, label = "beta[T0.5]"), size = 2, hjust = "left", parse = T) + 
  geom_text(data = est_ci95 %>% filter(temp_range == "T0.5"), 
            aes(x = 100, y = 135, label = label2), size = 2, hjust = "left") + 
  labs(x = "Body length (mm)",
       y = "Gonad weight (g)",
       color = NULL,
       fill = NULL) +
  scale_color_manual(values = c("#f8766d", "#00bfc4")) +
  scale_fill_manual(values = c("#f8766d", "#00bfc4")) +
  facet_grid(. ~ pop_name) +
  scale_color_manual(values = c("#00bfc4", "#f8766d")) +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  facet_grid(. ~ pop_name)  +
  theme(legend.position = "bottom", 
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.direction = "horizontal",
  ) 
p_main <- last_plot()

p_hist <- ggplot() +
  geom_histogram(data = data %>% left_join(df_pop_2line), aes(x = length), bins = 30) +
  facet_grid(~ pop_name) +
  theme_void() +
  theme(strip.text = element_blank())

p1 <- p_hist/p_main

# save file
(p1) + plot_layout(nrow = 2, heights = c(1, 3))

ggsave(last_plot(), file = file.path(dir_report, "figS_extrinsic-effect_temp_length.pdf"),
       device = cairo_pdf,
       width =  17, height = 11,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS_extrinsic-effect_temp_length.png"),
       width =  17, height = 11,
       units = "cm",
       dpi = 1200)