# BEFORE: 
#       2_analyze data
# AFTER:
#       WP1/report

# 1. SETUP --------------
library(tidyverse)
library(sf)
library(stars)
library(RColorBrewer)
library(patchwork)

# FUNCTION
# WORKING DIRECTORY
dir_report <- "./report" # not indicate dir_report 

# POPULATION NAME
# create df_pop
df_pop <- tibble(pop = c("4bc", "7a", "7fg", "7d"),
                 pop_name = factor(c("North Sea", "Irish Sea", "Bristol Channel, Celtic Sea North", "Eastern English Channel"),
                                   levels = c("North Sea", "Irish Sea", "Bristol Channel, Celtic Sea North", "Eastern English Channel")))

# THEME
theme_set(theme_bw())

# 2. LOAD DATA ------------
## 2.1. obs ----
### pre-processing obs ----
obs <- read_rds("./data/sol_gonads_2004_2023.rds")

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
         maturity_stage = factor(SPE.MaturityStageDescription, level = c("Immature", "Developing", "Spawning", "Spent")),
         maturity = if_else(SPE.MaturityStageDescription %in% c("Immature"), 0, 1),
         maturity_desc = if_else(SPE.MaturityStageDescription %in% c("Immature"), "Immature", "Mature"),
  ) %>%
  select(pop, year, cohort, week, age, length, body_weight, gonad_weight, gsi, maturity_stage, maturity, maturity_desc)

# data_strange
obs_strange <- obs_sub %>% 
  filter(gsi >= 100)

### filter weeks within gsi >= 50% peak gsi ----
# data
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

### remove very large/small values ----
quantile(obs_sub$gsi)
# upper and lower whisker
11.99+1.5*(11.99-5.85) # 21.2
5.85-1.5*(11.99-5.85) # -3.36

obs_sub <- obs_sub %>% mutate(id = row.names(obs_sub))
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
# View(data_obs %>% 
#        group_by(pop, year) %>%
#        summarize(n = n()) %>%
#        pivot_wider(names_from = pop, values_from = n) %>%
#        arrange(year)) 
# 

# sample size vs age
# View(data_obs2 %>%
#        group_by(pop, age) %>%
#        summarize(n = n()) %>%
#        arrange(age) %>%
#        pivot_wider(names_from = age, values_from = n))
# if only select age with at least 10 obs
# 4bc: 2-12
# 7a: 3-14
# 7d: 3-19
# 7fg: 3-16

data_obs <- data_obs %>%
  filter(age >= if_else(pop == "4bc", 2, 3), 
         age <= if_else(pop == "4bc", 12, 
                        if_else(pop == "7a", 14,
                                if_else(pop == "7d", 19, 16))))

table(data_obs$pop)

## 2.2. GIS data -----------

dir_gis <- "./data/admin"

# continents
continents <- read_sf(file.path(dir_gis, "esri_continent.gpkg"))

# countries
countries <- read_sf(file.path(dir_gis, "esri_countries.gpkg"))

# ices area 
ices_area <- read_sf(file.path(dir_gis, "ices_areas_sub_group_4abc_4326_new.gpkg")) 
ices_area <- ices_area %>% 
  filter(Area_27 %in% c("4abc", "7a", "7d", "7fg")) %>% 
  mutate(pop.name = if_else(Area_27 == "4abc", "North Sea",
                            if_else(Area_27 == "7a", "Irish Sea",
                                    if_else(Area_27 == "7d", "Eastern English Channel",
                                            "Celtic Sea")))) %>%
  st_simplify(dTolerance = 2000)

# distribution area - datras
dir_ices <- "./data/ices"
ices_datras <- read_sf(file.path(dir_ices, "hl_loc_4abc7adefg8ab.gpkg")) %>%
  filter(Area_27 %in% c("4abc", "7a", "7d", "7fg"))

# sbt - oras5 (2004-2022) 
dir_temp <- "./data/temp"
oras5 <- read_stars(file.path(dir_temp, 
                               "oras5.tif")) 
oras5_df <- oras5 %>% 
  st_crop(ices_datras) %>%
  slice(band, 553:780) %>%
  st_apply(1:2, mean) %>% 
  as.data.frame() %>%
  filter(is.na(mean) == F)

## 2.3. Temp autumn winter ----
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

# 3. PLOT DATA -------------

## main study area --------
 
# setup
# use quantile scale to aid visualization
## setup breaks and labels 
legend <- tibble(breaks = seq(0, 1, 0.1),
                 labels = quantile(oras5_df$mean, probs = seq(0, 1, 0.1))) %>%
  mutate(labels = round(labels, 1))
## create color scale
fill_colours = rev(brewer.pal(nrow(legend), "Spectral"))
## add quantile value to dataframe to plot
oras5_df$quant <- ecdf(oras5_df$mean)(oras5_df$mean)  

# pop name
sf_pop <- data.frame(name = c("North Sea", "Irish Sea", "Eastern English Channel", "Bristol Channel\nCeltic Sea North"),
                     x = c(3, -5, 0, -6.5),
                     y = c(57.75, 55.25, 49, 49.5)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) 

#x = c(3, -5, 0, -6.5),
#y = c(57.75, 55.25, 49.15, 49.8)) %>% 

# plot
st_bbox(ices_area)
ggplot() +
  geom_tile(data = oras5_df, aes(x = x, y = y, fill = quant)) +
  geom_sf(data = ices_datras, fill = NA, color = "black", linewidth = 0.5) +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.25) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 0.5) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_sf_text(data = sf_pop, aes(label = name), size = 6/.pt) +
  coord_sf(xlim = c(-9.5, 9.5), ylim = c(48.5, 62.5), expand = FALSE) +
  theme_bw() +
  scale_fill_gradientn(
    colours = fill_colours,
    breaks = slice(legend, c(1,3,5,7,9,11))$breaks,
    labels = sprintf("%.1f", slice(legend, c(1,3,5,7,9,11))$labels),
    limits = c(0,1),
    guide = guide_colorbar(title = "Annual temperature (°C) \n (2004-2022)", 
                           title.position = "top", #left
                           #title.vjust = 0.7,
                           #barwidth = 6,
                           title.hjust = 0.5)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme(legend.position = c(0.999, 0.03), #(0.99, 0.001)
        legend.justification = c(1, 0),
        legend.direction="horizontal",
        legend.key.height = unit(6, "pt"), #10
        legend.key.width = unit(12, "pt"), #15
        legend.title = element_text(size = 5),
        legend.text = element_text(size =4),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", color = NA)
        ) 
p1 <- last_plot()

## temperature trend ----
# over study period
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "4bc", year >= 2004))) #no
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "7a", year >= 2004))) #no
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "7fg", year >= 2004))) #no
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "7d", year >= 2004))) #yes

# over whole temp data range
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "4bc")))
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "7a")))
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "7fg")))
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "7fg", year >= 1970)))
summary(lm(c.temp ~ year, data_temp_aut_win %>% filter(pop == "7d")))

data_temp_aut_win <- data_temp_aut_win %>% 
  left_join(df_pop) %>%
  filter(pop %in% c("4bc", "7a", "7fg", "7d"))

# get 7fg only from 1970
data_temp_sub <- data_temp_aut_win %>% filter((pop != "7fg") | (pop == "7fg" & year >= 1970))

ggplot() +
  geom_rect(aes(xmin = 2004, xmax = 2022, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.2) +
  geom_line(data = data_temp_aut_win, 
            aes(x = year, y = temp, color = pop_name),
            linewidth = 0.5) +
  geom_smooth(data = data_temp_sub, aes(x = year, y = temp, color = pop_name), 
              method = "lm", se = F, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c")) +
  labs(x = "Year",
       y = "Autumn-winter temperature (°C)",
       color = "Population") +
  theme(#legend.position = "bottom",
        #legend.title.position = "top",
        #legend.title.align = 0.5,
        legend.position = c(0.25, 0.92),
        #legend.justification = c(1, 0),
        #legend.direction="horizontal",
        legend.key.height = unit(1, "pt"),
        legend.key.width = unit(7, "pt"),
        #legend.key.size = unit(3, "pt"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", color = NA)
  ) 
p2 <- last_plot()

## save plot -----
layout <- "
AB
"

p3 <- (p1 + p2) + 
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = "A")

ggsave(p3, file = file.path(dir_report, "fig1_sampling-site.png"),
       width = 17, height = 12, 
       units = "cm",  
       dpi = 1200)

ggsave(p3, file = file.path(dir_report, "fig1_sampling-site.pdf"),
       device = cairo_pdf,
       width = 17, height = 12, 
       units = "cm")

