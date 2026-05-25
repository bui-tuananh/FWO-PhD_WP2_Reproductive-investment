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
