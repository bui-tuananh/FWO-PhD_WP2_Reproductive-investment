## temperature effect 1oC ----
pred <- tibble()
for(pop_name in unique(data$pop)) {
  
  print(paste("processing", pop_name))
  
  model <- m3_w
  
  data_sub <- data %>% filter(pop == pop_name)
  var_value <- sort(unique(data_sub$log.body_weight))
  var_value2 <- c(0,1)
  
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
  mutate(temp_desc = if_else(c.temp == 0, "Average temperature (T0)", "Average temperature + 1°C (T1)"))

## plot full no line ----
ggplot() +
  geom_histogram(data = data %>% left_join(df_pop_2line) %>% 
                   filter(pop == "4bc", body_weight >= 100, body_weight <= 1000), 
                 aes(x = body_weight), bins = 200) +
  facet_grid(~ pop_name) +
  theme_void() +
  theme(strip.text = element_blank()) 
p_hist <- last_plot()

ggplot(data = pred %>% filter(pop == "4bc", exp(log.body_weight) >= 100, exp(log.body_weight) <= 1000)) +
  geom_line(aes(x = exp(log.body_weight), 
                y = exp(fit), 
                color = temp_desc),
            linewidth = 1) +
  geom_ribbon(aes(x = exp(log.body_weight),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = temp_desc),
              alpha = 0.05) +
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)",
       color = NULL,
       fill = NULL) +
  scale_color_manual(values = c("#00bfc4", "#f8766d")) +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  facet_grid(. ~ pop_name)  +
  theme(legend.position = "none") +
  #geom_vline(xintercept = 450, linetype = "dashed") +
  scale_x_continuous(breaks = c(100, 200, 300, 450, 600, 700, 800, 900, 1000) )
p_main <- last_plot()

(p_hist/p_main) + plot_layout(nrow = 2, heights = c(1, 3))
ggsave(last_plot(), file = file.path(dir_presentation, "vmsd_extrinsic-effect_temp_full_noline.png"),
       width =  11, height = 13,
       units = "cm",
       dpi = 1200)

## plot full ----
ggplot() +
  geom_histogram(data = data %>% left_join(df_pop_2line) %>% 
                   filter(pop == "4bc", body_weight >= 100, body_weight <= 1000), 
                 aes(x = body_weight), bins = 200) +
  facet_grid(~ pop_name) +
  theme_void() +
  theme(strip.text = element_blank()) 
p_hist <- last_plot()

ggplot(data = pred %>% filter(pop == "4bc", exp(log.body_weight) >= 100, exp(log.body_weight) <= 1000)) +
  geom_line(aes(x = exp(log.body_weight), 
                y = exp(fit), 
                color = temp_desc),
            linewidth = 1) +
  geom_ribbon(aes(x = exp(log.body_weight),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = temp_desc),
              alpha = 0.05) +
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)",
       color = NULL,
       fill = NULL) +
  scale_color_manual(values = c("#00bfc4", "#f8766d")) +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  facet_grid(. ~ pop_name)  +
  theme(legend.position = "none") +
  geom_vline(xintercept = 450, linetype = "dashed") +
  scale_x_continuous(breaks = c(100, 200, 300, 450, 600, 700, 800, 900, 1000) )
p_main <- last_plot()

(p_hist/p_main) + plot_layout(nrow = 2, heights = c(1, 3))
ggsave(last_plot(), file = file.path(dir_presentation, "vmsd_extrinsic-effect_temp_full.png"),
       width =  11, height = 13,
       units = "cm",
       dpi = 1200)

### plot small fish ----
ggplot(data = pred %>% filter(pop == "4bc", exp(log.body_weight) >= 100, exp(log.body_weight) < 450)) +
  geom_line(aes(x = exp(log.body_weight), 
                y = exp(fit), 
                color = temp_desc),
            linewidth = 1) +
  geom_ribbon(aes(x = exp(log.body_weight),
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = temp_desc),
              alpha = 0.05) +
  labs(x = "Body weight (g)",
       y = "Gonad weight (g)",
       color = NULL,
       fill = NULL) +
  scale_color_manual(values = c("#00bfc4", "#f8766d")) +
  scale_fill_manual(values = c("#00bfc4", "#f8766d")) +
  facet_grid(. ~ pop_name)  +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(100, 200, 300, 450) )

ggsave(last_plot(), file = file.path(dir_presentation, "vmsd_extrinsic-effect_temp_small.png"),
       width =  11, height = 11,
       units = "cm",
       dpi = 1200)

### diff ----
# pred_sum <- tibble(body_weight = seq(100, 2000, 50)) %>%
#   mutate(diff = exp(fixef(m3_w)["c.temp"]*0.5 + fixef(m3_w)["log.body_weight:c.temp"]*0.5*log(body_weight))*100 - 100) %>%
#   mutate(diff_name = sprintf("%.2f", diff))

pred_sum <- tibble(body_weight = seq(100, 2000, 50)) %>%
  mutate(diff = exp(fixef(m3_w)["c.temp"]*1 + fixef(m3_w)["log.body_weight:c.temp"]*1*log(body_weight))*100 - 100) %>%
  mutate(diff_name = sprintf("%.2f", diff)) %>%
  mutate(temp = "t0")
base <- tibble(body_weight = seq(100, 2000, 50),
               diff = 0,
               temp = "t1")

pred_sum <- bind_rows(pred_sum, base)

ggplot() +
  geom_linerange(data = pred_sum %>% filter(body_weight <= 1000), aes(x = body_weight, y = diff, ymin = 0, ymax = diff,  group = body_weight)) + 
  geom_point(data = pred_sum %>% filter(body_weight <= 1000), aes(x = body_weight, y = diff, color = temp)) + 
  #geom_hline(yintercept = 0) +
  #scale_x_continuous(breaks = seq(100, 1000, 50)) +
  scale_x_continuous(breaks = c(100, 200, 300, 450, 600, 700, 800, 900, 1000) ) +
  geom_vline(xintercept = 450, linetype = "dashed") +
  labs(x = "Weight (g)",
       y = "Relative change of gonad weight (%)") #+
#theme(legend.position = "none")

ggsave(last_plot(), file = file.path(dir_presentation, "vmsd_extrinsic-effect_temp_diff.png"),
       width =  11, height = 11,
       units = "cm",
       dpi = 1200)

### size range ----
ggplot(data = data %>% filter(pop %in% c("4bc") ), aes(x = body_weight, y = length, color = pop)) + 
  geom_point() + 
  scale_x_continuous(breaks = seq(100,1000,100)) +
  geom_hline(yintercept = 240) +
  geom_vline(xintercept = 450) +
  scale_y_continuous(breaks = seq(100, 600, 50))
