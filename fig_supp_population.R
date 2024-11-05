df_plot <- df_plot |>
  mutate(population = fct_collapse(population, 
                                   "ANC attendees" = c("ANC attendees", "ANC and FP attendees",
                                                       "ANC and GYN attendees"),
                                   "FP attendees" = c("FP attendees", "FP and GYN attendees")),
         population = factor(population,
                             levels = c("ANC attendees", "FP attendees", 
                                        "GYN attendees", "PHC/OPD attendees",
                                        "Students","Community members", 
                                        "HIV/STI prevention trial participants",
                                        "Population-representative survey participants")))


df2 |>
  ggplot() +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = sex), alpha = 0.15) +
  geom_line(aes(x = year, y = prev, colour = sex, linetype = sex), size = 0.6) +
  scale_fill_manual(values = colour_2) +
  scale_colour_manual(values = colour_2) +
  guides(colour=guide_legend(order=1, title = "")) +
  # lines on new colour scale
  ggnewscale::new_scale_color() +
  geom_linerange(data = df_plot, aes(x = year_mid, ymin = adj_prev_lwr, ymax = adj_prev_upr, colour = population),
                 size = 0.4) + #show.legend = FALSE, 
  geom_point(data = df_plot, aes(x = year_mid, y = adj_prev, colour = population, shape = sex), size = 1, fill = "white") +
  #scale_colour_manual(values = colour_3) +
  scale_shape_manual(values = c(19, 21)) +
  guides(colour=guide_legend(order=2, title = "", nrow=2),
         fill=guide_legend(order=1, title = "", nrow=1),
         shape=guide_legend(order=1, title = "", nrow=1),
         linetype=guide_legend(order=1, title = "", nrow=1)) +
  ggh4x::facet_wrap2(dplyr::vars(sti,region),
                     strip = ggh4x::strip_nested(text_x = ggh4x::elem_list_text(face = c("bold", "plain"),
                                                                                size = c(7*1.3,7*1.18)), # align with theme
                                                 by_layer_x = TRUE),
                     scales = "free_y",
                     nrow = 3,
                     axes = "x") +
  ggh4x::facetted_pos_scales(
    y = list(
      sti == "Chlamydia" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = scales::label_percent()),
      sti == "Gonorrhoea" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.3), labels = scales::label_percent()),
      sti == "Trichomoniasis" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = scales::label_percent()),
      sti == "Chlamydia" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = NULL),
      sti == "Gonorrhoea" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.3), labels = NULL),
      sti == "Trichomoniasis" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = NULL))) +
  mytheme +
  theme(panel.spacing.y = unit(0.8, "lines"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        legend.box = "vertical",
        legend.margin = margin(unit(c(t=0,r=0,b=0,l=0), "cm"))) +
  labs(x = "", y = "", tag = "B")

ggsave("./plots/fig_2_population_v2.png", width = 18, height = 20, unit = "cm", dpi = 700)  


# compare pops

df_gen <- read.csv("./results/prevalence_alldata_adjusted_general.csv") |>
  mutate(population = "Population-representative survey participants")

df_preventtrial <- read.csv("./results/prevalence_alldata_adjusted_preventtrial.csv") |>
  mutate(population = "HIV/STI prevention trial participants")

df_phcopd <- read.csv("./results/prevalence_alldata_adjusted_phcopd.csv") |>
  mutate(population = "PHC/OPD attendees")

df_anc <- read.csv("./results/prevalence_alldata_adjusted.csv") |>
  mutate(population = "ANC attendees")

df_models <- rbind(df_gen, df_preventtrial, df_phcopd, df_anc)  |>
  mutate(region = factor(region, levels = c("WCA", "EA", "SA"),
                         labels = c("Western and Central Africa", "Eastern Africa", "Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"),
                      labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         sex = fct_relevel(sex, "Female"),
         population = factor(population,
                             levels = c("ANC attendees", "PHC/OPD attendees",
                                        "HIV/STI prevention trial participants",
                                        "Population-representative survey participants","Other")))

df_plot <- df_plot |>
  mutate(population = fct_collapse(population, 
                                   "ANC attendees" = c("ANC attendees", "ANC and FP attendees",
                                                       "ANC and GYN attendees"),
                                   "Other" = c("FP attendees", "FP and GYN attendees","GYN attendees","Students","Community members")),
         population = factor(population,
                             levels = c("ANC attendees", "PHC/OPD attendees",
                                        "HIV/STI prevention trial participants",
                                        "Population-representative survey participants","Other")))

colour_3 <- c("#af4f2f","#df8d71","#1e5a46","#75884b","#5b859e")

df_models |>
  filter(sex == "Female") |>
  ggplot() +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = population), alpha = 0.15) +
  geom_line(aes(x = year, y = prev, colour = population), size = 0.6) +
  scale_fill_manual(values = colour_3) +
  scale_colour_manual(values = colour_3) +
  guides(colour=guide_legend(order=1, title = "")) +
  # lines on new colour scale
  ggnewscale::new_scale_color() +
  geom_linerange(data = df_plot, aes(x = year_mid, ymin = adj_prev_lwr, ymax = adj_prev_upr, colour = population),
                 size = 0.4) + #show.legend = FALSE, 
  geom_point(data = df_plot, aes(x = year_mid, y = adj_prev, colour = population), size = 1, fill = "white") +
  scale_colour_manual(values = colour_3) +
  scale_shape_manual(values = c(19, 21)) +
  guides(colour=guide_legend(order=2, title = "", nrow=2),
         fill=guide_legend(order=1, title = "", nrow=1),
         shape=guide_legend(order=1, title = "", nrow=1),
         linetype=guide_legend(order=1, title = "", nrow=1)) +
  ggh4x::facet_wrap2(dplyr::vars(sti,region),
                     strip = ggh4x::strip_nested(text_x = ggh4x::elem_list_text(face = c("bold", "plain"),
                                                                                size = c(7*1.3,7*1.18)), # align with theme
                                                 by_layer_x = TRUE),
                     scales = "free_y",
                     nrow = 3,
                     axes = "x") +
  ggh4x::facetted_pos_scales(
    y = list(
      sti == "Chlamydia" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = scales::label_percent()),
      sti == "Gonorrhoea" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.3), labels = scales::label_percent()),
      sti == "Trichomoniasis" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = scales::label_percent()),
      sti == "Chlamydia" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = NULL),
      sti == "Gonorrhoea" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.3), labels = NULL),
      sti == "Trichomoniasis" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = NULL))) +
  mytheme +
  theme(panel.spacing.y = unit(0.8, "lines"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        legend.box = "vertical",
        legend.margin = margin(unit(c(t=0,r=0,b=0,l=0), "cm"))) +
  labs(x = "", y = "", tag = "B")

ggsave("./plots/fig_2_population_v3.png", width = 18, height = 20, unit = "cm", dpi = 700)  
