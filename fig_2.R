library(patchwork)
library(ggh4x)

mytheme <- theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.2, "cm"),
        legend.position = "bottom",
        legend.margin = margin(unit(c(t=-5,r=5,b=0,l=5), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.2), face="bold",
                                  margin = margin(3,3,3,3),
                                  hjust = 0.5),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.25), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

df1 <- read.csv("./results/prevalence_2020_adjusted.csv") |>
  mutate(region = factor(region, levels = c("Western and Central", "Eastern", "Southern", "SSA"),
                         labels = c("WCA", "EA", "SA", "SSA")),
         sti = factor(sti, levels = c("CT", "NG", "TV"),
                      labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         sex = fct_relevel(sex, "Female"))

df2 <- read.csv("./results/prevalence_alldata_adjusted.csv") |>
  mutate(region = factor(region, levels = c("Western and Central", "Eastern", "Southern"),
                         labels = c("Western and Central Africa", "Eastern Africa", "Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"),
                      labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         sex = fct_relevel(sex, "Female")) |>
  filter(year >= 2005)

df_plot <- readxl::read_xlsx("./data/study_data.xlsx") |>
  mutate(region = fct_collapse(region_analysis, 
                               "Western and Central" = c("Western", "Central")),
         region = factor(region, levels = c("Western and Central", "Eastern", "Southern"),
                         labels = c("Western and Central Africa","Eastern Africa","Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"),
                      labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         sex = factor(sex, levels = c("Female", "Male", "Both sexes")),
         test = fct_collapse(test, "Wet mount" = c("Wet mount", "Wet mount and NAAT")),
         test = factor(test, levels = c("NAAT","Culture", "DFA", "ELISA", "Rapid antigen test", "Wet mount"))) |>
  filter(!sex == "Both sexes")

colour_1 <- c("grey30","#CAA633")

p1 <- df1 |>
  select(sti, region, sex, prev, lwr, upr) |>
  mutate(xpos = case_when(region == "WCA" ~ 1, region == "EA" ~ 2, 
                          region == "SA" ~ 3, region == "SSA" ~ 4.4)) |>
  ggplot(aes(x = xpos, y = prev, ymin = lwr, ymax = upr, group = sex)) +
  geom_col(aes(fill = sex), position = position_dodge(0.9), width = 0.8) +
  geom_linerange(aes(group = sex), position = position_dodge(0.9), size = 0.2) +
  geom_label(aes(label = paste0(sprintf("%.1f",prev*100)), y = prev + 0.01),
             position = position_dodge(0.9),
             size = 7/.pt, label.size = NA, 
             label.padding = unit(0.04,"lines"),
             fill = "white",
             colour = "grey30",
             alpha = 0.85,
             hjust = 0.5, vjust = 0.5) +
  ggh4x::facet_wrap2(~ sti,
                     scales = "free_y",
                     nrow = 1,
                     axes = "x") +
  mytheme +
  scale_x_continuous(breaks = c(1, 2, 3, 4.4),
                     labels = c("WCA", "EA", "SA", "SSA")) +
  # scale_x_continuous(breaks = c(1, 2, 3, 4.4),
  #                    labels = c("Western\nand\nCentral", "Eastern", "Southern", 
  #                               "Sub-Saharan\nAfrica")) +
  ggh4x::facetted_pos_scales(
    y = list(
      sti == "Chlamydia" ~ scale_y_continuous(limits = c(0, 0.25), 
                                              breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
                                              labels = scales::label_percent()),
      sti == "Gonorrhoea"  ~ scale_y_continuous(limits = c(0, 0.25), 
                                                breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
                                                labels = NULL),
      sti == "Trichomoniasis" ~ scale_y_continuous(limits = c(0, 0.25), 
                                                   breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
                                                   labels = NULL))) +
  scale_fill_manual(values = colour_1) +
  labs(x = "", y = "", fill = "", tag = "A")

colour_2 <- c("grey30","grey30")
colour_3 <- c("#af4f2f","#df8d71","#1e5a46","#75884b","#5b859e")
# MetBrewer::met.brewer("Redon", 6, direction=-1) |> as.vector()

p2 <- df2 |>
  ggplot() +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = sex), alpha = 0.15) +
  geom_line(aes(x = year, y = prev, colour = sex, linetype = sex), size = 0.6) +
  scale_fill_manual(values = colour_2) +
  scale_colour_manual(values = colour_2) +
  guides(colour=guide_legend(order=1, title = "")) +
  # lines on new colour scale
  ggnewscale::new_scale_color() +
  geom_linerange(data = df_plot, aes(x = year_mid, ymin = adj_prev_lwr, ymax = adj_prev_upr, colour = test),
                 size = 0.4) + #show.legend = FALSE, 
  geom_point(data = df_plot, aes(x = year_mid, y = adj_prev, colour = test, shape = sex), size = 1, fill = "white") +
  scale_colour_manual(values = colour_3) +
  scale_shape_manual(values = c(19, 21)) +
  guides(colour=guide_legend(order=2, title = ""),
         fill=guide_legend(order=1, title = ""),
         shape=guide_legend(order=1, title = ""),
         linetype=guide_legend(order=1, title = "")) +
  ggh4x::facet_wrap2(dplyr::vars(sti,region),
                     strip = strip_nested(text_x = elem_list_text(face = c("bold", "plain")),
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

(p1 / p2) + plot_layout(height = c(1,3.9))

ggsave("./plots/fig_2.png", width = 18, height = 25, unit = "cm", dpi = 700)  

