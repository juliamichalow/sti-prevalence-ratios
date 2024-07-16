## SENSITIVITY ANALYSIS
# Check whether trends vary for different population groups and/or time periods

library(tidyverse)
library(ggh4x)

## 1. Compare population groups ----
# Trends in prevalence among females, using different population groups to inform estimates 

# study observation data
df_obs <- readxl::read_xlsx("./data/study_data.xlsx") |>
  mutate(region = fct_collapse(region_analysis, 
                               "Western and Central" = c("Western", "Central")),
         region = factor(region, levels = c("Western and Central","Eastern","Southern"),
                         labels = c("Western and Central Africa","Eastern Africa","Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         sex = factor(sex, levels = c("Female", "Male", "Both sexes")))

df <- df_obs |> filter(sex == "Female") |>
  mutate(group = case_when(population == "ANC attendees" ~ "ANC attendees",
                           TRUE ~ "Other female populations")) |>
  mutate(group = factor(group, levels = c("ANC attendees", "Other female populations"))))

         
# model estimates
df_all <- read.csv("./results/prevalence_alldata_adjusted.csv") 
df_fm <- read.csv("./results/prevalence_alldata_fm.csv")
df_anc <- read.csv("./results/prevalence_alldata_anc.csv")

df_est <- rbind(df_all |> filter(sex == "Female") |> mutate(type = "Male and female (all populations)"),
                df_fm |> mutate(type = "Female (all populations)"),
                df_anc |> mutate(type = "ANC attendees")) |>
  mutate(type = factor(type, levels = c("ANC attendees", "Female (all populations)","Male and female (all populations)")),
         region = factor(region, levels = c("Western and Central","Eastern","Southern"),
                         labels = c("Western and Central Africa","Eastern Africa","Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")))

# plot theme
mytheme <- theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        panel.spacing.y = unit(0.8, "lines"),
        legend.position = "bottom",
        legend.spacing.x = unit(0.2, "cm"),
        legend.margin = margin(unit(c(t=0,r=15,b=0,l=0), "cm")),
        legend.box = 'vertical',
        #legend.margin = margin(unit(c(t=-10,r=5,b=0,l=5), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.2), face="bold",
                                  margin = margin(3,3,3,3)),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.25), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))
# plot

colour_1 <- c("#CC544D","#7577D1","#767676") 
colour_2 <- c("#CC544D","#7577D1")

df_est |> 
  ggplot() +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = type), alpha = 0.15) +
  geom_line(aes(x = year, y = prev, colour = type), size = 0.6) +
  scale_fill_manual(values = colour_1) +
  scale_colour_manual(values = colour_1) +
  guides(colour = guide_legend(order=2, title = "Estimate population"),
         fill = guide_legend(order=2, title = "Estimate population")) +
  # points on new colour scale
  ggnewscale::new_scale_color() +
  geom_point(data = df, aes(x = year_mid, y = adj_prev, colour = group), size = 1) +
  scale_colour_manual(values = colour_2) +
  guides(colour=guide_legend(order=1, title = "Study population")) +
  ggh4x::facet_wrap2(dplyr::vars(sti,region),
                     strip = strip_nested(text_x = elem_list_text(face = c("bold", "plain")),
                                          by_layer_x = TRUE),
                     scales = "free_y",
                     nrow = 3) +
  ggh4x::facetted_pos_scales(
    y = list(
      sti == "Chlamydia" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.4), labels = scales::label_percent()),
      sti == "Gonorrhoea" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.2), labels = scales::label_percent()),
      sti == "Trichomoniasis" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = scales::label_percent()),
      sti == "Chlamydia" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.4), labels = NULL),
      sti == "Gonorrhoea" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.2), labels = NULL),
      sti == "Trichomoniasis" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = NULL))) +
  mytheme +
  labs(x = "", y = "")

ggsave("./plots/fig_supp_sensitivity_1.png", width = 18, height = 20, unit = "cm", dpi = 700)  

df_est |> 
  filter(year == 2020) |>
  select(!c(lwr, upr)) |>
  pivot_wider(names_from = type, values_from = prev)

df_est |> 
  filter(year == 2005) |>
  select(!c(lwr, upr)) |>
  pivot_wider(names_from = type, values_from = prev)

## 2. Compare time periods ----
# Compare trends in prevalence among females, using different time periods to inform estimates

# study observation data
df <- df_obs |> filter(sex == "Female") |>
  mutate(group = case_when(year_mid >= 2005 ~ "2005-2023",
                           TRUE ~ "2000-2004")) |>
  mutate(group = factor(group, levels = c("2000-2004", "2005-2023")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")))


# model estimates
df_all <- read.csv("./results/prevalence_alldata_adjusted.csv") 
df_2005 <- read.csv("./results/prevalence_alldata_2005onwards.csv")

df_est <- rbind(df_all |> filter(sex == "Female") |> mutate(type = "2000-2023"),
                df_2005 |> filter(sex == "Female", year >= 2005) |> mutate(type = "2005-2023")) |>
  mutate(type = factor(type, levels = c("2000-2023", "2005-2023")),
         region = factor(region, levels = c("Western and Central","Eastern","Southern"),
                         labels = c("Western and Central Africa","Eastern Africa","Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")))

# plot

colour_1 <- c("#58D0BF","coral1")
colour_2 <- c("#58D0BF","coral1")

df_est |> 
  ggplot() +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = type), alpha = 0.15) +
  geom_line(aes(x = year, y = prev, colour = type), size = 0.6) +
  scale_fill_manual(values = colour_1) +
  scale_colour_manual(values = colour_1) +
  guides(colour = guide_legend(order=2, title = "Estimate time period"),
         fill = guide_legend(order=2, title = "Estimate time period")) +
  # points on new colour scale
  ggnewscale::new_scale_color() +
  geom_point(data = df, aes(x = year_mid, y = adj_prev, colour = group), size = 1) +
  scale_colour_manual(values = colour_2) +
  guides(colour=guide_legend(order=1, title = "Study mid-year")) +
  ggh4x::facet_wrap2(dplyr::vars(sti,region),
                     strip = strip_nested(text_x = elem_list_text(face = c("bold", "plain")),
                                          by_layer_x = TRUE),
                     scales = "free_y",
                     nrow = 3) +
  ggh4x::facetted_pos_scales(
    y = list(
      sti == "Chlamydia" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.4), labels = scales::label_percent()),
      sti == "Gonorrhoea" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.2), labels = scales::label_percent()),
      sti == "Trichomoniasis" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = scales::label_percent()),
      sti == "Chlamydia" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.4), labels = NULL),
      sti == "Gonorrhoea" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.2), labels = NULL),
      sti == "Trichomoniasis" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = NULL))) +
  mytheme +
  theme(panel.spacing.y = unit(0.8, "lines"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        legend.margin = margin(unit(c(t=0,r=15,b=0,l=0), "cm")),
        legend.direction = "horizontal") +
  labs(x = "", y = "")

ggsave("./plots/fig_supp_sensitivity_2.png", width = 18, height = 20, unit = "cm", dpi = 700)  

## 3. Compare test performance adjustment ----

df_adj <- read.csv("./results/prevalence_alldata_adjusted.csv") 
df_unadj <- read.csv("./results/prevalence_alldata_unadjusted.csv") 

colour_4 <- MetBrewer::met.brewer("Egypt", 2, direction = -1) 

# Female prevalence
rbind(df_adj |> mutate(type = "Adjusted"),
      df_unadj |> mutate(type = "Unadjusted")) |>  
  filter(!region == "SSA") |>
  mutate(type = fct_relevel(type, "Unadjusted"),
         region = factor(region, levels = c("Western and Central","Eastern","Southern"),
                         labels = c("Western and Central Africa","Eastern Africa","Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis"))) |>
  filter(sex == "Female") |>
  ggplot() +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = type), alpha = 0.15) +
  geom_line(aes(x = year, y = prev, colour = type), size = 0.6) +
  geom_point(data = df_obs |> filter(sex == "Female") |> mutate(type = "Adjusted"), 
             aes(x = year_mid, y = adj_prev, colour = type), size = 1) +
  geom_point(data = df_obs |> filter(sex == "Female") |> mutate(type = "Unadjusted"), 
             aes(x = year_mid, y = prev, colour = type), size = 1) +
  scale_fill_manual(values = colour_4) +
  scale_colour_manual(values = colour_4) +
  guides(colour = guide_legend(order=1, title = ""),
         fill = guide_legend(order=1, title = "")) +
  ggh4x::facet_wrap2(dplyr::vars(sti,region),
                     strip = strip_nested(text_x = elem_list_text(face = c("bold", "plain")),
                                          by_layer_x = TRUE),
                     scales = "free_y",
                     nrow = 3) +
  ggh4x::facetted_pos_scales(
    y = list(
      sti == "Chlamydia" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.4), labels = scales::label_percent()),
      sti == "Gonorrhoea" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.2), labels = scales::label_percent()),
      sti == "Trichomoniasis" & region == "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = scales::label_percent()),
      sti == "Chlamydia" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.4), labels = NULL),
      sti == "Gonorrhoea" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.2), labels = NULL),
      sti == "Trichomoniasis" & region != "Western and Central Africa" ~ scale_y_continuous(limits = c(0, 0.6), labels = NULL))) +
  mytheme +
  theme(panel.spacing.y = unit(0.8, "lines"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        legend.margin = margin(unit(c(t=0,r=15,b=0,l=0), "cm")),
        legend.direction = "horizontal") +
  labs(x = "", y = "")

ggsave("./plots/fig_supp_sensitivity_3.png", width = 18, height = 20, unit = "cm", dpi = 700)  

## 4 Compare test performance adjustment in 2020 ----

df_adj <- read.csv("./results/prevalence_2020_adjusted.csv") 
df_unadj <- read.csv("./results/prevalence_2020_unadjusted.csv") 

colour_4 <- MetBrewer::met.brewer("Egypt", 2, direction = -1) 

rbind(df_adj |> mutate(type = "Adjusted"),
      df_unadj |> mutate(type = "Unadjusted")) |>  
  mutate(type = fct_relevel(type, "Unadjusted"),
         region = factor(region, levels = c("Western and Central","Eastern","Southern","SSA"),
                         labels = c("WCA", "EA", "SA", "SSA")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         sex = fct_relevel(sex, "Female")) |>
  mutate(xpos = case_when(region == "WCA" ~ 1, region == "EA" ~ 2, 
                          region == "SA" ~ 3, region == "SSA" ~ 4.4)) |>
  ggplot(aes(x = xpos, y = prev, ymin = lwr, ymax = upr, group = type)) +
  geom_col(aes(fill = type), position = position_dodge(0.9), width = 0.8) +
  geom_linerange(aes(group = type), position = position_dodge(0.9), size = 0.2) +
  geom_label(aes(label = paste0(sprintf("%.1f",prev*100)), y = prev + 0.008),
             position = position_dodge(0.9),
             size = 7/.pt, label.size = NA, 
             label.padding = unit(0.04,"lines"),
             fill = "white",
             colour = "grey30",
             alpha = 0.85,
             hjust = 0.5, vjust = 0.5) +
  facet_grid(sex~sti) +
  mytheme +
  theme(legend.margin = margin(unit(c(t=-10,r=0,b=10,l=0), "cm")),
        strip.text.y = element_text(color="black", size = rel(1.2), face="bold",
                                    margin = margin(t=4,r=4,b=4,l=4))) +
  scale_x_continuous(breaks = c(1, 2, 3, 4.4),
                     labels = c("WCA", "EA", "SA","SSA")) +
  scale_y_continuous(labels = scales::label_percent(), 
                     breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
                     limits = c(0, 0.25)) +
  scale_fill_manual(values = colours) +
  #coord_cartesian(ylim=c(0, 0.25)) +
  labs(x = "", y = "", fill = "", tag = "")

ggsave("./plots/fig_supp_sensitivity_4.png", width = 18, height = 12, unit = "cm", dpi = 700)  
