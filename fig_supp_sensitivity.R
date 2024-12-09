## SENSITIVITY ANALYSES

library(tidyverse)
library(ggh4x)

# plot theme
mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(0.2, "cm"),
        panel.spacing.y = unit(0.8, "lines"),
        legend.position = "bottom",
        legend.spacing.x = unit(0.6, "cm"),
        legend.margin = margin(unit(c(t=0,r=15,b=0,l=0), "cm")),
        legend.box = 'vertical',
        #legend.margin = margin(unit(c(t=-10,r=5,b=0,l=5), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.2), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.2), face="bold",
                                  margin = margin(3,3,3,3),
                                  hjust = 0, vjust = 0.5),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.25), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

# study observation data
df_obs <- read.csv("./data/final-dataset-v8-adjusted.csv") |>
  filter(!sex == "Both sexes") |>
  filter(denom >= 15) |>
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Youth")) |>
  mutate(region = fct_collapse(region_analysis, "WCA" = c("WA", "CA")),
         region = factor(region, levels = c("WCA","EA","SA"), 
                         labels = c("Western and Central Africa", "Eastern Africa", "Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"),
                      labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         sex = factor(sex, levels = c("Female", "Male")),
         prev = num/denom) |>
  rowwise() |>
  mutate(prev_lwr = binom.test(num, denom)$conf.int[1],
         prev_upr = binom.test(num, denom)$conf.int[2])

# plot colours

colour_1 <- c("#df8d71","#732f30","grey50") 
colour_2 <- c("#df8d71","#732f30")
colour_3 <- c("grey50","aquamarine3","#008080")

# < Perf adjustment: region trends ----

df_adj <- read.csv("./results/prevalence_alldata_adjusted.csv") 
df_unadj <- read.csv("./results/prevalence_alldata_unadjusted.csv") 
df_naat_adj <- read.csv("./results/prevalence_alldata_naat_adjusted.csv")
df_naat_unadj <- read.csv("./results/prevalence_alldata_naat_unadjusted.csv")

# colour_4 <- c("#70C5A8""#008080","grey70","grey35")
colour_4 <- c("#70C5A8", "#008080", "#B9A8C7", "#6B4B8B")

# Female prevalence = all tests + NAAT
rbind(df_adj |> mutate(type = "All tests, adjusted"),
      df_unadj |> mutate(type = "All tests, unadjusted"),
      df_naat_adj |> mutate(type = "NAAT, adjusted"),
      df_naat_unadj |> mutate(type = "NAAT, unadjusted")) |>  
  filter(!region == "SSA") |>
  mutate(type = factor(type, levels = c("All tests, unadjusted", "All tests, adjusted","NAAT, unadjusted", "NAAT, adjusted")),
         region = factor(region, levels = c("WCA","EA","SA"),
                         labels = c("Western and Central Africa","Eastern Africa","Southern Africa")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis"))) |>
  filter(sex == "Female") |>
  ggplot() +
  # RIBBON
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = type), alpha = 0.1) +
  scale_fill_manual(values = colour_4) +
  guides(fill = guide_legend(order=1, "Estimate")) +
  # POINTS
  geom_point(data = df_obs |> filter(sex == "Female", !testcat == "NAAT") |> mutate(type = "Other tests, unadjusted"),
             aes(x = year_mid, y = prev, colour = type), size = 1.4, fill = "white", alpha = 0.7, stroke=NA) +
  geom_point(data = df_obs |> filter(sex == "Female", !testcat == "NAAT") |> mutate(type = "Other tests, adjusted"),
             aes(x = year_mid, y = adj_prev, colour = type), size = 1.4, fill = "white", alpha = 0.7, stroke=NA) +
  geom_point(data = df_obs |> filter(sex == "Female", testcat == "NAAT") |> mutate(type = "NAAT, unadjusted"),
             aes(x = year_mid, y = prev, colour = type), size = 1.4, fill = "white", alpha = 0.7, stroke=NA) +
  geom_point(data = df_obs |> filter(sex == "Female", testcat == "NAAT") |> mutate(type = "NAAT, adjusted"),
             aes(x = year_mid, y = adj_prev, colour = type), size = 1.4, fill = "white", alpha = 0.7, stroke=NA) +
  scale_colour_manual(values = colour_4,
                      breaks = c("Other tests, unadjusted", "Other tests, adjusted","NAAT, unadjusted", "NAAT, adjusted")) +
  guides(colour=guide_legend(order=2, "Observation")) +
  # LINES
  ggnewscale::new_scale_color() +
  geom_line(aes(x = year, y = prev, colour = type), size = 0.6) +
  scale_colour_manual(values = colour_4) +
  guides(colour=guide_legend(order=1, "Estimate")) +
  ggh4x::facet_wrap2(dplyr::vars(sti,region),
                     strip = strip_nested(text_x = elem_list_text(face = c("bold", "plain"),
                                                                  size = c(7.5*1.3,7.5*1.18)), # align with theme
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
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015, 2020)) +
  mytheme +
  theme(panel.spacing.y = unit(0.8, "lines"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        legend.margin = margin(unit(c(t=0,r=15,b=0,l=0), "cm")),
        legend.direction = "horizontal") +
  labs(x = "", y = "")

ggsave("./plots/fig_supp_sensitivity_1.png", width = 18, height = 20, unit = "cm", dpi = 700)  

# < Perf adjustment: SSA in 2020 ----

df_adj <- read.csv("./results/prevalence_2020_adjusted.csv") 
df_unadj <- read.csv("./results/prevalence_2020_unadjusted.csv") 
df_naat_adj <- read.csv("./results/prevalence_2020_naat_adjusted.csv")
df_naat_unadj <- read.csv("./results/prevalence_2020_naat_unadjusted.csv")

rbind(df_unadj |> mutate(type = "All tests, unadjusted"),
      df_adj |> mutate(type = "All tests, adjusted"),
      df_naat_unadj |> mutate(type = "NAAT, unadjusted"),
      df_naat_adj |> mutate(type = "NAAT, adjusted")) |> 
  mutate(type = factor(type, levels = c("All tests, unadjusted", "All tests, adjusted","NAAT, unadjusted", "NAAT, adjusted")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis"))) |>
  filter(region == "SSA") |>
  ggplot() +
  geom_col(aes(x = sex, y = prev, fill = type), 
           width = 0.8, position= position_dodge(0.9)) +
  geom_linerange(aes(x = sex, ymin = lwr, ymax = upr, group = type), position= position_dodge(0.9),
                 size = 0.2,
                 show.legend = FALSE) +
  geom_label(aes(label = paste0(sprintf("%.1f",prev*100)), x = sex, y = prev + 0.005, group = type),
             position = position_dodge(0.9),
             size = 7.5/.pt, label.size = NA, 
             label.padding = unit(0.035,"lines"),
             fill = "white",
             colour = "grey30",
             alpha = 0.85,
             hjust = 0.5, vjust = 0.5) +
  facet_wrap(~sti) +
  scale_fill_manual(values = colour_4) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0,0.108),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1)) +
  mytheme +
  labs(x = "", y = "", fill = "")

ggsave("./plots/fig_supp_sensitivity_2.png", width = 18, height = 7, unit = "cm", dpi = 700)  


# < Perf-adjustment: Ratio ----
df_ratio <- read.csv("./results/ratios.csv") |>
  filter(!type == "Pooled") |>
  mutate(diag = factor(diag, levels = c("All tests, unadjusted", "All tests, adjusted","NAAT, unadjusted", "NAAT, adjusted")),
         sti = factor(sti, levels = c("CT", "NG", "TV"), labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")),
         label = paste0(sprintf("%.2f",ratio)," (n=",n,")"))

df_ratio |>
  ggplot(aes(y = type, x = ratio, xmin = lwr, xmax = upr, colour = diag)) +
  geom_vline(xintercept = 1, size = 0.2, linetype = "dotted") +
  geom_linerange(size = 0.5, position = position_dodge(0.8)) +
  geom_point(size = 1.5, stroke = NA, position = position_dodge(0.8)) +
  geom_label(aes(label = label, y = type, x = 4.1, group = diag),
             position = position_dodge(0.8),
             size = 7.5/.pt, label.size = NA, 
             label.padding = unit(0.035,"lines"),
             fill = "white",
             colour = "grey30",
             alpha = 0.65,
             hjust = 1, vjust = 0.5) +
  facet_wrap(~sti) +
  mytheme +
  scale_x_log10(breaks = c(0.1, 1),
                limits = c(0.1, 4.1),
                labels = scales::number_format(accuracy = 0.1),
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(limits = c("Within study", "Between study"), 
                   labels = c("Within\nstudy","Between\nstudy")) +
  scale_colour_manual(values = colour_4) +
  labs(x = "Prevalence ratio (log scale)", y = "", colour = "")

ggsave("./plots/fig_supp_sensitivity_3.png", width = 18, height = 7, unit = "cm", dpi = 700)  


df_ratio |>
  filter(sti == "Chlamydia") |>
  ggplot(aes(y = type, x = ratio, xmin = lwr, xmax = upr, colour = diag)) +
  geom_vline(xintercept = 1, size = 0.2, linetype = "dotted") +
  geom_linerange(size = 0.5, position = position_dodge(0.8)) +
  geom_point(size = 1.5, stroke = NA, position = position_dodge(0.8)) +
  mytheme +
  scale_x_log10(breaks = c(0.5, 1),
                limits = c(0.4, 1),
                labels = scales::number_format(accuracy = 0.1),
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(limits = c("Within study", "Between study"), 
                   labels = c("Within\nstudy","Between\nstudy")) +
  scale_colour_manual(values = colour_4) +
  labs(x = "Prevalence ratio (log scale)", y = "", colour = "")

df_ratio |>
  filter(sti == "Gonorrhoea") |>
  ggplot(aes(y = type, x = ratio, xmin = lwr, xmax = upr, colour = diag)) +
  geom_vline(xintercept = 1, size = 0.2, linetype = "dotted") +
  geom_linerange(size = 0.5, position = position_dodge(0.8)) +
  geom_point(size = 1.5, stroke = NA, position = position_dodge(0.8)) +
  mytheme +
  scale_x_log10(breaks = c(0.5, 1),
                limits = c(0.4, 1.3),
                labels = scales::number_format(accuracy = 0.1),
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(limits = c("Within study", "Between study"), 
                   labels = c("Within\nstudy","Between\nstudy")) +
  scale_colour_manual(values = colour_4) +
  labs(x = "Prevalence ratio (log scale)", y = "", colour = "")

df_ratio |>
  filter(sti == "Trichomoniasis") |>
  ggplot(aes(y = type, x = ratio, xmin = lwr, xmax = upr, colour = diag)) +
  geom_vline(xintercept = 1, size = 0.2, linetype = "dotted") +
  geom_linerange(size = 0.5, position = position_dodge(0.8)) +
  geom_point(size = 1.5, stroke = NA, position = position_dodge(0.8)) +
  mytheme +
  scale_x_log10(breaks = c(0.5, 1),
                limits = c(0.1, 0.5),
                labels = scales::number_format(accuracy = 0.1),
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(limits = c("Within study", "Between study"), 
                   labels = c("Within\nstudy","Between\nstudy")) +
  scale_colour_manual(values = colour_4) +
  labs(x = "Prevalence ratio (log scale)", y = "", colour = "")
