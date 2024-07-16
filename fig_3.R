mytheme <- theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(),
        # plot.margin = unit(c(t=0.5,r=0,b=0.5,l=0), "cm"),
        panel.spacing = unit(0.2, "cm"),
        legend.position = "bottom",
        legend.margin = margin(unit(c(t=-5,r=0,b=0,l=0), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.0)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1)),
        strip.text = element_text(color="black", size = rel(1.2), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.25), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

df1 <- read.csv("./Tables/ratio_sex.csv") |>
  mutate(region = factor(region_analysis, levels = c("west and centre", "east", 
                                                     "south", "SSA"),
                         labels = c("CWA", "EA", "SA", "SSA")),
         type = factor(type, levels = c("between", "within"),
                       labels = c("Between study", "Within study")),
         sti = factor(ratio, levels = c("CT M-F", "NG M-F", "TV M-F"),
                      labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")))

MetBrewer::met.brewer("Tiepolo", 8)
as.vector(MetBrewer::met.brewer("Redon", 30))

colour <- c("#3B6947","#af4f2f")
scales::show_col(colour)


df1 |> filter(region == "SSA") |>
  ggplot(aes(x = sti, y = est, ymin = lwr, ymax = upr, group = type)) +
  geom_hline(yintercept = 1, size = 0.2, linetype = "dotted") +
  geom_pointrange(aes(colour = type), position = position_dodge(0.5), 
                  size = 0.1) +
  geom_label(aes(label = sprintf("%.2f",est), y = est*1.3),
             position = position_dodge(0.5),
             size = 6/.pt, label.size = NA,
             label.padding = unit(0.08,"lines"),
             fill = "white",
             colour = "grey30",
             alpha = 0.9,
             hjust = 0.5, vjust = 0.5) +
  mytheme +
  scale_y_continuous(breaks = c(0.01, 0.1, 1, 10), limits = c(0.01, 10)) +
  scale_colour_manual(values = colour) +
  coord_trans(y = "log10") +
  labs(x = "", y = "Ratio (log scale)", colour = "", tag = "")

ggsave("./Figures/fig_sexratio.tiff", width = 8, height = 8, unit = "cm", dpi = 700)


df1 |> 
  ggplot(aes(x = region, y = est, ymin = lwr, ymax = upr, group = type)) +
  geom_hline(yintercept = 1, size = 0.2, linetype = "dotted") +
  geom_pointrange(aes(colour = type), position = position_dodge(0.5), 
                  size = 0.1) +
  geom_label(aes(label = sprintf("%.2f",est), y = est*1.35),
             position = position_dodge(0.5),
             size = 6/.pt, label.size = NA,
             label.padding = unit(0.08,"lines"),
             fill = "white",
             colour = "grey30",
             alpha = 0.9,
             hjust = 0.5, vjust = 0.5) +
  facet_grid(~sti) +
  mytheme +
  scale_y_continuous(breaks = c(0.01, 0.1, 1, 10)) +
  scale_colour_manual(values = colour) +
  coord_trans(y = "log10") +
  labs(x = "", y = "Ratio (log scale)", colour = "", tag = "")

ggsave("./Figures/fig_sexratio_region.tiff", width = 20, height = 8, unit = "cm", dpi = 700)

df1 |> 
  filter(type == "Between study") |>
  ggplot(aes(x = region, y = est, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 1, size = 0.2, linetype = "dotted") +
  geom_pointrange(size = 0.08) +
  geom_label(aes(label = sprintf("%.2f",est), y = est*1.35),
             position = position_dodge(0.5),
             size = 6/.pt, label.size = NA,
             label.padding = unit(0.08,"lines"),
             fill = "white",
             colour = "grey30",
             alpha = 0.9,
             hjust = 0.5, vjust = 0.5) +
  facet_grid(~sti) +
  mytheme +
  scale_y_continuous(breaks = c(0.01, 0.1, 1, 10)) +
  scale_colour_manual(values = colour) +
  coord_trans(y = "log10") +
  labs(x = "", y = "Ratio (log scale)", colour = "", tag = "")

ggsave("./Figures/fig_sexratio_region.tiff", width = 18, height = 7, unit = "cm", dpi = 700)

