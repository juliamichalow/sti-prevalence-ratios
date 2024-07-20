# Forest plots (manual)
library(patchwork)

# Study data for plot ----

dat_plot <- read.csv("./data/data_withinstudy_adjusted.csv") |>
  pivot_wider(names_from = sex, values_from = c(adj_prev, adj_se)) |>
  mutate(ratio = adj_prev_Male/adj_prev_Female,
         ratio_log = log(ratio),
         ratio_se_log = sqrt((adj_se_Male / adj_prev_Male)^2 + (adj_se_Female / adj_prev_Female)^2),
         lwr = exp(ratio_log - 1.96 * ratio_se_log),
         upr = exp(ratio_log + 1.96 * ratio_se_log)) |>
  select(!c(starts_with("adj_"), ratio_log, ratio_se_log)) |>
  mutate(est = paste0(sprintf(ratio,fmt = '%#.2f')," (", sprintf(lwr,fmt = '%#.2f'), 
                      "-", sprintf(upr,fmt = '%#.2f'), ")"),
         est = case_when(est == "1.96 (0.38-10.07)" ~ "1.96 (0.38-10.1)", TRUE ~ est),
         year_mid = as.character(year_mid)) |>
  # Headings
  bind_rows(data.frame(study_name = NA, study_id = "Study",
                       year_mid = "Year", region = "Region",
                       country = "Country", pop_analysis = "Population",
                       sti = NA, ratio = NA, lwr = NA,
                       upr = NA, est = "Ratio (95% CI)"))

dat_ct <- dat_plot |> filter(sti == "CT" | is.na(sti)) |>
  mutate(region = factor(region, levels = c("Region","WCA", "EA", "SA"))) |>
  arrange(region, year_mid, country) |>
  mutate(study_id = forcats::fct_inorder(study_id)) |>
  mutate(study_id = fct_relevel(study_id, "Study"))

dat_ng <- dat_plot |> filter(sti == "NG" | is.na(sti)) |>
  mutate(region = factor(region, levels = c("Region","WCA", "EA", "SA"))) |>
  arrange(region, year_mid, country) |>
  mutate(study_id = forcats::fct_inorder(study_id)) |>
  mutate(study_id = fct_relevel(study_id, "Study"))

dat_tv <- dat_plot |> filter(sti == "TV" | is.na(sti)) |>
  mutate(region = factor(region, levels = c("Region","WCA", "EA", "SA"))) |>
  arrange(region, year_mid, country) |>
  mutate(study_id = forcats::fct_inorder(study_id)) |>
  mutate(study_id = fct_relevel(study_id, "Study"))

# Weighted summary ratios ----

dat_summary <- read.csv("./results/ratios.csv") |>
  relocate(type) |>
  mutate(est = paste0(sprintf(ratio,fmt = '%#.2f'), " (", sprintf(lwr,fmt = '%#.2f'), 
                      "-", sprintf(upr,fmt = '%#.2f'), ")"),
         type = factor(type, levels = c("Unadjusted","Within study", "Between study"),
                       labels = c("Pooled within-study ratio","Adjusted within-study ratio", "Adjusted between-study ratio")),
         sum = paste0("n = ",n, ", \u03C4","\u00B2"," = ", sprintf(tau_study,fmt = '%#.2f'))) |>
  arrange(sti, type)

dat_summary_ct <- dat_summary |> filter(sti == "CT")
dat_summary_ng <- dat_summary |> filter(sti == "NG")
dat_summary_tv <- dat_summary |> filter(sti == "TV")

# Plot ----  

# dat = study data
# dat_summary = weighted ratios
# dim = dimensions for top and bottom of forest
plot_forest <- function(dat, dat_summary, dim, who_ratio) {
  
  textsize <- 8
  
  # Studies at top of plot 
  
  p1_left <- dat |>
    ggplot(aes(y = fct_rev(study_id))) +
    geom_text(aes(x = 0, label = study_id), hjust = 0, 
              fontface = ifelse(dat$study_id == "Study", "bold", "plain"),
              size = textsize / .pt) +
    geom_text(aes(x = 1.05, label = region, hjust = 0),
              fontface = ifelse(dat$region == "Region", "bold", "plain"),
              size = textsize / .pt) +
    geom_text(aes(x = 1.6, label = country, hjust = 0),
              fontface = ifelse(dat$country == "Country", "bold", "plain"),
              size = textsize / .pt) +
    geom_text(aes(x = 2.5, label = year_mid, hjust = 0),
              fontface = ifelse(dat$year_mid == "Year", "bold", "plain"),
              size = textsize / .pt) +
    theme_void(base_size = textsize) +
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    coord_cartesian(xlim = c(0,2.8))
  
  p1_mid <- dat |>    
    ggplot(aes(y = fct_rev(study_id))) + 
    geom_segment(aes(x = 1, xend = 1), y = 0, yend = "Study", linetype = "solid",
                 size = 0.2, colour = "gray45") +
    geom_segment(aes(x = who_ratio, xend = who_ratio), y = 0, yend = "Study", linetype = "dashed",
                 size = 0.3, colour = "mediumblue") +
    geom_point(aes(x = ratio), shape=15, size=1.5, colour = "black") +
    geom_linerange(aes(xmin = lwr, xmax = upr), size = 0.4, colour = "black") +
    geom_label(aes(x = who_ratio, y = "Study", label = "Global ratio"),
               size = textsize*0.85 / .pt, label.size = NA,
               label.padding = unit(0.08, "lines"),
               colour = "mediumblue", fill = "white", fontface = "plain",
               hjust = 0.5, vjust = 0.5) +
    theme_classic(base_size = textsize) +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          # remove bottom axis and numbers
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    scale_x_log10() +
    coord_cartesian(xlim = c(0.039, 1.98))
  
  p1_right <- dat |>
    ggplot(aes(y = fct_rev(study_id))) +
    geom_text(aes(x = 0.05, label = est), hjust = 0,
              fontface = ifelse(dat$est == "Ratio (95% CI)", "bold", "plain"),
              size = textsize / .pt) +
    theme_void(base_size = textsize) +
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    coord_cartesian(xlim = c(0,0.3))
  
  # Summary measures at bottom of plot
  
  p2_left <- dat_summary |>
    ggplot(aes(y = fct_rev(type))) +
    geom_text(aes(x = 0, label = type), hjust = 0, fontface = "bold",
              size = textsize / .pt) +
    geom_text(aes(x = 1.6, label = sum), hjust = 0, fontface = "plain",
              colour = "grey30",
              size = textsize*0.85 / .pt) +
    theme_void(base_size = textsize) +
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    coord_cartesian(xlim = c(0,2.8))
  
  # p2_mid = summary measures using diamonds, which need 4 coordinates per ratio:
  # (x,y) corresponding to (lwr,y_mid), (ratio,y_up), (upr, y_mid), (ratio, y_low)
  # where y_mid, y_up and y_low define height of diamond
  
  height <- 0.5
  y_positions <- c(1, 2)  # Fixed y positions for each type
  
  diamond_data <- dat_summary |> 
    select(!est) |>
    # repeat ratio value so that two coordinates at centre
    mutate(ratio2 = ratio) |>
    pivot_longer(cols = c(ratio, ratio2, lwr, upr)) |>
    mutate(x = value,
           # Fixed y position for each ratio (reverse order)
           y_pos = c(rep(3,4), rep(2,4), rep(1,4)),
           y = case_when(name %in% c("lwr", "upr") ~ y_pos,
                         name %in% "ratio" ~ y_pos + height/2,
                         name %in% "ratio2" ~ y_pos - height/2)) |>
    # reorder coordinates
    mutate(type = factor(type, 
                         levels = c(" ", "Pooled within-study ratio", "Adjusted within-study ratio", "Adjusted between-study ratio")),
           name = factor(name, levels = c("lwr", "ratio", "upr", "ratio2"))) |>
    arrange(type,name)
  
  # Define x-axis breaks for p2_mid
  # Define the breaks based on the value of `who_ratio`
  x_breaks <- if (who_ratio < 0.5) {
    c(who_ratio, 0.5, 1.5)
  } else {
    c(0.1, 0.5, who_ratio, 1.5)
  }
  
  p2_mid <- ggplot(data = diamond_data) +
    geom_vline(xintercept = 1, linetype = "solid", size = 0.22, colour = "gray45") +
    geom_vline(xintercept = who_ratio, linetype="dashed", size = 0.3, colour = "mediumblue") +
    geom_polygon(aes(x = x, y = y, group = type),
                 fill = "grey", color = "black", size = 0.4) +
    theme_classic(base_size = textsize) +
    # NOTE ylim is exactly lined up with y-axes for p2_left and p2_right
    scale_x_log10(breaks = x_breaks,
                  labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(xlim = c(0.039, 1.98), ylim = c(0.55,3.45)) +
    scale_y_continuous(breaks = c(1,2,3,4)) +
    # make who_ratio blue + change other theme
    theme(axis.text.x = element_text(colour = ifelse(x_breaks == who_ratio, "mediumblue", "grey30"),
                                     face = ifelse(x_breaks == who_ratio, "bold", "plain"),
                                     size = rel(1)),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = rel(0.9), colour = "grey30"),
          axis.line.x = element_line(colour = "gray45", size = 0.3),
          axis.ticks.x = element_line(colour = "gray45", size = 0.3),
          plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    labs(x = "Ratio (log scale)")
  
  p2_right <- dat_summary |>
    ggplot(aes(y = fct_rev(type))) +
    geom_text(aes(x = 0.05, label = est), hjust = 0, fontface = "bold",
              size = textsize / .pt) +
    theme_void(base_size = textsize) +
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    coord_cartesian(xlim = c(0,0.3))
  
  # dotted line in between the two plots
  
  line <- ggplot() +
    geom_segment(aes(x = 0, xend = 5, y = 1.5, yend = 1.5), 
                 linetype = "solid", colour = "gray45", size = 0.15) +
    theme_void() +
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    coord_cartesian(ylim = c(0, 2), xlim = c(0.13, 4.93))
  
  # full plot
  
  p_top <- (p1_left + p1_mid + p1_right) +
    plot_layout(widths = c(2, 1, 0.8))  
  
  p_bot <- (p2_left + p2_mid + p2_right) +
    plot_layout(widths = c(1.999, 1, 0.8)) 
  
  (p_top / line / p_bot) +
    plot_layout(heights = dim) &
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"))
  
  # FOR PLOT WITHOUT DOTTED LINE 
  # (p_top / plot_spacer() / p_bot) +
  #   plot_layout(heights = dim) 
} 

p_ct <- plot_forest(dat_ct, dat_summary_ct, c(4, 0.1, 0.95), 0.8)
p_ng <- plot_forest(dat_ng, dat_summary_ng, c(4.1, 0.1, 0.94), 0.86)
p_tv <- plot_forest(dat_tv, dat_summary_tv, c(2, 0.1, 0.95), 0.1)

# FOR PLOT WITHOUT DOTTED LINE
# p_ct <- plot_forest(dat_ct, dat_summary_ct, c(3.9, -0.36, 0.6), 0.8)
# p_ng <- plot_forest(dat_ng, dat_summary_ng, c(4, -0.36, 0.6), 0.86)
# p_tv <- plot_forest(dat_tv, dat_summary_tv, c(2, -0.36, 0.6), 0.1)

p <- wrap_elements(p_ct) / wrap_elements(p_ng) / wrap_elements(p_tv) +
  plot_layout(height = c(1.9, 2, 1.25)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 9, face="bold"))

p

ggsave("./plots/fig_3.png", p, width = 16, height = 23, unit = "cm", dpi = 700)
