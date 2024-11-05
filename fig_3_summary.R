# Forest plots (manual)
library(patchwork)

# Weighted summary ratios ----

dat_summary <- read.csv("./results/ratios.csv") |>
  relocate(type) |>
  mutate(est = paste0(sprintf(ratio,fmt = '%#.2f'), " (", sprintf(lwr,fmt = '%#.2f'), 
                      "-", sprintf(upr,fmt = '%#.2f'), ")"),
         type = factor(type, levels = c("Unadjusted","Within study", "Between study"),
                       labels = c("Pooled within-study ratio","Adjusted within-study ratio", "Adjusted between-study ratio")),
         sum = paste0("n = ",n, ", \u03C4","\u00B2"," = ", sprintf(tau_study,fmt = '%#.2f'))) |>
  arrange(sti, type) |>
  # Headings
  bind_rows(data.frame(type = "Ratio type", sti = NA, year = NA, ratio = NA, lwr = NA,
                       upr = NA, n = NA, est = "Ratio (95% CI)", sum = NA))

dat_summary_ct <- dat_summary |> filter(sti == "CT" | is.na(sti)) |>
  mutate(type = forcats::fct_inorder(type)) |>
  mutate(type = fct_relevel(type, "Ratio type"))

dat_summary_ng <- dat_summary |> filter(sti == "NG"| is.na(sti))|>
  mutate(type = forcats::fct_inorder(type)) |>
  mutate(type = fct_relevel(type, "Ratio type"))

dat_summary_tv <- dat_summary |> filter(sti == "TV"| is.na(sti)) |>
  mutate(type = forcats::fct_inorder(type)) |>
  mutate(type = fct_relevel(type, "Ratio type"))

# Plot ----  

# dat = study data
# dat_summary = weighted ratios
# dim = dimensions for top and bottom of forest
plot_forest <- function(dat, dat_summary, dim, who_ratio, title) {
  
  textsize <- 12
  x_axis_min <- 0.08
  x_axis_max <- 1.5
  
  # Summary measures at bottom of plot
  
  p2_left <- dat_summary |>
    ggplot(aes(y = fct_rev(type))) +
    geom_text(aes(x = 0, label = type), hjust = 0, 
              fontface = ifelse(dat_summary$type == "Ratio type", "bold", "plain"),
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
    filter(!ratio == "Ratio type") |>
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
    # add heading line
    rbind(expand.grid(type = "Ratio type", sti = NA, year = NA, tau_study = NA, n = NA, sum = NA,
                      name = c("lwr", "ratio","ratio2","upr"), value = 0, x = 5, y_pos = 4, y = 4) |>
            mutate(y = case_when(name == "ratio" ~ 4.25,
                                 name == "ratio2" ~ 4.75, 
                                 TRUE ~ y))) |>
    # reorder coordinates
    mutate(type = factor(type, 
                         levels = c(" ", "Ratio type", "Pooled within-study ratio", "Adjusted within-study ratio", "Adjusted between-study ratio")),
           name = factor(name, levels = c("lwr", "ratio", "upr", "ratio2"))) |>
    arrange(type,name)

  
  # Define x-axis breaks for p2_mid
  # Define the breaks based on the value of `who_ratio`
  x_breaks <- if (who_ratio < 0.5) {
    c(who_ratio, 0.5, 1)
  } else {
    c(0.1, 0.5, who_ratio)
  }
  
  p2_mid <- ggplot(data = diamond_data) +
    geom_vline(xintercept = 1, linetype = "solid", size = 0.3, colour = "gray45") +
    geom_vline(xintercept = who_ratio, linetype="dashed", size = 0.4, colour = "mediumblue") +
    geom_polygon(aes(x = x, y = y, group = type),
                 fill = "grey", color = "black", size = 0.45) +
    theme_classic(base_size = textsize) +
    # NOTE ylim is exactly lined up with y-axes for p2_left and p2_right
    scale_x_log10(breaks = x_breaks,
                  labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(xlim = c(x_axis_min, x_axis_max), ylim = c(0.55,5)) +
    scale_y_continuous(breaks = c(1,2,3,4,5)) +
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
    geom_text(aes(x = 0.05, label = est), hjust = 0, 
              fontface = ifelse(dat_summary$est == "Ratio (95% CI)","bold","plain"),
              size = textsize / .pt) +
    theme_void(base_size = textsize) +
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    coord_cartesian(xlim = c(0,0.3))
  
  # full plot

  p_bot <- (p2_left + p2_mid + p2_right) +
    plot_layout(widths = c(1.999, 1, 0.8)) +
    plot_annotation(title = title) &
      theme(plot.margin = unit(c(t=0,r=0,b=0,l=5), "pt"),
            plot.title = element_text(size = textsize*1.1, face="bold", hjust = -0.02, #-0.014
                                      margin = margin(t = 0, r = 0, b = 4, l = 0)),
            plot.title.position = "plot")
  
} 

p_ct <- plot_forest(dat_ct, dat_summary_ct, c(26/3, 0.1, 1), 0.8, "A    Chlamydia")
p_ng <- plot_forest(dat_ng, dat_summary_ng, c(26/3, 0.1, 1), 0.86, "B    Gonorrhoea")
p_tv <- plot_forest(dat_tv, dat_summary_tv, c(12/3, 0.1, 1), 0.1, "C    Trichomoniasis")


p <- wrap_elements(p_ct) / wrap_elements(p_ng) / wrap_elements(p_tv) +
  plot_layout(height = c(1,1,1)) &
  theme(plot.margin = unit(c(t = 0.05, r = 0, b = 0, l = 0), "cm"))

p

ggsave("./plots/fig_3_summary.png", p, width = 22, height = 14, unit = "cm", dpi = 700)
