# Number of studies per country by symptom

# Load packages
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(patchwork)

africa <- ne_countries(scale = "small", type = "countries", continent = "Africa", returnclass = "sf") |>
  # filter(region_wb == "Sub-Saharan Africa") |>
  # Rename countries to match
  mutate(admin = case_when(admin == "Ivory Coast" ~ "Côte d'Ivoire",
                           admin == "Democratic Republic of the Congo" ~ "Democratic Republic of Congo",
                           admin == "Republic of the Congo" ~ "Congo",
                           admin == "Gambia" ~ "The Gambia",
                           admin == "eSwatini" ~ "Eswatini",
                           admin == "United Republic of Tanzania" ~ "Tanzania",
                           admin == "São Tomé and Principe" ~ "Sao Tome and Principe",
                           TRUE ~ admin)) |>
  select(admin, geometry)

notssa <- ne_countries(scale = "small", type = "countries", continent = "Africa", returnclass = "sf") |>
  filter(region_wb != "Sub-Saharan Africa")

# Calculate number of unique studies per country
dat <- readxl::read_xlsx("./data/study_data.xlsx") |>
  left_join(read.csv("./data/study_name.csv"),
            by = c("study_id")) |>
  filter(year_mid >= 2005) |>
  group_by(sti,sex,study_name,country) |>
  summarise(n_id = n_distinct(study_id)) |>
  mutate(n_study = case_when(!is.na(study_name) ~ 1, is.na(study_name) ~ n_id)) |>
  # Separate studies across multiple countries (double count)
  mutate(country = strsplit(country, ", ")) |>
  unnest(country) |>
  group_by(sti, sex, country) |>
  summarise(n_study = sum(n_study)) |>
  mutate(sti = factor(sti, levels = c("CT", "NG", "TV"),
                      labels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis")))

# Combine with map data
df_country <- africa |>
  crossing(expand.grid(sti = c("Chlamydia", "Gonorrhoea", "Trichomoniasis"),
                       sex = c("Female", "Male"))) |>
  mutate(sti = factor(sti, levels = c("Chlamydia", "Gonorrhoea", "Trichomoniasis"))) |>
  left_join(dat, by = c("admin" = "country", "sti", "sex")) |>
  mutate(n_study = case_when(admin %in% notssa$admin ~ 0, TRUE ~ n_study))

my_theme <- theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # move legend closer to figure
        legend.margin = margin(10,10,10,30),
        legend.key.width = unit(0.5, 'cm'),
        legend.key.height = unit(0.35, "cm"),
        #legend.box = 'horizontal',
        legend.position = "bottom",
        #axis.line = element_line(colour="black"),
        panel.spacing = unit(0.2,"cm"),
        # text size
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = rel(0.8)),
        legend.title = element_text(size = rel(0.9), face="bold", vjust=0.88),
        #change facet labels
        strip.text = element_text(color="black", size = rel(1.1), face="bold",vjust=1.5), 
        # axis ticks
        axis.ticks = element_blank(),
        # change facet label background and border
        strip.background = element_blank())

custom_palette <- viridis::viridis_pal(direction = -1)(100)
custom_palette[1] <- "white"

plot_country <- ggplot(df_country) +
  geom_sf(aes(fill = n_study, geometry = geometry)) +
  facet_grid(sex~sti) +
  scale_fill_gradientn(colors = custom_palette, breaks = c(1,10,20,30), na.value="lightgrey") +
  theme_minimal(base_size=7) +
  my_theme +
  labs(fill="Number of studies")

plot_country

ggsave("./plots/fig_supp_map.png", plot_country, width = 18, units="cm", dpi=700)
