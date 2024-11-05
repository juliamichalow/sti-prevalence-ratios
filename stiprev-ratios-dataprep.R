## STI PREVALENCE RATIOS
## PREP DATA AND ADJUST FOR DIAGNOSTIC ACCURACY

## Load packages
library(tidyverse)
library(rstan)
library(readxl)
library(writexl)

## Load data 
# df_wide <- read_xlsx("./../Data/final-dataset-v5.xlsx") |>
#   separate(study_id, into=c("author","year_publish"), sep = "(?<=\\D)\\s(?=\\d+[a-zA-Z]*$)", remove = FALSE) |>
#   select(!author) |>
#   mutate(sex = str_to_lower(sex),
#          year_publish = case_when(year_publish %in% c("2022a", "2022b") ~ 2022,
#                                   TRUE ~ as.numeric(year_publish)),
#          country = str_to_title(country)) |>
#   # exclude studies marked for exclusion (duplicate articles, etc)
#   filter(is.na(exclude))

df_wide <- read_xlsx("./../../Data/final-dataset-v7.xlsx")

# Put into long format
df_long <- df_wide |> 
  mutate(across(c(starts_with("ct_"), starts_with("ng_"), starts_with("tr_")), as.character)) |> 
  pivot_longer(cols=c(starts_with("ct_"),starts_with("ng_"),starts_with("tr_"))) |>
  separate(name, into=c("sti","variable")) |> 
  pivot_wider(names_from = "variable", values_from="value") |> 
  # drop rows with NA values across all selected columns (due to pivots)
  filter(!if_all(tested:testcat,is.na)) |> 
  mutate(sti = str_to_upper(sti)) |>
  mutate(across(c(positive,tested),as.integer),
         prevalence = as.numeric(prevalence)) |> 
  # remove non-essential columns for the analysis
  # select(!c(month_start, month_end, 
  #           starts_with("geography_"),starts_with("site_"),
  #           starts_with("name_"),starts_with("ownership_"),
  #           study_design, study_sample, n_participants, exclude, exclude_reason,
  #           specimen, test, datatype)) |>
  # reorder columns
  relocate(c(notes, search_v), .after = last_col()) 

# Join performance data

df_perform <- read_xlsx("./../../Data/diagnostic-test-performance.xlsx") |>
  select(sti,testcat,specimen,sex,sens,spec) |>
  mutate(sex = str_to_sentence(sex))

df <- df_long |>
  # drop rows with NR for testcat 
  filter(!(testcat == "NR")) |> 
  # drop results with chlamydia serology
  filter(!(sti=="CT" & specimencat == "blood")) |>
  # drop results with gram stain for women - not recommended
  filter(!(sti=="NG" & sex == "female" & testcat == "gram stain")) |>
  # join performance data
  left_join(df_perform, by = c("sti", "specimencat"="specimen", "testcat", "sex"))

# Confirm there are no empty joins
df |>
  filter(is.na(sens)) |>
  select(cov_id,study_id,sex,sti,specimencat,testcat) |>
  print(n=50)

# Adjust for sensitivity and specificity
# stan model
model <- "
data {
  int<lower=0> N;
  int<lower=0> sample_size[N];
  int<lower=0> x[N];
  vector<lower = 0, upper=1>[N] sens;
  vector<lower = 0, upper=1>[N] spec;
}
parameters {
  vector<lower=0, upper=1>[N] prev_true;
}
transformed parameters {
  vector<lower=0, upper=1>[N] prev_obs;
  prev_obs = prev_true .* sens + (1.0 - prev_true) .* (1.0 - spec);
}
model {
  x ~ binomial(sample_size, prev_obs);
}
"
stan_data <- list(N = nrow(df),
                  x = df$positive,
                  sample_size = df$tested,
                  sens = df$sens / 100.0,
                  spec = df$spec / 100.0)

fit <- stan(model_code = model, data = stan_data)

prev_true_out <- summary(fit, pars = "prev_true")$summary

prev_true_out <- prev_true_out[ , c("mean", "sd", "2.5%", "97.5%")] |>
  as.data.frame()  

df_out <- df |> bind_cols(prev_true_out) |>
  rename(adj_prev = mean, adj_se = sd, adj_prev_lwr = `2.5%`, adj_prev_upr = `97.5%`) |>
  # calculate new numerator and new denominator as a function of adj_prev and adj_se
  mutate(adj_denom = adj_prev*(1-adj_prev)/(adj_se^2),
         adj_num = adj_prev*adj_denom)

str(df_out)

# Reorder
df_final <- df_out |>
  relocate(c(sens, spec, starts_with("adj_")), .after = "testcat")

# Difference between year_publish and year_mid
df_final |>
  mutate(year_mid = as.numeric(year_mid)) |>
  filter(!is.na(year_mid)) |>
  mutate(year_diff = year_publish - year_mid) |>
  select(study_id, country, year_publish, year_mid, year_diff) |>
  unique() |>
  group_by() |>
  summarise(min = min(year_diff),
            max = max(year_diff),
            mean = mean(year_diff),
            median = median(year_diff))

df_final |>
  mutate(year_mid = as.numeric(year_mid)) |>
  filter(!is.na(year_mid)) |>
  mutate(year_diff = year_publish - year_mid) |>
  select(study_id, year_publish, year_mid, year_diff) |>
  unique() |>
  group_by() |>
  summarise(min = min(year_diff),
            max = max(year_diff),
            mean = mean(year_diff),
            median = median(year_diff))

# Replace year_mid for studies with missing data
df_save <- df_final |>
  mutate(year_mid = as.numeric(year_mid)) |>
  mutate(year_mid = case_when(is.na(year_mid) ~ year_publish - 3,
                              TRUE ~ year_mid)) |>
  # Remove studies with year_mid < 2000
  filter(!year_mid < 2000) |>
  # Rename variables
  mutate(sti = case_when(sti == "TR" ~ "TV", TRUE ~ sti)) |>
  rename(num = positive, denom = tested)

# save
write.csv(df_save,"./../../Data/final-dataset-v7-adjusted.csv", row.names = FALSE)
write.csv(df_save,"./data/final-dataset-v7-adjusted.csv", row.names = FALSE)

# visualise prev vs adj_prev

df_out |>
  ggplot(aes(prevalence, adj_prev, color = testcat)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~sti) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(labels=scales::percent_format()) +
  scale_y_continuous(labels=scales::percent_format()) +
  theme_bw() +
  labs(x = "Prevalence", y = "Adjusted prevalence", colour = "Test category")

