## STI PREVALENCE RATIOS

library(tidyverse)
library(glmmTMB)
library(sjPlot)

# PREP DATA ----

# Population by region, 15-49 years in 2020
df_pop <- read.csv("./data/unpopulation_dataportal_20241203105103.csv") |>
  rename_all(tolower) |>
  filter(time == 2020,
         sex %in% c("Female", "Male")) |>
  select(time, location, iso3, sex, age, value) |>
  mutate(region = case_when(
    iso3 %in% c("BWA","SWZ","LSO","NAM","ZAF") ~ "SA",
    iso3 %in% c("AGO","CMR","CAF","TCD","COG","COD","GNQ","GAB","STP") ~ "CA",
    iso3 %in% c("BEN","BFA","CPV","CIV","GMB","GHA","GIN","GNB",
                "LBR","MLI","MRT","NER","NGA","SEN","SLE","TGO") ~ "WA",
    iso3 %in% c("BDI","COM","DJI","ERI","ETH","KEN","MDG","MWI",
                "MUS","MOZ","RWA","SYC","SOM","SSD",
                "UGA","TZA","TZA","ZMB","ZWE") ~ "EA")) |>
  mutate(region = fct_collapse(region, 
                               "WCA" = c("CA", "WA")),
         region = factor(region,
                         levels = c("WCA", "EA", "SA"))) |>
  filter(!is.na(region)) |>
  group_by(region, sex) |>
  summarise(population = sum(value)) |>
  ungroup()

# Study data
df <- read.csv("./data/final-dataset-v8-adjusted.csv") |>
  # drop Both sexes
  filter(!sex == "Both sexes") |>
  # minimum sample size of 15
  filter(denom >= 15) |>
  mutate(year_regression = year_mid - 2012,
         sti = factor(sti, levels = c("CT", "NG", "TV")),
         year_grp = case_when(year_mid %in% c(2000:2004) ~ "2000-2004",
                              year_mid %in% c(2005:2009) ~ "2005-2009",
                              year_mid %in% c(2010:2014) ~ "2010-2014",
                              year_mid %in% c(2015:2019) ~ "2015-2019",
                              year_mid %in% c(2020:2022) ~ "2020-2024"),
         year_grp = factor(year_grp, 
                           levels = c("2000-2004", "2005-2009", "2010-2014", 
                                      "2015-2019", "2020-2024")),
         region = fct_collapse(region_analysis, "WCA" = c("WA", "CA")),
         region = factor(region, levels = c("SA", "EA", "WCA")),
         region_actual = factor(region_actual, levels = c("WA", "CA", "EA", "SA", "Multiple")),
         sex = factor(sex, levels = c("Female", "Male")),
         hiv_status = fct_relevel(hiv_status, "Mixed"),
         age_group = factor(age_group, levels = c("Adult", "Youth")),
         population = fct_collapse(population, 
                                   "ANC attendees" = c("ANC attendees", "ANC and FP attendees",
                                                       "ANC and GYN attendees"),
                                   "FP attendees" = c("FP attendees", "FP and GYN attendees")),
         population = factor(population,
                             levels = c("ANC attendees", "FP attendees", 
                                        "GYN attendees", "PHC/OPD attendees",
                                        "Students","Community members", 
                                        "HIV/STI prevention trial participants",
                                        "Population-representative survey participants")),
         rob_study_sample = factor(case_when(study_sample %in% c("Consecutive", "Convenience", "Respondent-driven sampling", "NR") ~ "Higher",
                                             study_sample %in% c("Random") ~ "Lower"),
                                   levels = c("Lower", "Higher")),
         rob_participant = factor(case_when(age_range == "NR" | (is.na(age_mean) & is.na(age_median)) | hiv_prevalence == "NR" | is.na(location) ~ "Higher",
                                            TRUE ~ "Lower"),
                                  levels = c("Lower", "Higher")),
         rob_measurement = factor(case_when(test %in% c("DFA and NAAT", "culture or NAAT", "NAAT and WM", "culture and WM") |
                                              specimen %in% c("genital fluid and urine", "genital fluid or urine") ~ "Higher",
                                            TRUE ~ "Lower"),
                                  levels = c("Lower", "Higher")),
         rob_precision = factor(case_when(denom < 100 ~ "Higher", denom >= 100 ~ "Lower"),
                                levels = c("Lower", "Higher")),
         test = fct_collapse(test, 
                             "DFA" = c("DFA", "DFA and NAAT"),
                             "Culture" = c("culture", "culture or NAAT"),
                             "Wet mount" = c("WM", "NAAT and WM", "culture and WM"),
                             "Rapid antigen test" = c("rapid antigen test")),
         test = factor(test, levels = c("NAAT","Culture", "DFA", "Rapid antigen test", "ELISA", "Wet mount")))


# FUNCTION ---- 

# Function that predicts prevalence for:
# 1. Region by year and sex
# 2. SSA mean by year and sex
# 3. SSA mean male-to-female prevalence ratio

prediction <- function(mod, sti, pop_pred = "ANC attendees") {
  
  df_pred <- crossing(sex = c("Female", "Male"),
                      region = c("WCA", "EA", "SA"),
                      year = 2000:2024) %>%
    mutate(year_regression = year - 2012,
           population = pop_pred,
           age_group = "Adult",
           hiv_status = "Mixed",
           test = "NAAT",
           study_id = NA)
  
  ## 1. Region prevalence
  
  pred <- predict(mod, newdata = df_pred, type="link", se.fit = TRUE, allow.new.levels=TRUE)
  df_region_prev <- df_pred |>
    mutate(log_prev = pred$fit, log_se = pred$se.fit,
           prev = exp(log_prev),
           lwr = exp(log_prev - 1.96 * log_se),
           upr = exp(log_prev + 1.96 * log_se),
           sti = sti) |>
    select(sti, sex, region, year, prev, lwr, upr)
  
  ## 2. SSA mean prevalence
  
  df_pred$log_prev <- predict(mod, newdata = df_pred)
  # covariance matrix for predictions (log scale)
  var_log_prev <- predict(mod, newdata = df_pred, se.fit = TRUE, cov.fit = TRUE)$cov.fit
  
  # Convert log_prev to prev
  df_pred$prev <- exp(df_pred$log_prev)
  
  # Delta method for variance of predictions (exponentiated scale)
  var_prev <- diag(df_pred$prev) %*% var_log_prev %*% diag(df_pred$prev)
  
  # Calculate weight per region and sex
  wt <- df_pop |>
    group_by(sex) |>
    mutate(wt = population/sum(population)) |>
    select(!population)
  
  # Create matrix of aggregated output wanted
  df_out <- distinct(df_pred, sex, year) %>%
    mutate(out_idx = row_number())
  
  # Assign row indices to the prediction data frame and
  # link them to the indices in the output data frame ->
  # Create the matrix of weights to aggregate
  df_join <- df_pred %>%
    select(sex, region, year) %>%
    mutate(mod_idx = row_number()) %>%
    left_join(wt, by = join_by(sex, region)) %>%
    inner_join(df_out, by = join_by(sex, year))
  
  # Note: do this with a sparse matrix if the problem becomes large
  Amat <- matrix(0, nrow = nrow(df_out), ncol = nrow(df_pred))
  Amat[cbind(df_join$out_idx, df_join$mod_idx)] <- df_join$wt
  
  # Aggregate prevalence
  df_out$prev <- as.numeric(Amat %*% df_pred$prev)
  var_prev_out <- Amat %*% var_prev %*% t(Amat) 
  
  # Convert to log prevalence
  df_out$log_prev <- log(df_out$prev)
  var_log_prev_out <- diag(1 / df_out$prev) %*% var_prev_out %*% diag(1 / df_out$prev)
  
  # Standard error
  df_out$se <- sqrt(diag(var_prev_out))
  
  # Final aggregate prediction and CI
  df_ssa_prev <- df_out |>
    mutate(lwr = prev - qnorm(0.975) * se,
           upr = prev + qnorm(0.975) * se,
           sti = sti) |>
    select(sti, sex, year, prev, lwr, upr)
  
  ## 3. SSA mean ratio
  
  # Calculate ratio
  df_ratio <- df_out %>%
    pivot_wider(id_cols = year, names_from = sex, values_from = out_idx) %>%
    mutate(ratio_idx = row_number())
  
  df_ratio$log_ratio <- df_out$log_prev[df_ratio$Male] - df_out$log_prev[df_ratio$Female]
  
  Aratio_mat <- matrix(0, nrow = nrow(df_ratio), ncol = nrow(df_out))
  Aratio_mat[cbind(df_ratio$ratio_idx, df_ratio$Male)] <- 1
  Aratio_mat[cbind(df_ratio$ratio_idx, df_ratio$Female)] <- -1
  
  df_ratio$log_ratio_check <- as.numeric(Aratio_mat %*% df_out$log_prev)
  
  var_log_ratio <- Aratio_mat %*% var_log_prev_out %*% t(Aratio_mat)
  
  df_ratio <- df_ratio %>%
    mutate(se_log_ratio = sqrt(diag(var_log_ratio)),
           lwr_log_ratio = log_ratio - qnorm(0.975) * se_log_ratio,
           upr_log_ratio = log_ratio + qnorm(0.975) * se_log_ratio,
           ratio = exp(log_ratio),
           se_ratio = ratio * se_log_ratio,
           lwr_ratio = exp(lwr_log_ratio),
           upr_ratio = exp(upr_log_ratio),
           sti = sti)
  
  df_ssa_ratio <- df_ratio |>
    select(sti, year, ratio, lwr_ratio, upr_ratio) |>
    rename(lwr = lwr_ratio, upr = upr_ratio) |>
    mutate(tau_study = insight::get_variance(mod)$var.random)
  
  return(list(region_prev = df_region_prev, ssa_prev = df_ssa_prev, ssa_ratio = df_ssa_ratio))
  
}

# Function that predicts:
# SSA mean prevalence ratio for time trend (mean aPR per year)

prediction_2 <- function(mod, sti, y1, y2, pop_pred = "ANC attendees") {
  
  df_pred <- crossing(sex = c("Female", "Male"),
                      region = c("WCA", "EA", "SA"),
                      year = 2000:2024) %>%
    mutate(year_regression = year - 2012,
           population = pop_pred,
           age_group = "Adult",
           hiv_status = "Mixed",
           test = "NAAT",
           study_id = NA)
  
  # predict prevalence
  df_pred$log_prev <- predict(mod, newdata = df_pred)
  df_pred$prev <- exp(df_pred$log_prev)
  
  # covariance matrix for predictions (log scale)
  var_log_prev <- predict(mod, newdata = df_pred, se.fit = TRUE, cov.fit = TRUE)$cov.fit
  
  # Delta method for variance of predictions (exponentiated scale)
  var_prev <- diag(df_pred$prev) %*% var_log_prev %*% diag(df_pred$prev)
  
  # Calculate weight per region and sex
  wt <- df_pop |>
    group_by(sex) |>
    mutate(wt = population/sum(population)) |>
    select(!population)
  
  # Create matrix of aggregated output wanted
  df_out <- distinct(df_pred, year) %>%
    mutate(out_idx = row_number())
  
  # Assign row indices to the prediction data frame and
  # link them to the indices in the output data frame ->
  # Create the matrix of weights to aggregate
  df_join <- df_pred %>%
    select(sex, region, year) %>%
    mutate(mod_idx = row_number()) %>%
    left_join(wt, by = join_by(sex, region)) %>%
    inner_join(df_out, by = join_by(year))
  
  # Note: do this with a sparse matrix if the problem becomes large
  Amat <- matrix(0, nrow = nrow(df_out), ncol = nrow(df_pred))
  Amat[cbind(df_join$out_idx, df_join$mod_idx)] <- df_join$wt
  
  # Aggregate prevalence
  df_out$prev <- as.numeric(Amat %*% df_pred$prev)
  var_prev_out <- Amat %*% var_prev %*% t(Amat) 
  
  # Convert to log prevalence
  df_out$log_prev <- log(df_out$prev)
  var_log_prev_out <- diag(1 / df_out$prev) %*% var_prev_out %*% diag(1 / df_out$prev)
  
  # Standard error
  df_out$se <- sqrt(diag(var_prev_out))
  
  # Calculate ratio
  df_ratio <- df_out %>%
    select(year, out_idx) |>
    pivot_wider(names_from = year, values_from = out_idx) %>%
    mutate(ratio_idx = row_number())
  
  # Use y2 and y1 as specified for the ratio
  df_ratio$log_ratio <- df_out$log_prev[df_ratio[[y2]]] - df_out$log_prev[df_ratio[[y1]]]
  
  Aratio_mat <- matrix(0, nrow = nrow(df_ratio), ncol = nrow(df_out))
  Aratio_mat[cbind(df_ratio$ratio_idx, df_ratio[[y2]])] <- 1
  Aratio_mat[cbind(df_ratio$ratio_idx, df_ratio[[y1]])] <- -1
  
  df_ratio$log_ratio_check <- as.numeric(Aratio_mat %*% df_out$log_prev)
  
  var_log_ratio <- Aratio_mat %*% var_log_prev_out %*% t(Aratio_mat)
  
  df_ratio <- df_ratio %>%
    select(c(log_ratio)) |>
    mutate(se_log_ratio = sqrt(diag(var_log_ratio)),
           lwr_log_ratio = log_ratio - qnorm(0.975) * se_log_ratio,
           upr_log_ratio = log_ratio + qnorm(0.975) * se_log_ratio,
           ratio = exp(log_ratio),
           se_ratio = ratio * se_log_ratio,
           lwr_ratio = exp(lwr_log_ratio),
           upr_ratio = exp(upr_log_ratio),
           sti = sti,
           years_for_ratio = paste0(y1," & ",y2)) |>
    mutate(est = paste0(sprintf(ratio,fmt = '%#3f')," (", sprintf(lwr_ratio,fmt = '%#.3f'), 
                        "-", sprintf(upr_ratio,fmt = '%#.3f'), ")"))
  
  df_ratio |>
    select(sti, years_for_ratio, est)
  
}

# BETWEEN STUDY ANALYSIS ----

# Exclude Mullick 2023 youth data point from meta-reg (keep only for within study ratios)
df_ct <- df |> filter(sti == "CT") |>
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Youth"))
df_ng <- df |> filter(sti == "NG") |>
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Youth"))
df_tv <- df |> filter(sti == "TV") |>
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Youth"))

# < Main model, with test adjustment ----

# Formula
form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  sex:region + population + age_group + hiv_status + test + (1 | study_id)

# Logit models
modct1 <- glmmTMB(form, data = df_ct, family = binomial(link="logit"))
modng1 <- glmmTMB(form, data = df_ng, family = binomial(link="logit"))
modtv1 <- glmmTMB(form, data = df_tv, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct1)
theta_ct <- as.vector(modct1$fit$par[names(modct1$fit$par) == "theta"])
fixed_ct <- as.vector(modct1$fit$par[names(modct1$fit$par) == "beta"])
#fixed_ct <- append(fixed_ct, 0, length(fixed_ct))

summary(modng1)
theta_ng <- as.vector(modng1$fit$par[names(modng1$fit$par) == "theta"])
fixed_ng <- as.vector(modng1$fit$par[names(modng1$fit$par) == "beta"])
#fixed_ng <- append(fixed_ng, 0, length(fixed_ng))

summary(modtv1)
theta_tv <- as.vector(modtv1$fit$par[names(modtv1$fit$par) == "theta"])
fixed_tv <- as.vector(modtv1$fit$par[names(modtv1$fit$par) == "beta"])
#fixed_tv <- append(fixed_tv, 0, length(fixed_tv))

# Log models
modct2 <- glmmTMB(form, data = df_ct, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng2 <- glmmTMB(form, data = df_ng, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv2 <- glmmTMB(form, data = df_tv, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

sjPlot::tab_model(modct1, modct2, modng1, modng2, modtv1, modtv2, p.style = "stars", wrap.labels = 100,
                  dv.labels = c("CT logit", "CT log", "NG logit", "NG log", "TV logit", "TV log"))
sjPlot::tab_model(modct2, modng2, modtv2, p.style = "stars", wrap.labels = 100)

# Predict by region
dat_region_annual <- rbind(
  prediction(modct2, "CT")$region_prev,
  prediction(modng2, "NG")$region_prev,
  prediction(modtv2, "TV")$region_prev)

dat_region_annual |> write.csv("./results/prevalence_alldata_adjusted.csv", row.names = FALSE)

# Predict weighted mean for SSA 
dat_ssa_annual <- rbind(
  prediction(modct2, "CT")$ssa_prev,
  prediction(modng2, "NG")$ssa_prev,
  prediction(modtv2, "TV")$ssa_prev) |>
  mutate(region = "SSA") |>
  relocate(sti, sex, region)

# Save 2020 predictions
rbind(dat_region_annual |> filter(year == 2020),
      dat_ssa_annual |> filter(year == 2020)) |>
  write.csv("./results/prevalence_2020_adjusted.csv", row.names = FALSE)

# Prediction for 2010
dat_ssa_annual |> filter(year == 2010)

# Compare 2005 and 2020
dat_ssa_annual |> filter(year %in% c(2005, 2020)) |> 
  select(!c(lwr,upr)) |> #
  pivot_wider(names_from = year, values_from = prev) |> 
  mutate(ratio = `2020`/`2005`)

# Compare 2010 and 2020
dat_ssa_annual |> filter(year %in% c(2010, 2020)) |> 
  select(!c(lwr,upr)) |> #
  pivot_wider(names_from = year, values_from = prev) |> 
  mutate(ratio = `2020`/`2010`)


# Prediction for 2012, 2016, 2020 (supp. table for comparison with WHO)
rbind(dat_region_annual, dat_ssa_annual) |>
  filter(year %in% c(2012, 2016, 2020)) |>
  mutate(prev = paste0(sprintf(prev*100,fmt = '%#.1f')," (", sprintf(lwr*100,fmt = '%#.1f'), 
                       "-", sprintf(upr*100,fmt = '%#.1f'), ")")) |>
  select(sti, region, sex, year, prev) |>
  pivot_wider(names_from = c(sti,sex), values_from = prev) |>
  mutate(region = factor(region, levels = c("WCA", "EA", "SA", "SSA"))) |>
  arrange(year, region) |>
  write.csv("./results/prevalence_comparewho.csv", row.names = FALSE)

# Ratio between study
ratio_btwn_adj <- rbind(
  prediction(modct2, "CT")$ssa_ratio,
  prediction(modng2, "NG")$ssa_ratio,
  prediction(modtv2, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Between study",
         diag = "All tests, adjusted",
         n = c(modct2$modelInfo$nobs, modng2$modelInfo$nobs, modtv2$modelInfo$nobs))

# Mean aPR per year
prediction_2(modct2, "CT", y1 = "2015", y2 = "2016")
prediction_2(modng2, "NG", y1 = "2015", y2 = "2016")
prediction_2(modtv2, "TV", y1 = "2015", y2 = "2016")

# Mean aPR between 2005 and 2020
prediction_2(modct2, "CT", y1 = "2005", y2 = "2020")
prediction_2(modng2, "NG", y1 = "2005", y2 = "2020")
prediction_2(modtv2, "TV", y1 = "2005", y2 = "2020")

# Mean aPR between 2010 and 2020
prediction_2(modct2, "CT", y1 = "2010", y2 = "2020")
prediction_2(modng2, "NG", y1 = "2010", y2 = "2020")
prediction_2(modtv2, "TV", y1 = "2010", y2 = "2020")

# < Main model, no test adjustment ----

# Formula
form <- cbind(num,(denom-num)) ~ region + year_regression:region + sex:region +
  population + age_group + hiv_status + test + (1 | study_id)

# Logit models
modct3 <- glmmTMB(form, data = df_ct, family = binomial(link="logit"))
modng3 <- glmmTMB(form, data = df_ng, family = binomial(link="logit"))
modtv3 <- glmmTMB(form, data = df_tv, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct3)
theta_ct <- as.vector(modct3$fit$par[names(modct3$fit$par) == "theta"])
fixed_ct <- as.vector(modct3$fit$par[names(modct3$fit$par) == "beta"])
#fixed_ct <- append(fixed_ct, 0, length(fixed_ct))

summary(modng3)
theta_ng <- as.vector(modng3$fit$par[names(modng3$fit$par) == "theta"])
fixed_ng <- as.vector(modng3$fit$par[names(modng3$fit$par) == "beta"])
#fixed_ng <- append(fixed_ng, 0, length(fixed_ng))

summary(modtv3)
theta_tv <- as.vector(modtv3$fit$par[names(modtv3$fit$par) == "theta"])
fixed_tv <- as.vector(modtv3$fit$par[names(modtv3$fit$par) == "beta"])
#fixed_tv <- append(fixed_tv, 0, length(fixed_tv)-1)

# Log models

modct4 <- glmmTMB(form, data = df_ct, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng4 <- glmmTMB(form, data = df_ng, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv4 <- glmmTMB(form, data = df_tv, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

sjPlot::tab_model(modct4, modng4, modtv4, p.style = "stars", wrap.labels = 100)

# Predict by region
dat_region_annual <- rbind(
  prediction(modct4, "CT")$region_prev,
  prediction(modng4, "NG")$region_prev,
  prediction(modtv4, "TV")$region_prev)

dat_region_annual |> write.csv("./results/prevalence_alldata_unadjusted.csv", row.names = FALSE)

# Predict weighted mean for SSA 
dat_ssa_annual <- rbind(
  prediction(modct4, "CT")$ssa_prev,
  prediction(modng4, "NG")$ssa_prev,
  prediction(modtv4, "TV")$ssa_prev) |>
  mutate(region = "SSA") |>
  relocate(sti, sex, region)

# Save 2020 predictions
rbind(dat_region_annual |> filter(year == 2020),
      dat_ssa_annual |> filter(year == 2020)) |>
  write.csv("./results/prevalence_2020_unadjusted.csv", row.names = FALSE)

dat_ssa_annual |> filter(year == 2020)

# Ratios
ratio_btwn_unadj <- rbind(
  prediction(modct4, "CT")$ssa_ratio,
  prediction(modng4, "NG")$ssa_ratio,
  prediction(modtv4, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Between study",
         diag = "All tests, unadjusted",
         n = c(modct4$modelInfo$nobs, modng4$modelInfo$nobs, modtv4$modelInfo$nobs))

# < Women only ----

df_ct_fm <- df |> filter(sti == "CT", sex == "Female")
df_ng_fm <- df |> filter(sti == "NG", sex == "Female")
df_tv_fm <- df |> filter(sti == "TV", sex == "Female")

# Formula
form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  population + age_group + hiv_status + test + (1 | study_id)

# Logit models
modct5 <- glmmTMB(form, data = df_ct_fm, family = binomial(link="logit"))
modng5 <- glmmTMB(form, data = df_ng_fm, family = binomial(link="logit"))
modtv5 <- glmmTMB(form, data = df_tv_fm, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct5)
theta_ct <- as.vector(modct5$fit$par[names(modct5$fit$par) == "theta"])
fixed_ct <- as.vector(modct5$fit$par[names(modct5$fit$par) == "beta"])

summary(modng5)
theta_ng <- as.vector(modng5$fit$par[names(modng5$fit$par) == "theta"])
fixed_ng <- as.vector(modng5$fit$par[names(modng5$fit$par) == "beta"])

summary(modtv5)
theta_tv <- as.vector(modtv5$fit$par[names(modtv5$fit$par) == "theta"])
fixed_tv <- as.vector(modtv5$fit$par[names(modtv5$fit$par) == "beta"])

# Log models
modct6 <- glmmTMB(form, data = df_ct_fm, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng6 <- glmmTMB(form, data = df_ng_fm, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv6 <- glmmTMB(form, data = df_tv_fm, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

sjPlot::tab_model(modct6, modng6, modtv6, p.style = "stars", wrap.labels = 100)

# Predict by region
dat_region_annual <- rbind(
  prediction(modct6, "CT")$region_prev,
  prediction(modng6, "NG")$region_prev,
  prediction(modtv6, "TV")$region_prev) |>
  filter(sex == "Female")

dat_region_annual |> write.csv("./results/prevalence_alldata_fm.csv", row.names = FALSE)

# Predict weighted mean for SSA 
dat_ssa_annual <- rbind(
  prediction(modct6, "CT")$ssa_prev,
  prediction(modng6, "NG")$ssa_prev,
  prediction(modtv6, "TV")$ssa_prev) |>
  mutate(region = "SSA") |>
  relocate(sti, sex, region) |>
  filter(sex == "Female")

# Save 2020 predictions
rbind(dat_region_annual |> filter(year == 2020),
      dat_ssa_annual |> filter(year == 2020)) |>
  write.csv("./results/prevalence_2020_fm.csv", row.names = FALSE)

# < ANC only ----

df_ct_anc <- df |> filter(sti == "CT", sex == "Female", population == "ANC attendees")
df_ng_anc <- df |> filter(sti == "NG", sex == "Female", population == "ANC attendees")
df_tv_anc <- df |> filter(sti == "TV", sex == "Female", population == "ANC attendees")

# Formula
form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  age_group + hiv_status + test + (1 | study_id)

# Logit models
modct7 <- glmmTMB(form, data = df_ct_anc, family = binomial(link="logit"))
modng7 <- glmmTMB(form, data = df_ng_anc, family = binomial(link="logit"))
modtv7 <- glmmTMB(form, data = df_tv_anc, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct7)
theta_ct <- as.vector(modct7$fit$par[names(modct7$fit$par) == "theta"])
fixed_ct <- as.vector(modct7$fit$par[names(modct7$fit$par) == "beta"])

summary(modng7)
theta_ng <- as.vector(modng7$fit$par[names(modng7$fit$par) == "theta"])
fixed_ng <- as.vector(modng7$fit$par[names(modng7$fit$par) == "beta"])

summary(modtv7)
theta_tv <- as.vector(modtv7$fit$par[names(modtv7$fit$par) == "theta"])
fixed_tv <- as.vector(modtv7$fit$par[names(modtv7$fit$par) == "beta"])

# Log models
modct8 <- glmmTMB(form, data = df_ct_anc, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng8 <- glmmTMB(form, data = df_ng_anc, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv8 <- glmmTMB(form, data = df_tv_anc, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

sjPlot::tab_model(modct8, modng8, modtv8, p.style = "stars", wrap.labels = 100)

# Predict by region
dat_region_annual <- rbind(
  prediction(modct8, "CT")$region_prev,
  prediction(modng8, "NG")$region_prev,
  prediction(modtv8, "TV")$region_prev) |>
  filter(sex == "Female")

dat_region_annual |> write.csv("./results/prevalence_alldata_anc.csv", row.names = FALSE)

# Predict weighted mean for SSA 
dat_ssa_annual <- rbind(
  prediction(modct8, "CT")$ssa_prev,
  prediction(modng8, "NG")$ssa_prev,
  prediction(modtv8, "TV")$ssa_prev) |>
  mutate(region = "SSA") |>
  relocate(sti, sex, region) |>
  filter(sex == "Female")

# Save 2020 predictions
rbind(dat_region_annual |> filter(year == 2020),
      dat_ssa_annual |> filter(year == 2020)) |>
  write.csv("./results/prevalence_2020_anc.csv", row.names = FALSE)

# < NAAT only, with test adjustment ----

df_ct_naat <- df |> filter(sti == "CT", testcat == "NAAT")
df_ng_naat <- df |> filter(sti == "NG", testcat == "NAAT")
df_tv_naat <- df |> filter(sti == "TV", testcat == "NAAT")

# Formula
form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  sex:region + population + age_group + hiv_status + (1 | study_id)

# Logit models
modct9 <- glmmTMB(form, data = df_ct_naat, family = binomial(link="logit"))
modng9 <- glmmTMB(form, data = df_ng_naat, family = binomial(link="logit"))
modtv9 <- glmmTMB(form, data = df_tv_naat, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct9)
theta_ct <- as.vector(modct9$fit$par[names(modct9$fit$par) == "theta"])
fixed_ct <- as.vector(modct9$fit$par[names(modct9$fit$par) == "beta"])
#fixed_ct <- append(fixed_ct, 0, length(fixed_ct))

summary(modng9)
theta_ng <- as.vector(modng9$fit$par[names(modng9$fit$par) == "theta"])
fixed_ng <- as.vector(modng9$fit$par[names(modng9$fit$par) == "beta"])
#fixed_ng <- append(fixed_ng, 0, length(fixed_ng))

summary(modtv9)
theta_tv <- as.vector(modtv9$fit$par[names(modtv9$fit$par) == "theta"])
fixed_tv <- as.vector(modtv9$fit$par[names(modtv9$fit$par) == "beta"])
#fixed_tv <- append(fixed_tv, 0, length(fixed_tv))

# Log models
modct10 <- glmmTMB(form, data = df_ct_naat, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng10 <- glmmTMB(form, data = df_ng_naat, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv10 <- glmmTMB(form, data = df_tv_naat, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

sjPlot::tab_model(modct10, modng10, modtv10, p.style = "stars", wrap.labels = 100)

# Predict by region
dat_region_annual <- rbind(
  prediction(modct10, "CT")$region_prev,
  prediction(modng10, "NG")$region_prev,
  prediction(modtv10, "TV")$region_prev)

dat_region_annual |> write.csv("./results/prevalence_alldata_naat_adjusted.csv", row.names = FALSE)

# Predict weighted mean for SSA 
dat_ssa_annual <- rbind(
  prediction(modct10, "CT")$ssa_prev,
  prediction(modng10, "NG")$ssa_prev,
  prediction(modtv10, "TV")$ssa_prev) |>
  mutate(region = "SSA") |>
  relocate(sti, sex, region)

# Save 2020 predictions
rbind(dat_region_annual |> filter(year == 2020),
      dat_ssa_annual |> filter(year == 2020)) |>
  write.csv("./results/prevalence_2020_naat_adjusted.csv", row.names = FALSE)

# Ratios
ratio_btwn_naat_adj <- rbind(
  prediction(modct10, "CT")$ssa_ratio,
  prediction(modng10, "NG")$ssa_ratio,
  prediction(modtv10, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Between study",
         diag = "NAAT, adjusted",
         n = c(modct10$modelInfo$nobs, modng10$modelInfo$nobs, modtv10$modelInfo$nobs))

# < NAAT only, no test adjustment ----

df_ct_naat <- df |> filter(sti == "CT", testcat == "NAAT")
df_ng_naat <- df |> filter(sti == "NG", testcat == "NAAT")
df_tv_naat <- df |> filter(sti == "TV", testcat == "NAAT")

# Formula
form <- cbind(num,(denom-num)) ~ region + year_regression:region + sex:region +
  population + age_group + hiv_status + (1 | study_id)

# Logit models
modct11 <- glmmTMB(form, data = df_ct_naat, family = binomial(link="logit"))
modng11 <- glmmTMB(form, data = df_ng_naat, family = binomial(link="logit"))
modtv11 <- glmmTMB(form, data = df_tv_naat, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct11)
theta_ct <- as.vector(modct11$fit$par[names(modct11$fit$par) == "theta"])
fixed_ct <- as.vector(modct11$fit$par[names(modct11$fit$par) == "beta"])
#fixed_ct <- append(fixed_ct, 0, length(fixed_ct))

summary(modng11)
theta_ng <- as.vector(modng11$fit$par[names(modng11$fit$par) == "theta"])
fixed_ng <- as.vector(modng11$fit$par[names(modng11$fit$par) == "beta"])
#fixed_ng <- append(fixed_ng, 0, length(fixed_ng))

summary(modtv11)
theta_tv <- as.vector(modtv11$fit$par[names(modtv11$fit$par) == "theta"])
fixed_tv <- as.vector(modtv11$fit$par[names(modtv11$fit$par) == "beta"])
#fixed_tv <- append(fixed_tv, 0, length(fixed_tv))

# Log models
modct12 <- glmmTMB(form, data = df_ct_naat, family = binomial(link="log"),
                   start = list(theta = theta_ct, beta = fixed_ct))

modng12 <- glmmTMB(form, data = df_ng_naat, family = binomial(link="log"),
                   start = list(theta = theta_ng, beta = fixed_ng))

modtv12 <- glmmTMB(form, data = df_tv_naat, family = binomial(link="log"),
                   start = list(theta = theta_tv, beta = fixed_tv))

sjPlot::tab_model(modct12, modng12, modtv12, p.style = "stars", wrap.labels = 100)

# Predict by region
dat_region_annual <- rbind(
  prediction(modct12, "CT")$region_prev,
  prediction(modng12, "NG")$region_prev,
  prediction(modtv12, "TV")$region_prev)

dat_region_annual |> write.csv("./results/prevalence_alldata_naat_unadjusted.csv", row.names = FALSE)

# Predict weighted mean for SSA 
dat_ssa_annual <- rbind(
  prediction(modct12, "CT")$ssa_prev,
  prediction(modng12, "NG")$ssa_prev,
  prediction(modtv12, "TV")$ssa_prev) |>
  mutate(region = "SSA") |>
  relocate(sti, sex, region)

# Save 2020 predictions
rbind(dat_region_annual |> filter(year == 2020),
      dat_ssa_annual |> filter(year == 2020)) |>
  write.csv("./results/prevalence_2020_naat_unadjusted.csv", row.names = FALSE)

# Ratios
ratio_btwn_naat_unadj <- rbind(
  prediction(modct12, "CT")$ssa_ratio,
  prediction(modng12, "NG")$ssa_ratio,
  prediction(modtv12, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Between study",
         diag = "NAAT, unadjusted",
         n = c(modct12$modelInfo$nobs, modng12$modelInfo$nobs, modtv12$modelInfo$nobs))

# < 2010 onwards----
df_ct_2010 <- df |> filter(sti == "CT") |>
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Youth")) |>
  filter(year_mid >= 2010)

df_ng_2010 <- df |> filter(sti == "NG") |>
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Youth")) |>
  filter(year_mid >= 2010)

df_tv_2010 <- df |> filter(sti == "TV") |>
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Youth")) |>
  filter(year_mid >= 2010)

# Formula
form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  sex:region + population + age_group + hiv_status + (1 | study_id)

# Logit models
modct1_2010 <- glmmTMB(form, data = df_ct_2010, family = binomial(link="logit"))
modng1_2010 <- glmmTMB(form, data = df_ng_2010, family = binomial(link="logit"))
modtv1_2010 <- glmmTMB(form, data = df_tv_2010, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct1_2010)
theta_ct <- as.vector(modct1_2010$fit$par[names(modct1_2010$fit$par) == "theta"])
fixed_ct <- as.vector(modct1_2010$fit$par[names(modct1_2010$fit$par) == "beta"])
#fixed_ct <- append(fixed_ct, 0, length(fixed_ct))

summary(modng1_2010)
theta_ng <- as.vector(modng1_2010$fit$par[names(modng1_2010$fit$par) == "theta"])
fixed_ng <- as.vector(modng1_2010$fit$par[names(modng1_2010$fit$par) == "beta"])
#fixed_ng <- append(fixed_ng, 0, length(fixed_ng))

summary(modtv1_2010)
theta_tv <- as.vector(modtv1_2010$fit$par[names(modtv1_2010$fit$par) == "theta"])
fixed_tv <- as.vector(modtv1_2010$fit$par[names(modtv1_2010$fit$par) == "beta"])
#fixed_tv <- append(fixed_tv, 0, length(fixed_tv))

# Log models
modct2_2010 <- glmmTMB(form, data = df_ct_2010, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng2_2010 <- glmmTMB(form, data = df_ng_2010, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv2_2010 <- glmmTMB(form, data = df_tv_2010, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

# Predict by region
dat_region_annual <- rbind(
  prediction(modct2_2010, "CT")$region_prev,
  prediction(modng2_2010, "NG")$region_prev,
  prediction(modtv2_2010, "TV")$region_prev)

dat_region_annual |> write.csv("./results/prevalence_alldata_2010onwards.csv", row.names = FALSE)

# Predict weighted mean for SSA 
dat_ssa_annual <- rbind(
  prediction(modct2_2010, "CT")$ssa_prev,
  prediction(modng2_2010, "NG")$ssa_prev,
  prediction(modtv2_2010, "TV")$ssa_prev) |>
  mutate(region = "SSA") |>
  relocate(sti, sex, region)

# Save 2020 predictions
rbind(dat_region_annual |> filter(year == 2020),
      dat_ssa_annual |> filter(year == 2020)) |>
  write.csv("./results/prevalence_2020_2010onwards.csv", row.names = FALSE)


# WITHIN STUDY ANALYSIS ----

sex_within <- df |> 
  # Filter out Mullick 2023 adult data point (focus on same pop for ratios)
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
  # Consider whether to filter out Paz-Soldan 2012 as different age groups
  # filter(!cov_id == "#23347") |>
  group_by(study_name, study_id, country, location, sti) |> 
  summarise(n=n_distinct(sex)) |> 
  filter(n>1) |>
  mutate(keep = "y")

df_sex_within <- df |>
  # Filter out Mullick 2023 adult data point (focus on same pop for ratios)
  filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
  # Consider whether to filter out Paz-Soldan 2012 as different age groups
  # filter(!cov_id == "#23347") |>
  left_join(sex_within) |>
  filter(!is.na(keep))

df_ct_within <- df_sex_within |> filter(sti=="CT")
df_ng_within <- df_sex_within |> filter(sti=="NG")
df_tv_within <- df_sex_within |> filter(sti=="TV")

# save
df_sex_within |>
  select(study_name, study_id, year_mid, region, country, population, age_group, sex, sti, adj_prev, adj_se) |>
  write.csv("./data/data_withinstudy_adjusted.csv", row.names = FALSE)

# < Main model, with test adjustment ----

# original formula
# form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
#  sex:region + population + age_group + hiv_status + (1 | study_id)

form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region +
  sex + age_group + hiv_status + test + (1 | study_id)
# exclude population as studies relatively distributed across other groups
# no ANC pop as reference group (obviously)
# including population makes time trend non-estimable for TV in EA  

# exclude region:sex interaction as too few studies in WCA. Results in very wide CIs for overall estimate. 

# Previous version used a simpler model
# form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
#   sex:region + (1 | study_id)
# Removed HIV-status 
  # No HIV positive pops, and predominantly mixed HIV status
  # TV dataset only mixed HIV status, so model does not fit with only one level of the variable
# Removed population - insufficient data
# Removed age_group - unable to calculate SE for TV with this included

# Logit models
modct1_w <- glmmTMB(form, data = df_ct_within, family = binomial(link="logit"))
modng1_w <- glmmTMB(form, data = df_ng_within, family = binomial(link="logit"))
modtv1_w <- glmmTMB(form, data = df_tv_within, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct1_w)
theta_ct_w <- as.vector(modct1_w$fit$par[names(modct1_w$fit$par) == "theta"])
fixed_ct_w <- as.vector(modct1_w$fit$par[names(modct1_w$fit$par) == "beta"])

summary(modng1_w)
theta_ng_w <- as.vector(modng1_w$fit$par[names(modng1_w$fit$par) == "theta"])
fixed_ng_w <- as.vector(modng1_w$fit$par[names(modng1_w$fit$par) == "beta"])

summary(modtv1_w)
theta_tv_w <- as.vector(modtv1_w$fit$par[names(modtv1_w$fit$par) == "theta"])
fixed_tv_w <- as.vector(modtv1_w$fit$par[names(modtv1_w$fit$par) == "beta"])

# Log models
modct2_w <- glmmTMB(form, data = df_ct_within, family = binomial(link="log"),
                    start = list(theta = theta_ct_w, beta = fixed_ct_w))

modng2_w <- glmmTMB(form, data = df_ng_within, family = binomial(link="log"),
                    start = list(theta = theta_ng_w, beta = fixed_ng_w))

modtv2_w <- glmmTMB(form, data = df_tv_within, family = binomial(link="log"),
                    start = list(theta = theta_tv_w, beta = fixed_tv_w))

sjPlot::tab_model(modct2_w, modng2_w, modtv2_w, p.style = "stars", wrap.labels = 100)

# Ratio within study
ratio_within_adj <- rbind(
  prediction(modct2_w, "CT")$ssa_ratio,
  prediction(modng2_w, "NG")$ssa_ratio,
  prediction(modtv2_w, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Within study",
         diag = "All tests, adjusted",
         n = c(modct2_w$modelInfo$nobs, modng2_w$modelInfo$nobs, modtv2_w$modelInfo$nobs))

# < Main model, no test adjustment ----

form <- cbind(num,(denom-num)) ~ region + year_regression:region +
  sex + age_group + hiv_status + test + (1 | study_id)

# Logit models
modct3_w <- glmmTMB(form, data = df_ct_within, family = binomial(link="logit"))
modng3_w <- glmmTMB(form, data = df_ng_within, family = binomial(link="logit"))
modtv3_w <- glmmTMB(form, data = df_tv_within, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct3_w)
theta_ct_w <- as.vector(modct3_w$fit$par[names(modct3_w$fit$par) == "theta"])
fixed_ct_w <- as.vector(modct3_w$fit$par[names(modct3_w$fit$par) == "beta"])

summary(modng3_w)
theta_ng_w <- as.vector(modng3_w$fit$par[names(modng3_w$fit$par) == "theta"])
fixed_ng_w <- as.vector(modng3_w$fit$par[names(modng3_w$fit$par) == "beta"])

summary(modtv3_w)
theta_tv_w <- as.vector(modtv3_w$fit$par[names(modtv3_w$fit$par) == "theta"])
fixed_tv_w <- as.vector(modtv3_w$fit$par[names(modtv3_w$fit$par) == "beta"])

# Log models
modct4_w <- glmmTMB(form, data = df_ct_within, family = binomial(link="log"),
                    start = list(theta = theta_ct_w, beta = fixed_ct_w))

modng4_w <- glmmTMB(form, data = df_ng_within, family = binomial(link="log"),
                    start = list(theta = theta_ng_w, beta = fixed_ng_w))

modtv4_w <- glmmTMB(form, data = df_tv_within, family = binomial(link="log"),
                    start = list(theta = theta_tv_w, beta = fixed_tv_w))

sjPlot::tab_model(modct4_w, modng4_w, modtv4_w, p.style = "stars", wrap.labels = 100)

# Ratio within study
ratio_within_unadj <- rbind(
  prediction(modct4_w, "CT")$ssa_ratio,
  prediction(modng4_w, "NG")$ssa_ratio,
  prediction(modtv4_w, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Within study",
         diag = "All tests, unadjusted",
         n = c(modct4_w$modelInfo$nobs, modng4_w$modelInfo$nobs, modtv4_w$modelInfo$nobs))

# < NAAT, with test adjustment ----

df_ct_within_naat <- df_sex_within |> filter(sti=="CT", testcat == "NAAT")
df_ng_within_naat <- df_sex_within |> filter(sti=="NG", testcat == "NAAT")
df_tv_within_naat <- df_sex_within |> filter(sti=="TV", testcat == "NAAT")

form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region +
  sex + age_group + hiv_status + (1 | study_id)

# Logit models
modct1_w_naat <- glmmTMB(form, data = df_ct_within_naat, family = binomial(link="logit"))
modng1_w_naat <- glmmTMB(form, data = df_ng_within_naat, family = binomial(link="logit"))
modtv1_w_naat <- glmmTMB(form, data = df_tv_within_naat, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct1_w_naat)
theta_ct_w <- as.vector(modct1_w_naat$fit$par[names(modct1_w_naat$fit$par) == "theta"])
fixed_ct_w <- as.vector(modct1_w_naat$fit$par[names(modct1_w_naat$fit$par) == "beta"])

summary(modng1_w_naat)
theta_ng_w <- as.vector(modng1_w_naat$fit$par[names(modng1_w_naat$fit$par) == "theta"])
fixed_ng_w <- as.vector(modng1_w_naat$fit$par[names(modng1_w_naat$fit$par) == "beta"])

summary(modtv1_w_naat)
theta_tv_w <- as.vector(modtv1_w_naat$fit$par[names(modtv1_w_naat$fit$par) == "theta"])
fixed_tv_w <- as.vector(modtv1_w_naat$fit$par[names(modtv1_w_naat$fit$par) == "beta"])
fixed_tv_w <- append(fixed_tv_w, 0, length(fixed_tv_w))

# Log models
modct2_w_naat <- glmmTMB(form, data = df_ct_within_naat, family = binomial(link="log"),
                    start = list(theta = theta_ct_w, beta = fixed_ct_w))

modng2_w_naat <- glmmTMB(form, data = df_ng_within_naat, family = binomial(link="log"),
                    start = list(theta = theta_ng_w, beta = fixed_ng_w))

modtv2_w_naat <- glmmTMB(form, data = df_tv_within_naat, family = binomial(link="log"),
                    start = list(theta = theta_tv_w, beta = fixed_tv_w))

sjPlot::tab_model(modct1_w_naat, modng1_w_naat, modtv1_w_naat, p.style = "stars", wrap.labels = 100)
sjPlot::tab_model(modct2_w_naat, modng2_w_naat, modtv2_w_naat, p.style = "stars", wrap.labels = 100)

# Ratio within study
ratio_within_naat_adj <- rbind(
  prediction(modct2_w_naat, "CT")$ssa_ratio,
  prediction(modng2_w_naat, "NG")$ssa_ratio,
  prediction(modtv2_w_naat, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Within study",
         diag = "NAAT, adjusted",
         n = c(modct2_w_naat$modelInfo$nobs, modng2_w_naat$modelInfo$nobs, modtv2_w_naat$modelInfo$nobs))

# < NAAT, no test adjustment ----
form <- cbind(num,(denom-num)) ~ region + year_regression:region +
  sex + age_group + hiv_status + (1 | study_id)

# Logit models
modct3_w_naat <- glmmTMB(form, data = df_ct_within_naat, family = binomial(link="logit"))
modng3_w_naat <- glmmTMB(form, data = df_ng_within_naat, family = binomial(link="logit"))
modtv3_w_naat <- glmmTMB(form, data = df_tv_within_naat, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct3_w_naat)
theta_ct_w <- as.vector(modct3_w_naat$fit$par[names(modct3_w_naat$fit$par) == "theta"])
fixed_ct_w <- as.vector(modct3_w_naat$fit$par[names(modct3_w_naat$fit$par) == "beta"])

summary(modng3_w_naat)
theta_ng_w <- as.vector(modng3_w_naat$fit$par[names(modng3_w_naat$fit$par) == "theta"])
fixed_ng_w <- as.vector(modng3_w_naat$fit$par[names(modng3_w_naat$fit$par) == "beta"])

summary(modtv3_w_naat)
theta_tv_w <- as.vector(modtv3_w_naat$fit$par[names(modtv3_w_naat$fit$par) == "theta"])
fixed_tv_w <- as.vector(modtv3_w_naat$fit$par[names(modtv3_w_naat$fit$par) == "beta"])
fixed_tv_w <- append(fixed_tv_w, 0, length(fixed_tv_w))

# Log models
modct4_w_naat <- glmmTMB(form, data = df_ct_within_naat, family = binomial(link="log"),
                         start = list(theta = theta_ct_w, beta = fixed_ct_w))

modng4_w_naat <- glmmTMB(form, data = df_ng_within_naat, family = binomial(link="log"),
                         start = list(theta = theta_ng_w, beta = fixed_ng_w))

modtv4_w_naat <- glmmTMB(form, data = df_tv_within_naat, family = binomial(link="log"),
                         start = list(theta = theta_tv_w, beta = fixed_tv_w))

sjPlot::tab_model(modct4_w_naat, modng4_w_naat, modtv4_w_naat, p.style = "stars", wrap.labels = 100)

# Ratio within study
ratio_within_naat_unadj <- rbind(
  prediction(modct4_w_naat, "CT")$ssa_ratio,
  prediction(modng4_w_naat, "NG")$ssa_ratio,
  prediction(modtv4_w_naat, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Within study",
         diag = "NAAT, unadjusted",
         n = c(modct4_w_naat$modelInfo$nobs, modng4_w_naat$modelInfo$nobs, modtv4_w_naat$modelInfo$nobs))

# < Pooled model (no covariates), with test adjustment ----

form <- cbind(adj_num,(adj_denom-adj_num)) ~ sex + (1 | study_id)

# Logit models

modct1_wn <- glmmTMB(form, data = df_ct_within, family = binomial(link="logit"))
modng1_wn <- glmmTMB(form, data = df_ng_within, family = binomial(link="logit"))
modtv1_wn <- glmmTMB(form, data = df_tv_within, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
theta_ct_wn <- as.vector(modct1_wn$fit$par[names(modct1_wn$fit$par) == "theta"])
fixed_ct_wn <- as.vector(modct1_wn$fit$par[names(modct1_wn$fit$par) == "beta"])

theta_ng_wn <- as.vector(modng1_wn$fit$par[names(modng1_wn$fit$par) == "theta"])
fixed_ng_wn <- as.vector(modng1_wn$fit$par[names(modng1_wn$fit$par) == "beta"])

theta_tv_wn <- as.vector(modtv1_wn$fit$par[names(modtv1_wn$fit$par) == "theta"])
fixed_tv_wn <- as.vector(modtv1_wn$fit$par[names(modtv1_wn$fit$par) == "beta"])

# Log models

modct2_wn <- glmmTMB(form, data = df_ct_within, family = binomial(link="log"),
                     start = list(theta = theta_ct_wn, beta = fixed_ct_wn))

modng2_wn <- glmmTMB(form, data = df_ng_within, family = binomial(link="log"),
                     start = list(theta = theta_ng_wn, beta = fixed_ng_wn))

modtv2_wn <- glmmTMB(form, data = df_tv_within, family = binomial(link="log"),
                     start = list(theta = theta_tv_wn, beta = fixed_tv_wn))

tab_model(modct2_wn, modng2_wn, modtv2_wn, p.style = "stars", wrap.labels = 100)

# Ratio within study unadjusted (no covariates) 
extract_ratio <- function(mod) {
  coef(summary(mod))$cond |>
    as.data.frame() |>
    rownames_to_column("variable") |>
    filter(variable == "sexMale") |>
    mutate(year = 2012,
           ratio = exp(Estimate),
           lwr = exp(Estimate - 1.96*`Std. Error`),
           upr = exp(Estimate + 1.96*`Std. Error`),
           tau_study = insight::get_variance(mod)$var.random) |>
    select(year, ratio, lwr, upr, tau_study)
}

ratio_pool_adj <- rbind(
  extract_ratio(modct2_wn) |> mutate(sti = "CT"),
  extract_ratio(modng2_wn) |> mutate(sti = "NG"),
  extract_ratio(modtv2_wn) |> mutate(sti = "TV")) |> 
  relocate(sti) |>
  mutate(type = "Pooled",
         diag = "All tests, adjusted",
         n = c(modct2_wn$modelInfo$nobs, modng2_wn$modelInfo$nobs, modtv2_wn$modelInfo$nobs))

# < Pooled model (no covariates), no test adjustment ----

form <- cbind(num,(denom-num)) ~ sex + (1 | study_id)

# Logit models

modct3_wn <- glmmTMB(form, data = df_ct_within, family = binomial(link="logit"))
modng3_wn <- glmmTMB(form, data = df_ng_within, family = binomial(link="logit"))
modtv3_wn <- glmmTMB(form, data = df_tv_within, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
theta_ct_wn <- as.vector(modct3_wn$fit$par[names(modct3_wn$fit$par) == "theta"])
fixed_ct_wn <- as.vector(modct3_wn$fit$par[names(modct3_wn$fit$par) == "beta"])

theta_ng_wn <- as.vector(modng3_wn$fit$par[names(modng3_wn$fit$par) == "theta"])
fixed_ng_wn <- as.vector(modng3_wn$fit$par[names(modng3_wn$fit$par) == "beta"])

theta_tv_wn <- as.vector(modtv3_wn$fit$par[names(modtv3_wn$fit$par) == "theta"])
fixed_tv_wn <- as.vector(modtv3_wn$fit$par[names(modtv3_wn$fit$par) == "beta"])

# Log models

modct4_wn <- glmmTMB(form, data = df_ct_within, family = binomial(link="log"),
                     start = list(theta = theta_ct_wn, beta = fixed_ct_wn))

modng4_wn <- glmmTMB(form, data = df_ng_within, family = binomial(link="log"),
                     start = list(theta = theta_ng_wn, beta = fixed_ng_wn))

modtv4_wn <- glmmTMB(form, data = df_tv_within, family = binomial(link="log"),
                     start = list(theta = theta_tv_wn, beta = fixed_tv_wn))

tab_model(modct4_wn, modng4_wn, modtv4_wn, p.style = "stars", wrap.labels = 100)

ratio_pool_unadj <- rbind(
  extract_ratio(modct4_wn) |> mutate(sti = "CT"),
  extract_ratio(modng4_wn) |> mutate(sti = "NG"),
  extract_ratio(modtv4_wn) |> mutate(sti = "TV")) |> 
  relocate(sti) |>
  mutate(type = "Pooled",
         diag = "All tests, unadjusted",
         n = c(modct4_wn$modelInfo$nobs, modng4_wn$modelInfo$nobs, modtv4_wn$modelInfo$nobs))

# ROB SENSITIVITY ANALYSES ----

form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  sex:region + population + age_group + hiv_status + test + 
  rob_study_sample + rob_participant + rob_measurement + rob_precision + (1 | study_id)

# Logit models
modct1_rob <- glmmTMB(form, data = df_ct, family = binomial(link="logit"))
modng1_rob <- glmmTMB(form, data = df_ng, family = binomial(link="logit"))
modtv1_rob <- glmmTMB(form, data = df_tv, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct1_rob)
theta_ct <- as.vector(modct1_rob$fit$par[names(modct1_rob$fit$par) == "theta"])
fixed_ct <- as.vector(modct1_rob$fit$par[names(modct1_rob$fit$par) == "beta"])

summary(modng1_rob)
theta_ng <- as.vector(modng1_rob$fit$par[names(modng1_rob$fit$par) == "theta"])
fixed_ng <- as.vector(modng1_rob$fit$par[names(modng1_rob$fit$par) == "beta"])

summary(modtv1_rob)
theta_tv <- as.vector(modtv1_rob$fit$par[names(modtv1_rob$fit$par) == "theta"])
fixed_tv <- as.vector(modtv1_rob$fit$par[names(modtv1_rob$fit$par) == "beta"])

# Log models
modct2_rob <- glmmTMB(form, data = df_ct, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng2_rob <- glmmTMB(form, data = df_ng, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv2_rob <- glmmTMB(form, data = df_tv, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

sjPlot::tab_model(modct2_rob, modng2_rob, modtv2_rob, p.style = "stars", wrap.labels = 100)

# ALL RATIOS ---- 

df_ratio <- rbind(ratio_btwn_adj,
      ratio_btwn_unadj,
      ratio_btwn_naat_adj,
      ratio_btwn_naat_unadj,
      ratio_within_adj,
      ratio_within_unadj,
      ratio_within_naat_adj,
      ratio_within_naat_unadj,
      ratio_pool_adj,
      ratio_pool_unadj) 

df_ratio |>
  arrange(sti) |>
  print(n=50)

write.csv(df_ratio, "./results/ratios.csv", row.names = FALSE)

# SAVE REGRESSION TABLES ----

# variable order
desired_order <- c("(Intercept)","WCA","EA", "SA",
                   "WCA:Male", "EA:Male","SA:Male", 
                   "WCA:Year", "EA:Year","SA:Year", 
                   "ANC attendees", "FP attendees", "GYN attendees", 
                   "PHC/OPD attendees", "Students", "Community members", 
                   "HIV/STI prevention trial participants", "Population-representative survey participants", 
                   "Female","Male", "Adult", "Youth",
                   "Mixed", "HIV negative",
                   "NAAT","Culture","DFA","ELISA","Rapid antigen test","Wet mount",
                   "Participant sampling", "Participant characterisation",
                   "Diagnostic method consistency", "Sample size adequacy")

reg_table <- function(model){
  
  coef(summary(model))$cond %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    mutate(est = sprintf(exp(Estimate),fmt = '%#.2f'),
           lwr = sprintf(exp(Estimate - 1.96*`Std. Error`),fmt = '%#.2f'),
           upr = sprintf(exp(Estimate + 1.96*`Std. Error`),fmt = '%#.2f'),           
           estimate = case_when(is.na(Estimate) ~ "-",
                                TRUE ~ paste0(est," (",lwr,"-",upr,")"))) %>%
    select(variable,estimate) %>%
    mutate(variable = gsub("region","", variable),
           variable = gsub("year_regression","Year",variable),
           variable = gsub("population","", variable),
           variable = gsub("sexMale","Male",variable),
           variable = gsub("age_group","",variable),
           variable = gsub("test","",variable),
           variable = gsub("Rapid antigen ","Rapid antigen test",variable),
           variable = gsub("hiv_status","",variable),
           variable = gsub("rob_","",variable),
           variable = gsub("participantHigher","Participant characterisation",variable),
           variable = gsub("Higher","",variable),
           variable = gsub("study_sample","Participant sampling",variable),
           variable = gsub("measurement","Diagnostic method consistency",variable),
           variable = gsub("precision","Sample size adequacy",variable)) |>
    rbind(data.frame(variable = c("ANC attendees","Adult", "Mixed", "SA", "NAAT"),
                     estimate = rep("Ref",5)))
}

# < Between study, adjusted ----
full_join(
  reg_table(modct2) |> rename(CT = estimate),
  reg_table(modng2) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv2) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_bwtn_all_adjusted.csv", row.names = FALSE)

insight::get_variance(modct2)
insight::get_variance(modng2)
insight::get_variance(modtv2)

sjPlot::tab_model(modct2, modng2, modtv2, p.style = "stars", wrap.labels = 100)

# < Between study, unadjusted ----
full_join(
  reg_table(modct4) |> rename(CT = estimate),
  reg_table(modng4) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv4) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_btwn_all_unadjusted.csv", row.names = FALSE)

insight::get_variance(modct4)
insight::get_variance(modng4)
insight::get_variance(modtv4)

sjPlot::tab_model(modct4, modng4, modtv4, p.style = "stars", wrap.labels = 100)

# < Between study, NAAT adjusted ----
full_join(
  reg_table(modct10) |> rename(CT = estimate),
  reg_table(modng10) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv10) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_btwn_NAAT_adjusted.csv", row.names = FALSE)

insight::get_variance(modct10)
insight::get_variance(modng10)
insight::get_variance(modtv10)

sjPlot::tab_model(modct10, modng10, modtv10, p.style = "stars", wrap.labels = 100)

# < Between study, NAAT unadjusted ----
full_join(
  reg_table(modct12) |> rename(CT = estimate),
  reg_table(modng12) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv12) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_bwtn_NAAT_unadjusted.csv", row.names = FALSE)

insight::get_variance(modct12)
insight::get_variance(modng12)
insight::get_variance(modtv12)

sjPlot::tab_model(modct12, modng12, modtv12, p.style = "stars", wrap.labels = 100)

# < Within study, adjusted ----
full_join(
  reg_table(modct2_w) |> rename(CT = estimate),
  reg_table(modng2_w) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv2_w) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |>
  write.csv("./tables/regression_within_all_adjusted.csv", row.names = FALSE)

insight::get_variance(modct2_w)
insight::get_variance(modng2_w)
insight::get_variance(modtv2_w)

sjPlot::tab_model(modct2_w, modng2_w, modtv2_w, p.style = "stars", wrap.labels = 100)

# < Within study, unadjusted ----
full_join(
  reg_table(modct4_w) |> rename(CT = estimate),
  reg_table(modng4_w) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv4_w) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |>
  write.csv("./tables/regression_within_all_unadjusted.csv", row.names = FALSE)

insight::get_variance(modct4_w)
insight::get_variance(modng4_w)
insight::get_variance(modtv4_w)

sjPlot::tab_model(modct4_w, modng4_w, modtv4_w, p.style = "stars", wrap.labels = 100)

# < Within study, NAAT adjusted ----
full_join(
  reg_table(modct2_w_naat) |> rename(CT = estimate),
  reg_table(modng2_w_naat) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv2_w_naat) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |>
  write.csv("./tables/regression_within_naat_adjusted.csv", row.names = FALSE)

insight::get_variance(modct2_w_naat)
insight::get_variance(modng2_w_naat)
insight::get_variance(modtv2_w_naat)

sjPlot::tab_model(modct2_w_naat, modng2_w_naat, modtv2_w_naat, p.style = "stars", wrap.labels = 100)

# < Within study, NAAT unadjusted ----
full_join(
  reg_table(modct4_w_naat) |> rename(CT = estimate),
  reg_table(modng4_w_naat) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv4_w_naat) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |>
  write.csv("./tables/regression_within_naat_unadjusted.csv", row.names = FALSE)

insight::get_variance(modct4_w_naat)
insight::get_variance(modng4_w_naat)
insight::get_variance(modtv4_w_naat)

sjPlot::tab_model(modct4_w_naat, modng4_w_naat, modtv4_w_naat, p.style = "stars", wrap.labels = 100)

# < ROB ----

full_join(
  reg_table(modct2_rob) |> rename(CT = estimate),
  reg_table(modng2_rob) |> rename(NG = estimate)) |>
  full_join(reg_table(modtv2_rob) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |>
  write.csv("./tables/regression_rob.csv", row.names = FALSE)

insight::get_variance(modct2_rob)
insight::get_variance(modng2_rob)
insight::get_variance(modtv2_rob)

sjPlot::tab_model(modct2_rob, modng2_rob, modtv2_rob, p.style = "stars", wrap.labels = 100)

# STUDY CHARACTERISTICS ----

# Number articles
df |> 
  group_by() |>
  summarise(n = n_distinct(study_id))

df |> 
  group_by(sti) |>
  summarise(n = n_distinct(study_id))

# Number studies
cross_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by() |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by() |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article)


tab_characteristics <- function(df){
  
  calc <- function(df, var){
    
    full_join(
      df |> 
        filter(!is.na(study_name)) |>
        group_by(!!sym(var)) |>
        summarise(n_study = n_distinct(study_name)),
      df |> 
        filter(is.na(study_name)) |>
        group_by(!!sym(var)) |>
        summarise(n_article = n_distinct(study_id))) |>
      mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
             n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
             n_tot = n_study + n_article)
  }
  
  rbind(calc(df, "region_actual") |> mutate(var = "Region") |> rename(group = region_actual),
        calc(df, "year_grp")|> mutate(var = "Study midpoint year")|> rename(group = year_grp),
        calc(df, "population") |> mutate(var = "Population group")|> rename(group = population),
        calc(df, "sex") |> mutate(var = "Sex")|> rename(group = sex),
        calc(df, "age_group") |> mutate(var = "Age group")|> rename(group = age_group),
        calc(df, "hiv_status")|> mutate(var = "HIV status")|> rename(group = hiv_status)) |>
    mutate(tot = cbind(
      df |> 
        filter(!is.na(study_name)) |>
        group_by() |>
        summarise(n_study = n_distinct(study_name)),
      df |> 
        filter(is.na(study_name)) |>
        group_by() |>
        summarise(n_article = n_distinct(study_id))) |>
        mutate(n_tot = n_study + n_article) |> pull(n_tot)) |>
    mutate(prop = n_tot/tot,
           col = paste0(n_tot, " (", sprintf("%.1f", prop*100), "%)")) |>
    select(var, group, col, tot)
  
}

# Between study 
t_characteristics_all <- full_join(
  tab_characteristics(df |> filter(sti == "CT")) |> rename(CT = col, CT_tot = tot),
  tab_characteristics(df |> filter(sti == "NG")) |> rename(NG = col, NG_tot = tot)) |>
  full_join(tab_characteristics(df |> filter(sti == "TV")) |> rename(TV = col, TV_tot = tot)) |>
  full_join(tab_characteristics(df) |> rename(ALL = col, ALL_tot = tot)) |>
  relocate(var, group, CT, NG, TV)

t_characteristics_all |> write.csv("./tables/study_characteristics_alldata.csv", row.names = FALSE)

# Within study
t_characteristics_within <- full_join(
  tab_characteristics(df_ct_within) |> rename(CT = col, CT_tot = tot),
  tab_characteristics(df_ng_within) |> rename(NG = col, NG_tot = tot)) |>
  full_join(tab_characteristics(df_tv_within) |> rename(TV = col, TV_tot = tot)) |>
  full_join(tab_characteristics(df_sex_within) |> rename(ALL = col, ALL_tot = tot)) |>
  relocate(var, group, CT, NG, TV)

t_characteristics_within |> write.csv("./tables/study_characteristics_within.csv", row.names = FALSE)

# < Countries ----
# Number countries

df |> 
  group_by(country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by() |>
  summarise(n=n_distinct(country_list))

df |> 
  group_by(country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  distinct(country_list) |> print(n=40)

df |> 
  group_by(sti, country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by(sti) |>
  summarise(n=n_distinct(country_list))

# Number observations per country
df |> 
  filter(!region_actual == "Multiple") |>
  group_by(region,country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by(region,country_list) |>
  summarise(n=sum(n)) |>
  arrange(region,-n)
  
df |> 
  filter(!region_actual == "Multiple") |>
  group_by(region,country) |>
  summarise(n = n()) |> arrange(region, -n) |> write.csv("temp2.csv")

df |> 
  filter(!region_actual == "Multiple") |>
  group_by(region) |>
  summarise(n = n()) 

# Number articles per country
df |> 
  group_by(country) |>
  summarise(n = n_distinct(study_id)) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by(country_list) |>
  summarise(n=sum(n)) |>
  arrange(-n)

# Number countries per region

df |> 
  group_by(region_actual, country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by(region_actual) |>
  summarise(n=n_distinct(country_list))

df |> 
  group_by(region_actual, country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by(region_actual) |>
  distinct(country_list) |>
  print(n=100)

# Number countries for men

df |> 
  group_by(sex, country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by(sex) |>
  distinct(country_list) |>
  print(n=100)

# Number studies per country - overall

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_id)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_art=sum(n))
) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_art = case_when(is.na(n_art) ~ 0, TRUE ~ n_art),
         n_tot = n_study + n_art) |>
  arrange(-n_tot) |>
  mutate(prop = n_tot/212)

# Number studies per region and country - overall

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(region, country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(region, country_list) |>
    summarise(n_study=sum(n)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(region, country) |>
    summarise(n = n_distinct(study_id)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(region, country_list) |>
    summarise(n_art=sum(n))
) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_art = case_when(is.na(n_art) ~ 0, TRUE ~ n_art),
         n_tot = n_study + n_art) |>
  arrange(region, -n_tot) |>
  mutate(prop = n_tot/211*100) |>
  print(n=100)

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(region_actual, country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(region_actual, country_list) |>
    summarise(n_study=sum(n)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(region_actual, country) |>
    summarise(n = n_distinct(study_id)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(region_actual, country_list) |>
    summarise(n_art=sum(n))
) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_art = case_when(is.na(n_art) ~ 0, TRUE ~ n_art),
         n_tot = n_study + n_art) |>
  arrange(region_actual, -n_tot) |>
  mutate(prop = n_tot/211*100) |>
  print(n=100)

# Number studies per country - CT
full_join(
  df |>
    filter(sti == "CT") |>
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df |>
    filter(sti == "CT") |>
    filter(is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_id)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_art=sum(n))
) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_art = case_when(is.na(n_art) ~ 0, TRUE ~ n_art),
         n_tot = n_study + n_art) |>
  arrange(-n_tot)|> print(n=100)

# Number studies per country - NG
full_join(
  df |>
    filter(sti == "NG") |>
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df |>
    filter(sti == "NG") |>
    filter(is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_id)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_art=sum(n))
) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_art = case_when(is.na(n_art) ~ 0, TRUE ~ n_art),
         n_tot = n_study + n_art) |>
  arrange(-n_tot) |> print(n=100)

# Number studies per country - TV
full_join(
  df |>
    filter(sti == "TV") |> 
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df |>
    filter(sti == "TV") |>
    filter(is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_id)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_art=sum(n))
) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_art = case_when(is.na(n_art) ~ 0, TRUE ~ n_art),
         n_tot = n_study + n_art) |>
  arrange(-n_tot) |> print(n=100)

# < Sex ----

# Number studies among men only
full_join(
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(!is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & is.na(Female)) |>
    group_by(sti) |>
    summarise(n_study=n_distinct(study_name)),
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & is.na(Female)) |>
    group_by(sti) |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# Number studies among women only
full_join(
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(!is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_study=n_distinct(study_name)),
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# Number studies among both sexes
full_join(
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(!is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_study=n_distinct(study_name)),
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# < Sex total ----

# Number studies among men only
cross_join(
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(!is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & is.na(Female)) |>
    group_by() |>
    summarise(n_study=n_distinct(study_name)),
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & is.na(Female)) |>
    group_by() |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# Number studies among women only
cross_join(
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(!is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_study=n_distinct(study_name)),
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# Number studies among both sexes
cross_join(
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(!is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_study=n_distinct(study_name)),
  df |> 
    filter(!(cov_id == "#23088" & sex == "Male" & age_group == "Adult")) |>
    filter(is.na(study_name)) |>
    select(study_name, study_id, year_mid, country, location, population, age_group, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# < Age range ----
# Classify age range groups
df_age_count <- df |>
  mutate(age_range_class = case_when(age_range == "NR" ~ "NR",
                                     age_min >= 15 & age_max <= 49 & !age_max=="NR" & !age_min=="NR" ~ "15-49",
                                     age_min < 15 & age_max > 49 & !age_max=="NR" & !age_min=="NR" ~ "<15 & >49",
                                     age_min < 15 & !age_min=="NR"~ "<15",
                                     age_max > 49 & !age_max=="NR" ~ ">49",
                                     age_min != "NR" & age_max == "NR" ~ "no upper"))

df_age_count |> count(age_range_class, age_range)

df_age_count |> count(cov_id,age_range_class) |>
  group_by(cov_id) |>
  summarise(n = n_distinct(age_range_class)) |>
  filter(n>1) |>
  pull(cov_id)

df_age_count |> filter(cov_id %in% c("#21877", "#22843", "#23088", "#31017")) |>
  select(cov_id, age_range, age_range_class, population, country, sex, sti)

df_age_count <- df_age_count |>
  filter(age_group == "Adult") |>
  filter(!(cov_id == "#23088" & age_range_class == "no upper"),
         !(cov_id == "#31017" & age_range_class == ">49"),
         !(cov_id == "#22843" & age_range_class == "NR"),
         !(cov_id == "#21877" & age_range_class == "15-49"))

full_join(
  df_age_count |> 
    filter(!is.na(study_name)) |>
    group_by(age_range_class) |>
    summarise(n_study=n_distinct(study_name)),
  df_age_count |> 
    filter(is.na(study_name)) |>
    group_by(age_range_class) |>
    summarise(n_art=n_distinct(study_id))) |>
  rowwise() |>
  mutate(n_tot = sum(n_study,n_art, na.rm=TRUE))

# < ROB ----

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(sti, rob_study_sample) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(sti, rob_study_sample) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_study_sample, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(sti, rob_participant) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(sti, rob_participant) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_participant, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)


full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(sti, rob_measurement) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(sti, rob_measurement) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_measurement, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)


full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(sti, rob_precision) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(sti, rob_precision) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_precision, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)

# Overall
full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(rob_study_sample) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(rob_study_sample) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_study_sample, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(rob_participant) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(rob_participant) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_participant, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(rob_measurement) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(rob_measurement) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_measurement, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)

full_join(
  df |> 
    filter(!is.na(study_name)) |>
    group_by(rob_precision) |>
    summarise(n_study = n_distinct(study_name)),
  df |> 
    filter(is.na(study_name)) |>
    group_by(rob_precision) |>
    summarise(n_article = n_distinct(study_id))) |>
  mutate(n_study = case_when(is.na(n_study) ~ 0, TRUE ~ n_study),
         n_article = case_when(is.na(n_article) ~ 0, TRUE ~ n_article),
         n_tot = n_study + n_article) |>
  select(!c(n_study, n_article)) |>
  pivot_wider(names_from = rob_precision, values_from = n_tot) |>
  mutate(n_tot = Lower + Higher,
         p_low = Lower/n_tot*100,
         p_high = Higher/n_tot*100)

