## STI PREVALENCE RATIOS

library(tidyverse)
library(glmmTMB)
library(sjPlot)

# PREP DATA ----

# Population by region, 15-49 years in 2020
df_pop <- read.csv("./data/unpopulation_dataportal_20240426093533.csv") |>
  rename_all(tolower) |>
  filter(time == 2020,
         sex %in% c("Female", "Male")) |>
  select(time, location, iso3, sex, age, value) |>
  mutate(region = case_when(
    iso3 %in% c("BWA","SWZ","LSO","NAM","ZAF") ~ "Southern",
    iso3 %in% c("AGO","CMR","CAF","TCD","COG","COD","GNQ","GAB","STP") ~ "Central",
    iso3 %in% c("BEN","BFA","CPV","CIV","GMB","GHA","GIN","GNB",
                "LBR","MLI","MRT","NER","NGA","SEN","SLE","TGO") ~ "Western",
    iso3 %in% c("BDI","COM","DJI","ERI","ETH","KEN","MDG","MWI",
                "MUS","MOZ","RWA","SYC","SOM","SSD",
                "UGA","TZA","TZA","ZMB","ZWE") ~ "Eastern")) |>
  mutate(region = fct_collapse(region, 
                               "Western and Central" = c("Central", "Western")),
         region = factor(region,
                         levels = c("Western and Central", "Eastern", "Southern"))) |>
  filter(!is.na(region)) |>
  group_by(region, sex) |>
  summarise(population = sum(value)) |>
  ungroup()

# Study data
df <- readxl::read_xlsx("./data/study_data.xlsx") |>
  mutate(year_regression = year_mid - 2015,
         sti = factor(sti, levels = c("CT", "NG", "TV")),
         year_grp = case_when(year_mid %in% c(2005:2009) ~ "2005-2009",
                              year_mid %in% c(2010:2014) ~ "2010-2014",
                              year_mid %in% c(2015:2019) ~ "2015-2019",
                              year_mid %in% c(2020:2022) ~ "2020-2024"),
         year_grp = factor(year_grp, 
                           levels = c("2005-2009", "2010-2014", 
                                      "2015-2019", "2020-2024")),
         region = fct_collapse(region_analysis, 
                                        "Western and Central" = c("Western", "Central")),
         region = factor(region, levels = c("Southern", "Eastern", "Western and Central")),
         region_actual = factor(region_actual, 
                                levels = c("Western", "Central", "Eastern", "Southern", "Multiple")),
         sex = factor(sex, levels = c("Female", "Male", "Both sexes")),
         hiv_status = fct_relevel(hiv_status, "Mixed"),
         age_group = factor(age_group, levels = c("Adult", "Youth")),
         population = factor(population, 
                             levels = c("ANC attendees", "FP attendees", 
                                        "GYN attendees", "PHC/OPD attendees",
                                        "Students","Community members", 
                                        "HIV/STI prevention trial participants",
                                        "Population-representative survey participants")))

# FUNCTION ---- 

# Function that predicts prevalence for:
# 1. Region by year and sex
# 2. SSA mean by year and sex
# 3. SSA mean male-to-female prevalence ratio

prediction <- function(mod, sti) {
  
  df_pred <- crossing(sex = c("Female", "Male"),
                      region = c("Western and Central", "Eastern", "Southern"),
                      year = 2005:2023) %>%
    mutate(year_regression = year - 2015,
           population = "ANC attendees",
           age_group = "Adult",
           hiv_status = "Mixed",
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

prediction_2 <- function(mod, sti, y1, y2) {
  
  df_pred <- crossing(sex = c("Female", "Male"),
                      region = c("Western and Central", "Eastern", "Southern"),
                      year = 2005:2023) %>%
    mutate(year_regression = year - 2015,
           population = "ANC attendees",
           age_group = "Adult",
           hiv_status = "Mixed",
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
    mutate(est = paste0(sprintf(ratio,fmt = '%#.1f')," (", sprintf(lwr_ratio,fmt = '%#.1f'), 
                        "-", sprintf(upr_ratio,fmt = '%#.1f'), ")"))
  
  df_ratio |>
    select(sti, years_for_ratio, est)
  
}

# BETWEEN STUDY ANALYSIS ----

df_ct <- df |> filter(sti == "CT")
df_ng <- df |> filter(sti == "NG")
df_tv <- df |> filter(sti == "TV")

# < Main model ----

# Formula
form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  sex:region + population + age_group + hiv_status + (1 | study_id)

# Logit models
modct1 <- glmmTMB(form, data = df_ct, family = binomial(link="logit"))
modng1 <- glmmTMB(form, data = df_ng, family = binomial(link="logit"))
modtv1 <- glmmTMB(form, data = df_tv, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct1)
theta_ct <- as.vector(modct1$fit$par[names(modct1$fit$par) == "theta"])
fixed_ct <- as.vector(modct1$fit$par[names(modct1$fit$par) == "beta"])
fixed_ct <- append(fixed_ct, 0, length(fixed_ct))

summary(modng1)
theta_ng <- as.vector(modng1$fit$par[names(modng1$fit$par) == "theta"])
fixed_ng <- as.vector(modng1$fit$par[names(modng1$fit$par) == "beta"])
fixed_ng <- append(fixed_ng, 0, length(fixed_ng))

summary(modtv1)
theta_tv <- as.vector(modtv1$fit$par[names(modtv1$fit$par) == "theta"])
fixed_tv <- as.vector(modtv1$fit$par[names(modtv1$fit$par) == "beta"])
fixed_tv <- append(fixed_tv, 0, length(fixed_tv) - 1)
fixed_tv <- append(fixed_tv, 0, length(fixed_tv))

# Log models
modct2 <- glmmTMB(form, data = df_ct, family = binomial(link="log"),
                  start = list(theta = theta_ct, beta = fixed_ct))

modng2 <- glmmTMB(form, data = df_ng, family = binomial(link="log"),
                  start = list(theta = theta_ng, beta = fixed_ng))

modtv2 <- glmmTMB(form, data = df_tv, family = binomial(link="log"),
                  start = list(theta = theta_tv, beta = fixed_tv))

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

# Compare 2020 and 2010
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
  mutate(region = factor(region, 
                                  levels = c("Western and Central", "Eastern", 
                                             "Southern", "SSA"))) |>
  arrange(year, region) |>
  write.csv("./results/prevalence_comparewho.csv", row.names = FALSE)

# Ratio between study
ratio_btwn <- rbind(
  prediction(modct2, "CT")$ssa_ratio,
  prediction(modng2, "NG")$ssa_ratio,
  prediction(modtv2, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Between study",
         n = c(modct2$modelInfo$nobs, modng2$modelInfo$nobs, modtv2$modelInfo$nobs))

# Mean aPR per year
prediction_2(modct2, "CT", y1 = "2015", y2 = "2016")
prediction_2(modng2, "NG", y1 = "2015", y2 = "2016")
prediction_2(modtv2, "TV", y1 = "2015", y2 = "2016")

# Mean aPR between 2008 and 2020
prediction_2(modct2, "CT", y1 = "2010", y2 = "2020")
prediction_2(modng2, "NG", y1 = "2010", y2 = "2020")
prediction_2(modtv2, "TV", y1 = "2010", y2 = "2020")


# < No diagnostic test adjustment ----

# Formula
form <- cbind(num,(denom-num)) ~ region + year_regression:region + sex:region +
  population + age_group + hiv_status + (1 | study_id)

# Logit models
modct3 <- glmmTMB(form, data = df_ct, family = binomial(link="logit"))
modng3 <- glmmTMB(form, data = df_ng, family = binomial(link="logit"))
modtv3 <- glmmTMB(form, data = df_tv, family = binomial(link="logit"))

# Use  parameter values from logit model as starting values for log model
summary(modct3)
theta_ct <- as.vector(modct3$fit$par[names(modct3$fit$par) == "theta"])
fixed_ct <- as.vector(modct3$fit$par[names(modct3$fit$par) == "beta"])
fixed_ct <- append(fixed_ct, 0, length(fixed_ct))

summary(modng3)
theta_ng <- as.vector(modng3$fit$par[names(modng3$fit$par) == "theta"])
fixed_ng <- as.vector(modng3$fit$par[names(modng3$fit$par) == "beta"])
fixed_ng <- append(fixed_ng, 0, length(fixed_ng))

summary(modtv3)
theta_tv <- as.vector(modtv3$fit$par[names(modtv3$fit$par) == "theta"])
fixed_tv <- as.vector(modtv3$fit$par[names(modtv3$fit$par) == "beta"])
fixed_tv <- append(fixed_tv, 0, length(fixed_tv)-1)
fixed_tv <- append(fixed_tv, 0, length(fixed_tv))

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

# < Women only ----

df_ct_fm <- df |> filter(sti == "CT", sex == "Female")
df_ng_fm <- df |> filter(sti == "NG", sex == "Female")
df_tv_fm <- df |> filter(sti == "TV", sex == "Female")

# Formula
form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  population + age_group + hiv_status + (1 | study_id)

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
  age_group + hiv_status + (1 | study_id)

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

# WITHIN STUDY ANALYSIS ----

sex_within <- df |> 
  left_join(read.csv("./data/study_name.csv"),
            by = c("study_id")) |>
  group_by(study_name, study_id, country, location, sti) |> 
  summarise(n=n_distinct(sex)) |> 
  filter(n>1) |>
  mutate(keep = "y")

df_sex_within <- df |>
  left_join(sex_within) |>
  filter(!is.na(keep))

df_ct_within <- df_sex_within |> filter(sti=="CT")
df_ng_within <- df_sex_within |> filter(sti=="NG")
df_tv_within <- df_sex_within |> filter(sti=="TV")

# < Main model ----

# original formula
# form <- region + year_regression:region + 
#   sex:region + population + age_group + hiv_status + (1 | study_id)

# Removed HIV-status 
# No HIV positive pops, and predominantly mixed HIV status
# TV dataset only mixed HIV status, so model does not fit with only one level of the variable

# Removed population - insufficient data

# Removed age_group - unable to calculate SE for TV with this included

form <- cbind(adj_num,(adj_denom-adj_num)) ~ region + year_regression:region + 
  sex:region + (1 | study_id)

# Logit models
modct1_w <- glmmTMB(form, data = df_ct_within, family = binomial(link="logit"))
modng1_w <- glmmTMB(form, data = df_ng_within, family = binomial(link="logit"))
modtv1_w <- glmmTMB(form, data = df_tv_within, family = binomial(link="logit"))

tab_model(modct1_w, modng1_w, modtv1_w, p.style = "stars", wrap.labels = 100)

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
fixed_tv_w <- append(fixed_tv_w, 0, length(fixed_tv_w)-3)
fixed_tv_w <- append(fixed_tv_w, 0, length(fixed_tv_w)-3)

# Log models
modct2_w <- glmmTMB(form, data = df_ct_within, family = binomial(link="log"),
                    start = list(theta = theta_ct_w, beta = fixed_ct_w))

modng2_w <- glmmTMB(form, data = df_ng_within, family = binomial(link="log"),
                    start = list(theta = theta_ng_w, beta = fixed_ng_w))

modtv2_w <- glmmTMB(form, data = df_tv_within, family = binomial(link="log"),
                    start = list(theta = theta_tv_w, beta = fixed_tv_w))

sjPlot::tab_model(modct2_w, modng2_w, modtv2_w, p.style = "stars", wrap.labels = 100)

# Ratio within study
ratio_within <- rbind(
  prediction(modct2_w, "CT")$ssa_ratio,
  prediction(modng2_w, "NG")$ssa_ratio,
  prediction(modtv2_w, "TV")$ssa_ratio) |>
  filter(year == 2020) |>
  mutate(type = "Within study",
         n = c(modct2_w$modelInfo$nobs, modng2_w$modelInfo$nobs, modtv2_w$modelInfo$nobs))

# < Pooled model (no covariates) ----

form <- cbind(adj_num,(adj_denom-adj_num)) ~ sex + (1 | study_id)

# Logit models

modct1_wn <- glmmTMB(form, data = df_ct_within, family = binomial(link="logit"))
modng1_wn <- glmmTMB(form, data = df_ng_within, family = binomial(link="logit"))
modtv1_wn <- glmmTMB(form, data = df_tv_within, family = binomial(link="logit"))

tab_model(modct1_wn, modng1_wn, modtv1_wn, p.style = "stars", wrap.labels = 100)

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
    mutate(year = 2015,
           ratio = exp(Estimate),
           lwr = exp(Estimate - 1.96*`Std. Error`),
           upr = exp(Estimate + 1.96*`Std. Error`),
           tau_study = insight::get_variance(mod)$var.random) |>
    select(year, ratio, lwr, upr, tau_study)
}

ratio_within_unadjusted <- rbind(
  extract_ratio(modct2_wn) |> mutate(sti = "CT"),
  extract_ratio(modng2_wn) |> mutate(sti = "NG"),
  extract_ratio(modtv2_wn) |> mutate(sti = "TV")) |> 
  relocate(sti) |>
  mutate(type = "Unadjusted",
         n = c(modct2_wn$modelInfo$nobs, modng2_wn$modelInfo$nobs, modtv2_wn$modelInfo$nobs))

# ALL RATIOS ---- 

rbind(ratio_btwn,
      ratio_within,
      ratio_within_unadjusted) |>
  write.csv("./results/ratios.csv", row.names = FALSE)

# SAVE REGRESSION TABLES ----

# variable order
desired_order <- c("(Intercept)","Western and Central","Eastern", "Southern",
                   "Western and Central:Male", "Eastern:Male","Southern:Male", 
                   "Western and Central:Both sexes", "Eastern:Both sexes","Southern:Both sexes", 
                   "Western and Central:Year", "Eastern:Year","Southern:Year", 
                   "ANC attendees", "FP attendees", "GYN attendees", 
                   "PHC/OPD attendees", "Students", "Community members", 
                   "HIV/STI prevention trial participants", "Population-representative survey participants", 
                   "Female","Male","Both sexes", "Adult", "Youth",
                   "Mixed", "HIV negative", "HIV positive")

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
           variable = gsub("sexBoth sexes","Both sexes",variable),
           variable = gsub("age_group","",variable),
           variable = gsub("hiv_status","",variable)) |>
    # add reference categories 
    rbind(data.frame(variable = c("ANC attendees","Adult", "Mixed", "Southern"),
                     estimate = rep("Ref",4)))
}

# between study
left_join(
  reg_table(modct2) |> rename(CT = estimate),
  reg_table(modng2) |> rename(NG = estimate)) |>
  left_join(reg_table(modtv2) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_alldata.csv", row.names = FALSE)

insight::get_variance(modct2)
insight::get_variance(modng2)
insight::get_variance(modtv2)

sjPlot::tab_model(modct2, modng2, modtv2, p.style = "stars", wrap.labels = 100)

# between study - unadjusted
left_join(
  reg_table(modct4) |> rename(CT = estimate),
  reg_table(modng4) |> rename(NG = estimate)) |>
  left_join(reg_table(modtv4) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_alldata_unadjusted.csv", row.names = FALSE)

insight::get_variance(modct4)
insight::get_variance(modng4)
insight::get_variance(modtv4)

sjPlot::tab_model(modct4, modng4, modtv4, p.style = "stars", wrap.labels = 100)

# between study - female pops
left_join(
  reg_table(modct6) |> rename(CT = estimate),
  reg_table(modng6) |> rename(NG = estimate)) |>
  left_join(reg_table(modtv6) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_alldata_female.csv", row.names = FALSE)

insight::get_variance(modct6)
insight::get_variance(modng6)
insight::get_variance(modtv6)

sjPlot::tab_model(modct6, modng6, modtv6, p.style = "stars", wrap.labels = 100)


# between study - ANC
left_join(
  reg_table(modct8) |> rename(CT = estimate),
  reg_table(modng8) |> rename(NG = estimate)) |>
  left_join(reg_table(modtv8) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |> 
  write.csv("./tables/regression_alldata_ANC.csv", row.names = FALSE)

insight::get_variance(modct8)
insight::get_variance(modng8)
insight::get_variance(modtv8)

sjPlot::tab_model(modct8, modng8, modtv8, p.style = "stars", wrap.labels = 100)


# within study
left_join(
  reg_table(modct2_w) |> rename(CT = estimate),
  reg_table(modng2_w) |> rename(NG = estimate)) |>
  left_join(reg_table(modtv2_w) |> rename(TV = estimate)) |>
  # reorder variables
  mutate(variable = factor(variable, levels = desired_order)) |>
  arrange(variable)  |>
  write.csv("./tables/regression_within.csv", row.names = FALSE)

insight::get_variance(modct2_w)
insight::get_variance(modng2_w)
insight::get_variance(modtv2_w)

sjPlot::tab_model(modct2_w, modng2_w, modtv2_w, p.style = "stars", wrap.labels = 100)

# STUDY CHARACTERISTICS ----

# Number articles
df |> 
  group_by() |>
  summarise(n = n_distinct(study_id))

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
df_study <- df |> 
  left_join(read.csv("./data/study_name.csv"),
            by = c("study_id"))

t_characteristics_all <- full_join(
  tab_characteristics(df_study |> filter(sti == "CT")) |> rename(CT = col, CT_tot = tot),
  tab_characteristics(df_study |> filter(sti == "NG")) |> rename(NG = col, NG_tot = tot)) |>
  full_join(tab_characteristics(df_study |> filter(sti == "TV")) |> rename(TV = col, TV_tot = tot)) |>
  full_join(tab_characteristics(df_study) |> rename(ALL = col, ALL_tot = tot)) |>
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
  group_by(country) |>
  summarise(n = n()) |>
  mutate(country_list = strsplit(country, ", ")) |>
  unnest(country_list) |>
  group_by(country_list) |>
  summarise(n=sum(n)) |>
  arrange(-n)

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
  df_study |> 
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df_study |> 
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
  mutate(prop = n_tot/124)

# Number studies per country - CT
full_join(
  df_study |>
    filter(sti == "CT") |>
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df_study |>
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
  df_study |>
    filter(sti == "NG") |>
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df_study |>
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
  df_study |>
    filter(sti == "TV") |> 
    filter(!is.na(study_name)) |>
    group_by(country) |>
    summarise(n = n_distinct(study_name)) |>
    mutate(country_list = strsplit(country, ", ")) |>
    unnest(country_list) |>
    group_by(country_list) |>
    summarise(n_study=sum(n)),
  df_study |>
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

# Number studies among both men and women
full_join(
  df_study |> 
    filter(!is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_study=n_distinct(study_name)),
  df_study |> 
    filter(is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# Number studies among both men and women - total
cbind(
  df_study |> 
    filter(!is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_study=n_distinct(study_name)),
  df_study |> 
    filter(is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(!is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# Number studies among both sexes (unstratified)
df_study |> 
  select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
  pivot_wider(names_from = sex, values_from = adj_prev) |>
  filter(!is.na(`Both sexes`)) |>
  group_by(sti) |>
  summarise(n=n_distinct(study_id))

# Number studies among both sexes (unstratified) - total
df_study |> 
  select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
  pivot_wider(names_from = sex, values_from = adj_prev) |>
  filter(!is.na(`Both sexes`)) |>
  group_by() |>
  summarise(n=n_distinct(study_id))

# Number studies among men only
df_study |> 
  select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
  pivot_wider(names_from = sex, values_from = adj_prev) |>
  filter(!is.na(Male) & is.na(Female)) |>
  group_by() |>
  summarise(n=n_distinct(study_id))

# Number studies among women only
full_join(
  df_study |> 
    filter(!is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_study=n_distinct(study_name)),
  df_study |> 
    filter(is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by(sti) |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)

# Number studies among women only - total
cbind(
  df_study |> 
    filter(!is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_study=n_distinct(study_name)),
  df_study |> 
    filter(is.na(study_name)) |>
    select(study_name, study_id, country, location, population, sex, sti, adj_prev) |> 
    pivot_wider(names_from = sex, values_from = adj_prev) |>
    filter(is.na(Male) & !is.na(Female)) |>
    group_by() |>
    summarise(n_art=n_distinct(study_id))) |>
  mutate(n_tot = n_study + n_art)
