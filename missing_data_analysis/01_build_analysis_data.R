library(tidyverse)
library(here)
library(janitor)

# Build harmonized analysis datasets for unadjusted paper-model workflows,
# deriving DASH nutrient score from recall data using Mellen/Stata logic.

output_dir <- here::here("Paper_Model", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(here::here("Data", "thrive_data.rds"))

first_non_missing_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

first_non_missing_num <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  x[[1]]
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

score_mellen <- function(v, full, mid, direction = c("low", "high")) {
  direction <- match.arg(direction)
  if (direction == "low") {
    ifelse(is.na(v), NA_real_, ifelse(v <= full, 1, ifelse(v <= mid, 0.5, 0)))
  } else {
    ifelse(is.na(v), NA_real_, ifelse(v >= full, 1, ifelse(v >= mid, 0.5, 0)))
  }
}

# Randomized ITT cohort
rand <- dat$`THRIVE Randomization Data`$Sheet1 %>%
  clean_names() %>%
  transmute(
    case = as.character(caselabel),
    arm_raw = as.character(randomization),
    arm = case_when(
      stringr::str_detect(arm_raw, "Intervention") ~ "Intervention",
      TRUE ~ "Control"
    )
  ) %>%
  mutate(arm = factor(arm, levels = c("Control", "Intervention"))) %>%
  select(case, arm) %>%
  distinct()

stopifnot(n_distinct(rand$case) == 80L)

# Baseline covariates from sociodemographic form
sd <- dat$`Sociodemographic Questionnaire` %>%
  clean_names() %>%
  transmute(
    case = as.character(case),
    demo_age,
    demo_sexatbirth,
    demo_ishispaniclatino
  ) %>%
  group_by(case) %>%
  summarise(
    age = first_non_missing_num(demo_age),
    sex = first_non_missing_chr(demo_sexatbirth),
    hispanic = first_non_missing_chr(demo_ishispaniclatino),
    .groups = "drop"
  ) %>%
  mutate(
    sex = factor(if_else(sex %in% c("Female", "Male"), sex, NA_character_), levels = c("Female", "Male")),
    hispanic = factor(if_else(hispanic %in% c("No", "Yes"), hispanic, NA_character_), levels = c("No", "Yes"))
  )

# SBP outcomes
bp_long <- dat$`Blood Pressure Measurement Form` %>%
  clean_names() %>%
  transmute(
    case = as.character(case),
    visit = as.character(bp_event),
    sbp = rowMeans(cbind(bp_sbp1, bp_sbp2, bp_sbp3), na.rm = TRUE)
  ) %>%
  mutate(sbp = if_else(is.nan(sbp), NA_real_, sbp)) %>%
  filter(case %in% rand$case, visit %in% c("Baseline", "12W", "24W")) %>%
  group_by(case, visit) %>%
  summarise(sbp = first_non_missing_num(sbp), .groups = "drop")

sbp_wide <- bp_long %>%
  pivot_wider(names_from = visit, values_from = sbp) %>%
  rename(
    sbp_bl = Baseline,
    sbp_12 = `12W`,
    sbp_24 = `24W`
  )

# DASH nutrient score (Mellen): recall-level from FIMrecalls04
visit_map <- c("0" = "Baseline", "12" = "12W", "24" = "24W")

fim04 <- dat$FIMrecalls04 %>%
  clean_names() %>%
  transmute(
    case = as.character(participant_id),
    visit_num = as.character(visit_number),
    kcal = as.numeric(energy_kcal),
    caloriesfromsfa = as.numeric(percent_calories_from_sfa),
    caloriesfromfat = as.numeric(percent_calories_from_fat),
    caloriesfromprotein = as.numeric(percent_calories_from_protein),
    cholesterolmg = as.numeric(cholesterol_mg),
    totaldietaryfiberg = as.numeric(total_dietary_fiber_g),
    magnesiummg = as.numeric(magnesium_mg),
    calciummg = as.numeric(calcium_mg),
    potassiummg = as.numeric(potassium_mg),
    sodiummg = as.numeric(sodium_mg)
  ) %>%
  filter(case %in% rand$case, visit_num %in% names(visit_map), !is.na(kcal), kcal > 0)

dash_recall <- fim04 %>%
  mutate(
    chol_per1000kcal = cholesterolmg / kcal * 1000,
    fiber_per1000kcal = totaldietaryfiberg / kcal * 1000,
    magnesium_per1000kcal = magnesiummg / kcal * 1000,
    calcium_per1000kcal = calciummg / kcal * 1000,
    potassium_per1000kcal = potassiummg / kcal * 1000,
    sodium_per1000kcal = sodiummg / kcal * 1000,
    sfa_dash = score_mellen(caloriesfromsfa, full = 6, mid = 11, direction = "low"),
    fat_dash = score_mellen(caloriesfromfat, full = 27, mid = 32, direction = "low"),
    protein_dash = score_mellen(caloriesfromprotein, full = 18, mid = 16.5, direction = "high"),
    cholesterol_dash = score_mellen(chol_per1000kcal, full = 71.4, mid = 107.1, direction = "low"),
    fiber_dash = score_mellen(fiber_per1000kcal, full = 14.8, mid = 9.5, direction = "high"),
    magnesium_dash = score_mellen(magnesium_per1000kcal, full = 238, mid = 158, direction = "high"),
    calcium_dash = score_mellen(calcium_per1000kcal, full = 590, mid = 402, direction = "high"),
    potassium_dash = score_mellen(potassium_per1000kcal, full = 2238, mid = 1534, direction = "high"),
    sodium_dash = score_mellen(sodium_per1000kcal, full = 1143, mid = 1286, direction = "low"),
    dash_nutrient = sfa_dash + fat_dash + protein_dash + cholesterol_dash + fiber_dash +
      magnesium_dash + calcium_dash + potassium_dash + sodium_dash,
    visit = recode(visit_num, !!!visit_map)
  )

dash_visit <- dash_recall %>%
  group_by(case, visit_num, visit) %>%
  summarise(
    dash_score = mean(dash_nutrient, na.rm = TRUE),
    n_recalls_visit = n(),
    dash_category = if_else(dash_score < 4.5, "Low", "High"),
    .groups = "drop"
  )

dash_wide <- dash_visit %>%
  select(case, visit, dash_score) %>%
  pivot_wider(names_from = visit, values_from = dash_score) %>%
  rename(
    dash_bl = Baseline,
    dash_12 = `12W`,
    dash_24 = `24W`
  )


ndsr_diff <- NA_real_
ndsr_n_common <- 0L
if ("NDSR_recalls04" %in% names(dat)) {
  ndsr04 <- dat$NDSR_recalls04 %>%
    clean_names() %>%
    transmute(
      case = as.character(participant_id),
      visit_num = as.character(visit_number),
      kcal = as.numeric(energy_kcal),
      caloriesfromsfa = as.numeric(percent_calories_from_sfa),
      caloriesfromfat = as.numeric(percent_calories_from_fat),
      caloriesfromprotein = as.numeric(percent_calories_from_protein),
      cholesterolmg = as.numeric(cholesterol_mg),
      totaldietaryfiberg = as.numeric(total_dietary_fiber_g),
      magnesiummg = as.numeric(magnesium_mg),
      calciummg = as.numeric(calcium_mg),
      potassiummg = as.numeric(potassium_mg),
      sodiummg = as.numeric(sodium_mg)
    ) %>%
    filter(case %in% rand$case, visit_num %in% names(visit_map), !is.na(kcal), kcal > 0) %>%
    mutate(
      chol_per1000kcal = cholesterolmg / kcal * 1000,
      fiber_per1000kcal = totaldietaryfiberg / kcal * 1000,
      magnesium_per1000kcal = magnesiummg / kcal * 1000,
      calcium_per1000kcal = calciummg / kcal * 1000,
      potassium_per1000kcal = potassiummg / kcal * 1000,
      sodium_per1000kcal = sodiummg / kcal * 1000,
      dash_nutrient =
        score_mellen(caloriesfromsfa, 6, 11, "low") +
        score_mellen(caloriesfromfat, 27, 32, "low") +
        score_mellen(caloriesfromprotein, 18, 16.5, "high") +
        score_mellen(chol_per1000kcal, 71.4, 107.1, "low") +
        score_mellen(fiber_per1000kcal, 14.8, 9.5, "high") +
        score_mellen(magnesium_per1000kcal, 238, 158, "high") +
        score_mellen(calcium_per1000kcal, 590, 402, "high") +
        score_mellen(potassium_per1000kcal, 2238, 1534, "high") +
        score_mellen(sodium_per1000kcal, 1143, 1286, "low")
    ) %>%
    group_by(case, visit_num) %>%
    summarise(dash_ndsr = mean(dash_nutrient, na.rm = TRUE), .groups = "drop")

  fim_visit_cmp <- dash_visit %>% select(case, visit_num, dash_fim = dash_score)
  cmp <- fim_visit_cmp %>% inner_join(ndsr04, by = c("case", "visit_num"))
  ndsr_n_common <- nrow(cmp)
  if (ndsr_n_common > 0) ndsr_diff <- max(abs(cmp$dash_fim - cmp$dash_ndsr), na.rm = TRUE)
}

# Unified participant-level spine
spine <- rand %>%
  left_join(sd, by = "case") %>%
  left_join(sbp_wide, by = "case") %>%
  left_join(dash_wide, by = "case")

stopifnot(nrow(spine) == 80L)
stopifnot(n_distinct(spine$case) == 80L)

# Follow-up long datasets (used for missingness detection at 12W/24W)
sbp_followup <- spine %>%
  transmute(
    case, arm, age, sex, hispanic,
    y_baseline = sbp_bl,
    y_12 = sbp_12,
    y_24 = sbp_24
  ) %>%
  pivot_longer(cols = c(y_12, y_24), names_to = "time", values_to = "y") %>%
  mutate(
    time = recode(time, y_12 = "12W", y_24 = "24W"),
    time = factor(time, levels = c("12W", "24W")),
    miss_y = as.integer(is.na(y)),
    outcome = "SBP"
  )

dash_followup <- spine %>%
  transmute(
    case, arm, age, sex, hispanic,
    y_baseline = dash_bl,
    y_12 = dash_12,
    y_24 = dash_24
  ) %>%
  pivot_longer(cols = c(y_12, y_24), names_to = "time", values_to = "y") %>%
  mutate(
    time = recode(time, y_12 = "12W", y_24 = "24W"),
    time = factor(time, levels = c("12W", "24W")),
    miss_y = as.integer(is.na(y)),
    outcome = "DASH"
  )

# Three-visit longitudinal datasets (Baseline/12W/24W) for unadjusted paper model
sbp_long3 <- spine %>%
  transmute(
    case, arm,
    y_Baseline = sbp_bl,
    y_12W = sbp_12,
    y_24W = sbp_24
  ) %>%
  pivot_longer(cols = starts_with("y_"), names_to = "time", values_to = "y") %>%
  mutate(
    time = recode(time, y_Baseline = "Baseline", y_12W = "12W", y_24W = "24W"),
    time = factor(time, levels = c("Baseline", "12W", "24W")),
    miss_y = as.integer(is.na(y)),
    outcome = "SBP"
  )

dash_long3 <- spine %>%
  transmute(
    case, arm,
    y_Baseline = dash_bl,
    y_12W = dash_12,
    y_24W = dash_24
  ) %>%
  pivot_longer(cols = starts_with("y_"), names_to = "time", values_to = "y") %>%
  mutate(
    time = recode(time, y_Baseline = "Baseline", y_12W = "12W", y_24W = "24W"),
    time = factor(time, levels = c("Baseline", "12W", "24W")),
    miss_y = as.integer(is.na(y)),
    outcome = "DASH"
  )

covariate_missing <- spine %>%
  summarise(
    n_rand = n(),
    miss_age = sum(is.na(age)),
    miss_sex = sum(is.na(sex)),
    miss_hispanic = sum(is.na(hispanic)),
    miss_sbp_bl = sum(is.na(sbp_bl)),
    miss_dash_bl = sum(is.na(dash_bl))
  )

followup_counts <- bind_rows(sbp_followup, dash_followup) %>%
  group_by(outcome, arm, time) %>%
  summarise(
    n_total = n(),
    n_observed = sum(!is.na(y)),
    pct_observed = 100 * mean(!is.na(y)),
    .groups = "drop"
  )

# DASH component/total summaries and QC checks
dash_component_summary <- dash_recall %>%
  summarise(
    n_recall_rows = n(),
    mean_dash_nutrient = mean(dash_nutrient, na.rm = TRUE),
    sd_dash_nutrient = sd(dash_nutrient, na.rm = TRUE),
    mean_sfa_dash = mean(sfa_dash, na.rm = TRUE),
    mean_fat_dash = mean(fat_dash, na.rm = TRUE),
    mean_protein_dash = mean(protein_dash, na.rm = TRUE),
    mean_cholesterol_dash = mean(cholesterol_dash, na.rm = TRUE),
    mean_fiber_dash = mean(fiber_dash, na.rm = TRUE),
    mean_magnesium_dash = mean(magnesium_dash, na.rm = TRUE),
    mean_calcium_dash = mean(calcium_dash, na.rm = TRUE),
    mean_potassium_dash = mean(potassium_dash, na.rm = TRUE),
    mean_sodium_dash = mean(sodium_dash, na.rm = TRUE)
  )

qc_checks <- bind_rows(
  tibble(check = "rand_n", value = n_distinct(rand$case), pass = n_distinct(rand$case) == 80),
  tibble(check = "dash_total_min", value = min(dash_recall$dash_nutrient, na.rm = TRUE), pass = min(dash_recall$dash_nutrient, na.rm = TRUE) >= 0),
  tibble(check = "dash_total_max", value = max(dash_recall$dash_nutrient, na.rm = TRUE), pass = max(dash_recall$dash_nutrient, na.rm = TRUE) <= 9),
  tibble(check = "dash_total_half_step", value = sum(abs(dash_recall$dash_nutrient * 2 - round(dash_recall$dash_nutrient * 2)) > 1e-8, na.rm = TRUE), pass = sum(abs(dash_recall$dash_nutrient * 2 - round(dash_recall$dash_nutrient * 2)) > 1e-8, na.rm = TRUE) == 0),
  tibble(check = "visit_n_baseline", value = sum(dash_visit$visit_num == "0"), pass = sum(dash_visit$visit_num == "0") == 80),
  tibble(check = "visit_n_12", value = sum(dash_visit$visit_num == "12"), pass = sum(dash_visit$visit_num == "12") == 51),
  tibble(check = "visit_n_24", value = sum(dash_visit$visit_num == "24"), pass = sum(dash_visit$visit_num == "24") == 55),
  tibble(check = "baseline_dash_mean", value = mean(dash_visit$dash_score[dash_visit$visit_num == "0"], na.rm = TRUE), pass = TRUE),
  tibble(check = "fim_ndsr_max_abs_diff", value = ndsr_diff, pass = if_else(is.na(ndsr_diff), TRUE, ndsr_diff < 1e-8)),
  tibble(check = "fim_ndsr_common_rows", value = ndsr_n_common, pass = TRUE)
)

# Save outputs
readr::write_csv(spine, file.path(output_dir, "01_analysis_spine_wide.csv"))
readr::write_csv(sbp_followup, file.path(output_dir, "01_sbp_followup_long.csv"))
readr::write_csv(dash_followup, file.path(output_dir, "01_dash_followup_long.csv"))
readr::write_csv(sbp_long3, file.path(output_dir, "01_sbp_long3.csv"))
readr::write_csv(dash_long3, file.path(output_dir, "01_dash_long3.csv"))
readr::write_csv(covariate_missing, file.path(output_dir, "01_covariate_missingness.csv"))
readr::write_csv(followup_counts, file.path(output_dir, "01_followup_counts.csv"))
readr::write_csv(dash_recall, file.path(output_dir, "01_dash_recall_level.csv"))
readr::write_csv(dash_visit, file.path(output_dir, "01_dash_visit_level.csv"))
readr::write_csv(dash_component_summary, file.path(output_dir, "01_dash_component_summary.csv"))
readr::write_csv(qc_checks, file.path(output_dir, "01_dash_qc_checks.csv"))

saveRDS(spine, file.path(output_dir, "01_analysis_spine_wide.rds"))
saveRDS(sbp_followup, file.path(output_dir, "01_sbp_followup_long.rds"))
saveRDS(dash_followup, file.path(output_dir, "01_dash_followup_long.rds"))
saveRDS(sbp_long3, file.path(output_dir, "01_sbp_long3.rds"))
saveRDS(dash_long3, file.path(output_dir, "01_dash_long3.rds"))

message("01_build_analysis_data completed. Outputs written to: ", output_dir)
