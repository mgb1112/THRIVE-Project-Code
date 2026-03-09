library(tidyverse)
library(here)
library(sandwich)
library(lmtest)

output_dir <- here::here("Paper_Model", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spine <- readRDS(file.path(output_dir, "01_analysis_spine_wide.rds"))
sbp_df <- readRDS(file.path(output_dir, "01_sbp_followup_long.rds"))
dash_df <- readRDS(file.path(output_dir, "01_dash_followup_long.rds"))

# Baseline SBP predictor (continuous) for both missingness outcomes.
sbp_df <- sbp_df %>%
  mutate(sbp_bl = as.numeric(y_baseline))

dash_df <- dash_df %>%
  left_join(
    spine %>% transmute(case = as.character(case), sbp_bl = as.numeric(sbp_bl)),
    by = "case"
  )

if (any(is.na(sbp_df$sbp_bl)) || any(is.na(dash_df$sbp_bl))) {
  stop("Baseline SBP is missing for one or more rows in missingness datasets.")
}

sbp_bl_center <- mean(spine$sbp_bl, na.rm = TRUE)

sbp_df <- sbp_df %>%
  mutate(sbp_bl_c = sbp_bl - sbp_bl_center)

dash_df <- dash_df %>%
  mutate(sbp_bl_c = sbp_bl - sbp_bl_center)

calc_visit_or <- function(coef_tbl, vcov_mat, outcome_label, n_participants, n_rows) {
  b <- coef_tbl
  V <- vcov_mat

  sbp_term <- "sbp_bl_c"
  int_term_a <- "time24W:sbp_bl_c"
  int_term_b <- "sbp_bl_c:time24W"
  int_term <- if (int_term_a %in% names(b)) int_term_a else int_term_b

  if (!(sbp_term %in% names(b))) {
    stop("Expected baseline SBP term not found in missingness model fixed effects.")
  }

  est12 <- b[[sbp_term]]
  se12 <- sqrt(V[sbp_term, sbp_term])

  int_est <- if (int_term %in% names(b)) b[[int_term]] else 0
  int_var <- if (int_term %in% colnames(V)) V[int_term, int_term] else 0
  covar <- if ((sbp_term %in% colnames(V)) && (int_term %in% colnames(V))) V[sbp_term, int_term] else 0

  est24 <- est12 + int_est
  se24 <- sqrt(V[sbp_term, sbp_term] + int_var + 2 * covar)

  z12 <- est12 / se12
  z24 <- est24 / se24

  tibble(
    outcome = outcome_label,
    visit = c("12W", "24W"),
    effect = c(
      "OR per +1 mmHg baseline SBP (centered) at 12W",
      "OR per +1 mmHg baseline SBP (centered) at 24W"
    ),
    log_odds = c(est12, est24),
    std.error = c(se12, se24),
    odds_ratio = exp(c(est12, est24)),
    conf.low = exp(c(est12 - 1.96 * se12, est24 - 1.96 * se24)),
    conf.high = exp(c(est12 + 1.96 * se12, est24 + 1.96 * se24)),
    p.value = c(2 * pnorm(-abs(z12)), 2 * pnorm(-abs(z24))),
    n_participants = n_participants,
    n_rows = n_rows
  )
}

fit_missingness_model <- function(df, outcome_label) {
  # Missingness detection with centered continuous baseline SBP:
  # miss_y ~ time * sbp_bl_c for follow-up visits (12W and 24W).
  # We use participant-clustered robust SE to account for repeated binary outcomes.
  model_df <- df
  form <- miss_y ~ time * sbp_bl_c

  fit <- glm(form, data = model_df, family = binomial())
  V_robust <- sandwich::vcovCL(fit, cluster = model_df$case, type = "HC0")
  ct <- lmtest::coeftest(fit, vcov. = V_robust)
  ct_df <- as.data.frame(unclass(ct))
  ct_df$term <- rownames(ct_df)
  rownames(ct_df) <- NULL

  names(ct_df) <- c("estimate", "std.error", "statistic", "p.value", "term")
  ct_df <- as_tibble(ct_df)

  # robust Wald CI
  ct_df <- ct_df %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    )

  terms_tbl <- ct_df %>%
    mutate(
      outcome = outcome_label,
      odds_ratio = exp(estimate),
      conf.low.or = exp(conf.low),
      conf.high.or = exp(conf.high)
    ) %>%
    select(
      outcome, term, estimate, std.error, statistic, p.value,
      conf.low, conf.high, odds_ratio, conf.low.or, conf.high.or
    )

  visit_tbl <- calc_visit_or(
    coef_tbl = stats::coef(fit),
    vcov_mat = V_robust,
    outcome_label = outcome_label,
    n_participants = n_distinct(model_df$case),
    n_rows = nrow(model_df)
  )

  fit_info <- tibble(
    outcome = outcome_label,
    n_participants = n_distinct(model_df$case),
    n_rows = nrow(model_df),
    formula = deparse(form),
    sbp_baseline_center_value = sbp_bl_center,
    model_type = "Binomial GLM (time*centered baseline SBP) with cluster-robust SE by participant"
  )

  list(fit = fit, terms = terms_tbl, visit = visit_tbl, info = fit_info)
}

sbp_fit <- fit_missingness_model(sbp_df, "SBP")
dash_fit <- fit_missingness_model(dash_df, "DASH")

terms_out <- bind_rows(sbp_fit$terms, dash_fit$terms)
visit_out <- bind_rows(sbp_fit$visit, dash_fit$visit)
fit_info <- bind_rows(sbp_fit$info, dash_fit$info)

readr::write_csv(terms_out, file.path(output_dir, "03_missingness_model_terms.csv"))
readr::write_csv(visit_out, file.path(output_dir, "03_missingness_visit_or.csv"))
readr::write_csv(fit_info, file.path(output_dir, "03_missingness_fit_info.csv"))

saveRDS(list(SBP = sbp_fit$fit, DASH = dash_fit$fit), file.path(output_dir, "03_missingness_model_fits.rds"))

message("03_missingness_detection completed. Outputs written to: ", output_dir)
