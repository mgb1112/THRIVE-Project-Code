library(tidyverse)
library(here)
library(lme4)
library(broom.mixed)

output_dir <- here::here("Paper_Model", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sbp_df <- readRDS(file.path(output_dir, "01_sbp_long3.rds"))
dash_df <- readRDS(file.path(output_dir, "01_dash_long3.rds"))

calc_visit_effects <- function(fit, outcome_label, n_participants, n_rows) {
  b <- lme4::fixef(fit)
  V <- as.matrix(stats::vcov(fit))

  int12_term <- "armIntervention:time12W"
  int24_term <- "armIntervention:time24W"

  if (!(int12_term %in% names(b)) || !(int24_term %in% names(b))) {
    stop("Expected interaction terms not found in fixed effects.")
  }

  est12 <- b[[int12_term]]
  est24 <- b[[int24_term]]

  se12 <- sqrt(V[int12_term, int12_term])
  se24 <- sqrt(V[int24_term, int24_term])

  z12 <- est12 / se12
  z24 <- est24 / se24

  tibble(
    outcome = outcome_label,
    visit = c("12W", "24W"),
    effect = c(
      "Interaction (Intervention change vs Control change) at 12W",
      "Interaction (Intervention change vs Control change) at 24W"
    ),
    estimate = c(est12, est24),
    std.error = c(se12, se24),
    conf.low = c(est12 - 1.96 * se12, est24 - 1.96 * se24),
    conf.high = c(est12 + 1.96 * se12, est24 + 1.96 * se24),
    p.value = c(2 * pnorm(-abs(z12)), 2 * pnorm(-abs(z24))),
    n_participants = n_participants,
    n_rows = n_rows
  )
}

fit_outcome_model <- function(df, outcome_label) {
  # Unadjusted paper model:
  # y ~ arm * time + (1|case), with Baseline as the reference visit.
  model_df <- df %>%
    filter(!is.na(y))

  form <- y ~ arm * time + (1 | case)

  fit <- lmer(form, data = model_df, REML = FALSE)

  terms_tbl <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE, conf.method = "Wald")
  if (!("p.value" %in% names(terms_tbl))) {
    terms_tbl <- terms_tbl %>%
      mutate(p.value = 2 * pnorm(-abs(statistic)))
  }
  terms_tbl <- terms_tbl %>%
    mutate(outcome = outcome_label) %>%
    select(outcome, term, estimate, std.error, statistic, p.value, conf.low, conf.high)

  visit_tbl <- calc_visit_effects(
    fit,
    outcome_label = outcome_label,
    n_participants = n_distinct(model_df$case),
    n_rows = nrow(model_df)
  )

  fit_info <- tibble(
    outcome = outcome_label,
    n_participants = n_distinct(model_df$case),
    n_rows = nrow(model_df),
    formula = deparse(form),
    model_type = "Unadjusted linear mixed-effects model"
  )

  list(fit = fit, terms = terms_tbl, visit = visit_tbl, info = fit_info)
}

sbp_fit <- fit_outcome_model(sbp_df, "SBP")
dash_fit <- fit_outcome_model(dash_df, "DASH")

terms_out <- bind_rows(sbp_fit$terms, dash_fit$terms)
visit_out <- bind_rows(sbp_fit$visit, dash_fit$visit)
fit_info <- bind_rows(sbp_fit$info, dash_fit$info)

readr::write_csv(terms_out, file.path(output_dir, "02_mixed_model_terms.csv"))
readr::write_csv(visit_out, file.path(output_dir, "02_mixed_model_visit_effects.csv"))
readr::write_csv(fit_info, file.path(output_dir, "02_mixed_model_fit_info.csv"))

saveRDS(list(SBP = sbp_fit$fit, DASH = dash_fit$fit), file.path(output_dir, "02_mixed_model_fits.rds"))

message("02_fit_paper_mixed_models completed. Outputs written to: ", output_dir)
