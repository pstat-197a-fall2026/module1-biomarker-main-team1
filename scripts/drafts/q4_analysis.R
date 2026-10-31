# Question 4: aligned class-baseline comparison
# Run from the project root: source("scripts/drafts/q4_analysis.R") or click Source.
# No teammate scripts need to run first; this reads the SAME clean data file.
# Keep scripts/drafts/q4_writeup.R alongside this file (our report helper).
#
# Alignment: same biomarker-clean.RData, ASD/TD forest labels, seed 101422,
# 80/20 split, logistic outcome ASD=TRUE, and fixed probability threshold >0.5.
# There is ONE final class benchmark, reproduced here with full-data selection.
# Q4 methods select/tune on training data. They reuse the team's global
# preprocessing for direct alignment; that retains preprocessing leakage.
# The test set was viewed before these corrections: all test results remain
# exploratory. Corrections fix alignment; they are not a fresh validation study.
#
# The 24 configurations and selection rules are unchanged. We do not expand
# the search or change the seed to chase a historical accuracy or a better score.
# Wait for all 15 validation folds, then inspect cv_summary, final_results,
# alternative_panel, simple_panel and inclass_panel. The report is regenerated too.

library(tidyverse)
library(randomForest)
library(rsample)
library(yardstick)
library(glmnet)

## 1. Settings and the team's supplied clean data -------------------------------------------------

project_root <- local({
  candidates <- c(".", "..", "../..")
  matches <- candidates[file.exists(file.path(candidates, "data", "biomarker-clean.RData"))]
  if (!length(matches)) stop("Open the team project before running this script.")
  normalizePath(matches[[1]], mustWork = TRUE)
})
output_dir <- file.path(project_root, "results", "sou-q4-final")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 101422
rf_trees <- 1000
cv_repeats <- 3
accuracy_margin <- 0.02  # Descriptive comparability margin, NOT equivalence proof.
coefficient_tolerance <- 1e-8 # Do not count floating-point residues as proteins.

# Specify the search BEFORE running it; do not change it to chase test scores.
# lambda_fraction is relative to the training fit's largest lambda. A larger
# penalty generally gives fewer nonzero protein coefficients.
configurations <- bind_rows(
  tibble(model = c("hard10", "hard25"), method = "glm",
         pool = c("hard10", "hard25"), k = NA_integer_, lambda_fraction = NA_real_),
  tibble(model = paste0("fuzzy", c(1, 2, 3, 5)), method = "glm",
         pool = "fuzzy", k = c(1L, 2L, 3L, 5L), lambda_fraction = NA_real_),
  crossing(pool = c("ttest", "rf"), k = 1:3) %>%
    mutate(model = paste0(pool, k), method = "glm", lambda_fraction = NA_real_),
  crossing(pool = c("fuzzy25", "all"),
           lambda_fraction = c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)) %>%
    mutate(model = paste0("lasso_", pool, "_", lambda_fraction),
           method = "lasso", k = NA_integer_)
)

clean_path <- file.path(project_root, "data", "biomarker-clean.RData")
input <- new.env(parent = emptyenv())
load(clean_path, envir = input)
analysis_data <- input$biomarker_clean %>%
  select(-ados) %>% mutate(subject_id = row_number(), .before = group)
protein_names <- setdiff(names(analysis_data), c("subject_id", "group"))
stopifnot(nrow(analysis_data) == 154L, length(protein_names) == 1317L,
          !anyNA(analysis_data[, protein_names]),
          all(as.matrix(analysis_data[, protein_names]) >= -3),
          all(as.matrix(analysis_data[, protein_names]) <= 3))

# Same sampling call, seed, row order and proportion as the team's Q3 scripts.
set.seed(seed)
outer_split <- initial_split(analysis_data, prop = 0.8)
train_data <- training(outer_split)
test_data <- testing(outer_split)
cat("Subjects:", nrow(train_data), "training;", nrow(test_data), "testing\n")

## 2. Extract model inputs -- no second transformation of clean values -------
# The team's preprocessing has ALREADY logged, scaled and clipped these values
# using all subjects. Reusing this file aligns Q4 with Q3 and the class model;
# it does not make preprocessing fully independent of the evaluation subjects.
prepare_data <- function(training_data, evaluation_data) {
  list(x_train = as.matrix(training_data[, protein_names]),
       x_evaluation = as.matrix(evaluation_data[, protein_names]),
       y_train = training_data$group == "ASD",
       y_evaluation = evaluation_data$group == "ASD",
       forest_group = factor(training_data$group, levels = c("ASD", "TD")))
}

## 3. Protein selection -- repeat inside each training-validation fold -------
rank_proteins <- function(data) {
  x <- data$x_train
  p_values <- apply(x, 2, function(z) {
    t.test(z[data$y_train], z[!data$y_train], var.equal = FALSE)$p.value
  })
  t_table <- tibble(protein = names(p_values), p_value = p_values) %>%
    arrange(p_value) %>%
    mutate(m = n(), hm = log(m) + 1 / (2 * m) - digamma(1),
           p.adj = m * hm * p_value / row_number())
  # Retain the group's adjusted-score formula for benchmark consistency.
  t_order <- t_table %>% arrange(p.adj, protein) %>% pull(protein)

  set.seed(seed)
  forest <- randomForest(x = x, y = data$forest_group,
                         ntree = rf_trees, importance = TRUE)
  rf_order <- names(sort(forest$importance[, "MeanDecreaseGini"], decreasing = TRUE))

  # Same fuzzy membership idea as Q3c, using the raw p-value ordering.
  fuzzy_table <- tibble(protein = t_table$protein, rank_t = seq_len(nrow(t_table))) %>%
    mutate(rank_rf = match(protein, rf_order),
           mem_t = 1 - (rank_t - 1) / (n() - 1),
           mem_rf = 1 - (rank_rf - 1) / (n() - 1),
           membership = pmin(mem_t, mem_rf)) %>%
    arrange(desc(membership), rank_t + rank_rf, protein)

  list(ttest = t_order, rf = rf_order, hard10 = intersect(head(t_order, 10), head(rf_order, 10)),
       hard25 = intersect(head(t_order, 25), head(rf_order, 25)),
       fuzzy = fuzzy_table$protein, fuzzy_table = fuzzy_table)
}

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)
score_predictions <- function(truth, probabilities) {
  tibble(truth = factor(truth, levels = c(FALSE, TRUE)),
         estimate = factor(probabilities > 0.5, levels = c(FALSE, TRUE)),
         probability = probabilities) %>%
    class_metrics(truth = truth, estimate = estimate, probability,
                  event_level = "second") %>%
    select(.metric, .estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)
}

## 4. Fit models: original logistic regression or sparse lasso ---------------

fit_configuration <- function(config, data, rankings) {
  panel <- switch(config$pool,
    hard10 = rankings$hard10, hard25 = rankings$hard25,
    fuzzy = head(rankings$fuzzy, config$k),
    ttest = head(rankings$ttest, config$k), rf = head(rankings$rf, config$k),
    fuzzy25 = head(rankings$fuzzy, 25), all = colnames(data$x_train))

  if (config$method == "glm") {
    model_data <- as.data.frame(data$x_train[, panel, drop = FALSE])
    model_data$class <- data$y_train
    fit <- glm(class ~ ., data = model_data, family = binomial())
    if (!fit$converged) stop("Logistic regression failed to converge: ", config$model)
    # Include a row-bearing data frame even for an empty, intercept-only panel.
    new_data <- as.data.frame(data$x_evaluation)
    probabilities <- as.numeric(predict(fit, newdata = new_data, type = "response"))
    coefficients <- tibble(protein = names(coef(fit)), coefficient = unname(coef(fit)))
    lambda <- NA_real_
  } else {
    x <- data$x_train[, panel, drop = FALSE]
    # One path per fit; specify relative penalty so folds use the same rule.
    fit <- glmnet(x, data$y_train, family = "binomial", alpha = 1,
                  standardize = FALSE, nlambda = 100, lambda.min.ratio = 0.001)
    lambda <- max(fit$lambda) * config$lambda_fraction
    probabilities <- as.numeric(predict(
      fit, newx = data$x_evaluation[, panel, drop = FALSE],
      s = lambda, type = "response"
    ))
    coefficients_matrix <- as.matrix(coef(fit, s = lambda))
    coefficients <- tibble(protein = rownames(coefficients_matrix),
                           coefficient = as.numeric(coefficients_matrix)) %>%
      filter(protein == "(Intercept)" | abs(coefficient) > coefficient_tolerance)
    panel <- setdiff(coefficients$protein, "(Intercept)")
  }
  list(metrics = score_predictions(data$y_evaluation, probabilities) %>%
         mutate(model = config$model, n_proteins = length(panel)),
       proteins = panel, coefficients = coefficients,
       probabilities = probabilities, lambda = lambda, fit = fit, config = config)
}

## 5. Repeated cross-validation: compare approaches without using test scores -

# Run this whole section, including the loop, before running the summaries.
# In RStudio select the entire for (...) { ... } block and press Cmd+Return.
set.seed(seed)
folds <- vfold_cv(train_data, v = 5, repeats = cv_repeats, strata = group)
cv_results <- vector("list", nrow(folds))
cv_panels <- vector("list", nrow(folds))
write_csv(configurations, file.path(output_dir, "search-configurations.csv"))

for (i in seq_len(nrow(folds))) {
  cat("Cross-validation fold", i, "of", nrow(folds), "\n")
  fold_data <- prepare_data(analysis(folds$splits[[i]]), assessment(folds$splits[[i]]))
  fold_rankings <- rank_proteins(fold_data)
  fitted <- map(seq_len(nrow(configurations)), function(j) {
    fit_configuration(configurations[j, ], fold_data, fold_rankings)
  })
  cv_results[[i]] <- map_dfr(fitted, "metrics") %>%
    mutate(repeat_id = folds$id[i], fold_id = folds$id2[i])
  cv_panels[[i]] <- map2_dfr(fitted, configurations$model, function(result, model) {
    tibble(model = model, protein = result$proteins)
  }) %>% mutate(repeat_id = folds$id[i], fold_id = folds$id2[i])
}
cv_results <- bind_rows(cv_results)
cv_panels <- bind_rows(cv_panels)

# Average within each repeat first. Reused subjects mean fold scores are not
# independent; variability below is descriptive, not a confidence interval.
cv_by_repeat <- cv_results %>%
  group_by(model, repeat_id) %>%
  summarise(across(c(accuracy, roc_auc, sensitivity, specificity, n_proteins), mean),
            .groups = "drop")
cv_summary <- cv_by_repeat %>% group_by(model) %>%
  summarise(mean_accuracy = mean(accuracy), sd_repeat_accuracy = sd(accuracy),
            mean_auc = mean(roc_auc), mean_sensitivity = mean(sensitivity),
            mean_specificity = mean(specificity), mean_panel_size = mean(n_proteins),
            .groups = "drop") %>%
  arrange(desc(mean_accuracy), mean_panel_size, desc(mean_auc))
print(cv_summary, n = Inf)

# Lock TWO approaches using training CV only. Do not inspect the final test
# table until these choices are made, or change the search to chase test scores.
# Primary: highest validation accuracy, then smaller mean panel, then higher AUC.
chosen_model <- cv_summary %>% filter(model != "hard10", mean_panel_size >= 1) %>%
  arrange(desc(mean_accuracy), mean_panel_size, desc(mean_auc)) %>%
  slice(1) %>% pull(model)
# Simpler alternative: limit the design to fixed panels of 1-3 proteins.
# Prefer the smallest size within 2 percentage points of the best SMALL model.
small_models <- configurations %>%
  filter(method == "glm", pool %in% c("fuzzy", "ttest", "rf"), k <= 3) %>%
  pull(model)
small_summary <- cv_summary %>% filter(model %in% small_models)
simple_model <- small_summary %>%
  filter(mean_accuracy >= max(mean_accuracy) - accuracy_margin) %>%
  arrange(mean_panel_size, desc(mean_accuracy), desc(mean_auc)) %>%
  slice(1) %>% pull(model)
chosen_config <- configurations %>% filter(model == chosen_model)
simple_config <- configurations %>% filter(model == simple_model)
selection_decisions <- tibble(role = c("primary", "simpler"),
                              model = c(chosen_model, simple_model))
write_csv(selection_decisions, file.path(output_dir, "cv-locked-choices.csv"))
cat("CV-selected approach:", chosen_model, "\n")
cat("CV-selected simpler approach:", simple_model, "\n")

## 6. Fit Q4 choices using training participants ----------------------------
final_data <- prepare_data(train_data, test_data)
final_rankings <- rank_proteins(final_data)
alternative_fit <- fit_configuration(chosen_config, final_data, final_rankings)
simple_fit <- fit_configuration(simple_config, final_data, final_rankings)
alternative_panel <- alternative_fit$proteins
simple_panel <- simple_fit$proteins

## 7. ONE class benchmark -- copy the original procedure into our own file ---
# Full-data t-tests/RF selection is intentionally retained ONLY for the class
# benchmark. Q4 choices above were locked using training cross-validation.
# The supplied clean values are used directly, not reconstructed from raw CSV.
# This reproduces the team's q3_base.R without sourcing or changing that file.
all_data <- prepare_data(analysis_data, analysis_data)
all_rankings <- rank_proteins(all_data)
inclass_panel <- all_rankings$hard10
baseline_fit <- fit_configuration(configurations %>% filter(model == "hard10"),
                                  final_data, list(hard10 = inclass_panel))
final_results <- bind_rows(
  baseline_fit$metrics %>% mutate(model = "In-class baseline (rerun)"),
  alternative_fit$metrics %>% mutate(model = paste0("Q4 alternative: ", chosen_model)),
  simple_fit$metrics %>% mutate(model = paste0("Q4 simpler: ", simple_model))
) %>% select(model, n_proteins, sensitivity, specificity, accuracy, roc_auc)
# Historical comments are kept OUT of the current-run comparison table.
historical_baseline <- tibble(source = "q3_base.R historical comments; not reproduced",
                              sensitivity = 0.875, specificity = 0.800,
                              accuracy = 0.839, roc_auc = 0.908)
benchmark_comparison <- final_results
print(final_results, width = Inf)
print(list(inclass_panel = inclass_panel, alternative_panel = alternative_panel,
           simple_panel = simple_panel))

## 8. Save results and verify the analysis -----------------------------------
# These checks validate the implementation, not the scientific hypothesis.
stopifnot(nrow(cv_results) == nrow(folds) * nrow(configurations),
          !anyNA(cv_results), length(intersect(train_data$subject_id,
                                              test_data$subject_id)) == 0,
          length(simple_panel) <= 3, !"ados" %in% protein_names,
          !"subject_id" %in% protein_names)
all_fits <- list(inclass = baseline_fit, primary = alternative_fit, simpler = simple_fit)
confusion_counts <- imap_dfr(all_fits, function(result, role) {
  truth <- test_data$group == "ASD"
  predicted <- result$probabilities > 0.5
  stopifnot(length(predicted) == nrow(test_data),
            all(is.finite(result$probabilities)),
            all(result$probabilities >= 0 & result$probabilities <= 1))
  tibble(role = role, model = if (role == "inclass") "inclass" else result$config$model,
         true_ASD = sum(truth & predicted), missed_ASD = sum(truth & !predicted),
         true_TD = sum(!truth & !predicted), false_ASD = sum(!truth & predicted),
         correct = sum(truth == predicted), total = length(truth))
})
protein_stability <- cv_panels %>%
  filter(model %in% c(chosen_model, simple_model, "hard10")) %>%
  count(model, protein, name = "folds_selected") %>%
  mutate(fraction_selected = folds_selected / nrow(folds)) %>%
  arrange(model, desc(folds_selected), protein)
final_panels <- imap_dfr(all_fits, function(result, role) {
  tibble(role = role, protein = result$proteins)
})
write_csv(cv_results, file.path(output_dir, "cv-fold-results.csv"))
write_csv(cv_summary, file.path(output_dir, "cv-summary.csv"))
write_csv(cv_panels, file.path(output_dir, "cv-selected-proteins.csv"))
write_csv(protein_stability, file.path(output_dir, "protein-selection-frequency.csv"))
write_csv(final_results, file.path(output_dir, "exploratory-test-results.csv"))
write_csv(benchmark_comparison, file.path(output_dir, "benchmark-comparison.csv"))
write_csv(confusion_counts, file.path(output_dir, "confusion-counts.csv"))
write_csv(final_panels, file.path(output_dir, "final-panels.csv"))
write_csv(alternative_fit$coefficients, file.path(output_dir, "selected-coefficients.csv"))
write_csv(simple_fit$coefficients, file.path(output_dir, "simpler-coefficients.csv"))
write_csv(final_rankings$fuzzy_table, file.path(output_dir, "training-fuzzy-ranks.csv"))
write_csv(historical_baseline, file.path(output_dir, "historical-baseline.csv"))
write_csv(tibble(subject_id = analysis_data$subject_id,
                 partition = if_else(analysis_data$subject_id %in% train_data$subject_id,
                                     "training", "testing")),
          file.path(output_dir, "subject-partitions.csv"))
write_csv(tibble(subject_id = test_data$subject_id, group = test_data$group,
                 inclass_probability = baseline_fit$probabilities,
                 alternative_probability = alternative_fit$probabilities,
                 simpler_probability = simple_fit$probabilities),
          file.path(output_dir, "exploratory-test-predictions.csv"))
write_csv(tibble(subject_id = test_data$subject_id, group = test_data$group,
                 probability_ASD = baseline_fit$probabilities,
                 predicted_group = if_else(baseline_fit$probabilities > 0.5,
                                           "ASD", "TD")),
          file.path(output_dir, "baseline-test-predictions.csv"))

cv_plot <- ggplot(cv_summary,
                 aes(x = mean_panel_size, y = mean_accuracy, label = model)) +
  geom_point(color = "#286696", size = 2.5, alpha = 0.75) +
  geom_text(data = cv_summary %>%
              filter(model %in% c("hard10", chosen_model, simple_model)),
            hjust = 0, nudge_x = 0.7, size = 3.2, check_overlap = TRUE) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = "Mean number of proteins selected across validation folds",
       y = "Mean training cross-validation accuracy",
       title = "Q4: predictive performance and panel size",
       subtitle = "Training-validation scores used for selection; not independent test estimates") +
  theme_minimal(base_size = 11)
ggsave(file.path(output_dir, "cv-comparison.png"), cv_plot, width = 9, height = 5)

# Save fitted models and inputs. Values are already processed: do not log or
# standardize them again when inspecting or reproducing these saved predictions.
save(analysis_data, train_data, test_data, protein_names, configurations, folds,
     seed, accuracy_margin, coefficient_tolerance,
     cv_results, cv_summary, cv_panels, protein_stability, chosen_model,
     chosen_config, final_rankings, alternative_panel,
     simple_model, simple_config, simple_panel, simple_fit,
     final_results, alternative_fit, inclass_panel, baseline_fit,
     historical_baseline, benchmark_comparison, confusion_counts, final_panels,
     final_data, cv_plot, file = file.path(output_dir, "sou-q4-results.RData"))
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))
write_csv(tibble(file = "data/biomarker-clean.RData", md5 = unname(tools::md5sum(clean_path))),
          file.path(output_dir, "input-fingerprint.csv"))

## 9. Generate the explanation and report -----------------------------------
# Generate the report using the Q4 helper in drafts/.
source(file.path(project_root, "scripts", "drafts", "q4_writeup.R"), local = TRUE)
cat("\nFinished. Results and report: ", output_dir, "\n", sep = "")
cat("Inspect cv_summary, final_results, inclass_panel, simple_panel, and alternative_panel.\n")
