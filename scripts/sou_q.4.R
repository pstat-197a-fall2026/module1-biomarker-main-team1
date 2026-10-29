# Question 4 -- Sou
# Question: can a smaller or alternative protein panel classify ASD/TD well?
#
# Q3a motivates training-only selection AND preprocessing.
# Q3b motivates testing a larger candidate pool, without assuming it helps.
# Q3c motivates fuzzy ranking instead of insisting on a top-10 intersection.
# New method: lasso logistic regression shrinks some coefficients to exactly
# zero, simultaneously selecting proteins and fitting the classifier.
#
# FIRST RUN IN RSTUDIO:
# 1. Open module1-biomarker.Rproj.
# 2. Run this in the Console: source("scripts/sou_q.4.R")
#    Alternatively, open this file and click Source (not Run on a single line).
# 3. Wait for cross-validation to finish, then inspect cv_summary/final_results.
#
# To work section by section, run the libraries and sections 1 through 9 IN ORDER.
# Later sections depend on earlier ones. In particular, cv_results is created
# by the entire loop in section 5; do not jump straight to cv_by_repeat.
# Other Q1-Q3 scripts do not need to be run first.
# Loading saved results lets you inspect tables, but does not load the libraries
# or define this script's helper functions. Run from the top to edit/recompute.
# The old starter in scripts/drafts/ is historical only; do not run it first.
# Current outputs go to results/sou-q4-final/; older outputs are not overwritten.
# IMPORTANT: this test set was viewed during the starter analysis. Subsequent
# test results are exploratory, even though choices below use training CV only.

library(tidyverse)
library(randomForest)
library(rsample)
library(yardstick)
library(glmnet)

## 1. Settings and raw data -------------------------------------------------

project_root <- if (file.exists("data/biomarker-raw.csv")) "." else ".."
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

raw_path <- file.path(project_root, "data", "biomarker-raw.csv")
headers <- read_csv(raw_path, col_names = FALSE, n_max = 2,
                    col_select = -(1:2), show_col_types = FALSE)
var_names <- tibble(
  name = as.character(unlist(headers[1, ], use.names = FALSE)),
  abbreviation = as.character(unlist(headers[2, ], use.names = FALSE))
) %>% drop_na()

raw_data <- read_csv(
  raw_path, skip = 2, col_select = -2L,
  col_names = c("group", "empty", var_names$abbreviation, "ados"),
  na = c("-", ""), show_col_types = FALSE
) %>%
  filter(!is.na(group)) %>% select(-ados) %>%
  mutate(subject_id = row_number(), .before = group)

protein_names <- setdiff(names(raw_data), c("subject_id", "group"))
stopifnot(all(as.matrix(raw_data[, protein_names]) > 0),
          !anyNA(raw_data[, protein_names]))

# Preserve Q3a's outer split rather than searching for a favorable partition.
set.seed(seed)
outer_split <- initial_split(raw_data, prop = 0.8)
train_data <- training(outer_split)
test_data <- testing(outer_split)
cat("Subjects:", nrow(train_data), "training;", nrow(test_data), "testing\n")

## 2. Preprocessing -- learn parameters using training subjects only ---------

prepare_data <- function(training_data, evaluation_data) {
  log_train <- log10(as.matrix(training_data[, protein_names]))
  log_evaluation <- log10(as.matrix(evaluation_data[, protein_names]))
  centers <- colMeans(log_train)
  scales <- apply(log_train, 2, sd)
  keep <- is.finite(scales) & scales > 0

  transform_values <- function(x) {
    z <- sweep(sweep(x[, keep, drop = FALSE], 2, centers[keep], "-"),
               2, scales[keep], "/")
    z[z > 3] <- 3
    z[z < -3] <- -3
    z
  }
  list(x_train = transform_values(log_train),
       x_evaluation = transform_values(log_evaluation),
       y_train = training_data$group == "ASD",
       y_evaluation = evaluation_data$group == "ASD")
}

## 3. Protein selection -- repeat inside every cross-validation fold ---------

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
  forest <- randomForest(x = x, y = factor(data$y_train),
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

## 6. Fit the locked choice and reference using all training subjects --------

final_data <- prepare_data(train_data, test_data)
final_rankings <- rank_proteins(final_data)
reference_fit <- fit_configuration(configurations %>% filter(model == "hard10"),
                                    final_data, final_rankings)
sou_fit <- fit_configuration(chosen_config, final_data, final_rankings)
simple_fit <- fit_configuration(simple_config, final_data, final_rankings)
reference_panel <- reference_fit$proteins
sou_panel <- sou_fit$proteins
simple_panel <- simple_fit$proteins
final_results <- bind_rows(reference_fit$metrics, sou_fit$metrics, simple_fit$metrics) %>%
  distinct(model, .keep_all = TRUE) %>%
  select(model, n_proteins, everything())
print(list(reference_panel = reference_panel, sou_panel = sou_panel))
print(final_results)

## 7. Reproduce the in-class benchmark without sourcing or editing team files --
# This copies the baseline procedure into OUR script: global preprocessing and
# full-data selection, followed by logistic fitting on the same 80% subjects.
# Its leakage is retained ONLY to benchmark the supplied method as requested.
all_log <- log10(as.matrix(raw_data[, protein_names]))
all_z <- scale(all_log)
all_z[all_z > 3] <- 3
all_z[all_z < -3] <- -3
baseline_p <- apply(all_z, 2, function(z) {
  t.test(z[raw_data$group == "ASD"], z[raw_data$group == "TD"],
         var.equal = FALSE)$p.value
})
baseline_t <- tibble(protein = names(baseline_p), p_value = baseline_p) %>%
  arrange(p_value) %>%
  mutate(m = n(), hm = log(m) + 1 / (2 * m) - digamma(1),
         p.adj = m * hm * p_value / row_number()) %>%
  slice_min(p.adj, n = 10) %>% pull(protein)
set.seed(seed)
baseline_forest <- randomForest(x = all_z, y = factor(raw_data$group),
                                 ntree = rf_trees, importance = TRUE)
baseline_rf <- tibble(protein = rownames(baseline_forest$importance),
                      gini = baseline_forest$importance[, "MeanDecreaseGini"]) %>%
  slice_max(gini, n = 10) %>% pull(protein)
inclass_panel <- intersect(baseline_t, baseline_rf)
baseline_data <- list(x_train = all_z[train_data$subject_id, , drop = FALSE],
                      x_evaluation = all_z[test_data$subject_id, , drop = FALSE],
                      y_train = train_data$group == "ASD",
                      y_evaluation = test_data$group == "ASD")
baseline_fit <- fit_configuration(configurations %>% filter(model == "hard10"),
                                  baseline_data, list(hard10 = inclass_panel))
inclass_benchmarks <- bind_rows(
  baseline_fit$metrics %>% mutate(model = "In-class rerun (full-data selection)"),
  tibble(model = "In-class historical comments (not reproduced)",
         n_proteins = NA_integer_, sensitivity = 0.875, specificity = 0.800,
         accuracy = 0.839, roc_auc = 0.908)
) %>% select(model, n_proteins, sensitivity, specificity, accuracy, roc_auc)
benchmark_comparison <- bind_rows(inclass_benchmarks, final_results)
print(benchmark_comparison)

## 8. Save results and verify the analysis -----------------------------------
# These checks validate the implementation, not the scientific hypothesis.
stopifnot(nrow(cv_results) == nrow(folds) * nrow(configurations),
          !anyNA(cv_results), length(intersect(train_data$subject_id,
                                              test_data$subject_id)) == 0,
          length(simple_panel) <= 3, !"ados" %in% protein_names,
          !"subject_id" %in% protein_names)
all_fits <- list(reference = reference_fit, primary = sou_fit,
                 simpler = simple_fit, inclass = baseline_fit)
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
write_csv(sou_fit$coefficients, file.path(output_dir, "selected-coefficients.csv"))
write_csv(simple_fit$coefficients, file.path(output_dir, "simpler-coefficients.csv"))
write_csv(final_rankings$fuzzy_table, file.path(output_dir, "training-fuzzy-ranks.csv"))
write_csv(tibble(subject_id = raw_data$subject_id,
                 partition = if_else(raw_data$subject_id %in% train_data$subject_id,
                                     "training", "testing")),
          file.path(output_dir, "subject-partitions.csv"))
write_csv(tibble(subject_id = test_data$subject_id, group = test_data$group,
                 inclass_probability = baseline_fit$probabilities,
                 reference_probability = reference_fit$probabilities,
                 sou_probability = sou_fit$probabilities,
                 simpler_probability = simple_fit$probabilities),
          file.path(output_dir, "exploratory-test-predictions.csv"))

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

# Save fitted models as well as their outputs. Each preprocessing step can be
# reproduced from train_data and prepare_data(); do not apply global scaling.
save(raw_data, train_data, test_data, protein_names, configurations, folds,
     seed, accuracy_margin, coefficient_tolerance,
     cv_results, cv_summary, cv_panels, protein_stability, chosen_model,
     chosen_config, final_rankings, reference_panel, sou_panel,
     simple_model, simple_config, simple_panel, simple_fit,
     final_results, sou_fit, reference_fit, inclass_panel, baseline_fit,
     inclass_benchmarks, benchmark_comparison, confusion_counts, final_panels,
     final_data, cv_plot, file = file.path(output_dir, "sou-q4-results.RData"))
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))
write_csv(tibble(file = raw_path, md5 = unname(tools::md5sum(raw_path))),
          file.path(output_dir, "input-fingerprint.csv"))

## 9. Generate the explanation and report -----------------------------------
# This helper is Sou's own new file, not a teammate's script.
source(file.path(project_root, "scripts", "sou_q4_writeup.R"), local = TRUE)
cat("\nFinished. Results and report: ", output_dir, "\n", sep = "")
cat("Inspect cv_summary, benchmark_comparison, simple_panel, and sou_panel.\n")
