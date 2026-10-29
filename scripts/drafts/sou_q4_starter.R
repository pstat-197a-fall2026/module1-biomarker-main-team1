# Question 4 -- Sou's starter analysis
# Goal: find a smaller protein panel with comparable classification accuracy.
# This is a starting approach, not a claim that a smaller panel will succeed.
#
# Insights from Q3:
# 1. Split before selecting proteins; keep the test set for final evaluation.
# 2. More proteins did not improve threshold-based accuracy in Q3b, so try fewer.
# 3. Use Q3c's fuzzy ranking to find proteins supported by BOTH methods without
#    requiring them to appear in both top-10 lists.
#
# Run from the project root or the scripts directory.
# Required packages: tidyverse, randomForest, rsample, yardstick.

library(tidyverse)
library(randomForest)
library(rsample)
library(yardstick)

## 1. Load data and define choices BEFORE looking at test results ------------

data_path <- if (file.exists("data/biomarker-clean.RData")) {
  "data/biomarker-clean.RData"
} else {
  "../data/biomarker-clean.RData"
}
load(data_path)

# Starting simplification: reuse the group's processed data for consistency.
# LIMITATION: its scaling/trimming used all subjects. Feature selection below
# is training-only, but preprocessing is not fully isolated. A further extension
# is to start from raw data and estimate scaling inside each training fold.
biomarker_data <- biomarker_clean %>% select(-ados)

seed <- 101422
panel_sizes <- c(1, 2, 3, 5) # TODO: decide candidate sizes before final testing.
accuracy_margin <- 0.02     # Working definition: within 2 percentage points
                           # of the best mean CV accuracy; not an equivalence test.
rf_trees <- 1000

# Same outer split seed/proportion as Q3a; do not change it to seek better scores.
set.seed(seed)
outer_split <- initial_split(biomarker_data, prop = 0.8)
train_data <- training(outer_split)
test_data <- testing(outer_split)

## 2. Rank proteins using only the data passed to this function --------------

rank_proteins <- function(training_data) {
  predictors <- training_data %>% select(-group)

  # Same Welch t-tests and adjusted-score formula used in the group scripts.
  t_ranking <- map_dfr(names(predictors), function(protein) {
    tibble(
      protein = protein,
      p_value = t.test(
        training_data[[protein]] ~ training_data$group,
        var.equal = FALSE
      )$p.value
    )
  }) %>%
    arrange(p_value) %>%
    mutate(
      m = n(),
      hm = log(m) + 1 / (2 * m) - digamma(1),
      p.adj = m * hm * p_value / row_number()
    )

  set.seed(seed)
  rf_fit <- randomForest(
    x = predictors, y = factor(training_data$group),
    ntree = rf_trees, importance = TRUE
  )
  rf_ranking <- tibble(
    protein = rownames(rf_fit$importance),
    gini = rf_fit$importance[, "MeanDecreaseGini"]
  ) %>% arrange(desc(gini))

  # Reference: original top-10 hard intersection, now selected on training data.
  hard_panel <- intersect(
    t_ranking %>% slice_min(p.adj, n = 10) %>% pull(protein),
    rf_ranking %>% slice_max(gini, n = 10) %>% pull(protein)
  )

  # Follow Q3c's membership construction, but take the best k proteins instead
  # of applying its 0.95 cutoff. This lets us directly control panel size.
  fuzzy_ranking <- t_ranking %>%
    select(protein) %>%
    mutate(rank_t = row_number()) %>%
    inner_join(
      rf_ranking %>% mutate(rank_rf = row_number()) %>%
        select(protein, rank_rf),
      by = "protein"
    ) %>%
    mutate(
      mem_t = 1 - (rank_t - 1) / (n() - 1),
      mem_rf = 1 - (rank_rf - 1) / (n() - 1),
      mem_fuzzy = pmin(mem_t, mem_rf)
    ) %>%
    arrange(desc(mem_fuzzy), rank_t + rank_rf, protein)

  list(hard = hard_panel, fuzzy = fuzzy_ranking$protein)
}

## 3. Helper: fit logistic regression and evaluate one panel -----------------

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

evaluate_panel <- function(training_data, evaluation_data, proteins) {
  model_data <- training_data %>%
    select(all_of(proteins)) %>%
    mutate(class = training_data$group == "ASD")
  # An empty hard intersection gives an intercept-only reference model.
  fit <- glm(class ~ ., data = model_data, family = binomial())
  probabilities <- predict(fit, newdata = evaluation_data, type = "response")

  tibble(
    truth = factor(evaluation_data$group == "ASD", levels = c(FALSE, TRUE)),
    estimate = factor(probabilities > 0.5, levels = c(FALSE, TRUE)),
    probability = probabilities
  ) %>%
    class_metrics(
      truth = truth, estimate = estimate, probability,
      event_level = "second"
    ) %>%
    select(.metric, .estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate) %>%
    mutate(n_proteins = length(proteins))
}

## 4. Choose panel size using cross-validation within TRAINING data only -----

set.seed(seed)
folds <- vfold_cv(train_data, v = 5, strata = group)

cv_results <- map2_dfr(folds$splits, folds$id, function(split, fold_id) {
  fold_train <- analysis(split)
  fold_validation <- assessment(split)
  # Repeat selection inside EVERY fold, rather than ranking once before CV.
  rankings <- rank_proteins(fold_train)

  map_dfr(panel_sizes, function(k) {
    evaluate_panel(fold_train, fold_validation, head(rankings$fuzzy, k)) %>%
      mutate(k = k, fold = fold_id)
  })
})

cv_summary <- cv_results %>%
  group_by(k) %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy = sd(accuracy),
    mean_auc = mean(roc_auc),
    .groups = "drop"
  ) %>% arrange(k)
print(cv_summary)

# Prefer the smallest panel within the prespecified margin of the best score.
chosen_k <- cv_summary %>%
  filter(mean_accuracy >= max(mean_accuracy) - accuracy_margin) %>%
  slice_min(k, n = 1) %>% pull(k)

## 5. Final comparison -- evaluate only AFTER fixing the approach above ------

# Fit the selection procedure again using all outer-training subjects.
final_rankings <- rank_proteins(train_data)
reference_panel <- final_rankings$hard
sou_panel <- head(final_rankings$fuzzy, chosen_k)

print(list(reference_panel = reference_panel, sou_panel = sou_panel))

final_results <- bind_rows(
  evaluate_panel(train_data, test_data, reference_panel) %>%
    mutate(model = "Training-only hard intersection"),
  evaluate_panel(train_data, test_data, sou_panel) %>%
    mutate(model = paste0("Sou fuzzy-ranked top ", chosen_k))
) %>% select(model, n_proteins, everything())
print(final_results)

# Historical benchmark: copied from q3_base.R, NOT rerun by this script.
# Its feature selection used all subjects, so it is an optimistic comparison.
inclass_recorded <- tibble(
  model = "Recorded in-class baseline (full-data selection)",
  sensitivity = 0.875, specificity = 0.800,
  accuracy = 0.839, roc_auc = 0.908
)
print(inclass_recorded)

## 6. TODO: interpret results in your own words ------------------------------
# - Which proteins did the reference and proposed panels contain?
# - Is the proposed panel actually smaller? Fewer than 10 does not necessarily
#   mean fewer than the original INTERSECTION; compare n_proteins explicitly.
# - How do test accuracy, sensitivity, specificity, and AUROC compare?
# - Does this support a simpler panel with comparable accuracy? If not, say so.
# - Explain how each Q3 insight motivated this approach.
# - Note the small test set, preprocessing limitation, and any glm warnings
#   (e.g. separation/convergence). Do not suppress or ignore fitting problems.
# - After inspecting final test results, do not retune against that same test
#   set and describe the resulting score as an untouched test evaluation.
