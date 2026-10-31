# Run from the project root with:
# source("scripts/analysis-main.R", chdir = TRUE)
# The existing Q1-Q3 code expects scripts/ as its working directory.

### QUESTION 1

library(tidyverse)
library(ggplot2)
library(tidyr)


# get names
var_names <- read_csv('../data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean <- read_csv('../data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group))

# Sample 10 proteins from DataFrame
set.seed(980)
protein_sample <- sample(names(biomarker_clean), 10)

# Plot density distributions
biomarker_clean |> 
  select(all_of(protein_sample)) |> 
  pivot_longer(everything(), names_to = "variable", values_to = "value") |> 
  ggplot(aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of Randomly Sampled Proteins",
       x = "Value", y = "Density")

# Findings
# Many of the proteins in the sample are heavily right skewed, which implies 
# most of the participants in the study have protein
# counts for these proteins in the very low end of the study's range of protein
# counts. A log transformation helps normalize these protein values so that
# the data looks more bell-shaped. Normalizing the data allows us to apply
# statistical methods which rely on the data being roughly normal.

# QUESTION 2

library(tidyverse)
library(kableExtra)

# get names
var_names <- read_csv('../data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# read in data
biomarker_clean <- read_csv('../data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # Create subject_id column to keep track of each subject
  mutate(subject_id = row_number()) %>%
  # Log transform, center, and scale, but do not trim
  mutate(across(
    -c(subject_id, group, ados),
    ~ scale(log10(.x))[, 1]
  )) |> 
  select(subject_id, group, ados, everything())

# Pivot proteins into long format
protein_long <- biomarker_clean |> 
  pivot_longer(
    cols = -c(subject_id, group, ados),
    names_to = "protein",
    values_to = "z_score"
  ) |> 
  # New column to determine if value is outlier
  mutate(is_outlier = abs(z_score) > 3)

# Count how many outliers occur for each observation, proportions
subject_outliers <- protein_long |> 
  group_by(subject_id, group) |> 
  summarise(
    n_measured = sum(!is.na(z_score)),
    n_outliers = sum(is_outlier, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    prop_outliers = if_else(
      n_measured > 0,
      n_outliers / n_measured,
      NA_real_
    )
  ) |> 
  arrange(desc(n_outliers))

# Subjects with the most outlying protein measurements
print(subject_outliers, n = 30)

# Looking at the amount of outliers for each subject, we have one person who
# has outlier values in over 10% of the proteins, with 157 out of 1317 measured
# proteins having values 3 standard deviations away from the mean.
# There are 7 total subjects with five percent or more of the proteins being
# outliers, with 6 of those subjects having over 100 protein outlier values.
# After that, the pure count drops sharply to below 50 outlier proteins per
# subject, and besides the top 23, the rest only have less than 19 proteins
# for which they have outlier values.

# See the counts of ASD vs TD grouped by number of outlier proteins
subject_outliers |> 
  count(group, n_outliers, name = "n_subjects") |> 
  pivot_wider(
    names_from = group,
    values_from = n_subjects,
    values_fill = 0
  ) |> 
  arrange(n_outliers) |> 
  kable()

# Based on initial checks, there does not seem to be a significant pattern
# to separate ASD vs TD for each individual count. In general there's
# Very few subjects that have more than 20 outlier values.

# Short summary of the groups' outlier counts
subject_outliers |> 
  group_by(group) |> 
  summarise(
    n_subjects = n(),
    mean_outliers = mean(n_outliers),
    median_outliers = median(n_outliers),
    max_outliers = max(n_outliers),
    pct_subjects_with_outliers = 100 * mean(n_outliers > 0),
    pct_measurements_outlying =
      100 * sum(n_outliers) / sum(n_measured),
    .groups = "drop"
  ) |> 
  kable()

# In this summary aggregating by groups, it looks as if there's actually less
# outliers in the ASD group. Intuitively this does not make sense as we would
# not expect outlier values to be associated with typically developing,
# so this might be due to small sample size / outlier values skewing the
# ASD group.

# Visualize particular subjects who stand out relative to the others
ggplot(
  subject_outliers,
  aes(x = reorder(factor(subject_id), n_outliers),
      y = n_outliers,
      fill = group)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Subject",
    y = "Number of outlying protein measurements",
    fill = "Group",
    title = "Outlying protein measurements per subject"
  ) +
  theme_minimal()

# As we can see, subject 154 has by far the most outlier protein measurements.
# Subjects 108, 9, 121, 52, and 77 have over 100 protein measurements with
# outlier values, separating themselves from the rest of the group having
# 77 or less outlier values. It appears that the typically developing group
# may have more outliers than the ASD group, though this may be due to small
# sample sizes.


# QUESTION 3
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

load('../data/biomarker-clean.RData')

# NEED FOR RENDERING
# Define classification evaluation metrics
class_metrics <- yardstick::metric_set(
  yardstick::sensitivity,
  yardstick::specificity,
  yardstick::accuracy,
  yardstick::roc_auc
)

# base analysis
## MULTIPLE TESTING 
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set (WHAT WAS CHANGED)

# eval metrics
class_metrics <- yardstick::metric_set(
  yardstick::sensitivity,
  yardstick::specificity,
  yardstick::accuracy,
  yardstick::roc_auc
)

# prepare predictions
test_predictions <- testing(biomarker_split) %>%
  mutate(
    pred = predict(
      fit,
      newdata = testing(biomarker_split),
      type = "response"
    ),
    truth = factor(class, levels = c(FALSE, TRUE)),
    estimate = factor(pred > 0.5, levels = c(FALSE, TRUE))
  )

# evaluate model
base <- class_metrics(
  test_predictions,
  truth = truth,
  estimate = estimate,
  pred,
  event_level = "second"
)

# MODIFICATION 1

# Modification 1: Training/Testing Split from the beginning

# using same seed
set.seed(101422)
# SPLIT DATA BEFORE!!!
biomarker_split <- initial_split(biomarker_clean, prop = 0.8)
biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)

## MULTIPLE TESTING 
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_train %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_train %>%
  select(-c(group, ados))

response <- biomarker_train %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar_train <- biomarker_train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_test <- biomarker_test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)


# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker_sstar_train, 
           family = 'binomial')

# Generate predictions on testing data
biomarker_sstar_test <- biomarker_sstar_test %>%
  mutate(pred = predict(fit, newdata = biomarker_sstar_test,
                        type = "response"))

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

# Evaluate predictions
mod1 <- biomarker_sstar_test %>%
  mutate(
    truth = factor(class, levels = c(FALSE, TRUE)),
    estimate = factor(pred > 0.5, levels = c(FALSE, TRUE))
  ) %>%
  class_metrics(
    truth = truth,
    estimate = estimate,
    pred,
    event_level = "second"
  )

# Modification 2:

# Modification 2: choose a larger number (more than ten) of top predictive proteins using each selection method
# only change from the base are choosing 25 each from the t test selection and the rf selection
# so only changed 2 numbers

## MULTIPLE TESTING 
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 25) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 25) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set (WHAT WAS CHANGED)

# eval metrics
class_metrics <- yardstick::metric_set(
  yardstick::sensitivity,
  yardstick::specificity,
  yardstick::accuracy,
  yardstick::roc_auc
)

# prepare predictions
test_predictions <- testing(biomarker_split) %>%
  mutate(
    pred = predict(
      fit,
      newdata = testing(biomarker_split),
      type = "response"
    ),
    truth = factor(class, levels = c(FALSE, TRUE)),
    estimate = factor(pred > 0.5, levels = c(FALSE, TRUE))
  )

# evaluate model
mod2 <- class_metrics(
  test_predictions,
  truth = truth,
  estimate = estimate,
  pred,
  event_level = "second"
)

# Modification 3

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)


# hard intersection
proteins_hard <- intersect(proteins_s1, proteins_s2)

# fuzzy intersection
# put t-test and RF scores for all proteins in one table
rank_ttest <- ttests_out %>%
  select(protein) %>%
  mutate(rank_ttest = row_number())

rank_rf <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(rank_rf = row_number()) %>%
  select(protein, rank_rf)

# change each rank into a membership grade (1=best, 0=worst)
n_prot <- nrow(rank_ttest)
fuzzy_tbl <- rank_ttest %>%
  inner_join(rank_rf, by = 'protein') %>%
  mutate(mem_ttest = 1 - (rank_ttest - 1) / (n_prot - 1),
         mem_rf = 1 - (rank_rf - 1) / (n_prot - 1),
         mem_fuzzy = pmin(mem_ttest, mem_rf)) %>%
  arrange(desc(mem_fuzzy))

# apply cutoff (0.95)
proteins_fuzzy <- fuzzy_tbl %>%
  filter(mem_fuzzy >= 0.95) %>%
  pull(protein)

# compare the two panels
proteins_hard
proteins_fuzzy

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- proteins_fuzzy

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
test_results <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(truth = factor(class, levels = c(FALSE, TRUE)),
         estimate = factor(pred > 0.5, levels = c(FALSE, TRUE)))

mod3 <- test_results %>%
  class_metrics(truth = truth,
                estimate = estimate,
                pred,
                event_level = 'second')

# complete table
comparison_table <- bind_rows(
  Base = base,
  `Modification 1` = mod1,
  `Modification 2` = mod2,
  `Modification 3` = mod3,
  .id = "Model"
) %>%
  select(Model, .metric, .estimate) %>%
  pivot_wider(names_from = Model, values_from = .estimate) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

knitr::kable(
  comparison_table,
  caption = "Comparison of Classifier Performance",
  col.names = c(
    "Metric", "Base", "Modification 1",
    "Modification 2", "Modification 3"
  )
)


# QUESTION 4

# Keep the Q3 results for the final comparison.
q3_results <- comparison_table
q3_baseline <- base
q3_model_results <- bind_rows(
  `In-class baseline (Q3)` = base,
  `Q3 modification 1: early split` = mod1,
  `Q3 modification 2: top 25` = mod2,
  `Q3 modification 3: fuzzy` = mod3,
  .id = "model"
) %>%
  select(model, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

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

# Check that Q3 and Q4 reproduce the same class baseline.
for (.q4_metric in q3_baseline$.metric) {
  stopifnot(abs(q3_baseline$.estimate[q3_baseline$.metric == .q4_metric] -
                  baseline_fit$metrics[[.q4_metric]]) < 1e-10)
}
rm(.q4_metric)

# Combine the Q3 and Q4 metrics with one class baseline.
combined_results <- bind_rows(
  q3_model_results,
  final_results %>%
    filter(model != "In-class baseline (rerun)") %>%
    select(model, sensitivity, specificity, accuracy, roc_auc)
)
print(combined_results, width = Inf)
write_csv(combined_results,
          file.path(output_dir, "combined-q3-q4-comparison.csv"))
