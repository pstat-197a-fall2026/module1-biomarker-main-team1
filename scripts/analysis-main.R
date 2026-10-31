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
