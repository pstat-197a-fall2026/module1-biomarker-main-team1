library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('../data/biomarker-clean.RData')

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
biomarker_sstar_test %>%
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


# results from modification 1

#.metric     .estimator .estimate
#<chr>       <chr>          <dbl>
#1 sensitivity binary         0.75 
#2 specificity binary         0.8  
#3 accuracy    binary         0.774
#4 roc_auc     binary         0.871


# The sensitivity, accuracy, and roc_auc from the modified analysis are all lower than the original metrics.
# However, this does not mean that our model performance worsened after the modification.
# Instead, the initial analysis likely had data leakage, leading it to overestimate the model performance
# By doing the training/testing split before any of the analysis, we do the feature extraction steps with only the training data.
# Thus, the metrics after modification are a better representation of our actual model performance.




