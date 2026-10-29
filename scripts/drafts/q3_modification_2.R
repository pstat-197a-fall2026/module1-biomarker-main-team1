library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('../../data/biomarker-clean.RData')


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
class_metrics(
  test_predictions,
  truth = truth,
  estimate = estimate,
  pred,
  event_level = "second"
)

# results from modification 2

## # A tibble: 4 × 3
##   .metric     .estimator .estimate
##   <chr>       <chr>          <dbl>
## 1 sensitivity binary         0.812
## 2 specificity binary         0.867
## 3 accuracy    binary         0.839
## 4 roc_auc     binary         0.946

# The number of proteins selected from both feature selection methods was changed from 10 to 25.
# From the results, the classification performance showed mixed changes for this modification. 
# Sensitivity decreased, accuracy remained the same, and specificity and roc_auc increased.
# This suggests that the modified model is a bit less sensitive to ASD while but a bit more specific to TD.
# The roc_auc increase may mean that the model got better at distinguishing between ASD and TD individuals from predicted probabilities.
# But, the classification accuracy is the same.
# Overall it seems increasing the number of selected proteins improves some aspects of the model performance but worsens otehrs.