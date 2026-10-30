library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

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

test_results %>%
  class_metrics(truth = truth,
                estimate = estimate,
                pred,
                event_level = 'second')