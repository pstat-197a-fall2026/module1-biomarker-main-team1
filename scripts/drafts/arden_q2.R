library(tidyverse)
library(kableExtra)

# get names
var_names <- read_csv('../../data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# read in data
biomarker_clean <- read_csv('../../data/biomarker-raw.csv', 
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
  theme_minimal() +
  theme(axis.text.y = element_blank())

# As we can see, subject 154 has by far the most outlier protein measurements.
# Subjects 108, 9, 121, 52, and 77 have over 100 protein measurements with
# outlier values, separating themselves from the rest of the group having
# 77 or less outlier values. It appears that the typically developing group
# may have more outliers than the ASD group, though this may be due to small
# sample sizes.
