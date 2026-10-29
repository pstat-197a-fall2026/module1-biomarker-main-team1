library(tidyverse)
library(ggplot2)
library(tidyr)


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

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean <- read_csv('../../data/biomarker-raw.csv', 
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