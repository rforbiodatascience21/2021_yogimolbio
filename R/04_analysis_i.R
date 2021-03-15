# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("dplyr")
library("purrr")
library("broom")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
gravier_clean_aug <- read_tsv(file = "data/03_gravier_clean_aug.tsv.gz")


# Wrangle data ------------------------------------------------------------
gravier_clean_aug_long <- 
  gravier_clean_aug %>%
  pivot_longer(cols = contains("g"),
               names_to = "gene",
               values_to = "log2_expr_level")


gravier_clean_aug_long_nested <- 
  gravier_clean_aug_long %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()

set.seed(7) 
gravier_clean_aug_long_nested <-
  gravier_clean_aug_long_nested %>% 
  sample_n(size = 100)

gravier_clean_aug_long_nested <- gravier_clean_aug_long_nested %>% 
  mutate(mdl = map(data, ~glm(outcome ~ log2_expr_level, 
                              data = .x,
                              family = binomial(link = "logit"))))

gravier_clean_aug_long_nested <- gravier_clean_aug_long_nested %>%
  mutate(mdl_tidy = map(mdl, tidy, conf.int = TRUE)) %>% 
  unnest(mdl_tidy)

gravier_clean_aug_long_nested <- 
  gravier_clean_aug_long_nested %>% 
  filter(term == "log2_expr_level" ) %>% 
  mutate( identified_as = case_when(p.value < 0.05 ~ "significant",
                                    p.value > 0.05 ~"non_significant"))

gravier_clean_aug_wide = gravier_clean_aug %>%
  select(outcome, pull(gravier_clean_aug_long_nested, gene))

# Model data
pca_fit <- gravier_clean_aug_wide %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale = TRUE)


# Visualise data ----------------------------------------------------------
pca_fit %>%
  augment(gravier_clean_aug_wide) %>% 
  mutate(outcome = factor(outcome)) %>% 
  ggplot(mapping = aes(x = .fittedPC1,
                       y = .fittedPC2,
                       color = outcome))+
  geom_point(size = 1.5) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title = "PCA plot")

# Write data --------------------------------------------------------------
write_tsv(x = gravier_clean_aug_wide,
          path = "data/04_gravier_clean_aug_wide.tsv.gz")
ggsave(filename = "results/04_PCA_plot.png", width = 16, height = 9, dpi = 72)
