# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
gravier_clean <- read_tsv(file = "data/02_gravier_clean.tsv.gz")


# Wrangle data ------------------------------------------------------------
gravier_clean_aug <- gravier_clean %>% 
  mutate(outcome = case_when(value == "good" ~ 0,
                                            value == "poor" ~ 1))


# Write data --------------------------------------------------------------
write_tsv(x = gravier_clean_aug,
          path = "data/03_gravier_clean_aug.tsv.gz")
