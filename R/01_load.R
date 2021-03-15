# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
load("data/_raw/gravier.RData")


# Wrangle data ------------------------------------------------------------
 
gravier_x <- gravier %>% 
  pluck(.,"x") %>% 
  as_tibble()

gravier_y <- gravier %>% 
  pluck(.,"y") %>% 
  as_tibble()



# Write data --------------------------------------------------------------
write_tsv(x = gravier_x,
          path = "data/01_gravier_x.tsv.gz")
write_tsv(x = gravier_y,
          path = "data/01_gravier_y.tsv.gz")
