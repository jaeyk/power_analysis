
if (!require(pacman)) install.packages("pacman")

pacman::p_load(tidyverse, here)

review <- read_csv(here("raw_data", "linked_fate_review.csv"))

sum(review$Asian == 1, na.rm = T)/nrow(review)

max(review$Pub.year)