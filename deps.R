# This script contains a list of packages used by scripts in this repository.
# You can safely ignore this. If you want to set up your machine to get all
# the R package dependencies, follow the instructions at this link:
# https://github.com/bodkan/cesky-krumlov-2025

# renv::init(bare = TRUE)

# install.packages(c("combinat", "cowplot", "dplyr", "tidyr", "ggplot2", "rmarkdown", "yaml", "slendr"))

# renv::snapshot()

library(combinat)
library(cowplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rmarkdown)
library(yaml)
library(slendr)
setup_env(agree = TRUE)
