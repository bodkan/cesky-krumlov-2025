# This script contains a list of packages used by scripts in this repository.
# You can safely ignore this. If you want to set up your machine to get all
# the R package dependencies, follow the instructions at this link:
# https://github.com/bodkan/cesky-krumlov-2025

install.packages(c("combinat", "cowplot", "dplyr", "ggplot2", "readr", "scales", "tidyr", "slendr"))
                   # "admixr", "ggrepel", "viridis", "smartsnp", "sf", "stars", "rnaturalearth",
                   # "devtools"))
# devtools::install_github("bodkan/slendr")

library(combinat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(tidyr)
library(rmarkdown)
library(yaml)
library(slendr)
setup_env(agree = TRUE)
