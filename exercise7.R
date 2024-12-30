library(slendr)
init_env()

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# African ancestral population
afr <- population("AFR", time = 65000, N = 5000)

# first migrants out of Africa
ooa <- population("OOA", parent = afr, time = 60000, N = 5000, remove = 27000)

# Eastern hunter-gatherers
ehg <- population("EHG", parent = ooa, time = 28000, N = 5000, remove = 6000)

# European population
eur <- population("EUR", parent = ehg, time = 25000, N = 5000)

# Anatolian farmers
ana <- population("ANA", time = 28000, N = 5000, parent = ooa, remove = 4000)

# Yamnaya steppe population
yam <- population("YAM", time = 8000, N = 5000, parent = ehg, remove = 2500)

# define gene-flow events
gf <- list(
  gene_flow(from = ana, to = yam, rate = 0.4, start = 7900, end = 7800),
  gene_flow(from = ana, to = eur, rate = 0.5, start = 6000, end = 5000),
  gene_flow(from = yam, to = eur, rate = 0.65, start = 4000, end = 3500),
  gene_flow(from = afr, to = eur, rate = 0.07, start = 3000, end = 0)
)

# purely neutral model ----------------------------------------------------

model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30
)

schedule <- rbind(
  schedule_sampling(model, times = 0, list(eur, 50), list(afr, 50)),
  schedule_sampling(model, times = 6000, list(ehg, 50)),
  schedule_sampling(model, times = 4000, list(ana, 50)),
  schedule_sampling(model, times = 2500, list(yam, 50))
)

plot_model(model, proportions = TRUE)
plot_model(model, proportions = TRUE, samples = schedule)

ts <- msprime(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = schedule) %>%
  ts_mutate(mutation_rate = 1e-8)

# get a vector of names of every recorded individual
samples <- ts_names(ts, split = "pop")
samples

# compute genome-wide Tajima's D for each population
ts_tajima(ts, sample_sets = samples)

# pre-compute genomic windows for window-based computation of Tajima's D
windows <- as.integer(seq(0, ts$sequence_length, length.out = 100))
windows

# compute genome-wide Tajima's D for each population in individual windows
tajima_wins <- ts_tajima(ts, sample_sets = samples, windows = windows)
tajima_wins
tajima_wins[1, ]$D

# the numeric vector format of the result above is a bit annoying to work with
# and visualize, so let's use some tidyverse "magic" to format it properly
tajima_df <- process_tajima(tajima_wins)
tajima_df

# now let's visualize the window-based Tajima's D along the simulated genome
# (hint: we don't actually expect anything interesting here because we ran
# a purely neutral simulation, so this is more like a control)
plot_tajima(tajima_df)

# model with selection ----------------------------------------------------

extension <- substitute_values(
  template = "exercise7_slim.txt",
  origin_pop = "EUR",
  s = 0.1,
  onset_time = 10000
)

model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30,
  extension = extension
)

tstart <- Sys.time()
unlink("exercise7")

slim(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = schedule, path = "exercise7")

tend <- Sys.time()
tend - tstart

dir("exercise7")

# allele frequency trajectories -------------------------------------------

traj_df <- read_trajectory(path = "exercise7", origin_pop = "EUR", s = 0.1, onset = 10000)
traj_df

plot_grid(
  plot_model(model, proportions = TRUE),
  plot_trajectory(traj_df),
  nrow = 1, rel_widths = c(0.7, 1)
)

# Tajima's D --------------------------------------------------------------

ts <- ts_read("exercise7/slim.trees", model)

ts_coalesced(ts, return_failed = TRUE) %>% sample(1)

ts_tree(ts, i = .Last.value) %>% ts_draw(labels = TRUE)

  ts_recapitate(Ne = 5000, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)

samples <- ts_names(ts, split = "pop")
samples

ts_tajima(ts, sample_sets = samples)

windows <- as.integer(seq(0, ts$sequence_length, length.out = 100))
windows

tajima_wins <- ts_tajima(ts, sample_sets = samples, windows = windows)

tajima_df <-
  tajima_wins %>%
  unnest %>%
  group_by(set) %>%
  mutate(window = row_number()) %>%
  ungroup
tajima_df

ggplot(tajima_df, aes(window, D, color = set)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 50, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(ylim = c(-4, 4)) +
  theme_minimal()
