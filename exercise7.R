devtools::load_all("~/Projects/slendr")
init_env()

# African ancestral population
afr <- population("AFR", time = 65000, N = 3000)

# first migrants out of Africa
ooa <- population("OOA", parent = afr, time = 60000, N = 500, remove = 27000) %>%
  resize(N = 2000, time = 40000, how = "step")

# Eastern hunter-gatherers
ehg <- population("EHG", parent = ooa, time = 28000, N = 1000, remove = 6000)

# European population
eur <- population("EUR", parent = ehg, time = 25000, N = 5000)

# Anatolian farmers
ana <- population("ANA", time = 28000, N = 3000, parent = ooa, remove = 4000)

# Yamnaya steppe population
yam <- population("YAM", time = 8000, N = 500, parent = ehg, remove = 2500)

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

# model with selection ----------------------------------------------------

extension <- substitute_values(
  template = "exercise7_slim.txt",
  s = 0.03, onset_time = 15000,
  origin_pop = "EHG"
)

model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30,
  extension = extension
)

ts <- slim(model, sequence_length = 1e6, recombination_rate = 1e-8, samples = schedule, verbose = TRUE) %>%
  ts_recapitate(recombination_rate = 1e-8, Ne = 3000)

tstart <- Sys.time()

slim(model, sequence_length = 1e6, recombination_rate = 1e-8, samples = schedule, path = "exercise7", verbose = T)

tend <- Sys.time()
tend - tstart

dir("exercise7")

# allele frequency trajectories -------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)

traj_df <- read_tsv("exercise7/trajectory_15000_EHG.tsv") %>%
  mutate(pop = factor(pop, levels = c("AFR", "OOA", "EHG", "ANA", "EUR", "YAM")))

p_traj <- ggplot(traj_df, aes(time, freq, color = pop)) +
  geom_line() +
  scale_x_reverse() +
  facet_wrap(~ pop) +
  coord_cartesian(ylim = c(0, 1))

plot_grid(plot_model(model, proportions = TRUE), p_traj, nrow = 1)

# Tajima's D --------------------------------------------------------------

ts <- ts_read("exercise7/slim.trees", model) %>%
  ts_recapitate(Ne = 3000, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)

samples <- ts_names(ts, split = "pop")
samples

ts_tajima(ts, sample_sets = samples)
