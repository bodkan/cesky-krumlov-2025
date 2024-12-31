library(slendr)
init_env()

source("utils.R")





# part 1 -- building a purely neutral model -------------------------------

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
  gene_flow(from = ana, to = yam, rate = 0.4, start = 7900, end = 7000),
  gene_flow(from = ana, to = eur, rate = 0.5, start = 6000, end = 5000),
  gene_flow(from = yam, to = eur, rate = 0.65, start = 4000, end = 3500)
)

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





# part 2 -- simulating a tree sequence and computing Tajima's D -----------

# because the focus of this exercise is selection, we'll be using slendr's
# SLiM simulation engine, not msprime (although this part of the exercise
# is still just a neutral simulation)
ts <- slim(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = schedule) %>%
  ts_recapitate(Ne = 5000, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)

# still, notice just how much faster an msprime simulation run is compared to SLiM!
msprime(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = schedule) %>%
  ts_mutate(mutation_rate = 1e-8)

# inspect the table of all individuals recorded in our tree sequence
ts_samples(ts)

# tskit functions in slendr generally operate on vectors (or lists) of individual
# names, like those produced by ts_samples() above -- let's get get a vector of
# such names using another helper function ts_names()
samples <- ts_names(ts, split = "pop")
samples

# compute genome-wide Tajima's D for each population -- note that we don't
# expect to see any significant differences because no population experienced
# natural selection (yet)
ts_tajima(ts, sample_sets = samples)





# part 3 -- computing Tajima's D in windows -------------------------------

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





# part 4 -- simulating positive selection ---------------------------------

# Our SLiM selection extension for slendr is in 'exercise7_slim.txt'. When
# you inspect that file, you will se pretty much standard SLiM code with one
# exception: the somewhat strange {{elements}} in curly brackets. Those
# parameters of the selection model, which must be first "instantiated" using
# the function `substitute_values()`.
extension <- substitute_values(
  template = "exercise7_slim.txt",
  origin_pop = "EHG",
  s = 0.1,
  onset_time = 12000
)

# When we take a look at the modified extension contents in the terminal using
# something like the unix less command, we no longer see the {{parameters}} but
# we see concrete values instead.
extension

# Using the SLiM extension is simple -- we simply provide it as an additional
# argument to the `compile_model()` function:
model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30,
  extension = extension   # <======== this is different to the compilation above!
)
# We can finally run our selection simulation!

# This time our model not only produces a tree sequence, but it also generates
# a table of allele frequencies in each population (see the file
# 'exercise7_slim.txt'). We need to be able to load both of these files after
# the simulation and thus need a path to a location we can find those files.
# We can do this by specifying `path = TRUE`.
path <- slim(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = schedule, path = TRUE)

# We can verify that the path not only contains a tree-sequence file but also
# the table of allele frequencies.
list.files(path)

# allele frequency trajectories -------------------------------------------

traj_df <- read_trajectory(path)
traj_df

plot_trajectory(traj_df)

# Comparing the trajectories side-by-side with the demographic model reveals
# some obvious patterns of both selection and demographic history.
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
