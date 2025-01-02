library(slendr)
init_env()

source("utils.R")





# Part 1 -- building a purely neutral model -------------------------------

# African ancestral population
afr <- population("AFR", time = 65000, N = 5000)

# First migrants out of Africa
ooa <- population("OOA", parent = afr, time = 60000, N = 5000, remove = 27000)

# Eastern hunter-gatherers
ehg <- population("EHG", parent = ooa, time = 28000, N = 5000, remove = 6000)

# European population
eur <- population("EUR", parent = ehg, time = 25000, N = 5000)

# Anatolian farmers
ana <- population("ANA", time = 28000, N = 5000, parent = ooa, remove = 4000)

# Yamnaya steppe population
yam <- population("YAM", time = 8000, N = 5000, parent = ehg, remove = 2500)

# Define gene-flow events
gf <- list(
  gene_flow(from = ana, to = yam, rate = 0.75, start = 7500, end = 6000),
  gene_flow(from = ana, to = eur, rate = 0.5, start = 6000, end = 5000),
  gene_flow(from = yam, to = eur, rate = 0.6, start = 4000, end = 3500)
)

# Compile all populations into a single slendr model object
model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30
)

# Schedule the sampling from four European populations roughly before their
# disappearance (or before the end of the simulation)
schedule <- rbind(
  schedule_sampling(model, times = 0, list(eur, 50)),
  schedule_sampling(model, times = 6000, list(ehg, 50)),
  schedule_sampling(model, times = 4000, list(ana, 50)),
  schedule_sampling(model, times = 2500, list(yam, 50))
)

# Verify visually that our model is correct
plot_model(model, samples = schedule)
plot_model(model, proportions = TRUE)





# Part 2 -- simulating a tree sequence and computing Tajima's D -----------

# Although the point of this exercise is to simulate selection, let's first
# simulate a normal neutral model using slendr's msprime engine as a sanity check
ts <- msprime(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = schedule) %>%
  ts_mutate(mutation_rate = 1e-8)

# Inspect the table of all individuals recorded in our tree sequence
ts_samples(ts)

# tskit functions in slendr generally operate on vectors (or lists) of individual
# names, like those produced by ts_samples() above -- let's get get a vector of
# such names using another helper function ts_names()
samples <- ts_names(ts, split = "pop")
samples

# Compute genome-wide Tajima's D for each population -- note that we don't
# expect to see any significant differences because no population experienced
# natural selection (yet)
ts_tajima(ts, sample_sets = samples)





# Part 3 -- computing Tajima's D in windows -------------------------------

# Pre-compute genomic windows for window-based computation of Tajima's D
windows <- round(seq(0, ts$sequence_length, length.out = 100))
windows

# Compute genome-wide Tajima's D for each population in individual windows
tajima_wins <- ts_tajima(ts, sample_sets = samples, windows = windows)
tajima_wins
tajima_wins[1, ]$D

# The numeric vector format of the result above is a bit annoying to work with
# and visualize, so let's use some tidyverse "magic" to format it properly
# using a helper function `process_tajima()` I wrote for you
tajima_df <- process_tajima(tajima_wins)
tajima_df

# Now let's visualize the window-based Tajima's D along the simulated genome
# using another helper function `plot_tajima()` (hint: we still don't expect
# anything interesting here because we ran a purely neutral simulation, so this
# is more like a control)
plot_tajima(tajima_df)





# Part 4 -- simulating positive selection ---------------------------------

# Our SLiM selection extension for slendr is in 'exercise5_slim.txt'. When
# you inspect that file, you will se pretty much standard SLiM code with one
# exception: the somewhat strange {{elements}} in curly brackets. Those
# parameters of the selection model, which must be first "instantiated" using
# the function `substitute_values()`.

# First, let's simulate selection happening only in the EUR lineage.
extension <- substitute_values(template = "exercise5_slim.txt",
  origin_pop = "EUR",
  s = 0.2,
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
  extension = extension   # <======== this is missing in the neutral example!
)
# We can finally run our selection simulation

# This time our model not only produces a tree sequence, but it also generates
# a table of allele frequencies in each population (see the file
# 'exercise5_slim.txt'). We need to be able to load both of these files after
# the simulation and thus need a path to a location we can find those files.
# We can do this by specifying `path = TRUE`.
path <- slim(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = schedule, path = TRUE, random_seed = 59879916)

# We can verify that the path not only contains a tree-sequence file but also
# the table of allele frequencies.
list.files(path)

# Allele frequency trajectories -------------------------------------------

traj_df <- read_trajectory(path)
traj_df

plot_trajectory(traj_df)

# Comparing the trajectories side-by-side with the demographic model reveals
# some obvious patterns of both selection and demographic history.
plot_grid(
  plot_model(model),
  plot_trajectory(traj_df),
  nrow = 1, rel_widths = c(0.7, 1)
)

# Tajima's D --------------------------------------------------------------

# We need an additional call to `ts_recapitate()` to ensure that the tree
# sequence is fully coalesced (and thus contains older polymorphisms)
ts <- ts_read(file.path(path, "slim.trees"), model) %>%
  ts_recapitate(Ne = 5000, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)

samples <- ts_names(ts, split = "pop")
samples

# Perhaps not unexpectedly, gneome-wide Tajima's D doesn't reveal any significant
# deviations even in this case.
ts_tajima(ts, sample_sets = samples)

windows <- as.integer(seq(0, ts$sequence_length, length.out = 100))

# compute genome-wide Tajima's D for each population in individual windows
tajima_wins <- ts_tajima(ts, sample_sets = samples, windows = windows)
tajima_df <- process_tajima(tajima_wins)
tajima_df

plot_tajima(tajima_df)





# Bonus exercises ---------------------------------------------------------





# Bonus 1 -----------------------------------------------------------------

# Vary the uniform recombination rate and observe what happens with Tajima's D
# in windows along the genome.

# Solution: just modify the value of the `recombination_rate =` argument provided
# to the `slim()` function above.





# Bonus 2 -----------------------------------------------------------------

# Simulate the origin of the beneficial allele in the EHG population -- what
# do the trajectories look like now? How does that change the Tajima's D
# distribution along the genome in our European populations?

# Solution -- use this extension in the `slim()` call, and repeat the rest.
extension <- substitute_values(
  template = "exercise5_slim.txt",
  origin_pop = "EHG",
  s = 0.1,
  onset_time = 12000
)
model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30,
  extension = extension
)





# Bonus 3 -- practice your tidyverse chops --------------------------------

# Re-implement my `process_tajima()` function on your own. Hard mode: do this
# using base R (no tidyverse allowed!).
