library(slendr)
init_env()






# Part 1: Building a demographic model in slendr --------------------------

# Chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 5000)

# Two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)

# Neanderthal population splitting at 600 ky ago from modern humans
# (becomes extinct by 40 ky ago)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# Neanderthal introgression event (3% admixture between 55-50 kya)
gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 50000)

# Compile the entire model into a single slendr R object
model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,
  generation_time = 30
)







# Part 2: Inspecting the model visually -----------------------------------

# Verify the correctness of our model visually (different arguments of
# the plotting function can make certain features more (or less) visible)
plot_model(model)
plot_model(model, sizes = FALSE)
plot_model(model, sizes = FALSE, log = TRUE)
plot_model(model, sizes = FALSE, log = TRUE, proportions = TRUE)






# Part 3: Simulating genomic data -----------------------------------------

ts <- msprime(model, sequence_length = 1e6, recombination_rate = 1e-8)

# `debug = TRUE` instructs slendr's built-in msprime script to print out
# msprime's own debugger information! This can be very useful for debugging,
# in addition to the visualization of the model as shown above.
ts <- msprime(model, sequence_length = 1e6, recombination_rate = 1e-8, debug = TRUE)

# For debugging of technical issues (with msprime, with slendr, or both), it is
# very useful to have the `msprime()` function dump the "raw" command-line to
# run the simulation on the terminal using plain Python interpreter
msprime(model, sequence_length = 1e6, recombination_rate = 1e-8, run = FALSE)






# Part 4: Inspecting the tree-sequence object -----------------------------

# Typing out the object with the result shows that it's a good old tskit
# tree-sequence object
ts

# slendr provides a helper function which allows access to all the low-level
# components of every tree-sequence object
ts_table(ts, "nodes")
ts_table(ts, "edges")
ts_table(ts, "individuals")
ts_table(ts, "mutations")
ts_table(ts, "sites")


# slendr provides a convenient function `ts_samples()` which allows us to
# inspect the contents of a simulated tree sequence in a more human-readable,
# simplified way. We can see that our tree sequence contains a massive number
# of individuals. Too many, in fact (we recorded every single individual alive
# at the end of our simulation -- something we're unlikely to be ever lucky
# enough to have, regardless of which species we study)
ts_samples(ts) %>% nrow()





# Part 5 -- scheduling sampling events ------------------------------------

# We can precisely define which individuals (from which populations, and at
# which times) should be recorded in a tree sequence using the slendr
# function `schedule_sampling()`.

# Here we scheduled the sampling of two Neanderthals at 70kya and 40kya
nea_samples <- schedule_sampling(model, times = c(70000, 40000), list(nea, 1))
nea_samples # (this function produces a plain old data frame!)

# Here we schedule one Chimpanzee sample, 5 African samples, and 10 European samples
present_samples <- schedule_sampling(model, times = 0, list(chimp, 1), list(afr, 5), list(eur, 10))

# We also schedule the recording of one European sample between 50kya and 2kya,
# every 2000 years
emh_samples <- schedule_sampling(model, times = seq(40000, 2000, by = -2000), list(eur, 1))

# Because those functions produce nothing but a data frame, we can bind
# individual sampling schedules together
schedule <- rbind(nea_samples, present_samples, emh_samples)
schedule

# The schedules can be also visualized on the graphical representation of
# a slendr model (although the result is often a bit wonky)
plot_model(model, sizes = FALSE, log = TRUE, samples = schedule)






# Part 6 -- simulating a defined set of individuals -----------------------

# Let's simulate a more realistic set of samples, using the schedule we've
# defined above

# The command below will likely take a few minutes to run, so feel free to go
# down from 100 Mb sequence_length to even 10Mb (it doesn't matter much)
# tstart <- Sys.time()
ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8, samples = schedule, random_seed = 1269258439) %>%
  ts_mutate(mutation_rate = 1e-8, random_seed = 1269258439)
# tend <- Sys.time()
# tend - tstart # Time difference of 2.141642 mins

# We can save a tree sequence object using a slendr function `ts_write()` (this
# can be useful if we want to save the results of a simulation for later use).
dir.create("data", showWarnings = FALSE)
ts_write(ts, "data/introgression.trees")

# A saved slendr tree sequence can be later read via function `ts_read()` (but
# note that this function needs slendr-specific metadata stored in a model
# object!)
ts <- ts_read("data/introgression.trees", model = model)

# Inspect the (tskit/Python-based) summary of the (now smaller) tree sequence
ts


# Get the table of all recorded samples in the tree sequence
ts_samples(ts)

# Compute the count of individuals in different time points
library(dplyr)

ts_samples(ts) %>% group_by(pop, time == 0) %>% tally %>% select(pop, n)
