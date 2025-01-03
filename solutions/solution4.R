library(slendr)
init_env()

source(here::here("utils.R"))





# Create a model, simulate tree sequence, get EIGENSTRAT data -------------

popZ <- population("popZ", time = 3000, N = 5000)
popX <- population("popX", time = 1500, N = 5000, parent = popZ)
popY <- population("popY", time = 1500, N = 5000, parent = popZ)

model <- compile_model(list(popZ, popX, popY), generation_time = 1)

schedule <- schedule_sampling(model, times = seq(3000, 0, by = -200), list(popZ, 10), list(popX, 10), list(popY, 10))

# Use the function plot_model() to make sure that the model and the sampling schedule
# are defined correctly (there's no such thing as too many sanity checks when doing research)
plot_model(model, proportions = TRUE, samples = schedule)

ts <- msprime(model, samples = schedule, sequence_length = 50e6, recombination_rate = 1e-8, random_seed = 1702182272) %>%
  ts_mutate(1e-8)

# Save the samples table to a data frame (and process it a bit for tidier plotting later)
samples <- ts_samples(ts) %>% mutate(pop = factor(pop, levels = c("popZ", "popX", "popY")))

samples %>% group_by(pop) %>% count()





# Part 1 -- create EIGENSTRAT for different subsets of samples ------------

# Our goal is to investigate how is PCA structure affected by different factors related
# to sampling (sample sizes, dates of samples, etc.). Right now, our tree sequence contains
# every single sample. What we will do now is create smaller subsets of the entire
# tree sequence to only defined sets of samples using the process called "simplification",
# then export the genomic data from the tree sequence into the EIGENSTRAT file format.

# EIGENSTRAT with only "present-day" X and Y individuals
subset <- filter(samples, pop %in% c("popX", "popY"), time == 0)
subset

ts_XY0 <- ts_simplify(ts, simplify_to = subset$name)
ts_eigenstrat(ts_XY0, "data/XY0")

# EIGENSTRAT with all X and Y individuals, present-day and ancient
subset <- filter(samples, pop %in% c("popX", "popY"))
subset

ts_XYall <- ts_simplify(ts, simplify_to = subset$name)
ts_eigenstrat(ts_XYall, "data/XYall")

# EIGENSTRAT with only "present-day" X,Y, and Z individuals
subset <- filter(samples, time == 0)
print(subset, n = Inf)

ts_XYZ0 <- ts_simplify(ts, simplify_to = subset$name)
ts_eigenstrat(ts_XYZ0, "data/XYZ0")

# EIGENSTRAT file with all individuals
ts_eigenstrat(ts, "data/XYZall")

# We can verify the contents of the EIGENSTRAT data by inspecting them in
# the bash terminal (cat, less, etc.). Alternatively, we can use the R
# package admixr to do this in R. For instance, we could do the following.

library(admixr)

eigen <- eigenstrat("data/XYZ0")
eigen

read_ind(eigen)
read_snp(eigen)
read_geno(eigen)



# Part 2 -- computing and visualizing PCA patterns ------------------------

plot_pca("data/XY0", ts_XY0, color_by = "pop")
plot_pca("data/XYall", ts_XYall, color_by = "time")
plot_pca("data/XYZ0", ts_XYZ0, color_by = "pop")

# The first two PCs can separate lineages from the split onwards, but can't
# distinguish the ancestral Z individuals
plot_pca("data/XYZall", ts, pc = c(1, 2), color_by = "time")
plot_pca("data/XYZall", ts, pc = c(1, 2), color_by = "pop")

# PC 2 vs 3 -- PC 2 doesn't really separate Z before and after the split, but
# PC 3 does reveal the tree structure!
plot_pca("data/XYZall", ts, pc = c(2, 3), color_by = "time")
plot_pca("data/XYZall", ts, pc = c(2, 3), color_by = "pop")

# PC 3 vs 4 -- PC 3 separates individuals before and after split, but doesn't
# really separate the populations. PC 4 does, kind of?
plot_pca("data/XYZall", ts, pc = c(3, 4), color_by = "time")
plot_pca("data/XYZall", ts, pc = c(3, 4), color_by = "pop")

# PC 4 vs 5 -- Hard to say what's going on, but it looks pretty. :)
plot_pca("data/XYZall", ts, pc = c(4, 5), color_by = "time")
plot_pca("data/XYZall", ts, pc = c(4, 5), color_by = "pop")

# ... same!
plot_pca("data/XYZall", ts, pc = c(5, 6), color_by = "time")
plot_pca("data/XYZall", ts, pc = c(5, 6), color_by = "pop")





# Bonus -- code up a PCA using another software you know ------------------

# In addition to `ts_eigenstrat()`, slendr also provides an interface to
# tskit's VCF functionality as its R function `ts_vcf()`. Use these tools
# to test whether your own PCA software of choice produces the same results!

