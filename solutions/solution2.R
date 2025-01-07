library(slendr)
init_env()

library(dplyr)
library(ggplot2)

# (the following model is copied from exercise1.R) ------------------------

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

# This can allow us to save a bit of computational time (and also enhances
# reproducibility because anyone with the .trees file can follow the downstream
# workflow)
ts <- ts_read("data/introgression.trees", model = model)

# It's always a good idea to check that the data we're working with is
# really what we think it is!
ts_samples(ts) %>% nrow
ts_samples(ts) %>% group_by(pop, time == 0) %>% tally %>% select(pop, n)






# Part 1: Computing nucleotide diversity ----------------------------------

# Let's first get a named list of individuals in each group we want to be
# working with (slendr tree-sequence statistic functions generally operate
# with this kind of structure)
sample_sets <- ts_names(ts, split = "pop")
sample_sets

# Compute nucleotide diversity (pi) in each population
pi_pop <- ts_diversity(ts, sample_sets = sample_sets)
arrange(pi_pop, diversity)

# Compute the same thing in each individual separately (so computing the
# individual-based heterozygosity). We can do this by passing the vector of
# individual names directory as the `sample_sets =` argument, rather than
# in a list of groups as we did above.

# For convenience, we first get a table of all individuals (which of course
# contains also their names) and just add their heterozygosities as a new column.
pi_df <- ts_samples(ts)
pi_df$name

pi_df$diversity <- ts_diversity(ts, sample_sets = pi_df$name)$diversity
pi_df

ggplot(pi_df, aes(pop, diversity, color = pop, group = pop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw()





# Part 2: Computing pairwise divergence -----------------------------------

# We will again use the `sample_sets` list of individual names defined above
sample_sets

div_df <- ts_divergence(ts, sample_sets)
arrange(div_df, divergence)






# Part 3: Detecting Neanderthal admixture in Europeans --------------------

#                                BABA - ABBA
# f4(AFR, Test; NEA, CHIMP) ~   -------------
#                                   #SNPs

# We will be comparing two modes of tskit tree-sequence-based computation
mode <- "branch"
# mode <- "site"

# Comparing two Africans vs Neanderthal should not reveal any deviation from
# the null hypothesis (this should be consistent with a tree with no admixture)
f4_null <- ts_f4(ts, W = "AFR_1", X = "AFR_2", Y = "NEA_1", Z = "CHIMP_1", mode = mode)
f4_null

# On the other hand, an African-European comparison should reveal an excess
# of sharing of Neanderthal alleles with Europeans (i.e. more ABBA sites)
f4_alt <- ts_f4(ts, W = "AFR_1", X = "EUR_1", Y = "NEA_1", Z = "CHIMP_1", mode = mode)
f4_alt

# We can see that the second test has ~50 times higher f3, although this is
# not a real test of significance (no Z-score or standard error as given by
# jackknife procedure in ADMIXTOOLS)
f4_alt$f4 / f4_null$f4


# Part 4: Detecting Neanderthal admixture in Europeans --------------------

# Let's compute the f4 statistic for all Africans and Europeans to see the
# f4 introgression patterns more clearly
f4_afr <- lapply(sample_sets$AFR, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1", mode = mode)) %>% bind_rows()
f4_afr
f4_eur <- lapply(sample_sets$EUR, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1", mode = mode)) %>% bind_rows()
f4_eur

# Let's add population columns to each of the two results, and bind them together
# for plotting
f4_afr$pop <- "AFR"
f4_eur$pop <- "EUR"

f4_results <- rbind(f4_afr, f4_eur)

f4_results %>%
  ggplot(aes(pop, f4, color = pop)) +
  geom_boxplot() +
  geom_jitter() +
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("f4(AFR, EUR; NEA, CHIMP)") +
  theme_bw()

# Why the difference between the "branch" and "site" modes?
# See this tutorial (and particularly directly the linked section):
# https://tskit.dev/tutorials/no_mutations.html#genealogy-based-measures-are-less-noisy




# Bonus exercises ---------------------------------------------------------






# Bonus 1 -- outgroup f3 statistic ----------------------------------------

# How do the outgroup f3 results compare to your expectation based on simple
# population relationships (and to the divergence computation above)?

# f3(A, B; C) = E[ (A - C) * (B - C) ]
# This means that in tskit, C is the outgroup (different from ADMIXTOOLS!)

# We can compute f3 for individuals...
ts_f3(ts, B = "AFR_1", C = "EUR_1", A = "CHIMP_1")

# ... but also whole populations (or population samples)
ts_f3(ts, B = sample_sets["AFR"], C = sample_sets["EUR"], A = "CHIMP_1")

ts_f3(ts, B = sample_sets["AFR"], C = sample_sets["NEA"], A = "CHIMP_1")

ts_f3(ts, B = sample_sets["EUR"], C = sample_sets["NEA"], A = "CHIMP_1")





# Bonus 2 -- outgroup f3 as a linear combination of f2 --------------------

# standard f3
ts_f3(ts, B = "AFR_1", C = "AFR_2", A = "CHIMP_1")

# a "homemade" f3 statistic as a linear combination of f2 statistics
# f3(A, B; C) = f2(A, C) + f2(B, C) - f2(A, B) / 2
homemade_f3 <- (
  ts_f2(ts, A = "AFR_1", B = "CHIMP_1")$f2 +
  ts_f2(ts, A = "AFR_2", B = "CHIMP_1")$f2 -
  ts_f2(ts, A = "AFR_1", B = "AFR_2")$f2
) / 2
homemade_f3





# Bonus 3 -- trajectory of Neanderthal ancestry in Europe over time -------

# Extract table with names and times of sampled Europeans (ancient and present day)
eur_inds <- ts_samples(ts) %>% filter(pop == "EUR")
eur_inds

# Compute f4-ration statistic (this will take ~30s) -- note that we can provide
# a vector of names for the X sample set to the `ts_f4ratio()` function
nea_ancestry <- ts_f4ratio(ts, X = eur_inds$name, A = "NEA_1", B = "NEA_2", C = "AFR_1", O = "CHIMP_1")

# Add the age of each sample to the table of proportions
nea_ancestry$time <- eur_inds$time
nea_ancestry

nea_ancestry %>%
  ggplot(aes(time, alpha)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = 2, color = "red", linewidth = 0.5) +
  xlim(40000, 0) +
  coord_cartesian(ylim = c(0, 0.1)) +
  labs(x = "time [years ago]", y = "Neanderthal ancestry proportion") +
  theme_bw() +
  ggtitle("Neanderthal ancestry proportion in Europeans over time")

# Let's test the significance of the decline over time using a linear model
summary(lm(alpha ~ time, data = nea_ancestry))

# Does this match observation in real data?
# See figure 1 in this paper: https://www.pnas.org/doi/full/10.1073/pnas.1814338116,
# and figure 2 in this paper: https://www.nature.com/articles/nature17993





# Bonus 4 -- how many unique f4 quartets are there? -----------------------

# # install a combinatorics R package
# install.packages("combinat")

library(combinat)

# These are the four samples we can create quartet combinations from
quartet <- c("AFR_1", "EUR_1", "NEA_1", "CHIMP_1")
quartets <- permn(quartet)
quartets

# How many permutations there are in total?
#   4! = 4 * 3 * 2 * 1 = 24
factorial(4)

# We should therefore have 24 different quartet combinations of samples
length(quartets)

# Loop across all quartets, computing the corresponding f4 statistic (we want
# to do this using branch lengths, not mutations, as the mutation-based computation
# would involve statistical noise)
all_f4s <- lapply(quartets, function(q) ts_f4(ts, q[1], q[2], q[3], q[4], mode = "branch"))

# Bind the list of f4 results into a single data frame and inspect the results
all_f4s <- bind_rows(all_f4s) %>% arrange(abs(f4))
print(all_f4s, n = Inf)

# Narrow down the results to only unique f4 values
distinct(all_f4s, f4, .keep_all = TRUE)
distinct(all_f4s, abs(f4), .keep_all = TRUE)
