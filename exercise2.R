library(slendr)
init_env()

library(dplyr)
library(ggplot2)

# (the following model is copied from exercise1.R) ------------------------

# chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 15000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)

# Neanderthal population splitting at 600 ky ago from modern humans
# (becomes extinct by 40 ky ago)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# Neanderthal introgression event (3% admixture between 55-50 kya)
gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 50000)

# compile the entire model into a single object
model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,
  generation_time = 30
)

plot_model(model, proportions = TRUE)
plot_model(model, proportions = TRUE, log = TRUE)

# solution to exercise #2 starts here -------------------------------------

# if the command below takes way too much time to run, feel free to go down from
# 100 Mb sequence_length to even 10Mb (it doesn't matter much)
tstart <- Sys.time()
ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)
tend <- Sys.time()
tend - tstart

# define sampling events
nea_samples <- schedule_sampling(model, times = c(70000, 40000), list(nea, 1))
present_samples <- schedule_sampling(model, times = 0, list(chimp, 1), list(afr, 10), list(eur, 10))
emh_samples <- schedule_sampling(model, times = seq(50000, 2000, by = -2000), list(eur, 1))

schedule <- rbind(nea_samples, present_samples, emh_samples)

tstart <- Sys.time()
ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8, samples = schedule) %>%
  ts_mutate(mutation_rate = 1e-8)
tend <- Sys.time()
tend - tstart

# ts_write(ts, "exercise2.trees")
# ts <- ts_read("exercise2.trees", model)

# make sure that the tree sequence contains exactly the samples that we've scheduled
# to be recorded in the output

samples <- ts_samples(ts)

samples %>% group_by(pop, present_day = time == 0) %>% tally

# Exercise #4 ------------------------------------------------------------------

# first get a (named!) list of individuals in each population

sample_sets <- ts_names(ts, split = "pop")
sample_sets

# compute nucleotide diversity (pi) in each population

pi_pop <- ts_diversity(ts, sample_sets)

arrange(pi_pop, diversity)

# heterozygosity in a single sample
pi_ind <- ts_diversity(ts, samples$name)
pi_ind <- inner_join(samples, pi_ind, by = join_by(name == set))

ggplot(pi_ind, aes(pop, diversity, color = pop, group = pop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()


# compute pairwise population divergence

div <- ts_divergence(ts, sample_sets)

arrange(div, divergence)


# outgroup f3

ts_f3(ts, B = "AFR_1", C = "EUR_1", A = "CHIMP_1")

ts_f3(ts, B = sample_sets["AFR"], C = sample_sets["EUR"], A = "CHIMP_1")

ts_f3(ts, B = sample_sets["AFR"], C = sample_sets["NEA"], A = "CHIMP_1")

ts_f3(ts, B = sample_sets["EUR"], C = sample_sets["NEA"], A = "CHIMP_1")


# outgroup f3 as a linear combination of f2 statistics

# branch-based version
ts_f3(ts, B = "AFR_1", C = "AFR_2", A = "CHIMP_1", mode = "branch")

homemade_f3 <- (
  ts_f2(ts, A = "AFR_1", B = "CHIMP_1", mode = "branch")$f2 +
  ts_f2(ts, A = "AFR_2", B = "CHIMP_1", mode = "branch")$f2 -
  ts_f2(ts, A = "AFR_1", B = "AFR_2", mode = "branch")$f2
) / 2
homemade_f3


# admixture detection -----------------------------------------------------

#                            BABA - ABBA
# D(AFR, EUR; NEA, CHIMP) = -------------
#                            BABA + ABBA

ts_f4(ts, W = "AFR_1", X = "AFR_2", Y = "NEA_1", Z = "CHIMP_1")

ts_f4(ts, W = "AFR_1", X = "EUR_1", Y = "NEA_1", Z = "CHIMP_1")

f4_afr <- lapply(sample_sets$AFR, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()
f4_eur <- lapply(sample_sets$EUR, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()

f4_afr$pop <- "AFR"
f4_eur$pop <- "EUR"

rbind(f4_afr, f4_eur) %>%
  ggplot(aes(pop, f4, color = pop)) +
  geom_boxplot() +
  geom_jitter() +
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("f4(AFR, EUR; NEA, CHIMP)")


# time-series

# extract table with names and times of sampled Europeans (ancient and present day)
eur_inds <- ts_samples(ts) %>% filter(pop == "EUR")
eur_inds

# compute f4-ration statistic (this will take ~30s)
nea_ancestry <- ts_f4ratio(ts, X = eur_inds$name, A = "NEA_1", B = "NEA_2", C = "AFR_1", O = "CHIMP_1")
nea_ancestry

# add the computed proportions to the data frame
eur_inds$ancestry <- nea_ancestry$alpha

eur_inds %>%
  ggplot(aes(time, ancestry)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = 2, color = "red", linewidth = 0.5) +
  xlim(40000, 0) +
  coord_cartesian(ylim = c(0, 0.1)) +
  labs(x = "time [years ago]", y = "Neanderthal ancestry proportion")


# unique quartets

# # install a combinatorics R package
# install.packages("combinat")

quartet <- c("AFR_1", "EUR_1", "NEA_1", "CHIMP_1")

quartets <- combinat::permn(quartet)

# how many permutations there are in total?
#   4! = 4 * 3 * 2 * 1 = 24

length(quartets)

# loop across all quartets, computing the corresponding f4 statistic
all_f4s <- lapply(quartets, function(q) ts_f4(ts, q[1], q[2], q[3], q[4], mode = "branch"))

# bind the list of f4 results into a single data frame
all_f4s <- bind_rows(all_f4s) %>% arrange(abs(f4))
print(all_f4s, n = Inf)

distinct(all_f4s, f4, .keep_all = TRUE)
distinct(all_f4s, abs(f4), .keep_all = TRUE)