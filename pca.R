devtools::load_all("~/Projects/slendr")
init_env()

source("utils.R")

library(smartsnp)

library(dplyr)
library(scales)
library(cowplot)
library(ggplot2)
library(viridis)

popZ <- population("popZ", time = 3000, N = 5000)
popX <- population("popX", time = 1500, N = 5000, parent = popZ)
popY <- population("popY", time = 1500, N = 5000, parent = popZ)

gf <- gene_flow(from = popX, to = popY, rate = 0.5, start = 200, end = 0)
model <- compile_model(list(popZ, popX, popY), generation_time = 1, gene_flow = NULL)

schedule <- schedule_sampling(model, times = seq(3000, 0, by = -200), list(popZ, 10), list(popX, 10), list(popY, 10))

# Use the function plot_model() to make sure that the model and the sampling schedule
# are defined correctly (there's no such thing as too many sanity checks when doing research)
plot_model(model, proportions = TRUE, samples = schedule)

ts <- msprime(model, samples = schedule, sequence_length = 50e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

# Verify that all the samples we've scheduled for recording really are in the tree sequence
ts_samples(ts)
ts_samples(ts) %>% filter(pop == "popX") %>% group_by(pop, time) %>% tally
ts_samples(ts) %>% filter(pop == "popY") %>% group_by(pop, time) %>% tally
ts_samples(ts) %>% filter(pop == "popZ") %>% group_by(pop, time) %>% tally %>% print(n = Inf)

# Save the samples table to a data frame (and process it a bit for tidier plotting later)
samples <- ts_samples(ts) %>% mutate(pop = factor(pop, levels = c("popZ", "popX", "popY")))
samples

# Our goal is to investigate how is PCA structure affected by different factors related
# to sampling (sample sizes, dates of samples, etc.). Right now, our tree sequence contains
# every single sample. Let's create smaller subsets of the data using process called
# "simplification".
samples_XYpresent <- samples %>% filter(pop %in% c("popX", "popY"), time == 0)
samples_XYpresent %>% group_by(pop, time) %>% tally
ts_simplify(ts, simplify_to = samples_XYpresent$name) %>% ts_eigenstrat("model_XYpresent")

samples_XYall <- samples %>% filter(pop %in% c("popX", "popY"))
samples_XYall %>% group_by(pop, time) %>% tally %>% print(n = Inf)
ts_simplify(ts, simplify_to = samples_XYall$name) %>% ts_eigenstrat("model_XYall")

samples_XYZpresent <- samples %>% filter(time == 0)
ts_simplify(ts, simplify_to = samples_XYZpresent$name) %>% ts_eigenstrat("model_XYZpresent")

samples_XYZall <- samples
ts_eigenstrat(ts, "model_XYZall")

samples_XYall_Zancient <- samples %>% filter(pop %in% c("popX", "popY") | (pop == "popZ" & time >= 1500))
ts_eigenstrat(ts, "model_XYall_Zancient")

p_model <- plot_model(model, proportions = TRUE, samples = schedule)
p_model

plot_grid(p_model, plot_pca("model_XYpresent", samples_XYpresent, color = "pop"))
plot_grid(p_model, plot_pca("model_XYall", samples_XYall, color = "time"))

plot_grid(p_model, plot_pca("model_XYZpresent", samples_XYZpresent, color = "pop"))

plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(1, 2)))
plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(1, 2), color = "pop"))
plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(2, 3)))
plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(2, 3), color = "pop"))
plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(3, 4)))
plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(3, 4), color = "pop"))
plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(4, 5)))
plot_grid(p_model, plot_pca("model_XYZall", samples_XYZall, pc = c(4, 5), color = "pop"))

# when the whole X and Y trajectories are in but only ancient Z lineage, X and Y are
# not actually separated from one another (only based on time)
# -- also note that this series seems to kind of replicate the previous series of
#    figures, except that X and Y are always mingled
plot_grid(p_model, plot_pca("model_XYall_Zancient", samples_XYall_Zancient, pc = c(1, 2), color = "time"))
plot_grid(p_model, plot_pca("model_XYall_Zancient", samples_XYall_Zancient, pc = c(1, 2), color = "pop"))
plot_grid(p_model, plot_pca("model_XYall_Zancient", samples_XYall_Zancient, pc = c(2, 3), color = "pop"))
plot_grid(p_model, plot_pca("model_XYall_Zancient", samples_XYall_Zancient, pc = c(3, 4), color = "pop"))
plot_grid(p_model, plot_pca("model_XYall_Zancient", samples_XYall_Zancient, pc = c(4, 5), color = "pop"))





