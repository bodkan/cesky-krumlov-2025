devtools::load_all("~/Projects/slendr")
init_env()

library(cowplot)

source("utils.R")

model <- landscape_model(rate = 0.1, Ne = 10000)

plot_model(model, proportions = TRUE)

plot_map(model, gene_flow = TRUE)

schedule <- landscape_sampling(model, n = 25)

plot_model(model, proportions = TRUE, samples = schedule)

ts <- msprime(model, samples = schedule, sequence_length = 50e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

samples <- ts_samples(ts)

ts_eigenstrat(ts, "exercise5")

plot_grid(
  plot_map(model, gene_flow = TRUE, labels = TRUE),
  plot_pca("exercise5", samples, color = "pop")
)
