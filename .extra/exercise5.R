devtools::load_all("~/Projects/slendr")
init_env()

library(cowplot)

source("utils.R")

model <- landscape_model(rate = 0.3, Ne = 10000)

pdf("exercise5.pdf", width = 8, height = 5)
for (n in c(50, 25, 10, 5, 2, 1)) {
  # n <- 10
  schedule <- landscape_sampling(model, n)

  ts <- msprime(model, samples = schedule, sequence_length = 20e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

  samples <- ts_samples(ts)

  ts_eigenstrat(ts, "exercise5")

  plot_grid(
    plot_map(model, gene_flow = TRUE, labels = TRUE),
    plot_pca("exercise5", samples, color = "pop") + ggtitle(paste("n = ", n))
  ) %>% print()
}
dev.off()

# per-population Ne values ------------------------------------------------

Ne <- list(
  p1 = 10000,
  p2 = 10000,
  p3 = 10000,
  p4 = 10000,
  p5 = 100,
  p6 = 10000,
  p7 = 10000,
  p8 = 10000,
  p9 = 10000,
  p10 = 10000
)

model <- landscape_model(rate = 0.3, Ne = Ne)

n <- list(
  p1 = 50,
  p2 = 50,
  p3 = 50,
  p4 = 50,
  p5 = 5,
  p6 = 50,
  p7 = 50,
  p8 = 50,
  p9 = 50,
  p10 = 50
)

schedule <- landscape_sampling(model, n)

ts <- msprime(model, samples = schedule, sequence_length = 20e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

samples <- ts_samples(ts)

ts_eigenstrat(ts, "exercise5")

plot_grid(
  plot_map(model, gene_flow = TRUE, labels = TRUE),
  plot_pca("exercise5", samples, color = "pop")
)
