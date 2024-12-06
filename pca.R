devtools::load_all("~/Projects/slendr")
init_env()

library(smartsnp)

library(dplyr)
library(scales)
library(cowplot)
library(ggplot2)
library(viridis)

popA <- population("popA", time = 50000, N = 10000)
popB <- population("popB", time = 30000, N = 10000, parent = popA)
popC <- population("popC", time = 20000, N = 10000, parent = popB)
popD <- population("popD", time = 10000, N = 10000, parent = popC)
popE <- population("popE", time = 5000, N = 10000, parent = popD)
popF <- population("popF", time = 2000, N = 10000, parent = popE)

plot_pca <- function(prefix, ts, pc = c(1, 2)) {
  if (length(pc) != 2)
    stop("The 'pc' argument of 'plot_pca' must be an integer vector of length two", call. = FALSE)

  groups <- ts_samples(ts)$pop
  suppressMessages(pca <- smart_pca(snp_data = paste0(prefix, ".geno"), sample_group = groups))

  pc_cols <- paste0("PC", pc)
  pca_df <- pca$pca.sample_coordinates[, pc_cols]
  pca_df$pop <- groups
  pca_df$time <- ts_samples(ts)$time

  variance_explained <- pca$pca.eigenvalues[2, ] %>% {. / sum(.) * 100} %>% round(1)

  ggplot(pca_df) +
    geom_point(aes(x = !!sym(pc_cols[1]), y = !!sym(pc_cols[2]), shape = pop, color = time)) +
    labs(x = sprintf("%s [%.1f %%]", pc_cols[1], variance_explained[pc[1]]),
         y = sprintf("%s [%.1f %%]", pc_cols[2], variance_explained[pc[2]])) +
    theme_bw() +
    scale_color_viridis_c()
}

# model without gene flow -------------------------------------------------

model_nogf <- compile_model(list(popA, popB, popC, popD, popE, popF), generation_time = 30)

plot_model(model_nogf)

# gf <- gene_flow(from = popA, to = popF, rate = 0.2, start = 1900, end = 2000)
# model_nogf <- compile_model(list(popA, popB, popC, popD, popE, popF), generation_time = 1, simulation_length = 2000, gene_flow = gf)

samples_nogf <- schedule_sampling(model_nogf, times = seq(50000, 0, by = -1000), list(popA, 2), list(popB, 2), list(popC, 2), list(popD, 2), list(popE, 2))
# samples_nogf <- schedule_sampling(model_nogf, times = 2000, list(popA, 25), list(popB, 25), list(popC, 25), list(popD, 25), list(popE, 25), list(popF, 25))

ts_nogf <- msprime(model_nogf, samples = samples_nogf, sequence_length = 10e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

ts_eigenstrat(ts_nogf, "pops_ABCDEF_nogf", chrom = "chr1", outgroup = "outgroup")

# comparison --------------------------------------------------------------

plot_grid(plot_model(model_nogf, proportions = TRUE), plot_pca("pops_ABCDEF_nogf", ts_nogf), nrow = 1)












# simple model ------------------------------------------------------------

popAnc <- population("popAnc", time = 2000, N = 10000, remove = 1499)
popX <- population("popX", time = 1000, N = 10000, parent = popAnc)
popY <- population("popY", time = 1000, N = 10000, parent = popAnc)

model_nogf <- compile_model(list(popAnc, popX, popY), generation_time = 1)

# gf <- gene_flow(from = popA, to = popF, rate = 0.2, start = 1900, end = 2000)
# model_nogf <- compile_model(list(popA, popB, popC, popD, popE, popF), generation_time = 1, simulation_length = 2000, gene_flow = gf)

samples_nogf <- schedule_sampling(model_nogf, times = seq(2000, 0, by = -100), list(popAnc, 10), list(popX, 10), list(popY, 10))
# samples_nogf <- schedule_sampling(model_nogf, times = 2000, list(popA, 25), list(popB, 25), list(popC, 25), list(popD, 25), list(popE, 25), list(popF, 25))

ts_nogf <- msprime(model_nogf, samples = samples_nogf, sequence_length = 10e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

ts_eigenstrat(ts_nogf, "pops_AncXY_nogf", chrom = "chr1", outgroup = "outgroup")

# comparison --------------------------------------------------------------

plot_grid(plot_model(model_nogf, proportions = TRUE), plot_pca("pops_AncXY_nogf", ts_nogf), nrow = 1)
