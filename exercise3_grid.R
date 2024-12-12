#
# This is a bonus advanced version. Ignore this unless you're an R expert. :)
#
# Otherwise, look at the solution in exercise3_simple.R!
#

library(slendr)
init_env()

simulate_afs <- function(Ne) {
  # create a slendr model with a given N size
  pop <- population("pop", N = Ne, time = 1)
  model <- compile_model(pop, generation_time = 1, simulation_length = 100000)

  # simulate a tree sequence
  ts <-
    msprime(model, sequence_length = 10e6, recombination_rate = 1e-8) %>%
    ts_mutate(mutation_rate = 1e-8)

  # get a random sample of names of 10 individuals
  samples <- ts_names(ts) %>% sample(10)

  # compute the AFS vector (dropping the '0-th' element added by tskit)
  afs <- ts_afs(ts, list(samples))[-1]

  afs
}

afs_observed <- c(2520, 1449, 855, 622, 530, 446, 365, 334, 349, 244,
                  264, 218,  133, 173, 159, 142, 167, 129, 125, 143)

# generate regularly spaced values of potential Ne values (our parameter grid)
Ne_grid <- seq(from = 1000, to = 30000, by = 1000)
Ne_grid

library(parallel)

# compute AFS (in parallel) across the entire grid of possible Ne values
afs_grid <- mclapply(Ne_grid, simulate_afs, mc.cores = detectCores())
names(afs_grid) <- Ne_grid

afs_grid

# plot the observed AFS...
plot(afs_observed, type = "b", col = "black", lwd = 3, xlab = "allele count bin", ylab = "count")
# ... and overlay the simulated AFS vectors on top of it
for (i in seq_along(Ne_grid)) {
  lines(afs_grid[[i]], lwd = 0.5)
}
legend("topright", legend = c("observed AFS", "simulated AFS"), fill = c("black", "gray"))

# compute mean-squared error of the AFS produced by each Ne value across the grid
errors <- sapply(afs_grid, function(sim_afs) {
  sum((sim_afs - afs_observed)^2) / length(sim_afs)
})

plot(Ne_grid, errors, ylab = "error")
abline(v = Ne_grid[which.min(errors)], col = "red")
legend("topright", legend = paste("minimum error Ne =", Ne_grid[which.min(errors)]), fill = "red")

# plot the AFS again, but this time highlight the most likely spectrum
# (i.e. the one which gave the lowest RMSE value)
plot(afs_observed, type = "b", col = "black", lwd = 3, xlab = "allele count bin", ylab = "count")
for (i in seq_along(Ne_grid)) {
  color <- if (i == which.min(errors)) "red" else "gray"
  width <- if (i == which.min(errors)) 2 else 0.75
  lines(afs_grid[[i]], lwd = width, col = color)
}
legend("topright", legend = c("observed AFS", paste("best fitting Ne =", Ne_grid[which.min(errors)])),
       fill = c("black", "red"))

#
# !!!!! SPOILER ALERT !!!!!
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# true Ne was 6543!