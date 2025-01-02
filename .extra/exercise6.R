library(slendr)
init_env()

library(cowplot)

source("utils.R")

# (the following model is copied from exerciser4.R) ------------------------

popZ <- population("popZ", time = 3000, N = 5000)
popX <- population("popX", time = 1500, N = 5000, parent = popZ)
popY <- population("popY", time = 1500, N = 5000, parent = popZ)

gf <- gene_flow(from = popX, to = popY, rate = 0.2, start = 800, end = 799)

model <- compile_model(list(popZ, popX, popY), generation_time = 1, gene_flow = gf)

schedule <- rbind(
  schedule_sampling(model, times = seq(1500, 0, by = -100), list(popY, 50)),
  schedule_sampling(model, times = 0, list(popX, 10), list(popZ, 10))
)

# Use the function plot_model() to make sure that the model and the sampling schedule
# are defined correctly (there's no such thing as too many sanity checks when doing research)
plot_model(model, proportions = TRUE, samples = schedule)

ts <- msprime(model, samples = schedule, sequence_length = 100e6, recombination_rate = 1e-8)

# solution to exercise #5 starts here -------------------------------------

library(dplyr)
library(ggplot2)

tracts <- ts_tracts(ts, census = 800)
tracts

# Select sampled individuals from different ages (remember, we recorded oldest individuals
# at 1000 generations ago, the youngest individuals at 0 generations ago). Use the function
# plot_tracts to visualize their ancestry tracts. What do you see in terms of the number of
# tracts and the length of tracts across these individuals? Can you eyeball some distinctive
# pattern?

subset_inds <- tracts %>% group_by(time) %>% distinct(name) %>% slice_sample(n = 1) %>% pull(name)
subset_inds

plot_tracts(tracts, subset_inds)

# It looks like the older individual is, the closer they lived to the start of the admixture event,
# and the longer the tracts they carry will be. Compute the average length of a ancestry tracts in each
# sample age group and visualize the length distribution of these tracts based on their age:

tracts %>%
  group_by(time) %>%
  summarise(mean(length))

ggplot(tracts, aes(length, color = factor(time))) +
  geom_density() +
  coord_cartesian(xlim = c(0, 3e6)) +
  theme_bw()


# Now, try to work backwards. Assuming you have the following distribution of tract lengths...

tracts <- filter(tracts, time == 0, length <= 1e6)
bins <- hist(tracts$length, breaks = 50, plot = FALSE)
length <- bins$mids
density <- bins$density

plot(length, density)

lambda_mle <- 1 / mean(tracts$length)
lambda_mle / 1e-8

y_mle <- dexp(length, rate = lambda_mle)
lines(length, y_mle, lty = 2, col = "darkgreen", lwd = 2)

nls_res <- nls(density ~ SSasymp(length, Asym, R0, lrc))
nls_res

lambda_nls <- exp(unname(coef(nls_res)["lrc"]))
lambda_nls / 1e-8

y_nls <- predict(nls_res, newdata = data.frame(length = length))
lines(length, y_nls, lty = 2, col = "purple", lwd = 2)

legend("topright", fill = c("darkgreen", "purple"),
       legend = c(paste("MLE, t =", round(lambda_mle / 1e-8, 1), "generations ago"),
                  paste("MLE, t =", round(lambda_nls / 1e-8, 1), "generations ago")))
