library(slendr)
init_env()

# chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 5000)

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

# verify the model visually
plot_model(model)
plot_model(model, sizes = FALSE)
plot_model(model, sizes = FALSE, log = TRUE)
plot_model(model, sizes = FALSE, log = TRUE, proportions = TRUE)

# schedule sampling of a defined set of individuals
nea_samples <- schedule_sampling(model, times = c(70000, 40000), list(nea, 1))
present_samples <- schedule_sampling(model, times = 0, list(chimp, 1), list(afr, 5), list(eur, 10))
emh_samples <- schedule_sampling(model, times = runif(n = 40, min = 10000, max = 40000), list(eur, 1))

schedule <- rbind(nea_samples, present_samples, emh_samples)
schedule

plot_model(model, sizes = FALSE, log = TRUE, samples = schedule)
