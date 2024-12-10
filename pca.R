devtools::load_all("~/Projects/slendr")
init_env()

library(smartsnp)

library(dplyr)
library(scales)
library(cowplot)
library(ggplot2)
library(viridis)

popAnc <- population("popAnc", time = 2000, N = 5000, remove = 1499)
popX <- population("popX", time = 1500, N = 5000, parent = popAnc)
popY <- population("popY", time = 1500, N = 5000, parent = popAnc)

model <- compile_model(list(popAnc, popX, popY), generation_time = 1)

schedule <- schedule_sampling(model, times = seq(2000, 0, by = -100), list(popAnc, 10), list(popX, 10), list(popY, 10))

ts <- msprime(model, samples = schedule, sequence_length = 50e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

samples_XY_present <- ts_samples(ts) %>% filter(pop %in% c("popX", "popY"), time == 0)
ts_simplify(ts, simplify_to = samples_XY_present$name) %>% ts_eigenstrat("model_XY_present")

samples_XY_all <- ts_samples(ts) %>% filter(pop %in% c("popX", "popY"))
ts_simplify(ts, simplify_to = samples_XY_all$name) %>% ts_eigenstrat("model_XY_all")

samples_AncXY_all <- ts_samples(ts)
ts_eigenstrat(ts, "model_AncXY_all")

plot_grid(plot_model(model), plot_pca("model_XY_present", samples_XY_present))
plot_grid(plot_model(model), plot_pca("model_XY_all", samples_XY_all))
plot_grid(plot_model(model), plot_pca("model_AncXY_all", samples_XY_all))



plot_pca <- function(prefix, samples, pc = c(1, 2), color = c("time", "pop")) {
  if (length(pc) != 2)
    stop("The 'pc' argument of 'plot_pca' must be an integer vector of length two", call. = FALSE)

  color <- match.arg(color)

  suppressMessages(pca <- smart_pca(snp_data = paste0(prefix, ".geno"), program_svd = "bootSVD", sample_group = samples$pop))

  pc_cols <- paste0("PC", pc)
  pca_df <- pca$pca.sample_coordinates[, pc_cols]
  pca_df$pop <- samples$pop
  pca_df$time <- samples$time

  variance_explained <- pca$pca.eigenvalues[2, ] %>% {. / sum(.) * 100} %>% round(1)

  if (color == "time") {
    gg_point <- geom_point(aes(x = !!sym(pc_cols[1]), y = !!sym(pc_cols[2]), shape = pop, color = time))
    gg_theme <- scale_color_viridis_c(option = "viridis")
  } else {
    gg_point <- geom_point(aes(x = !!sym(pc_cols[1]), y = !!sym(pc_cols[2]), shape = pop, color = pop))
    gg_theme <- scale_color_discrete(drop = FALSE)
  }

  ggplot(pca_df) +
    gg_point +
    scale_shape_discrete(drop = FALSE) +
    labs(x = sprintf("%s [%.1f %%]", pc_cols[1], variance_explained[pc[1]]),
         y = sprintf("%s [%.1f %%]", pc_cols[2], variance_explained[pc[2]])) +
    theme_bw() +
    gg_theme
}




# gene-flow model ---------------------------------------------------------

devtools::load_all("~/Projects/slendr")
init_env()

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

schedule <- schedule_sampling(model, times = seq(2000, 0, by = -100), list(popZ, 10), list(popX, 10), list(popY, 10))

ts <- msprime(model, samples = schedule, sequence_length = 50e6, recombination_rate = 1e-8) %>% ts_mutate(1e-8)

samples <- ts_samples(ts) %>% mutate(pop = factor(pop, levels = c("popZ", "popX", "popY")))

samples_XYpresent <- samples %>% filter(pop %in% c("popX", "popY"), time == 0)
ts_simplify(ts, simplify_to = samples_XYpresent$name) %>% ts_eigenstrat("model_XYpresent")

samples_XYall <- samples %>% filter(pop %in% c("popX", "popY"))
ts_simplify(ts, simplify_to = samples_XYall$name) %>% ts_eigenstrat("model_XYall")

samples_XYZpresent <- samples %>% filter(time == 0)
ts_simplify(ts, simplify_to = samples_XYZpresent$name) %>% ts_eigenstrat("model_XYZpresent")

samples_XYZall <- samples
ts_eigenstrat(ts, "model_XYZall")

samples_XYall_Zancient <- samples %>% filter(pop %in% c("popX", "popY") | (pop == "popZ" & time >= 1500))
ts_eigenstrat(ts, "model_XYall_Zancient")

p_model <- plot_model(model, proportions = TRUE, samples = schedule)

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





