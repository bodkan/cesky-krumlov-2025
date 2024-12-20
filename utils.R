# PCA ---------------------------------------------------------------------

library(smartsnp)
library(dplyr)
library(ggplot2)
library(ggrepel)

plot_pca <- function(prefix, samples, pc = c(1, 2), color = c("time", "pop"), return = c("plot", "pca", "both")) {
  if (length(pc) != 2)
    stop("The 'pc' argument of 'plot_pca' must be an integer vector of length two", call. = FALSE)

  return <- match.arg(return)
  color <- match.arg(color)

  suppressMessages(pca <- smart_pca(snp_data = paste0(prefix, ".geno"), program_svd = "bootSVD", sample_group = samples$pop))

  pc_cols <- paste0("PC", pc)
  pca_df <- pca$pca.sample_coordinates[, pc_cols]
  pca_df$pop <- factor(samples$pop, levels = unique(samples$pop))
  pca_df$time <- samples$time

  variance_explained <- pca$pca.eigenvalues[2, ] %>% {. / sum(.) * 100} %>% round(1)

  pop_df <- group_by(pca_df, pop, time) %>% summarise_all(mean)

  if (color == "time") {
    gg_point <- geom_point(aes(x = !!dplyr::sym(pc_cols[1]), y = !!dplyr::sym(pc_cols[2]), shape = pop, color = time))
    gg_label <- geom_label_repel(data = pop_df, aes(label = pop, x = !!dplyr::sym(pc_cols[1]), y = !!dplyr::sym(pc_cols[2]),
                                                    shape = pop, color = time), show.legend = FALSE)
    gg_theme <- scale_color_viridis_c(option = "viridis")
  } else {
    gg_label <- geom_label_repel(data = pop_df, aes(label = pop, x = !!dplyr::sym(pc_cols[1]), y = !!dplyr::sym(pc_cols[2]),
                                                    color = pop), show.legend = FALSE)
    if (length(unique(samples$pop)) > 6) {
      gg_point <- geom_point(aes(x = !!dplyr::sym(pc_cols[1]), y = !!dplyr::sym(pc_cols[2]), color = pop))
    } else {
      gg_point <- geom_point(aes(x = !!dplyr::sym(pc_cols[1]), y = !!dplyr::sym(pc_cols[2]), shape = pop, color = pop))
    }
    gg_theme <- scale_color_discrete(drop = FALSE)
  }

  plot <- ggplot(pca_df) +
    gg_point +
    gg_label +
    scale_shape_discrete(drop = FALSE) +
    labs(x = sprintf("%s [%.1f %%]", pc_cols[1], variance_explained[pc[1]]),
         y = sprintf("%s [%.1f %%]", pc_cols[2], variance_explained[pc[2]])) +
    theme_bw() +
    gg_theme

  if (return == "plot")
    return(plot)
  else if (return == "pca")
    return(pca)
  else
    return(list(plot = plot, pca = pca))
}

landscape_model <- function(rate, Ne) {
  start_pops <- 5000
  start_gf <- 8000
  simulation_length <- 20000

  xrange <- c(-90, -20)
  yrange <- c(-58, 15)

  map <- world(xrange = xrange, yrange = yrange, crs = "EPSG:31970")

  # non-spatial ancestral population
  p_anc <- population("p_ancestor", N = Ne, time = 1, remove = start_pops + 1)

  # spatial populations
  p1 <- population("p1", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-75, 0), radius = 200e3)
  p2 <- population("p2", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-60, 5), radius = 200e3)
  p3 <- population("p3", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-65, -5), radius = 200e3)
  p4 <- population("p4", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-60, -20), radius = 200e3)
  p5 <- population("p5", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-65, -35), radius = 200e3)
  p6 <- population("p6", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-69, -42), radius = 200e3)
  p7 <- population("p7", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-51, -10), radius = 200e3)
  p8 <- population("p8", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-45, -15), radius = 200e3)
  p9 <- population("p9", N = Ne, time = start_pops, parent = p_anc, map = map, center = c(-71, -12), radius = 200e3)

  gf <- list(
    gene_flow(p1, p2, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p2, p1, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p1, p3, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p3, p1, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p2, p7, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p7, p2, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p7, p8, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p8, p7, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p4, p7, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p7, p4, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p4, p5, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p5, p4, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p5, p6, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p6, p5, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p1, p9, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p9, p1, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p4, p9, rate, start = start_gf, end = simulation_length, overlap = FALSE),
    gene_flow(p9, p4, rate, start = start_gf, end = simulation_length, overlap = FALSE)
  )

  suppressWarnings(model <- compile_model(
    populations = list(p_anc, p1, p2, p3, p4, p5, p6, p7, p8, p9), gene_flow = gf,
    generation_time = 1, simulation_length = simulation_length,
    serialize = FALSE
  ))

  return(model)
}

landscape_sampling <- function(model, n) {
  schedule <- schedule_sampling(model, times = model$orig_length,
                                list(model$populations$p1, n),
                                list(model$populations$p2, n),
                                list(model$populations$p3, n),
                                list(model$populations$p4, n),
                                list(model$populations$p5, n),
                                list(model$populations$p6, n),
                                list(model$populations$p7, n),
                                list(model$populations$p8, n),
                                list(model$populations$p9, n))
  schedule
}

# ancestry tracts ---------------------------------------------------------

plot_tracts <- function(tracts, ind) {
  ind_tracts <- filter(tracts, name %in% ind) %>%
    mutate(haplotype = paste0(name, "\nhap. ", haplotype))
  ind_tracts$haplotype <- factor(ind_tracts$haplotype, levels = unique(ind_tracts$haplotype[order(ind_tracts$node_id)]))

  ggplot(ind_tracts) +
    geom_rect(aes(xmin = left, xmax = right, ymin = 1, ymax = 2, fill = name), linewidth = 1) +
    labs(x = "coordinate along a chromosome [bp]") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    ) +
    facet_grid(haplotype ~ .) +
    expand_limits(x = 0) +
    scale_x_continuous(labels = scales::comma)
}
