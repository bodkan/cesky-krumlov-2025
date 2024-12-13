# PCA ---------------------------------------------------------------------

library(smartsnp)

plot_pca <- function(prefix, samples, pc = c(1, 2), color = c("time", "pop"), return = c("plot", "pca", "both")) {
  if (length(pc) != 2)
    stop("The 'pc' argument of 'plot_pca' must be an integer vector of length two", call. = FALSE)

  return <- match.arg(return)
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

  plot <- ggplot(pca_df) +
    gg_point +
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
