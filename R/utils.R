#' Plot Marginal Component Effects of a \code{cjbart} Object
#'
#' @description Plots either the observation-level marginal component effect (OMCE) or the individual-level marginal component effect (IMCE). By default, all attribute-levels in the model are plotted.
#' @param x Object of class \code{cjbart}, the result of running [OMCE()]
#' @param covar Character string detailing the covariate over which to analyse heterogeneous effects
#' @param type Character string, either "imce" or "omce", to plot individual- and observation-level effects respectively.
#' @param plot_levels Optional vector of conjoint attribute names to plot. If not supplied, all attributes within the conjoint model will be plotted.
#' @return Plot of marginal component effects.
plot.cjbart <- function(x, covar, type = "imce", plot_levels = NULL) {

  if (type == "imce") {
    plot_data <- x$imce
  } else if (type == "omce") {
    plot_data <- x$omce
  } else {
    stop("type not recognised -- please check you have specified either 'imce' or 'omce'")
  }

  plot_data <- plot_data %>%
    pivot_longer(cols = x$att_levels, names_to = "att", values_to = "mce")

  if (!is.null(plot_levels)) {
    plot_data <- filter(att %in% plot_levels,
                        .data = plot_data)
  }

  if (nrow(plot_data) == 0) {
    stop("Filtering attribute levels renders data with no observations -- please check plot_levels argument.")
  }

  plot_data <- group_by(.data = plot_data,
                        att) %>%
    arrange(mce, by_group = TRUE) %>%
    mutate(x_order = 1:n())

  base_plot <- ggplot(plot_data,
                      aes_string(x = "x_order", y = "mce", color = covar)) +
    facet_wrap(~att, scales = "free") +
    geom_point(alpha = 0.7) +
    ylab(ifelse(type == "imce","IMCE","OMCE")) +
    xlab(ifelse(type == "imce","Individual","Observation")) +
    labs(color = str_to_sentence(covar)) +
    theme(legend.position = "bottom",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  if (typeof(plot_data[[covar]]) == "double") {

    final_plot <- base_plot +
      scale_color_gradient(low = "dodgerblue3", high = "goldenrod1")

  } else {

    final_plot <- base_plot +
      scale_colour_manual(values=c("#000000", # Colour-blind friendly pallete
                                   "#E69F00",
                                   "#56B4E9",
                                   "#009E73",
                                   "#F0E442",
                                   "#0072B2",
                                   "#D55E00",
                                   "#CC79A7"))
  }

  return(final_plot)
}

#' Summarizing \code{cjbart} Marginal Component Effect Estimates
#'
#' @description \code{summary} method for class "cjbart"
#' @param x Object of class \code{cjbart}, the result of running [OMCE()]
#' @return Tibble summarising the average marginal component effect, the minimum and maximum values, and standard deviations for each attribute-level.
summary.cjbart <- function(x) {

  # IMCE summary
  IMCE_only <- x$imce %>% select(all_of(x$att_levels))

  AMCE <- colMeans(IMCE_only)
  mins <- apply(IMCE_only,2,min)
  maxs <- apply(IMCE_only,2,max)
  sds <- apply(IMCE_only,2,sd)

  att_names <- x$att_levels %>%
    ifelse(nchar(.) > 30, paste0(substr(.,1,30),"..."), .)

  summary_tab <- tibble(Level = att_names,
                        AMCE = AMCE,
                        `Min.` = mins,
                        `Max.` = maxs,
                        `Std.Dev` = sds,
                        row.names = NULL)

  message("Summary table of individual marginal component effects (IMCEs)")
  print(summary_tab)

}
