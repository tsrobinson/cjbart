#' Plot Marginal Component Effects of a \code{cjbart} Object
#' @description Plots observation-level or individual-level marginal component effects (OMCE and IMCE respectively). By default, all attribute-levels in the model are plotted.
#' @param x Object of class \code{cjbart}, the result of running [OMCE()]
#' @param covar Character string detailing the covariate over which to analyse heterogeneous effects
#' @param type Character string, either "imce" or "omce", to plot individual- and observation-level effects respectively.
#' @param plot_levels Optional vector of conjoint attribute names to plot. If not supplied, all attributes within the conjoint model will be plotted.
#' @param ... Additional arguments for plotting the marginal component effects (see below).
#' @return Plot of marginal component effects.
#' @method plot cjbart
#' @export
plot.cjbart <- function(x, plot_levels = NULL, type = "imce", covar, ...) {

  if (type == "imce") {
    plot_data <- x$imce
  } else if (type == "omce") {
    plot_data <- x$omce
  } else {
    stop("type not recognised -- please check you have specified either 'imce' or 'omce'")
  }

  plot_data <- plot_data %>%
    tidyr::pivot_longer(cols = x$att_levels, names_to = "att", values_to = "mce")

  if (!is.null(plot_levels)) {
    plot_data <- dplyr::filter(rlang::.data$att %in% plot_levels,
                               .data = plot_data)
  }

  if (nrow(plot_data) == 0) {
    stop("Filtering attribute levels renders data with no observations -- please check plot_levels argument.")
  }

  plot_data <- dplyr::group_by(.data = plot_data,
                               rlang::.data$att) %>%
    dplyr::arrange(rlang::.data$mce, by_group = TRUE) %>%
    dplyr::mutate(x_order = 1:dplyr::n())

  base_plot <- ggplot2::ggplot(plot_data,
                               ggplot2::aes_string(x = "x_order",
                                                   y = "mce",
                                                   color = covar)
                               ) +
    ggplot2::facet_wrap(~rlang::.data$att, scales = "free") +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::ylab(ifelse(type == "imce","IMCE","OMCE")) +
    ggplot2::xlab(ifelse(type == "imce","Individual","Observation")) +
    ggplot2::labs(color = stringr::str_to_sentence(covar)) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())

  if (typeof(plot_data[[covar]]) == "double") {

    final_plot <- base_plot +
      ggplot2::scale_color_gradient(low = "dodgerblue3", high = "goldenrod1")

  } else {

    final_plot <- base_plot +
      ggplot2::scale_colour_manual(values=c("#000000", # Colour-blind friendly pallete
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
#' @description \code{summary} method for class "cjbart"
#' @param object Object of class \code{cjbart}, the result of running [OMCE()]
#' @param ... Further arguments (not currently used)
#' @return Tibble summarising the average marginal component effect, the minimum and maximum values, and standard deviations for each attribute-level.
#' @method summary cjbart
#' @export
summary.cjbart <- function(object, ...) {

  # IMCE summary
  IMCE_only <- object$imce %>% dplyr::select(tidyselect::all_of(object$att_levels))

  AMCE <- colMeans(IMCE_only)
  mins <- apply(IMCE_only,2,min)
  maxs <- apply(IMCE_only,2,max)
  sds <- apply(IMCE_only,2,stats::sd)

  att_names <- object$att_levels %>%
    ifelse(test = nchar(.) > 30,
           yes = paste0(substr(.,1,30),"..."),
           no = .)

  summary_tab <- dplyr::tibble(Level = att_names,
                        AMCE = AMCE,
                        `Min.` = mins,
                        `Max.` = maxs,
                        `Std.Dev` = sds,
                        row.names = NULL)

  message("Summary table of individual marginal component effects (IMCEs)")
  print(summary_tab)

}
