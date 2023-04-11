#' Plot Marginal Component Effects of a \code{cjbart} Object
#' @description Plots observation-level or individual-level marginal component effects (OMCE and IMCE respectively). By default, all attribute-levels in the model are plotted.
#' @param x Object of class \code{cjbart}, the result of running [cjbart::IMCE()]
#' @param covar Character string detailing the covariate over which to analyze heterogeneous effects
#' @param plot_levels Optional vector of conjoint attribute levels to plot. If not supplied, all attributes within the conjoint model will be plotted.
#' @param se Boolean determining whether to show an estimated 95% confidence interval
#' @param ... Additional arguments for plotting the marginal component effects (see below).
#' @return Plot of marginal component effects.
#' @importFrom rlang .data
#' @method plot cjbart
#' @export
plot.cjbart <- function(x, covar = NULL, plot_levels = NULL, se = TRUE,  ...) {

  exp_data <- x$imce

  plot_data <- tidyr::pivot_longer(cols = x$att_levels,
                                   names_to = "level",
                                   values_to = "imce",
                                   data = exp_data)

  plot_data <- merge(x = plot_data, y = x$att_lookup,
                     all.x = TRUE,
                     by = "level")

  if (se) {

    ci_lower <- tidyr::pivot_longer(cols = x$att_levels,
                                    names_to = "level",
                                    values_to = "imce_lower",
                                    data = x$imce_lower)

    ci_upper <- tidyr::pivot_longer(cols = x$att_levels,
                                    names_to = "level",
                                    values_to = "imce_upper",
                                    data = x$imce_upper)

    plot_data <- merge(x = plot_data,
                       y = ci_lower[,c(x$id,"level","imce_lower")],
                       all.x = TRUE,
                       by = c("level",x$id))

    plot_data <- merge(x = plot_data,
                       y = ci_upper[,c(x$id,"level","imce_upper")],
                       all.x = TRUE,
                       by = c("level",x$id))
  }

  if (!is.null(plot_levels)) {
    plot_data <- plot_data[plot_data$level %in% plot_levels,]

    if (nrow(plot_data) == 0) {
      stop("Failed to find attribute plot_levels in data. This may be because the conjoint has multiple attributes with the same level.\n If this is the case, try specifying the required levels as '<attribute>_level', e.g., '<age>_80+'.")
    }

  }

  plot_data <- do.call(

    rbind,

    by(
      plot_data,
      INDICES = plot_data$level,
      FUN = function(x) {
        x_ordered <- x[order(x$imce),]
        x_ordered$x_order <- 1:nrow(x_ordered)
        return(x_ordered)
        }
     )
  )

  plot_data$facet_labels <- paste0("atop(italic('",plot_data$Attribute,":'),bold('",plot_data$Level,"'))")

  base_plot <- ggplot2::ggplot(plot_data,
                                 ggplot2::aes_string(x = "x_order",
                                                     y = "imce")) +

    ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed") +

    {if (!is.null(covar)) {ggplot2::aes_string(color = covar)}} +
    {if (se) {ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "imce_lower", ymax = "imce_upper"), color = NA, fill = "grey60", alpha = 0.7)}} +

    ggplot2::facet_wrap(~.data$facet_labels, scales = "free", labeller = "label_parsed") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::ylab("IMCE") +
    ggplot2::xlab("Individual") +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())

  if (!is.null(covar)) {

    if (typeof(plot_data[[covar]]) == "double") {

      base_plot <- base_plot +
        ggplot2::scale_color_gradient(low = "#0072B2", high = "#E69F00")

    } else {

      base_plot <- base_plot +
        ggplot2::scale_colour_manual(values=c("#000000", # Colour-blind friendly pallete
                                              "#E69F00",
                                              "#56B4E9",
                                              "#009E73",
                                              "#F0E442",
                                              "#0072B2",
                                              "#D55E00",
                                              "#CC79A7"))
    }

  }


  return(base_plot)
}

#' Summarizing \code{cjbart} Marginal Component Effect Estimates
#' @description \code{summary} method for class "cjbart"
#' @param object Object of class \code{cjbart}, the result of running [cjbart::IMCE()]
#' @param ... Further arguments (not currently used)
#' @return Data frame summarizing the average marginal component effect (AMCE), the minimum and maximum values, and standard deviations for each attribute-level.
#' @note To calculate the AMCE with Bayesian credible intervals, please use the \code{AMCE()} function instead.
#' @method summary cjbart
#' @seealso [cjbart::AMCE()]
#' @export
summary.cjbart <- function(object, ...) {

  # IMCE summary
  IMCE_only <- subset(object$imce, select = object$att_levels)

  AMCE <- colMeans(IMCE_only)
  mins <- apply(IMCE_only,2,min)
  maxs <- apply(IMCE_only,2,max)
  sds <- apply(IMCE_only,2,stats::sd)

  summary_tab <- data.frame(level = object$att_levels,
                            AMCE = AMCE,
                            `Min.` = mins,
                            `Max.` = maxs,
                            `Std.Dev` = sds,
                            row.names = NULL)

  summary_tab <- merge(x = object$att_lookup,
                       y = summary_tab,
                       all.y = TRUE,
                       by = "level")

  summary_tab$level <- NULL

  message("Individual marginal component effects (IMCEs)")
  summary_tab <- tidyr::tibble(summary_tab)
  return(summary_tab)

}

#' Plot Variable Importance Matrix for Heterogeneity Analysis
#' @description Plots a heatmap of variable importance, across predicted IMCEs. By default, all attribute-levels and covariates in the model are plotted.
#' @param x Object of class \code{cjbart}, the result of running [cjbart::IMCE()]
#' @param covars Optional vector of covariate names to plot. By default, all included covariates are shown.
#' @param att_levels Optional vector of attribute-levels to plot. By default, all attribute-levels are shown.
#' @param ... Additional arguments (not currently used)
#' @return Plot of covariate importance scores
#' @method plot cjbart.vimp
#' @export
plot.cjbart.vimp <- function(x, covars = NULL, att_levels = NULL,  ...) {

  plot_data <- x$results

  if (!is.null(covars)) {

    plot_data <- plot_data[plot_data$covar %in% covars,]

  }

  if (!is.null(att_levels)) {

    plot_data <- merge(x = plot_data,
                       y = x$att_lookup,
                       by = c("Attribute","Level"),
                       all.x = TRUE)

    plot_data <- plot_data[plot_data$level %in% att_levels,]

  }

  if (is.null(att_levels) | length(att_levels) > 1) {

    ggplot2::ggplot(plot_data,
                    ggplot2::aes_string(x = "covar",
                                        y = "Level",
                                        fill = "importance")) +
      ggplot2::facet_grid(Attribute ~ ., space = "free", scales = "free") +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low="white", high="firebrick1") +
      ggplot2::labs(x = "Covariates", y = "Attribute-level", fill = "Importance") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust = 1, hjust = 1)) +
      ggplot2::theme_classic()

  } else {

    message("Plotting importance scores for single attribute-level (with 95 percent confidence intervals)")

    ggplot2::ggplot(plot_data,
                    ggplot2::aes_string(y = "covar",
                                        x = "importance",
                                        xmin = "lower2.5",
                                        xmax = "upper97.5")) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbarh(height = 0.2) +
      ggplot2::labs(x = "Standardised Importance",
                    y = "Covariates") +
      ggplot2::xlim(min(plot_data$lower2.5)-1, max(plot_data$upper97.5)+1) +
      ggplot2::labs(title = paste0(unique(plot_data$Attribute),": ", unique(plot_data$Level))) +
      ggplot2::theme_classic()

  }

}

