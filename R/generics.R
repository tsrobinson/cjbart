#' Plot Marginal Component Effects of a \code{cjbart} Object
#' @description Plots observation-level or individual-level marginal component effects (OMCE and IMCE respectively). By default, all attribute-levels in the model are plotted.
#' @param x Object of class \code{cjbart}, the result of running [cjbart::IMCE()]
#' @param covar Character string detailing the covariate over which to analyze heterogeneous effects
#' @param plot_levels Optional vector of conjoint attribute names to plot. If not supplied, all attributes within the conjoint model will be plotted.
#' @param se Boolean determining whether to show an estimated 95% confidence interval
#' @param ... Additional arguments for plotting the marginal component effects (see below).
#' @return Plot of marginal component effects.
#' @importFrom rlang .data
#' @method plot cjbart
#' @export
plot.cjbart <- function(x, covar = NULL, plot_levels = NULL, se = TRUE,  ...) {

  data <- x$imce

  plot_data <- tidyr::pivot_longer(cols = x$att_levels,
                                   names_to = "att",
                                   values_to = "imce",
                                   data = data)

  if (se) {

    ci_lower <- tidyr::pivot_longer(cols = x$att_levels,
                                    names_to = "att",
                                    values_to = "imce_lower",
                                    data = x$imce_lower)

    ci_upper <- tidyr::pivot_longer(cols = x$att_levels,
                                    names_to = "att",
                                    values_to = "imce_upper",
                                    data = x$imce_upper)

    plot_data <- cbind(plot_data, imce_lower = ci_lower$imce_lower, imce_upper =  ci_upper$imce_upper)

  }

  if (!is.null(plot_levels)) {
    plot_data <- plot_data[plot_data$att %in% plot_levels,]

    if (nrow(plot_data) == 0) {
      stop("Filtering attribute levels renders data with no observations -- please check plot_levels argument.")
    }

  }

  plot_data <- do.call(

    rbind,

    by(
      plot_data,
      INDICES = plot_data$att,
      FUN = function(x) {
        x_ordered <- x[order(x$imce),]
        x_ordered$x_order <- 1:nrow(x_ordered)
        return(x_ordered)
        }
     )
  )

  base_plot <- ggplot2::ggplot(plot_data,
                                 ggplot2::aes_string(x = "x_order",
                                                     y = "imce")) +

    ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed") +

    {if (!is.null(covar)) {ggplot2::aes_string(color = covar)}} +
    {if (se) {ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "imce_lower", ymax = "imce_upper"), color = NA, fill = "grey60", alpha = 0.7)}} +

    ggplot2::facet_wrap(~.data$att, scales = "free") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::ylab("IMCE") +
    ggplot2::xlab("Individual") +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())

  if (!is.null(covar)) {

    if (typeof(plot_data[[covar]]) == "double") {

      base_plot <- base_plot +
        ggplot2::scale_color_gradient(low = "dodgerblue3", high = "goldenrod1")

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
#' @return Data frame summarizing the average marginal component effect, the minimum and maximum values, and standard deviations for each attribute-level.
#' @method summary cjbart
#' @example inst/examples/basic_workflow.R
#' @export
summary.cjbart <- function(object, ...) {

  # IMCE summary
  IMCE_only <- subset(object$imce, select = object$att_levels)

  AMCE <- colMeans(IMCE_only)
  mins <- apply(IMCE_only,2,min)
  maxs <- apply(IMCE_only,2,max)
  sds <- apply(IMCE_only,2,stats::sd)

  att_names <-  ifelse(test = nchar(object$att_levels) > 30,
                       yes = paste0(substr(object$att_levels,1,30),"..."),
                       no = object$att_levels)

  summary_tab <- data.frame(Level = att_names,
                            AMCE = AMCE,
                            `Min.` = mins,
                            `Max.` = maxs,
                            `Std.Dev` = sds,
                            row.names = NULL)

  message("Summary table of individual marginal component effects (IMCEs)")
  print(summary_tab)

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

  plot_data <- x

  if (!is.null(covars)) {

    plot_data <- plot_data[plot_data$covar %in% covars,]

  }

  if (!is.null(att_levels)) {

    plot_data <- plot_data[plot_data$outcome %in% att_levels,]

  }

  if (is.null(att_levels) | length(att_levels) > 1) {

    ggplot2::ggplot(plot_data,
                    ggplot2::aes_string(x = "covar",
                                        y = "outcome",
                                        fill = "importance")) +
      ggplot2::facet_grid(attribute ~ ., space = "free", scales = "free", switch = "y") +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low="white", high="firebrick2") +
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
      ggplot2::theme_classic()

  }

}

