#' Estimate Variable Importance Metrics for \code{cjbart} Object
#' @description Estimates random forest variable importance scores for multiple attribute-levels of a conjoint experiment.
#' @param imces Object of class \code{cjbart}, the result of running [cjbart::IMCE()]
#' @param levels An optional vector of attribute-levels to generate importance metrics for. By default, all attribute-levels are analyzed.
#' @param covars An optional vector of covariates to include in the importance metric check. By default, all covariates are included in each importance model.
#' @param cores Number of CPU cores used during VIMP estimation. Each extra core will result in greater memory consumption. Assigning more cores than outcomes will not further boost performance.
#' @param ... Extra arguments (used to check for deprecated argument names)
#' @return A "long" data.frame of variable importance scores for each combination of covariates and attribute-levels, as well as the estimated 95% confidence intervals for each metric.
#' @details Having generated a schedule of individual-level marginal component effect estimates, this function fits a random forest model for each attribute-level using the supplied covariates as predictors. It then calculates a variable importance measure (VIMP) for each covariate. The VIMP method assesses how important each covariate is in terms of partitioning the predicted individual-level effects distribution, and can thus be used as an indicator of which variables drive heterogeneity in the IMCEs.
#'
#' To recover a VIMP measure, we used permutation-based importance metrics recovered from random forest models estimated using [randomForestSRC::rfsrc()]. To permute the data, this function uses random node assignment, whereby cases are randomly assigned to a daughter node whenever a tree splits on the target variable \insertCite{@see @ishwaran2008random}{cjbart}. Importance is defined in terms of how random node assignment degrades the performance of the forest. Higher degradation indicates a variable is more important to prediction.
#'
#' Variance estimates of each variable's importance are subsequently recovered using the delete-d jackknife estimator developed by \insertCite{ishwaran2019standard;textual}{cjbart}. The jackknife method has inherent bias correction properties, making it particularly effective for variable selection exercises such as identifying drivers of heterogeneity.
#' @references \insertAllCited{}
#' @seealso [randomForestSRC::rfsrc()] and [randomForestSRC::subsample()]
#' @importFrom Rdpack reprompt
#' @export
het_vimp <- function(imces, levels = NULL, covars = NULL, cores = 1, ...) {

  if (is.null(levels)) {
    levels <- imces$att_levels
  }

  extra_args = names(list(...))
  deprecated_args = c("model","outcomes")
  deprec_present <- deprecated_args[deprecated_args %in% extra_args]

  if (length(deprec_present) > 0) {
    warning("\nDetected deprecated argument names!\n The following argument names were changed in cjbart v0.3:",
            "\n * `model` (old) -> 'imces' (new) \n * `outcomes` -> `levels`",
            "\n Please check any output and update your scripts accordingly.")
  }

  if (cores == 1 | .Platform$OS.type != "unix") {

    if (cores > 1) {
      warning("Parallelization is not available on your operating system. Reverting to single-core computation.")
    }

    vimps <- lapply(levels,
                    function (x) rf_vimp(model = imces, outcome = x, covars = covars))
  } else {

    vimps <- parallel::mclapply(levels,
                      function (x) rf_vimp(model = imces, outcome = x, covars = covars),
                      mc.cores = cores)

  }

  full_results <- do.call("rbind", vimps)

  full_results <- merge(x = imces$att_lookup,
                        y = full_results,
                        all.y = TRUE,
                        by.x = "level", by.y = "outcome")

  full_results$level <- NULL

  out_obj <- list(results = full_results,
                  att_lookup = imces$att_lookup)

  class(out_obj) <- c("cjbart.vimp")

  return(out_obj)
}

#' Estimate a Single Variable Importance Metric for \code{cjbart} Object
#' @description Estimates random forest variable importance scores for a single attribute-level of a conjoint experiment. This function is for advanced use. Users should typically use the [cjbart::het_vimp()] function.
#' @param model Object of class \code{cjbart}, the result of running [cjbart::IMCE()]
#' @param outcome Character string detailing the covariate over which to analyze heterogeneous effects
#' @param covars An optional vector of covariates to include in the importance metric check. When \code{covars = NULL} (the default), all covariates are included in the importance model.
#' @return Data.frame of variable importance scores for each covariate in the model, as well as values for the estimated 95% confidence interval for each importance score.
#' @export
rf_vimp <- function(model, outcome, covars = NULL) {

  trial_data <- model$imce
  trial_data$outcome <- trial_data[[outcome]]

  if (is.null(covars)) {

    covars <- colnames(model$imce)[!(colnames(model$imce) %in% c(model$att_levels, model$id))]

  }

  trial_data <- trial_data[,c("outcome", covars)]

  trial_data <- stats::na.omit(trial_data)

  n_na <- nrow(model$imce) - nrow(trial_data)

  if (n_na != 0) {
    warning(paste0(n_na," observations removed due to missing covariate data in the model."))
  }

  for (i in 1:length(names(trial_data))) {
    if (typeof(trial_data[,i]) == "character") {
      trial_data[,i] <- as.factor(trial_data[,i])
    }
  }

  message(paste0("Calculating importance metrics for attribute-level: ",outcome))
  rf_mod <- randomForestSRC::rfsrc(outcome ~ ., data = trial_data, importance = "permute")
  vimp_ci <- randomForestSRC::extract.subsample(randomForestSRC::subsample(rf_mod), raw = TRUE)$ci.jk.Z

  att_results <- data.frame(outcome = outcome,
                            covar = colnames(vimp_ci),
                            importance = vimp_ci["50%",],
                            lower2.5 = vimp_ci["2.5%",],
                            upper97.5 = vimp_ci["97.5%",])

  class(att_results) <- c("cjbart.vimp","data.frame")

  return(att_results)

}
