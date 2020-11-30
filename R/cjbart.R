#' Generate Conjoint Model Using BART
#'
#' @description A wrapper for the \code{BART::pbart()} function.
#' @param data A data.frame or coercible, containing all attributes, covariates, the outcome and id variables to analyse.
#' @param Y_var Character string -- the outcome variable
#' @param id_var Character string -- the variable identifying individual respondents
#' @param ... Other arguments passed to pbart
#' @return Tibble summarising the average marginal component effect, the minimum and maximum values, and standard deviations for each attribute-level.
#' @details Please note, \code{cjbart} currently only works for a binary outcome.
#' @seealso [BART::pbart()]
#' @export
cjbart <- function(data, Y_var, id_var = NULL, ...) {

  if (missing(Y_var)) {
    stop("Please declare the output variable using the Y_var argument.")
  }

  if (!missing(id_var)) {
    if (!(id_var %in% names(data))) {
      stop("Supplied id_var not present in data.")
    }
  }

  train_X <- data %>%
    dplyr::select(-{{id_var}},-{{Y_var}}) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    as.data.frame(.data)

  train_Y <- data[[Y_var]]

  train_model <- BART::pbart(x.train = train_X, y.train = train_Y, ...)

  return(train_model)
}

#' Heterogeneous Effects Analysis of Conjoint Results
#'
#' @description \code{OMCE} calculates the marginal component effects at both the observation- and individual-level, from a BART-estimated conjoint model.
#' @param data A data.frame or coercible, containing all attributes, covariates, the outcome and id variables to analyse.
#' @param model A model object, the result of running \code{cjbart}
#' @param attribs Vector of attribute names
#' @param ref_levels Vector of reference levels, used to calculate marginal effects
#' @param Y_var Character string -- the outcome variable
#' @param id_var Character string -- the variable identifying individual respondents
#' @param cores Number of CPU cores used during prediction phase
#' @details The OMCE estimates are the result of subtracting the predicted value of each observation under the reference-level category from the predicted value of each observation under the given attribute level.
#' If an attribute has \code{k} levels, then this will yield \code{k-1} estimates per observation.
#' The IMCE is the average of the OMCEs for each individual within the data.
#' @return \code{OMCE} returns an object of type "cjbart", a list object.
#' \item{omce}{A data.frame containing the observation-level marginal effects}
#' \item{imce}{A data.frame containing the individual-level marginal effects}
#' \item{att_levels}{A vector containing the attribute levels}
#' @seealso [cjbart()]
#' @import BART
#' @importFrom stats predict
#' @export
OMCE <- function(data, model, attribs, ref_levels, Y_var, id_var, cores = 1) {

  data <- as.data.frame(data)

  # Data frame to store OMCEs
  results <- data %>%
    dplyr::select_if(!(names(.data) %in% attribs)) %>%
    dplyr::select(-{{Y_var}})

  # Vector to store attribute names (for future function calls)
  out_levels <- c()

  for (i in 1:length(attribs)) {

    message("Calculating OMCEs for attribute: ", attribs[i], " [",i,"/",length(attribs),"]")

    att_levels <- unique(data[[attribs[i]]][data[[attribs[i]]] != ref_levels[i]])

    out_levels <- c(out_levels, as.character(att_levels))

    X_pred0 <- data %>%
      dplyr::select(-{{Y_var}},
                               -{{id_var}}) %>%
      dplyr::mutate_if(is.character, as.factor)

    X_pred0[[attribs[i]]] <- factor(ref_levels[i],
                                    levels = levels(data[[attribs[i]]]))

    for (att_level in att_levels) {

      X_pred1 <- data %>%
        dplyr::select(-{{Y_var}},
                                 -{{id_var}}) %>%
        dplyr::mutate_if(is.character, as.factor)

      X_pred1[[attribs[i]]] <- factor(att_level,
                                      levels = levels(data[[attribs[i]]]))

      # Get predictions
      pred1 <- predict(model,
                       newdata = BART::bartModelMatrix(X_pred1),
                       mc.cores = cores)

      pred1_prob <- stats::pnorm(colMeans(pred1$yhat.test)) # Converts probit to predicted probabilities

      pred0 <- predict(model,
                       newdata = BART::bartModelMatrix(X_pred0),
                       mc.cores = cores)

      pred0_prob <- stats::pnorm(colMeans(pred0$yhat.test))

      # Get OMCE for single attribute-level comparison
      att_level_OMCE <- pred1_prob - pred0_prob

      # Store results in data.frame
      results[[as.character(att_level)]] <- att_level_OMCE

    }

  }

  ## IMCE

  covars <- results %>%
    dplyr::select(-all_of(out_levels)) %>%
    dplyr::distinct(.data)

  if(nrow(covars) != length(unique(results[[id_var]]))) {
    stop("Covariates vary within id.")
  }

  if(!(id_var %in% names(covars))) {
    stop("Could not find id variable in covariate matrix")
  }

  results_imce <- results %>%
    dplyr::group_by_at(id_var) %>%
    dplyr::summarise_at(out_levels, mean) %>%
    dplyr::left_join(covars, by = {{id_var}})

  out_obj <- list(omce = results,
                  imce = results_imce,
                  att_levels = out_levels)

  class(out_obj) <- "cjbart"

  return(out_obj)

}





