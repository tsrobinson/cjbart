#' Generate Conjoint Model Using BART
#'
#' @description A wrapper for the \code{BART::pbart()} function.
#' @param data A data.frame or coercible, containing all attributes, covariates, the outcome and id variables to analyse.
#' @param Y_var Character string -- the outcome variable
#' @param id_var Character string -- the variable identifying individual respondents
#' @param ... Other arguments passed to pbart
#' @return Tibble summarising the average marginal component effect, the minimum and maximum values, and standard deviations for each attribute-level.
#' @details Please note, \code{cjbart} currently only works for a binary outcome.
#' @importFrom BART pbart
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

  train_vars <- names(data)[!(names(data) %in% c(id_var,Y_var))]

  train_X <- data[, train_vars]

  train_X <- .char_to_fact(train_X)

  train_Y <- data[[Y_var]]

  train_model <- BART::mc.pbart(x.train = as.data.frame(train_X),
                             y.train = train_Y, ...)

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
#' @importFrom stats predict
#' @export
OMCE <- function(data, model, attribs, ref_levels, Y_var, id_var, cores = 1) {

  data <- as.data.frame(data)

  # Check incoming data

  test_data <- data[,!(names(data) %in% c(attribs,Y_var))]
  test_data <- test_data[!duplicated(test_data),]

  if(!(id_var %in% names(test_data))) {
    stop("Could not find id variable in covariate matrix")
  }

  if(nrow(test_data) != length(unique(test_data[[id_var]]))) {
    stop("Covariates vary within id: data must not contain covariates that vary across observations of the same subject")
  }

  rm(test_data)
  gc()

  # Data frame to store OMCEs
  results <- data[,!(names(data) %in% c(attribs, Y_var))]

  # Data frame for predicting outcomes
  train_vars <- names(data)[!(names(data) %in% c(id_var,Y_var))]
  data_predict <- data[,train_vars]

  # Vector to store attribute names (for future function calls)
  out_levels <- c()

  for (i in 1:length(attribs)) {

    message("Calculating OMCEs for attribute: ", attribs[i], " [",i,"/",length(attribs),"]")

    att_levels <- unique(data[[attribs[i]]][data[[attribs[i]]] != ref_levels[i]])

    out_levels <- c(out_levels, as.character(att_levels))

    X_pred0 <- data_predict

    X_pred0 <- .char_to_fact(X_pred0)

    X_pred0[[attribs[i]]] <- factor(ref_levels[i],
                                    levels = levels(X_pred0[[attribs[i]]]))

    phat_0 <- .quiet(

      predict(

        model,
        newdata = BART::bartModelMatrix(X_pred0),
        mc.cores = cores
      )
    )$prob.test.mean

    for (att_level in att_levels) {

      X_pred1 <- X_pred0

      X_pred1[[attribs[i]]] <- factor(att_level,
                                      levels = levels(X_pred0[[attribs[i]]]))

      # Get predictions

      phat_1 <- .quiet(

          predict(

            model,
            newdata = BART::bartModelMatrix(X_pred1),
            mc.cores = cores

            )
          )$prob.test.mean

      ## Note, prob.test.mean is equivalent to
      # stats::pnorm(colMeans(pred0$yhat.test))

      # Get OMCE for single attribute-level comparison
      att_level_OMCE <- phat_1 - phat_0

      # Store results in data.frame
      results[[as.character(att_level)]] <- att_level_OMCE

    }

  }

  ## IMCE

  message("Calculating IMCEs")

  covars <- results[,!(names(results) %in% out_levels)]

  covars <- covars[!duplicated(covars),]

  agg_formula <- stats::as.formula(
    paste0(
      "cbind(",
      paste0(paste0("`",out_levels,"`"), collapse = ", "),
      ") ~ ",
      id_var
      )
    )

  results_imce <- stats::aggregate(formula = agg_formula,
                            data = results,
                            FUN = mean)

  out_obj <- list(omce = results,
                  imce = results_imce,
                  att_levels = out_levels)

  class(out_obj) <- "cjbart"

  return(out_obj)

}

