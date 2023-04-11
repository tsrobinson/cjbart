#' Generate Conjoint Model Using BART
#'
#' @description A wrapper for the [BART::pbart()] function.
#' @param data A data.frame, containing all attributes, controls, the outcome and id variables to analyze.
#' @param Y Character string -- the outcome variable
#' @param type Type of conjoint experiment -- either "choice" (for forced-choice outcomes) or "rating" (for interval ratings). If NULL (default), the function will attempt to automatically detect the outcome type.
#' @param id Character string -- variable identifying individual respondents (optional)
#' @param round Character string -- variable identifying rounds of the conjoint experiment
#' @param use_round Boolean -- whether to include the round indicator column when training the BART model (default = \code{TRUE})
#' @param cores Integer -- number of CPU cores used in model training
#' @param ... Other arguments passed to [BART::pbart()]
#' @return A trained [BART::pbart()] model that can be passed to [cjbart::IMCE()]
#' @details Please note, \code{cjbart} currently only works for a binary outcome.
#' @importFrom BART pbart
#' @seealso [BART::pbart()]
#' @export
#' @example inst/examples/cjbart_model_example.R
cjbart <- function(data, Y, type = NULL, id = NULL, round = NULL, use_round = TRUE, cores = 1, ...) {

  if (missing(Y)) {
    stop("Please declare the output variable using the Y argument.")
  }

  if (!missing(id)) {
    if (!(id %in% names(data))) {
      stop("Supplied id column not present in data.")
    }
  } else {
    id <- NULL
  }

  if (!missing(round)) {
    if (!(round %in% names(data))) {
      stop("Supplied round column not present in data.")
    }
  } else {
    round <- NULL
  }

  if (use_round) {
    train_vars <- names(data)[!(names(data) %in% c(id,Y))]
  } else {
    train_vars <- names(data)[!(names(data) %in% c(id,Y,round))]
  }

  train_X <- data[, train_vars]

  train_X <- .char_to_fact(train_X)

  train_Y <- data[[Y]]

  if (!missing(type)) {
    if (!(type %in% c("choice","rating"))) {
      stop("Invalid type argument: must be 'choice' (forced-choice outcomes) or 'rating' (interval ratings outcomes), or NULL.")
    }
  } else {
    if (setequal(unique(train_Y),c(0,1))) {
      type = "choice"
      message("Detected forced-choice outcome")
    } else {
      type = "rating"
      message("Detected interval (rating) outcome")
    }
  }

  # Store factor levels
  factor_levels <- list()
  for (v in colnames(train_X)) {

    if (is.factor(train_X[[v]])) {

      factor_levels[[v]] <- levels(train_X[[v]])

    }

  }

  if (.Platform$OS.type=='unix') {

    if (type == "choice") {
      train_model <- BART::mc.pbart(x.train = as.data.frame(train_X),
                                    y.train = train_Y, mc.cores = cores, ...)
    } else {
      train_model <- BART::mc.wbart(x.train = as.data.frame(train_X),
                                    y.train = train_Y, mc.cores = cores, ...)
    }


  } else {

    if (type == "choice") {
      train_model <- BART::pbart(x.train = as.data.frame(train_X),
                                 y.train = train_Y, ...)
    } else {
      train_model <- BART::wbart(x.train = as.data.frame(train_X),
                                 y.train = train_Y, ...)
    }
  }

  gc()

  train_model$id_col <- id
  train_model$round_col <- round
  train_model$Y_col <- Y
  train_model$factor_levels <- factor_levels
  train_model$type <- type

  return(train_model)

}

#' Heterogeneous Effects Analysis of Conjoint Results
#'
#' @description \code{IMCE} calculates the individual-level marginal component effects from a BART-estimated conjoint model.
#' @param data A data.frame, containing all attributes, covariates, the outcome and id variables to analyze.
#' @param model A model object, the result of running \code{cjbart()}
#' @param attribs Vector of attribute names for which IMCEs will be predicted
#' @param ref_levels Vector of reference levels, used to calculate marginal effects
#' @param method Character string, setting the variance estimation method to use. When method is "parametric", a typical combined variance estimate is employed; when \code{method = "bayes"}, the 95% posterior interval is calculated; and when \code{method = "rubin"}, combination rules are used to combine the variance analogous to in multiple imputation analysis.
#' @param alpha Number between 0 and 1 -- the significance level used to compute confidence/posterior intervals. When \code{method = "bayes"}, the posterior interval is calculated by taking the alpha/2 and (1-alpha/2) quantiles of the posterior draws. When \code{method = "rubin"}, the confidence interval equals the IMCE +/- \code{qnorm(alpha/2)}. By default, alpha is 0.05 i.e. generating a 95% confidence/posterior interval.
#' @param keep_omce Boolean, indicating whether to keep the OMCE-level results (default = \code{FALSE})
#' @param cores Number of CPU cores used during prediction phase
#' @param skip_checks Boolean, indicating whether to check the structure of the data (default = \code{FALSE}). Only set this to \code{TRUE} if you are confident that the data is structured appropriately
#' @details The OMCE estimates are the result of subtracting the predicted value of each observation under the reference-level category from the predicted value of each observation under the given attribute level.
#' If an attribute has *k* levels, then this will yield *k-1* estimates per observation.
#' The IMCE is the average of the OMCEs for each individual within the data.
#' @return \code{IMCE} returns an object of type "cjbart", a list object.
#' \item{omce}{A data.frame containing the observation-level marginal effects}
#' \item{imce}{A data.frame containing the individual-level marginal effects}
#' \item{imce_upper}{A data.frame containing the upper bound of the IMCE confidence/credible interval}
#' \item{imce_lower}{A data.frame containing the lower bound of the IMCE confidence/credible interval}
#' \item{att_levels}{A vector containing the attribute levels}
#' @seealso [cjbart::cjbart()]
#' @importFrom stats predict
#' @example inst/examples/basic_workflow.R
#' @export
IMCE <- function(data,
                 model,
                 attribs,
                 ref_levels,
                 method = "bayes",
                 alpha = 0.05,
                 keep_omce = FALSE,
                 cores = 1,
                 skip_checks = FALSE) {

  data <- as.data.frame(data)

  # Get variables from trained model
  Y <- model$Y_col
  round <- model$round_col
  id <- model$id_col
  unq_ids <- unique(data[[id]])
  type <- model$type

  # Check optional args
  if (!(method %in% c("average","bayes","rubin"))) {
    stop("Variance estimation method must be in c('parametric','bayes','rubin'). See ?OMCE for more details.")
  }


  if (!skip_checks) {
    test_data <- data[,!(names(data) %in% c(attribs, Y, round))]

    if (is.data.frame(test_data)) { # Skip if only one variable left (assume it is id)

      test_data <- test_data[!duplicated(test_data),]

      if(!(id %in% names(test_data))) {
        stop("Could not find id variable in data")
      }

      if(nrow(test_data) != length(unique(test_data[[id]]))) {
        warning("Covariates vary within id: if this is not intentional, please check your data")
      }

    }

    if (!sum(sapply(attribs, function(x) class(data[[x]])) %in% c("character","factor")) == length(attribs)) {

      stop("Conjoint attribute columns must be character vectors")

    }

    rm(test_data)
    gc()
  }

  # Check attribute-level unique and rename if needed
  att_lev_list <- lapply(attribs, function (x) data.frame(Attribute = x, level = unique(data[[x]])))
  att_lookup <- do.call("rbind", att_lev_list)

  if (length(unique(att_lookup$level)) < nrow(att_lookup)) {

    ref_levels <- paste0("<",attribs,">_", ref_levels)
    for (att in attribs) {
      data[[att]] <- paste0("<",att,">_", data[[att]])
    }

    warning("Found duplicate level names across attributes: ",
            "[",unique(att_lookup$level[duplicated(att_lookup$level)]),"].",
            "\n*All* levels in data will be renamed to format <attribute>_level, ",
            "e.g. ", ref_levels[1],
            "\n You can use the `att_lookup` table in the return object to compare new and old names")

    for (i in 1:length(attribs)) {
      if (!(ref_levels[i] %in% data[[attribs[i]]])) {
        stop("Failed to re-label data. Try manually converting conjoint data so all attribute-level names are unique.")
      }
    }

    att_lev_list <- lapply(attribs, function (x) data.frame(Attribute = x, level = unique(data[[x]])))
    att_lookup <- do.call("rbind", att_lev_list)
  }

  # Add label to lookup (does nothing if no duplicates)
  att_lookup$Level <- sub("<.*>_", "", att_lookup$level)


  # Frame to store OMCEs
  results <- data[,!(names(data) %in% c(attribs, Y))]

  # If only attribs and id in data.frame, conditional to correct formatting
  if (!is.data.frame(results)) {

    if (identical(results,data[[id]])) {
      results <- as.data.frame(results)
      names(results) <- id
    } else {
      stop("Results frame creation failed")
    }

  }

  # Variance and interval frames
  var_omce <- data.frame(row.names = 1:nrow(data))
  imce_lower <- imce_upper <- data.frame(row.names = 1:length(unq_ids))

  # Data frame for predicting outcomes
  train_vars <- names(data)[!(names(data) %in% c(id,Y))]
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


    if (type == "choice") {
      phat_0 <- .quiet(

        predict(

          model,
          newdata = BART::bartModelMatrix(X_pred0),
          mc.cores = cores
        )
      )$prob.test
    } else {
      phat_0 <- .quiet(

        predict(

          model,
          newdata = BART::bartModelMatrix(X_pred0),
          mc.cores = cores
        )
      )
    }

    for (att_level in att_levels) {

      X_pred1 <- X_pred0

      X_pred1[[attribs[i]]] <- factor(att_level,
                                      levels = levels(X_pred0[[attribs[i]]]))

      # Get predictions
      if (type == "choice") {
        phat_1 <- .quiet(

          predict(

            model,
            newdata = BART::bartModelMatrix(X_pred1),
            mc.cores = cores

          )
        )$prob.test
      } else {
        phat_1 <- .quiet(

          predict(

            model,
            newdata = BART::bartModelMatrix(X_pred1),
            mc.cores = cores

          )
        )
      }


      ## Note, prob.test.mean is equivalent to
      # stats::pnorm(colMeans(pred_0$yhat.test))

      # Get OMCE for single attribute-level comparison and store
      results[[as.character(att_level)]] <- colMeans(phat_1) - colMeans(phat_0)
      var_omce[[as.character(att_level)]] <- apply(phat_1 - phat_0, 2, stats::var)

      if (method == "bayes") {

        # Save interval as vector to make code easier to read.
        intvl <- c(alpha/2, (1-alpha/2))

        # Calculate distribution of marginal effects
        var_z <- phat_1 - phat_0

        # Calculate IMCE interval at this point to avoid holding many frames in memory
        imce_ci <- sapply(unq_ids, function(s) stats::quantile(var_z[,data[[id]] == s], intvl))

        imce_lower[[as.character(att_level)]] <- imce_ci[1,]
        imce_upper[[as.character(att_level)]] <- imce_ci[2,]

        rm(var_z, imce_ci)

      }

      rm(X_pred1, phat_1)
      gc()

    }

  }

  ## IMCE

  message("Calculating IMCEs")

  covars <- results[,!(names(results) %in% c(out_levels))]

  # In case only id is supplied, make sure covariates stored as data.frame
  if (is.data.frame(covars)) {
    covars <- covars[!duplicated(covars),]
  } else {
    covars <- data.frame(covars[!duplicated(covars)])
    names(covars) <- id
  }

  agg_formula <- stats::as.formula(
    paste0(
      "cbind(",
      paste0(paste0("`",out_levels,"`"), collapse = ", "),
      ") ~ ",
      id
      )
    )

  results_imce <- stats::aggregate(agg_formula,
                                   data = results,
                                   FUN = mean)

  # Double check nrow now we have covariates recovered
  if (!(nrow(covars) == nrow(results_imce))) {

    warning("Number of unique covariate rows does not match number of ids -- attempting to merge data, but please check results.")

  }

  results_imce <- merge(results_imce, covars, by = id)

  if (method == "bayes") {

    # Note bounds are calculated above to prevent having to OMCE draws
    imce_upper <- as.data.frame(imce_upper)
    imce_lower <- as.data.frame(imce_lower)

    imce_lower[[id]] <- imce_upper[[id]] <- unq_ids

  } else if (method == "rubin") {

    results_var <- sapply(colnames(var_omce), function (x) {
      sapply(results_imce[[id]], function (y) {

        .combine(theta = results[results[[id]] == y, x],
                 var_theta = var_omce[results[[id]] == y,x])
        }
      )
    })

    imce_upper <- sapply(colnames(results_var), function (x) results_imce[[x]] + stats::qnorm(1-(alpha/2))*results_var[,x])
    imce_lower <- sapply(colnames(results_var), function (x) results_imce[[x]] + stats::qnorm(alpha/2)*results_var[,x])

    imce_upper <- as.data.frame(imce_upper)
    imce_lower <- as.data.frame(imce_lower)

    imce_lower[[id]] <- imce_upper[[id]] <- results_imce[[id]]

  } else if (method == "average") {

    results_var <- sapply(colnames(var_omce), function (x) {
      sapply(results_imce[[id]], function (y) {

        mean(var_omce[results[[id]] == y,x])
      }
      )
    })

    imce_upper <- sapply(colnames(results_var), function (x) results_imce[[x]] + stats::qnorm(1-(alpha/2))*sqrt(results_var[,x]))
    imce_lower <- sapply(colnames(results_var), function (x) results_imce[[x]] + stats::qnorm(alpha/2)*sqrt(results_var[,x]))

    imce_upper <- as.data.frame(imce_upper)
    imce_lower <- as.data.frame(imce_lower)

    imce_lower[[id]] <- imce_upper[[id]] <- results_imce[[id]]
  }

  out_obj <- list(imce = results_imce,
                  imce_lower = imce_lower,
                  imce_upper = imce_upper,
                  alpha = alpha,
                  att_levels = out_levels,
                  id = id,
                  att_lookup = att_lookup,
                  omce = NULL,
                  imce_var = NULL,
                  round = NULL)

  # Fill in optional elements

  if (keep_omce) {
    out_obj$omce <- results
  }

  if (method %in% c("rubin","average")) {
    out_obj$imce_var = results_var
  }

  if (!is.null(round)) {
    out_obj$round = round
  }

  class(out_obj) <- "cjbart"

  return(out_obj)

}


#' Inspect Round-Level Marginal Component Effect (RMCE)
#'
#' @description \code{RMCE} calculates the round-level marginal component effects from a cjbart model.
#' @param imces An object of class "cjbart", the result of calling the \code{cjbart::IMCE()} function
#' @details The RMCE estimates are the result of averaging the OMCEs within each round, for each subject in the experiment.
#' The RMCE is the intermediate causal quantity between OMCEs and IMCEs, and can be useful for inspecting whether there are any carryover or stability issues across rounds.
#' @return \code{IMCE} returns a data frame of RMCEs.
#' @seealso [cjbart::cjbart()] and [cjbart::IMCE()]
#' @export
RMCE <- function(imces) {

  if (!inherits(imces, "cjbart")) {
    stop("imces must be the result of calling IMCE()")
  }

  if (is.null(imces$omce)) {
    stop("To calculate the RMCE, keep_omce = TRUE when IMCE() is called")
  }

  agg_formula <- stats::as.formula(
    paste0(
      "cbind(",
      paste0(paste0("`",imces$att_levels,"`"), collapse = ", "),
      ") ~ ",
      imces$round
    )
  )

  results_rmce <- stats::aggregate(agg_formula,
                                   data = imces$omce,
                                   FUN = mean)

}

#' Average Marginal Component Effect Estimation with Credible Interval
#'
#' @description \code{AMCE} calculates the average marginal component effects from a BART-estimated conjoint model.
#' @param data A data.frame, containing all attributes, covariates, the outcome and id variables to analyze.
#' @param model A model object, the result of running \code{cjbart()}
#' @param attribs Vector of attribute names for which IMCEs will be predicted
#' @param ref_levels Vector of reference levels, used to calculate marginal effects
#' @param method Character string, setting the variance estimation method to use. When method is "parametric", a typical combined variance estimate is employed; when \code{method = "bayes"}, the 95% posterior interval is calculated; and when \code{method = "rubin"}, combination rules are used to combine the variance analogous to in multiple imputation analysis.
#' @param alpha Number between 0 and 1 -- the significance level used to compute confidence/posterior intervals. When \code{method = "bayes"}, the posterior interval is calculated by taking the alpha/2 and (1-alpha/2) quantiles of the posterior draws. When \code{method = "rubin"}, the confidence interval equals the IMCE +/- \code{qnorm(alpha/2)}. By default, alpha is 0.05 i.e. generating a 95% confidence/posterior interval.
#' @param cores Number of CPU cores used during prediction phase
#' @param skip_checks Boolean, indicating whether to check the structure of the data (default = \code{FALSE}). Only set this to \code{TRUE} if you are confident that the data is structured appropriately
#' @details The AMCE estimates are the average of all computed OMCEs.
#' @return \code{AMCE} returns an object of type "cjbart", a list object.
#' \item{amces}{A data.frame containing the average marginal component effects}
#' \item{alpha}{The significance level used to compute the credible interval}
#' @seealso [cjbart::cjbart()]
#' @importFrom stats predict
#' @export
AMCE <- function(data,
                 model,
                 attribs,
                 ref_levels,
                 method = "bayes",
                 alpha = 0.05,
                 cores = 1,
                 skip_checks = FALSE) {

  data <- as.data.frame(data)

  # Get variables from trained model
  Y <- model$Y_col
  round <- model$round_col
  id <- model$id_col
  type <- model$type

  # Check optional args
  if (method != "bayes") {
    warning("Only the Bayesian credible interval strategy is available at present. The method argument has been overridden to 'bayes'.")
  }

  if (!skip_checks) {
    test_data <- data[,!(names(data) %in% c(attribs, Y, round))]

    if (is.data.frame(test_data)) { # Skip if only one variable left (assume it is id)

      test_data <- test_data[!duplicated(test_data),]

      if(!(id %in% names(test_data))) {
        stop("Could not find id variable in data")
      }

      if(nrow(test_data) != length(unique(test_data[[id]]))) {
        warning("Covariates vary within id: if this is not intentional, please check your data")
      }

    }

    if (!sum(sapply(attribs, function(x) class(data[[x]])) %in% c("character","factor")) == length(attribs)) {

      stop("Conjoint attribute columns must be character vectors")

    }

    rm(test_data)
    gc()
  }

  # Frame to store AMCEs
  results <- data.frame(attribute = as.character(),
                        level = as.character(),
                        AMCE = as.numeric(),
                        AMCE_lower = as.numeric(),
                        AMCE_upper = as.numeric())

  # Data frame for predicting outcomes
  train_vars <- names(data)[!(names(data) %in% c(id,Y))]
  data_predict <- data[,train_vars]

  # Vector to store attribute names (for future function calls)
  out_levels <- c()

  for (i in 1:length(attribs)) {

    message("Calculating AMCEs for attribute: ", attribs[i], " [",i,"/",length(attribs),"]")

    att_levels <- unique(data[[attribs[i]]][data[[attribs[i]]] != ref_levels[i]])

    out_levels <- c(out_levels, as.character(att_levels))

    X_pred0 <- data_predict

    X_pred0 <- .char_to_fact(X_pred0)

    X_pred0[[attribs[i]]] <- factor(ref_levels[i],
                                    levels = levels(X_pred0[[attribs[i]]]))

    if (type == "choice") {
      phat_0 <- .quiet(

        predict(

          model,
          newdata = BART::bartModelMatrix(X_pred0),
          mc.cores = cores
        )
      )$prob.test
    } else {
      phat_0 <- .quiet(

        predict(

          model,
          newdata = BART::bartModelMatrix(X_pred0),
          mc.cores = cores
        )
      )$yhat.test
    }


    for (att_level in att_levels) {

      X_pred1 <- X_pred0

      X_pred1[[attribs[i]]] <- factor(att_level,
                                      levels = levels(X_pred0[[attribs[i]]]))

      # Get predictions
      if (type == "choice") {
        phat_1 <- .quiet(

          predict(

            model,
            newdata = BART::bartModelMatrix(X_pred1),
            mc.cores = cores

          )
        )$prob.test
      } else {
        phat_1 <- .quiet(

          predict(

            model,
            newdata = BART::bartModelMatrix(X_pred1),
            mc.cores = cores

          )
        )$yhat.test
      }

      # Get OMCE for single attribute-level comparison and store
      amce_est <- mean(colMeans(phat_1) - colMeans(phat_0))

      if (method == "bayes") {

        # Save interval as vector to make code easier to read.
        intvl <- c(alpha/2, (1-alpha/2))

        # Calculate distribution of marginal effects
        var_z <- phat_1 - phat_0

        # Calculate IMCE interval at this point to avoid holding many frames in memory
        imce_ci <- stats::quantile(var_z, intvl)

        amce_lower <- imce_ci[1]
        amce_upper  <- imce_ci[2]

        rm(var_z, imce_ci)

      }

      results[nrow(results)+1,] <- list(attribute = attribs[i],
                                     level = as.character(att_level),
                                     AMCE = amce_est,
                                     AMCE_lower = amce_lower,
                                     AMCE_upper = amce_upper)

      rm(X_pred1, phat_1)
      gc()

    }

  }

  out_obj <- list(amces = results,
                  alpha = alpha)

  class(out_obj) <- "cjbart"

  return(out_obj)

}
