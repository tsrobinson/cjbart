#' Generate Conjoint Model Using BART
#'
#' @description A wrapper for the \code{BART::pbart()} function.
#' @param data A data.frame or coercible, containing all attributes, covariates, the outcome and id variables to analyse.
#' @param Y Character string -- the outcome variable
#' @param id Character string -- the variable identifying individual respondents
#' @param ... Other arguments passed to pbart
#' @return Tibble summarising the average marginal component effect, the minimum and maximum values, and standard deviations for each attribute-level.
#' @details Please note, \code{cjbart} currently only works for a binary outcome.
#' @importFrom BART pbart
#' @seealso [BART::pbart()]
#' @export
cjbart <- function(data, Y, id = NULL, ...) {

  if (missing(Y)) {
    stop("Please declare the output variable using the Y argument.")
  }

  if (!missing(id)) {
    if (!(id %in% names(data))) {
      stop("Supplied id not present in data.")
    }
  }

  train_vars <- names(data)[!(names(data) %in% c(id,Y))]

  train_X <- data[, train_vars]

  train_X <- .char_to_fact(train_X)

  train_Y <- data[[Y]]

  train_model <- BART::mc.pbart(x.train = as.data.frame(train_X),
                             y.train = train_Y, ...)

  return(train_model)

}

#' Heterogeneous Effects Analysis of Conjoint Results
#'
#' @description \code{IMCE} calculates the individual-level marginal component effects from a BART-estimated conjoint model.
#' @param data A data.frame or coercible, containing all attributes, covariates, the outcome and id variables to analyse.
#' @param model A model object, the result of running \code{cjbart}
#' @param attribs Vector of attribute names
#' @param ref_levels Vector of reference levels, used to calculate marginal effects
#' @param Y Character string -- the outcome variable
#' @param id Character string -- the variable identifying individual respondents
#' @param method Character string -- the variance estimation method to use. When method = "parametric", a typical combined variance estimate is employed; when method = "bayes", the 95th posterior interval is calculated; and when method = "rubin", combination rules are used to combine the variance analogous to in multiple imputation analysis.
#' @param alpha Number between 0 and 1 -- the significance level used to compute confidence/posterior intervals. When method = "bayes", the posterior interval is calculated by taking the alpha/2 and (1-alpha/2) quantiles of the posterior draws. When method = "rubin", the confidence interval equals the IMCE +/- qnorm(alpha/2). By default, alpha = 0.05 i.e. generating a 95% confidence/posterior interval.
#' @param keep_omce Boolean, indicating whether to keep the OMCE-level results (default = \code{FALSE})
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
IMCE <- function(data, model, attribs, ref_levels, Y, id, method = "bayes", alpha = 0.05, keep_omce = FALSE, cores = 1) {

  data <- as.data.frame(data)

  # Check optional args
  if (!(method %in% c("parametric","bayes","rubin"))) {
    stop("Variance estimation method must be in c('parametric','bayes','rubin'). See ?OMCE for more details.")
  }

  # Check incoming data

  test_data <- data[,!(names(data) %in% c(attribs,Y))]
  test_data <- test_data[!duplicated(test_data),]

  if(!(id %in% names(test_data))) {
    stop("Could not find id variable in covariate matrix")
  }

  if(nrow(test_data) != length(unique(test_data[[id]]))) {
    stop("Covariates vary within id: data must not contain covariates that vary across observations of the same subject")
  }

  rm(test_data)
  gc()

  # Frames to store OMCEs, variances, and confidence intervals
  results <- data[,!(names(data) %in% c(attribs, Y))]
  var_omce <- data.frame(row.names = 1:nrow(data))
  imce_lower <- imce_upper <- data.frame(row.names = 1:length(unique(data[[id]])))

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

    phat_0 <- .quiet(

      predict(

        model,
        newdata = BART::bartModelMatrix(X_pred0),
        mc.cores = cores
      )
    )#$prob.test.mean

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
          )#$prob.test.mean

      ## Note, prob.test.mean is equivalent to
      # stats::pnorm(colMeans(pred_0$yhat.test))

      # Get OMCE for single attribute-level comparison and store
      results[[as.character(att_level)]] <- phat_1$prob.test.mean - phat_0$prob.test.mean
      var_omce[[as.character(att_level)]] <- apply(phat_1$prob.test - phat_0$prob.test, 2, var)

      if (method == "bayes") {

        # Save interval as vector to make code easier to read.
        intvl <- c(alpha/2, (1-alpha/2))

        # Calculate distribution of marginal effects
        var_z <- phat_1$prob.test - phat_0$prob.test

        # Calculate IMCE interval at this point to avoid holding many frames in memory
        imce_ci <- sapply(

          unique(data[[id]]), function(s) quantile(var_z[,data[[id]] == s], intvl)

        )

        imce_lower[[as.character(att_level)]] <- imce_ci[1,]
        imce_upper[[as.character(att_level)]] <- imce_ci[2,]

      }

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
      id
      )
    )

  results_imce <- stats::aggregate(formula = agg_formula,
                            data = results,
                            FUN = mean)

  # Double check nrow now we have covariates recovered
  if (!(nrow(covars) == nrow(results_imce))) {

    warning("Number of unique covariate rows does not match number of ids -- attempting to merge data, but please check results.")

  }

  results_imce <- merge(results_imce, covars, by = id)

  if (method == "rubin") {

    results_var <- sapply(colnames(var_omce), function (x) {
      sapply(results_imce[[id]], function (y) {

        .combine(theta = results[results[[id]] == y, x],
                 var_theta = var_omce[results[[id]] == y,x])
        }
      )
    })

    imce_upper <- sapply(colnames(results_var), function (x) results_imce[[x]] + qnorm(1-(alpha/2))*results_var[,x])
    imce_lower <- sapply(colnames(results_var), function (x) results_imce[[x]] + qnorm(alpha/2)*results_var[,x])
  }

  imce_upper <- cbind(results_imce[[id]], imce_upper)
  imce_lower <- cbind(results_imce[[id]], imce_lower)

  colnames(imce_upper)[1] <- id
  colnames(imce_lower)[1] <- id

  out_obj <- list(imce = results_imce,
                  imce_lower = as.data.frame(imce_lower),
                  imce_upper = as.data.frame(imce_upper),
                  alpha = alpha,
                  att_levels = out_levels,
                  id = id)

  if (keep_omce) {
    out_obj$omce <- results
  }

  if (method == "rubin") {
    out_obj$imce_var = results_var
  }

  class(out_obj) <- "cjbart"

  return(out_obj)

}

