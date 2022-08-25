#' Population-Weighted Heterogeneous Effects Analysis of Conjoint Results
#'
#' @description \code{pIMCE} calculates the population individual-level marginal component effects from a BART-estimated conjoint model, using marginal attribute distributions specified by the researcher.
#' @param model A model object, the result of running \code{cjbart()}
#' @param covar_data A data.frame of covariate information to predict pIMCEs over
#' @param attribs Vector of attribute names
#' @param l Name of the attribute of interest
#' @param l_1 Attribute-level of interest for attribute *l*
#' @param l_0 Reference level for attribute *l*
#' @param marginals A named list where every element is a named vector of marginal probabilities for each corresponding attribute-level.  For example, \code{marginals = list("A1" = c("q" = 0.4, "r" = 0.6), "A2" = c("x" = 0.7, "y" = 0.2, "z" = 0.1))}
#' @param method Character string, setting the variance estimation method to use. When method is "parametric", a typical combined variance estimate is employed; when \code{method = "bayes"}, the 95% posterior interval is calculated; and when \code{method = "rubin"}, combination rules are used to combine the variance analogous to in multiple imputation analysis.
#' @param alpha Number between 0 and 1 -- the significance level used to compute confidence/posterior intervals. When \code{method = "bayes"}, the posterior interval is calculated by taking the alpha/2 and (1-alpha/2) quantiles of the posterior draws. When \code{method = "rubin"}, the confidence interval equals the IMCE +/- \code{qnorm(alpha/2)}. By default, alpha is 0.05 i.e. generating a 95% confidence/posterior interval.
#' @param cores Number of CPU cores used during prediction phase
#' @param skip_checks Boolean, indicating whether to check the structure of the data (default = \code{FALSE}). Only set this to \code{TRUE} if you are confident that the data is structured appropriately
#' @param verbose Boolean, indicating whether to print progress (default = TRUE)
#' @details This function calculates the population-weighted IMCE, which takes into account the population distribution of profiles. Rather than average over the multiple OMCE estimates, this function generates estimated treatment effects for *all* possible potential outcomes along all attributes except the attribute of interest, and then marginalizes these over the supplied marginal distributions. Uncertainty estimates are recovered using credible intervals.
#' @return \code{pIMCE} returns a data.frame of population-weighted estimates, credible interval bounds, and the covariate information supplied
#' @seealso [cjbart::cjbart()]
#' @importFrom stats predict
#' @export
pIMCE <- function(model,
                  covar_data,
                  attribs,
                  l,
                  l_1,
                  l_0,
                  marginals,
                  method = "bayes",
                  alpha = 0.05,
                  cores = 1,
                  skip_checks = FALSE,
                  verbose = TRUE) {

  # Get variables from trained model
  Y <- model$Y_col
  round <- model$round_col
  id_var <- model$id_col
  type <- model$type

  # Check optional args
  if (!(method %in% c("average","bayes","rubin"))) {
    stop("Variance estimation method must be in c('parametric','bayes','rubin'). See ?OMCE for more details.")
  }

  if (!skip_checks) {

    # Check marginals
    for (attrib in attribs[attribs != l]) {
      if (!exists(attrib, where= marginals)) {

        unq_levels <- model$factor_levels[[attrib]]
        unq_levels_n <- length(unq_levels)
        marginals[[attrib]] <- rep(1/unq_levels_n,unq_levels_n)
        names(marginals[[attrib]]) <- unq_levels

        warning(paste0("Marginal probabilities not provided for attribute ", attrib,". Assuming a uniform distribution."))
      } else if (sum(marginals[[attrib]]) != 1) {
        stop(paste0("Marginal probabilities must sum to 1. First error found for attribute ", attrib))
      }
    }

  }

  for (v in names(model$factor_levels)) {
    if (v %in% colnames(covar_data)) {
      covar_data[[v]] <- factor(covar_data[[v]], levels = model$factor_levels[[v]])
    }
  }

  if (id_var %in% names(covar_data)) {
    id_vec <- covar_data[[id_var]]
    covar_data[[id_var]] <- NULL
  } else {
    id_vec <- NULL
  }

  # Generate potential outcomes for L-1 attributes
  attribs_res <- attribs[attribs != l]
  levels_res <- model$factor_levels[attribs_res]
  pot_outs_res <- expand.grid(levels_res)

  for (v in names(pot_outs_res)) {
    pot_outs_res[[v]] <- factor(pot_outs_res[[v]], levels = model$factor_levels[[v]])
  }

  tau_wgts <- apply(pot_outs_res, 1, function(x) {
    prod(sapply(attribs[attribs != l], function (y) marginals[[y]][[x[[y]]]]))
  })

  preds <- sapply(1:nrow(covar_data), function (i) {

    if (verbose) {
      cat(paste0("\rPredicting covariate profile ",i,"/",nrow(covar_data),""))
    }

    pred_data <- cbind(pot_outs_res,
                       covar_data[i,], row.names = NULL)

    phats <- list()
    for (l_d in c(l_1, l_0)) {

      pred_data[[l]] <- factor(l_d, levels = model$factor_levels[[l]])

      if (type == "choice") {
        phats[[l_d]] <- .quiet(predict(
          model,
          newdata = BART::bartModelMatrix(pred_data)[,names(model$varcount.mean)],
          mc.cores = cores
        )$prob.test)
      } else {
        phats[[l_d]] <- .quiet(predict(
          model,
          newdata = BART::bartModelMatrix(pred_data)[,names(model$varcount.mean)],
          mc.cores = cores
        )$yhat.test)
      }
    }

    omce_mat <- phats[[l_1]] - phats[[l_0]]
    pIMCE_z <- omce_mat%*%tau_wgts
    pIMCE_est <- mean(pIMCE_z)
    pIMCE_lower <- stats::quantile(pIMCE_z, alpha/2)[[1]]
    pIMCE_upper <- stats::quantile(pIMCE_z, 1-alpha/2)[[1]]

    return(c(pIMCE = pIMCE_est, pIMCE_lower = pIMCE_lower, pIMCE_upper = pIMCE_upper))

  })

  if (verbose) {
    cat("\n")
  }

  if (!is.null(id_vec)) {
    covar_data[[id_var]] <- id_vec
  }

  return(cbind(t(preds), covar_data))
}

