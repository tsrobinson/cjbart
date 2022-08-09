#' Population-Weighted Heterogeneous Effects Analysis of Conjoint Results
#'
#' @description \code{pIMCE} calculates the population individual-level marginal component effects from a BART-estimated conjoint model, using marginal attribute distributions specified by the researcher.
#' @param data A data.frame, containing all attributes, covariates, the outcome and id variables to analyze.
#' @param model A model object, the result of running \code{cjbart()}
#' @param attribs Vector of attribute names
#' @param l Name of the attribute of interest
#' @param l_1 Attribute-level of interest for attribute *l*
#' @param l_0 Reference level for attribute *l*
#' @param marginals A named list where every element is a named vector of marginal probabilities for each corresponding attribute-level.  For example, \code{marginals = list("A1" = c("q" = 0.4, "r" = 0.6), "A2" = c("x" = 0.7, "y" = 0.2, "z" = 0.1))}
#' @param covars A data.frame of covariate information to predict pIMCEs over
#' @param method Character string, setting the variance estimation method to use. When method is "parametric", a typical combined variance estimate is employed; when \code{method = "bayes"}, the 95% posterior interval is calculated; and when \code{method = "rubin"}, combination rules are used to combine the variance analogous to in multiple imputation analysis.
#' @param alpha Number between 0 and 1 -- the significance level used to compute confidence/posterior intervals. When \code{method = "bayes"}, the posterior interval is calculated by taking the alpha/2 and (1-alpha/2) quantiles of the posterior draws. When \code{method = "rubin"}, the confidence interval equals the IMCE +/- \code{qnorm(alpha/2)}. By default, alpha is 0.05 i.e. generating a 95% confidence/posterior interval.
#' @param cores Number of CPU cores used during prediction phase
#' @param skip_checks Boolean, indicating whether to check the structure of the data (default = \code{FALSE}). Only set this to \code{TRUE} if you are confident that the data is structured appropriately
#' @details This function calculates the population-weighted IMCE, which takes into account the population distribution of profiles. Rather than average over the multiple OMCE estimates, this function generates estimated treatment effects for *all* possible potential outcomes along all attributes except the attribute of interest, and then marginalizes these over the supplied marginal distributions. Uncertainty estimates are recovered using credible intervals.
#' @return \code{pIMCE} returns a data.frame of population-weighted estimates, credible interval bounds, and the covariate information supplied
#' @seealso [cjbart::cjbart()]
#' @importFrom stats predict
#' @export
pIMCE <- function(data,
                  model,
                  attribs,
                  l,
                  l_1,
                  l_0,
                  covars,
                  marginals,
                  method = "bayes",
                  alpha = 0.05,
                  cores = 1,
                  skip_checks = FALSE) {
  
  data <- as.data.frame(data)
  
  # Get variables from trained model
  Y <- model$Y_col
  round <- model$round_col
  id_var <- model$id_col
  
  # Check optional args
  if (!(method %in% c("average","bayes","rubin"))) {
    stop("Variance estimation method must be in c('parametric','bayes','rubin'). See ?OMCE for more details.")
  }
  
  if (!skip_checks) {
    
    # Check marginals
    for (attrib in attribs[attribs != l]) {
      if (!exists(attrib, where= marginals)) {
        stop(paste0("Could not find marginal probabilities for attribute ", attrib,". Please check the marginals argument."))
      } else if (sum(marginals[[attrib]]) != 1) {
        stop(paste0("Marginal probabilities must sum to 1. First error found for attribute ", attrib))
      }
    }
    
  }
  
  data <- .char_to_fact(data)
  attribute_levels <- unique(data[[l]])
  
  if (id_var %in% covars) {
    covars[[id_var]] <- NULL
  }
  
  # Generate potential outcomes for L-1 attributes
  attribs_res <- attribs[attribs != l]
  levels_res <- apply(data[,attribs_res],2, unique, simplify = FALSE)
  pot_outs_res <- expand.grid(levels_res)
  
  tau_wgts <- apply(pot_outs_res, 1, function(x) {
    prod(sapply(attribs[attribs != l], function (y) marginals[[y]][[x[[y]]]]))
  })
  
  preds <- sapply(1:nrow(covars), function (i) {
    
    pred_data <- cbind(pot_outs_res, covars[i,names(covars)[!(names(covars) == id_var)]], row.names = NULL)
    
    phats <- list()
    for (l_d in c(l_1, l_0)) {
      D_test <- .char_to_fact(pred_data)
      
      D_test[[l_d]] <- factor(l_d, levels = attribute_levels)
      
      phats[[l_d]] <- .quiet(predict(
        model,
        newdata = BART::bartModelMatrix(D_test),
        mc.cores = cores
      )$prob.test)
    }
    
    tau_mat <- (phats[[l_1]] - phats[[l_0]])%*%tau_wgts
    pIMCE_est <- mean(tau_mat)
    pIMCE_lower <- stats::quantile(tau_mat, alpha/2)[[1]]
    pIMCE_upper <- stats::quantile(tau_mat, 1-alpha/2)[[1]]
    
    return(list(pIMCE = pIMCE_est, pIMCE_lower = pIMCE_lower, pIMCE_upper = pIMCE_upper))
    
  })
  
  return(cbind(t(preds), covars))
}
  
 