% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{summary.cjbart}
\alias{summary.cjbart}
\title{Summarizing \code{cjbart} Marginal Component Effect Estimates}
\usage{
\method{summary}{cjbart}(object, ...)
}
\arguments{
\item{object}{Object of class \code{cjbart}, the result of running \code{\link[=IMCE]{IMCE()}}}

\item{...}{Further arguments (not currently used)}
}
\value{
Data frame summarizing the average marginal component effect, the minimum and maximum values, and standard deviations for each attribute-level.
}
\description{
\code{summary} method for class "cjbart"
}
\examples{
subjects <- 5
rounds <- 2
profiles <- 2
obs <- subjects*rounds*profiles

fake_data <- data.frame(A = sample(c("a1","a2"), obs, replace = TRUE),
                        B = sample(c("b1","b2"), obs, replace = TRUE),
                        id1 = rep(1:subjects, each=rounds),
                        stringsAsFactors = TRUE)

fake_data$Y <- sample(c(0,1), obs, replace = TRUE)

cj_model <- cjbart(data = fake_data,
                   Y = "Y",
                   id = "id1")

## Skip if not Unix due to longer CPU time
if (.Platform$OS.type=='unix') {

  het_effects <- IMCE(data = fake_data,
                      model = cj_model,
                      attribs = c("A","B"),
                      ref_levels = c("a1","b1"),
                      cores = 1)

  summary(het_effects)
}


}
