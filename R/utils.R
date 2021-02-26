# Convert character columns to factors with warning
.char_to_fact <- function(data) {

  missing_factors <- sapply(data, is.character)

  if (sum(missing_factors) > 0) {

    missing_vars <- names(missing_factors)[missing_factors]

    data[,missing_vars] <- Map(as.factor, subset(data, select = missing_vars))

    warning("The following variables were converted to factors: ",
            paste0(missing_vars, collapse = ", ")
    )

  }

  return(data)
}

# Suppress console output where verbose option not present
.quiet <- function(x) {

  # Code from Hadley Wickham
  # See: https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html

  sink(tempfile())
  on.exit(sink())
  invisible(force(x))

}

# Calculate combined variance
# See: https://www.emathzone.com/tutorials/basic-statistics/combined-variance.html
.combine_variance <- function(X, V, Dr = 1000) {

  if (length(X) != length(V)) {
    stop("Incorrect dimensions to combined variance calculation")
  }

  (Dr*(sum(V) + sum((X-mean(X))^2)))/(length(X)*Dr)

}

.var_to_se <- function(V_mat, N) {

  sqrt(V_mat)/sqrt(N)

}
