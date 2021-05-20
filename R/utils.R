# Convert character columns to factors with warning
.char_to_fact <- function(data) {

  missing_factors <- sapply(data, is.character)

  if (sum(missing_factors) > 0) {

    missing_vars <- names(missing_factors)[missing_factors]

    data[,missing_vars] <- Map(as.factor, subset(data, select = missing_vars))

    # message("The following variables were converted to factors: ",
    #         paste0(missing_vars, collapse = ", ")
    # )

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


# Implementation of Rubin's combination rules
.combine <- function(theta, var_theta) {

  m <- length(theta)

  Q_bar <- (1/m)*sum(theta)
  U_bar <- (1/m)*sum(var_theta)

  demean <- (theta-Q_bar)^2

  B <- (1/(m-1)) * sum(demean)

  Q_bar_var <- U_bar + (1 + (1/m))*B
  Q_bar_se <- sqrt(Q_bar_var)

  v_m <- (m-1)*(1+(U_bar/((1+m^-1)*B)))^2

  std_err = Q_bar_se

  return(std_err)

}
