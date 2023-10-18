#' @title nof1.normal.rjags
#'
#' @description Creates an Rjags model file code object based on the information
#' provided by the user that will later be used for individual analysis.
#'
#' @param nobs Total number of observations in the study
#' @param corr.y Indicator for whether the correlation among measurements
#' should be modeled. The default is \code{F}.
#' @param bs.trend Indicator for whether the model should adjust for trend
#' using splines. The default is \code{F}.
#' @param bs_df The b-spline basis matrix from the bs function in the splines
#' package
#' @param beta.prior Prior for the treatment-specific intercept. It should be a
#' list, where the first element is the distribution and the next two are the
#' parameters associated with the distribution. For example,
#' list("dnorm", 0, 1e-6) gives a normal prior with mean 0 and standard
#' deviation 1e-6. If truncation is desired, the last two parameters should be
#' the upper and lower limits for the truncation.
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform,
#' gamma, and half normal. It should also be a list of length 3, where first
#' element is the distribution (one of dunif, dgamma, or dhnorm) and the next
#' two are the parameters associated with the distribution.
#' @param rho.prior Prior for the correlated error model. It should be also be
#' a list of length 3, where the first element is the distribution and the next
#' two are the parameters associated with the distribution.
#' @param eta.prior Prior for modelling spline terms. It should also be a list
#' of length 3, where first element is the distribution and the next two are
#' the parameters associated with the distribution.
#' @return The Rjags model file code for individual analysis
#' @export

nof1.normal.rjags <- function(nobs,
                              corr.y,
                              Treat.name,
                              bs.trend,
                              bs_df,
                              beta.prior = list("dnorm", 0, 1e-6),
                              hy.prior = list("dgamma", 0.001, 0.001),
                              rho.prior = list("dunif", -1, 1),
                              eta.prior = list("dnorm", 0, 1e-6)){

  code <- paste0("model{")

  code <- paste0(code,
                 "\n\tfor (i in 1:", nobs, ") {")

  if (corr.y) {
    code <- paste0(code,
                   "\n\t\tY[i] ~ dnorm(m[i], prec*((1 - equals(i, 1)) * 1 +
                   equals(i, 1) * (1-rho^2)))",
                   "\n\t\tm[i] <- mu[i] + rho * e[i]")
  } else {
    code <- paste0(code,
                   "\n\t\tY[i] ~ dnorm(m[i], prec)",
                   "\n\t\tm[i] <- mu[i]")
  }

  code <- paste0(code,
                 "\n\t\tmu[i] <- beta_", Treat.name[1], "*Treat_",
                 Treat.name[1], "[i]")

  if (length(Treat.name) > 1){
    for(i in Treat.name[2:length(Treat.name)]){
      code <- paste0(code, " + beta_", i, "*Treat_", i, "[i]")
    }
  }

  if (bs.trend){
    for (i in 1:bs_df){
      code <- paste0(code, " + eta_", i, "*bs", i, "[i]")
    }
  }

  code <- paste0(code,
                 "\n\t}")

  if (corr.y) {
    code <- paste0(code,
                   "\n",
                   "\n\te[1] <- 0",
                   "\n\tfor (i in 2:", nobs, ") {",
                   "\n\t\te[i] <- (1 - equals(i, 1)) * (Y[i-1] - mu[i-1]) +
                   equals(i, 1) * 0",
                   "\n\t}")

    code <- paste0(code,
                   "\n\trho ~ ", rho.prior[[1]], "(",
                   rho.prior[[2]], ", ", rho.prior[[3]], ")")
  }

  for(i in Treat.name){
    code <- paste0(code,
                   "\n\tbeta_", i, " ~ ", beta.prior[[1]], "(",
                   beta.prior[[2]], ", ", beta.prior[[3]], ")")

    if (length(beta.prior) > 3) {

      if (!is.na(beta.prior[[4]])) {
        code <- paste0(code,
                       " T(", beta.prior[[4]], ",")
      } else {
        code <- paste0(code,
                       " T(,")
      }

      if(!is.na(beta.prior[[5]])) {
        code <- paste0(code,
                       beta.prior[[5]], ")")
      } else {
        code <- paste0(code,
                       ")")
      }
    }
  }
  code <- paste0(code, "\n")

  if (bs.trend){
    for (i in 1:bs_df){
      code <- paste0(code, "\n\teta_", i, " ~ ", eta.prior[[1]],
                     "(", eta.prior[[2]], ", ", eta.prior[[3]], ")")
    }
  }

  if (hy.prior[[1]] == "dunif") {
    code <- paste0(code,
                   "\n\tprec <- pow(sd, -2)",
                   "\n\tsd ~ ", hy.prior[[1]], "(", hy.prior[[2]],
                   ", ", hy.prior[[3]], ")")
  } else if(hy.prior[[1]] == "dgamma"){
    code <- paste0(code,
                   "\n\tsd <- pow(prec, -0.5)",
                   "\n\tprec ~ dgamma(", hy.prior[[2]], ", ",
                   hy.prior[[3]], ")",
                   "\n\tlogprec <- log(prec)")
  } else if(hy.prior[[1]] == "dhnorm"){
    code <- paste0(code,
                   "\n\tsd ~ dnorm(", hy.prior[[2]], ", ",
                   hy.prior[[3]], ")T(0,)",
                   "\n\tprec <- pow(sd, -2)")
  }

  code <- paste0(code, "\n}")
  return(code)
}
