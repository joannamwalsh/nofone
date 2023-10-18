#' @title nof1data
#'
#' @description Creates an N of 1 object containing the specified data, priors,
#' sample specifications, and jags model code to be passed to \code{jags.fit}
#' for individual analysis
#'
#' @param Y Outcome of the study. This should be a vector with \code{NA}'s
#' included in time order.
#' @param Treat Treatment indicator vector with same length as the outcome.
#' Can be character or numeric.
#' @param response Type of outcome. Right now, it can only be "normal" for a
#' continuous outcome
#' @param corr.y Indicator for whether the correlation among measurements
#' should be modeled. The default is \code{F}.
#' @param bs.trend Indicator for whether the model should adjust for trend
#' using splines. The default is \code{F}.
#' @param y.time Parameter used for modeling splines (time when the outcome
#' is measured). This should be a numeric vector with the same length as
#' the outcome.
#' @param knots.bt.block Parameter used for modeling splines. Indicator for
#' whether or not knots should be set at the end of each block (except for the
#' last block). If \code{TRUE}, user should then specify \code{block.no};
#' if \code{FALSE}, user should then specify \code{bs.df}.
#' @param block.no Block numbers used for modeling splines for the setting where
#' the knots are set at the end of each block. This should be a vector of block
#' numbers with the same length as the outcome.
#' @param bs.df Degrees of freedom for modeling splines when knots are not set
#' at the end of each block.
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
#' @param inits Initial values for the parameters being sampled. If left
#' unspecified, the program will generate reasonable initial values.
#' @param n.chains Number of chains to run
#' @param max.run Maximum number of iterations that the user is willing to run.
#' If the algorithm is not converging, it will run up to \code{max.run}
#' iterations before printing a message that it did not converge
#' @param setsize Number of iterations that are run between convergence checks.
#' If the algorithm converges quickly, a large \code{setsize} is not needed. The
#' number that is printed between each convergence check is the gelman-rubin
#' diagnostic, and the function checks whether this number is below the
#' user-specified convergence limit (\code{conv.limit}).
#' @param n.run Final number of independent iterations (i.e. minimum
#' effective sample size) that the user wants to store. If the user wants to
#' store fewer iterations after the algorithm has converged, the sequence will
#' be thinned. If the user wants to store more iterations, extra iterations will
#' be run to reach the specified number of runs.
#' @param conv.limit Convergence limit for Gelman and Rubin's convergence
#' diagnostic.
#' @return An object of class "nof1.data" that is a list containing the
#' following:
#' \item{nobs}{Total number of observations in the study}
#' \item{Treat.name}{Treatment names}
#' \item{n.Treat}{Number of treatments}
#' \item{response}{Type of outcome}
#' \item{corr.y}{Indicator for whether the correlation among measurements
#' will be modeled}
#' \item{bs.trend}{Indicator for whether the model will adjust for trend
#' using splines}
#' \item{priors}{Priors that the code will be using. Default priors are used if
#' priors were not specified}
#' \item{bs_df} {The b-spline basis matrix from the bs function in the splines
#' package}
#' \item{code}{Rjags model file code that is generated using information
#' provided by the user with \code{nof1.normal.rjags}. To view model file
#' inside R, use \code{cat(nof1$code).}}
#' \item{pars.save}{Parameters that will be saved}
#' \item{data}{A list containing Y (the outcome), vectors of 0s and 1s in
#' the model matrix for each \emph{Treat.name}, indicating whether or not the
#' response was recorded under the treatment, and base spline information}
#' \item{n.chains}{Number of chains to run}
#' \item{max.run}{Maximum number of iterations that the user is willing to run}
#' \item{setsize}{Number of iterations that are run between convergence checks}
#' \item{n.run}{Final number of independent iterations (i.e. minimum
#' effective sample size) that the user wants to store}
#' \item{conv.limit}{Convergence limit for Gelman and Rubin's convergence
#' diagnostic}
#' \item{inits}{Initial values for the parameters being sampled}
#' @export


nof1.data <- function(Y,
                      Treat,
                      response = "normal",
                      corr.y = F,
                      bs.trend = F,
                      y.time = NULL,
                      knots.bt.block = NULL,
                      block.no = NULL,
                      bs.df = NULL,
                      beta.prior = list("dnorm", 0, 1e-6),
                      hy.prior = list("dgamma", 0.001, 0.001),
                      rho.prior = list("dunif", -1, 1),
                      eta.prior = list("dnorm", 0, 1e-6),
                      inits = NULL,
                      n.chains = 3,
                      max.run = 100000,
                      setsize = 10000,
                      n.run = 50000,
                      conv.limit = 1.05){

  #create variable for the number of observations
  nobs <- length(Y)

  #replace blank spaces with "\\_" in the Treat variable
  Treat <- gsub(" ", "\\_", Treat)

  #sort to make the treatment order the same every time
  Treat.name <- sort(unique(Treat))

  #create a list of all information
  nof1 = list(nobs = nobs,
              Treat.name = Treat.name,
              n.Treat = length(Treat.name),
              response = response,
              corr.y = corr.y,
              bs.trend = bs.trend)

  #for splines:
  if (bs.trend){
    # center the y.time variable
    cent.y.time <- y.time - mean(y.time, na.rm=T)
    # default knots at end of each block
    if (knots.bt.block){
      knots <- cent.y.time[cumsum(rle(block.no)$lengths)]
      # remove knot at the end of last block
      knots <- knots[-length(knots)]
      bs.design.matrix <- bs(cent.y.time, knots = knots)
    }
    # otherwise set desired degrees of freedom
    else {
      bs.design.matrix <- bs(cent.y.time, df = bs.df)
    }

    nof1$bs_df <- ncol(bs.design.matrix)

    # save the columns in the bs.design.matrix in nof1
    for (i in 1:(nof1$bs_df)){
      nof1[[paste0("bs", i)]] <- bs.design.matrix[, i]
    }
  }

  #set up priors
  prior.data <- list(response = response,
                     beta.prior = beta.prior,
                     hy.prior = hy.prior,
                     rho.prior = rho.prior,
                     eta.prior = eta.prior)

  #add prior information to nof1 object
  nof1 <- c(nof1, prior.data)

  #set up code for model
  nof1$code <- nof1.normal.rjags(nof1$nobs,
                                 nof1$corr.y,
                                 nof1$Treat.name,
                                 nof1$bs.trend,
                                 nof1$bs_df,
                                 nof1$beta.prior,
                                 nof1$hy.prior,
                                 nof1$rho.prior,
                                 nof1$eta.prior)

  #check that setsize is smaller than max.run
  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }

  #set up parameters to save
  pars.save <- NULL
  for(i in nof1$Treat.name){
    pars.save <- c(pars.save, paste0("beta_", i))
  }

  pars.save <- c(pars.save, "sd")

  if (nof1$bs.trend) {
    for(i in 1:nof1$bs_df){
      pars.save <- c(pars.save, paste0("eta_", i))
    }
  }

  if (nof1$corr.y) {
    pars.save <- c(pars.save, "rho")
  }

  nof1$pars.save <- pars.save

  #create data object for analysis
  data <- list(Y = Y)

  for(i in Treat.name){
    nam <- paste("Treat_", i, sep = "")
    nam <- assign(nam, as.numeric(Treat == i))
    data[[paste0("Treat_", i)]] <- nam
    data[[paste0("Treat_", i)]][is.na(data[[paste0("Treat_", i)]])] <- 0
  }

  if (nof1$bs.trend) {
    for (i in 1:nof1$bs_df){
      data[[paste0("bs", i)]] <- nof1[[paste0("bs", i)]]
    }
  }

  nof1$data <- data

  nof1$n.chains = n.chains
  nof1$max.run = max.run
  nof1$setsize = setsize
  nof1$n.run = n.run
  nof1$conv.limit = conv.limit

  if (is.null(inits)) {
    inits <- nof1.inits.normal(nof1$data,
                               nof1$Treat.name,
                               nof1$bs.trend,
                               nof1$bs_df,
                               nof1$n.chains)
  }

  nof1$inits <- inits

  #define class of model
  class(nof1) <- "nof1.data"

  #return complete list
  return(nof1)

}
