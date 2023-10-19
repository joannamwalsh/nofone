#' @title nof1.data
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
#' \item{bs_df}{The b-spline basis matrix from the bs function in the splines
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
#' @import splines
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
      bs.design.matrix <- splines::bs(cent.y.time, knots = knots)
    }
    # otherwise set desired degrees of freedom
    else {
      bs.design.matrix <- splines::bs(cent.y.time, df = bs.df)
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
                               nof1$beta.prior,
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


#' @title nof1.normal.rjags
#'
#' @description Creates an Rjags model file code object based on the information
#' provided by the user that will later be used for individual analysis.
#'
#' @param nobs Total number of observations in the study
#' @param corr.y Indicator for whether the correlation among measurements
#' should be modeled. The default is \code{F}.
#' @param Treat.name Vector of unique treatment names
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

#' @title nof1.inits.normal
#'
#' @description The nof1.inits.normal function generates initial parameter
#' values to be passed to the \code{jags.fit} function.
#'
#' @param data A list containing Y (the outcome), vectors of 0s and 1s in
#' the model matrix for each \emph{Treat.name}, indicating whether or not the
#' response was recorded under the treatment, and base spline information
#' @param Treat.name Vector of unique treatment names
#' @param beta.prior Prior for the treatment-specific intercept. It should be a
#' list, where the first element is the distribution and the next two are the
#' parameters associated with the distribution. For example,
#' list("dnorm", 0, 1e-6) gives a normal prior with mean 0 and standard
#' deviation 1e-6. If truncation is desired, the last two parameters should be
#' the upper and lower limits for the truncation.
#' @param bs.trend Indicator for whether the model should adjust for trend
#' using splines. The default is \code{F}.
#' @param bs_df The b-spline basis matrix from the bs function in the splines
#' package
#' @param n.chains Number of chains to run
#' @return A list of the initial values for each treatment-specific intercept
#' and spline term to be used in the \code{jags.fit} function
#' @import stats
#' @export

nof1.inits.normal <- function(data,
                              Treat.name,
                              beta.prior,
                              bs.trend,
                              bs_df,
                              n.chains){
  Y <- data$Y
  inits <- tryCatch({
    Treat.matrix <- NULL
    for(i in Treat.name){
      Treat.matrix <- cbind(Treat.matrix, data[[paste0("Treat_", i)]])
    }

    if (bs.trend){
      for (i in 1:bs_df){
        Treat.matrix <- cbind(Treat.matrix, data[[paste0("bs", i)]])
      }
    }

    model <- lm(Y ~ Treat.matrix - 1)
    co <- coef(summary(model))

    # Generate initial values
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    for(i in 1:n.chains){

      for(j in 1:length(Treat.name)){
        initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- co[j,1] +
          rnorm(1) * co[j,2]

        # Take care of bounded beta values if uniform prior prespecified
        if (beta.prior[[1]] == "dunif") {
          initial.values[[i]][[paste0("beta_", Treat.name[j])]] <-
            max(initial.values[[i]][[paste0("beta_",
                                            Treat.name[j])]],
                beta.prior[[2]])
          initial.values[[i]][[paste0("beta_", Treat.name[j])]] <-
            min(initial.values[[i]][[paste0("beta_", Treat.name[j])]],
                beta.prior[[3]])
        }
      }

      if(bs.trend){
        for(j in 1:bs_df){
          initial.values[[i]][[paste0("eta_", j)]] <-
            co[length(Treat.name)+j, 1] +
            rnorm(1, 0, co[length(Treat.name)+j, 2])
        }
      }
    }

    if(any(is.nan(unlist(initial.values)))) initial.values <- NULL

    #if initial value generated is too big (i.e. model didn't work because
    # of sparse data), just set initial value to be NULL
    if(!is.null(initial.values)){
      if(any(abs(unlist(initial.values)) > 100)) initial.values <- NULL
    }
    initial.values
  }, error = function(err){
    print(paste("Initial value not working: ",err))
    return(NULL)
  }, warning = function(warning){
    print(paste("Warning: ", warning))
    return(NULL)
  })
  return(inits)
}

#' @title jags.fit
#'
#' @description The jags.fit function calls the \code{jags.model}, \code{adapt},
#' and \code{coda.samples} functions from the \code{rjags} package to run the
#' model and extract samples.
#'
#' @param nof1 An nof1 object of class "nof1.data" created from the
#' \code{\link{nof1.data}} function
#' @return A list with the following elements:
#' \item{burnin}{Half of the converged sequence is thrown out as a burnin}
#' \item{n.thin}{If the number of iterations that the user wants (\code{n.run})
#' is less than the number of converged sequences after burnin, the sequence is
#'  thinned and the thinning interval is stored}
#' \item{samples}{MCMC samples stored using jags. The returned samples have the
#' form of mcmc.list and can be directly applied to coda functions}
#' \item{max.gelman}{Maximum Gelman and Rubin's convergence diagnostic
#' calculated for the final sample}
#' \item{dic}{Deviance information criterion from \code{dic.samples}}
#' \item{eff}{Minimum effective sample size out of all saved parameters}
#' @import coda
#' @import rjags
#' @import stats
#' @export

jags.fit <- function(nof1){

  pars.save = nof1$pars.save
  n.chains = nof1$n.chains
  max.run = nof1$max.run
  setsize = nof1$setsize
  n.run = nof1$n.run
  conv.limit = nof1$conv.limit

  # set.seed(seed)
  mod = rjags::jags.model(textConnection(nof1$code),
                          data = nof1$data,
                          inits = nof1$inits,
                          n.chains = n.chains,
                          n.adapt = setsize,
                          quiet = TRUE)

  # adapt model
  adapted <- FALSE
  count <- 0
  while(!adapted){
    adapted <- rjags::adapt(mod, setsize, end.adaptation = FALSE)
    count <- count + 1
    if(count == 100){
      stop("algorithm has not adapted")
    }
  }

  print('model fit and adapted')

  # draw samples
  samples <- rjags::coda.samples(model = mod,
                                 variable.names = pars.save,
                                 n.iter  = setsize)

  # check convergence
  max.gelman <- find.max.gelman(samples)
  print(max.gelman)
  check <- max.gelman[[1]] > conv.limit

  # draw samples until convergence
  if(check) {
    count <- 1
    while (check & (count < max.run/setsize)) {
      samples2 <- rjags::coda.samples(mod,
                                      variable.names = pars.save,
                                      n.iter = setsize)
      samples <- add.mcmc(samples, samples2)

      count <- count + 1

      max.gelman <- find.max.gelman(samples)
      print(max.gelman)
      check <- max.gelman[[1]] > conv.limit
    }
  }

  if (check) {
    stop("code didn't converge according to gelman-rubin diagnostics.
         try increasing your max.run")
  }
  # create the mcmc object
  start <- coda::mcpar(samples[[1]])[1]
  end <- coda::mcpar(samples[[1]])[2]
  mid <- (end + start-1)/2
  burnin <- ceiling(end - mid)
  samples <- window(samples, start = mid+1, end = end, frequency = 1)
  #keep the last half of the converged sequence
  samples <- new.mcmc(samples)

  # check if the samples reach the number of n.run
  eff.conv <- coda::effectiveSize(samples)[pars.save]
  print('Min ESS after conv:')
  print(which(eff.conv==min(eff.conv)))
  print(min(eff.conv))

  #keep adding samples until desired effective sample size is reached
  while (min(eff.conv) < n.run) {

    extra.run <- max(setsize,ceiling((n.run/min(eff.conv)-1) *
                                       nrow(samples[[1]])))
    added_samples <- rjags::coda.samples(mod,
                                         variable.names = pars.save,
                                         n.iter = extra.run)
    samples <- add.mcmc(samples, added_samples)

    eff.conv <- coda::effectiveSize(samples)[pars.save]
    print('Current Min ESS:')
    print(which(eff.conv==min(eff.conv)))
    print(min(eff.conv))
  }
  n.thin <- ceiling(nrow(samples[[1]])*n.chains/n.run/2)
  #Divide by 2 to make sure there are enough samples
  samples.thinned <- window(samples, 1, dim(samples[[1]])[1], n.thin)
  eff.conv <- coda::effectiveSize(samples.thinned)[pars.save]
  print('Min ESS While Thinning:')
  print(which(eff.conv==min(eff.conv)))
  print(min(eff.conv))

  while(min(eff.conv) < n.run) {
    extra.run <- max(setsize,ceiling((n.run/min(eff.conv)-1) *
                                       nrow(samples[[1]])))
    added_samples <- rjags::coda.samples(mod, variable.names = pars.save,
                                         n.iter = extra.run)
    samples <- add.mcmc(samples, added_samples)

    n.thin <- ceiling(nrow(samples[[1]])*n.chains/n.run/2)
    samples.thinned <- window(samples, 1, dim(samples[[1]])[1], n.thin)
    eff.conv <- coda::effectiveSize(samples.thinned)[pars.save]
    print('Min ESS While Thinning:')
    print(which(eff.conv==min(eff.conv)))
    print(min(eff.conv))
  }

  # find DIC
  dic <- rjags::dic.samples(mod, n.iter = 10^4)

  out <-list(burnin = burnin,
             n.thin = n.thin,
             samples = samples,
             max.gelman = max.gelman,
             dic = dic,
             eff = min(coda::effectiveSize(samples)[pars.save]))
  return(out)
}

#' @title add.mcmc
#'
#' @description Combines two objects from \code{coda.samples} and converts the
#' result to an \code{mcmc.list} object complete with the names of the saved
#' parameters.
#'
#' @param x A list of samples created from the \code{coda.samples} function.
#' @param y A list of samples created from the \code{coda.samples} function.
#' @return An \code{mcmc.list} object containing both samples
#' @import coda
#' @export

add.mcmc <- function(x, y){

  n.chains <- length(x)
  n.var <- coda::nvar(x)
  newobjects <- vector("list", length = n.chains)

  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA,
                              nrow = 0,
                              ncol = n.var,
                              dimnames = list(NULL,
                                              dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- rbind(x[[i]], y[[i]])
    newobjects[[i]] <- coda::mcmc(newobjects[[i]])
  }
  coda::mcmc.list(newobjects)
}

#' @title find.max.gelman
#'
#' @description Calculates the maximum Gelman and Rubin's convergence diagnostic
#' for a sample generated from \code{coda.samples}.
#'
#' @param samples A list of samples created from the \code{coda.samples} function.
#' @return Returns the maximum Gelman and Rubin's convergence diagnostic for the
#' samples.
#' @import coda
#' @export

find.max.gelman <- function(samples){
  samples <- lapply(samples, function(x) { x[,colSums(abs(x)) != 0] })
  gr.uni.stats <- coda::gelman.diag(samples, multivariate = FALSE)$psrf[,1]
  largest.gr <- max(gr.uni.stats[!is.nan(gr.uni.stats)])
  out = list(largest.gr)
  names(out) <- names(which(gr.uni.stats == largest.gr))
  return(out)
}

#' @title new.mcmc
#'
#' @description Converts an object from \code{coda.samples} to an
#' \code{mcmc.list} object complete with the names of the saved parameters.
#'
#' @param x A list of samples created from the \code{coda.samples} function.
#' @return An \code{mcmc.list} object containing the samples
#' @import coda
#' @export

new.mcmc <- function(x){
  n.chains <- length(x)
  n.var <- coda::nvar(x)
  newobjects <- vector("list", length = n.chains)

  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var,
                              dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- x[[i]]
    newobjects[[i]] <- coda::mcmc(newobjects[[i]])
  }
  coda::mcmc.list(newobjects)
}

