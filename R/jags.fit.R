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
  start <- mcpar(samples[[1]])[1]
  end <- mcpar(samples[[1]])[2]
  mid <- (end + start-1)/2
  burnin <- ceiling(end - mid)
  samples <- window(samples, start = mid+1, end = end, frequency = 1)
  #keep the last half of the converged sequence
  samples <- new.mcmc(samples)

  # check if the samples reach the number of n.run
  eff.conv <- effectiveSize(samples)[pars.save]
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

    eff.conv <- effectiveSize(samples)[pars.save]
    print('Current Min ESS:')
    print(which(eff.conv==min(eff.conv)))
    print(min(eff.conv))
  }
  n.thin <- ceiling(nrow(samples[[1]])*n.chains/n.run/2)
  #Divide by 2 to make sure there are enough samples
  samples.thinned <- window(samples, 1, dim(samples[[1]])[1], n.thin)
  eff.conv <- effectiveSize(samples.thinned)[pars.save]
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
    eff.conv <- effectiveSize(samples.thinned)[pars.save]
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
             eff = min(effectiveSize(samples)[pars.save]))
  return(out)
}
