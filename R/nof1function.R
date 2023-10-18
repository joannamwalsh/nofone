#' @title nof1function
#'
#' @description This is the wrapper function that takes in a user specified
#' dataset, variables, priors, etc., runs the individual analyses on each
#' individual in the dataset, and outputs a summary table, a probability bar
#' plot, a forest plot, and density plots, dic information, and parameter
#' summary statistics for each individual.
#'
#' @param data A dataframe containing all of the N of 1 data
#' @param ID The name of the ID variable in the dataset
#' @param groupvars The names of grouping variables that the user would like to
#' be included in the summary table (not including the stratifying variables)
#' @param stratvars The names of variables that the user would like the plots to
#'  be faceted by. A separate plot will be created for each level of the first
#' stratifying variable, and these plots will be arranged next to each other. If
#' a second stratifying variable is specified, this variable will be the facet
#' on the x axis for each of the plots.
#' @param response Type of outcome. Right now, it can only be "normal" for a
#' continuous outcome
#' @param Y The name of the outcome variable in the dataset. This variable
#' should be a vector with NA's included in time order.
#' @param treatment The name of the treatment variable in the dataset. This
#' variable should have the same length as the outcome variable and can be
#' character or numeric.
#' @param corr.y Indicator for whether the correlation among measurements
#' should be modeled. The default is \code{F}.
#' @param bs.trend Indicator for whether the model should adjust for trend using
#' splines. The default is \code{F}.
#' @param y.time Parameter used for modeling splines (time when the outcome
#' is measured)
#' @param knots.bt.block Parameter used for modeling splines. Indicator for
#' whether or not knots should be set
#' at the end of each block (except for the last block). If \code{TRUE}, user
#' should then specify \code{block.no};
#' if \code{FALSE}, user should then specify \code{bs.df}.
#' @param block.no Block number used for modeling splines for the setting where
#' the knots are set at the end
#' of each block. Block number with the same length as the outcome.
#' @param bs.df Degrees of freedom for modeling splines when knots are not set
#' at the end of each block.
#' @param beta.prior Prior for the treatment-specific intercept. It should be a
#' list, where the first element is the distribution and the next two are the
#' parameters associated with the distribution. For example,
#' list("dnorm", 0, 1e-6) gives a normal prior with mean 0 and standard
#' deviation 1e-6. If truncation is desired, the last two parameters should be
#' the upper and lower limits for the truncation.
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform,
#' gamma, and half normal. It should also be a list of length 3, where the first
#' element is the distribution (one of dunif, dgamma, or dhnorm)
#' and the next two are the parameters associated with the distribution.
#' @param rho.prior Prior for the correlated error model. It should be also be
#' a list of length 3, where the first element is the distribution
#' and the next two are the parameters associated with the distribution.
#' @param eta.prior Prior for modelling spline terms. It should also be a list
#' of length 3, where first element is the distribution
#' and the next two are the parameters associated with the distribution.
#' @param n.chains Number of chains to run
#' @param max.run Maximum number of iterations that the user is willing to run.
#' If the algorithm is not
#' converging, it will run up to \code{max.run} iterations before printing a
#' message that it did not converge
#' @param setsize Number of iterations that are run between convergence checks.
#' If the algorithm converges
#' quickly, a large \code{setsize} is not needed. The number that is printed
#' between each convergence check is the
#' gelman-rubin diagnostic, and the function checks whether this number is below
#' the user-specified conv.limit.
#' @param n.run Final number of iterations that the user wants to store. If the
#' user wants to store fewer iterations
#' after the algorithm has converged, the sequence will be thinned. If the user
#' wants to store more iterations,
#' extra iterations will be run to reach the specified number of runs.
#' @param conv.limit Convergence limit for Gelman and Rubin's convergence
#' diagnostic.
#' @param clinicaldiff The lowest possible practical clinical difference between
#' treatment coefficients
#' @param alpha The desired alpha value for confidence intervals. The default is
#' 0.05.
#' @param a A vector of random seeds that is at least as long as the number of
#' individuals in the dataset
#' @return A list of the following elements:
#' \item{indivsummary}{A table of summary statistics for each individual,
#' including their ID, the specified \code{groupvars} and \code{stratvars},
#' the means and medians of the treatment specific intercepts, the posterior
#' quantiles of each treatment minus each other treatment, the posterior
#' probability that each treatment coefficient minus each other is greater than
#' 0, and the probability that each treatment is at least \code{clinicaldiff}
#' greater than and less than each other
#' treatment}
#' \item{indivresults}{A list of the following elements for each individual:
#' a list of kernel density plots for each parameter for that individual, a
#' list of kernel density plots for each parameter for that individual, and
#' the model dic for that individual's analysis}
#' \item{probplot}{A probability bar plot for each treatment comparison,
#' stratified by the specified \code{stratvars}}
#' \item{forestplot}{A winsorized forest error bar plot, stratified by the
#' specified \code{stratvars}}
#' @export

nof1function <- function(data,
                         ID,
                         groupvars,
                         stratvars,
                         response = "normal",
                         Y,
                         treatment,
                         corr.y,
                         bs.trend,
                         y.time,
                         knots.bt.block,
                         block.no,
                         bs.df,
                         beta.prior = list("dnorm", 0, 1e-6),
                         hy.prior = list("dgamma", 0.001, 0.001),
                         rho.prior = list("dunif", -1, 1),
                         eta.prior = list("dnorm", 0, 1e-6),
                         n.chains,
                         max.run,
                         setsize,
                         n.run,
                         conv.limit,
                         clinicaldiff,
                         alpha = 0.05,
                         a) {

  ## SET UP
  pt.ids <- unique(data[[ID]])
  res.ind <- c()
  indresults <- list()

  ## LOOP THROUGH INDIVIDUALS
  for (i in 1:length(pt.ids)) {

    set.seed(a[i])
    data.ind <- data[data[[ID]] == pt.ids[i], ]

    ##CREATE THE N OF 1 DATA OBJECT
    #if the response is the same for all treatments and times, or if there is
    #only one treatment observed, there are no results to summarize
    #and no comparisons to make

    if (length(unique(data.ind[[treatment]])) != 1 &
        length(unique(data.ind[[Y]][!is.na(data.ind[[Y]])])) != 1) {

      nof1.df <- nof1.data(Y              = data.ind[[Y]],
                           Treat          = data.ind[[treatment]],
                           response       = "normal",
                           corr.y         = corr.y,
                           bs.trend       = bs.trend,
                           y.time         = data.ind[[y.time]],
                           knots.bt.block = knots.bt.block,
                           block.no       = ifelse(is.null(block.no),
                                                   NULL,
                                                   data.ind[[block.no]]),
                           bs.df          = bs.df,
                           beta.prior     = beta.prior,
                           hy.prior       = hy.prior,
                           rho.prior      = rho.prior,
                           eta.prior      = eta.prior,
                           inits           = NULL,
                           n.chains        = n.chains,
                           max.run         = max.run,
                           setsize         = setsize,
                           n.run           = n.run,
                           conv.limit      = conv.limit)

      ## RUN THE ANALYSIS
      outcome <- jags.fit(nof1.df)

      ## EXTRACT THE SAMPLES
      samples <- do.call(rbind, outcome$samples)

      ## INDIVIDUAL RESULTS
      indiv.summ <- indiv.summary(samples)
      indresults[[as.character(pt.ids[i])]] <- list(indiv.summ[1],
                                                    indiv.summ[2],
                                                    outcome$dic)
      ## SUMMARIZE TREATMENT COMPARISONS
      res <- results.table(indiv = pt.ids[i],
                           treatment = treatment,
                           groupvars = groupvars,
                           data.ind = data.ind,
                           data = data,
                           stratvars = stratvars,
                           nof1 = nof1.df,
                           Y = Y,
                           samples = samples,
                           clinicaldiff = clinicaldiff,
                           alpha = alpha)

      ## COMBINE ALL INDIVIDUALS
      res.ind <- rbind(res.ind, res)
    }
  }

  ## CREATE PLOTS
  plots <- plotarrange(res.ind = res.ind,
                        groupvars = groupvars,
                        stratvars = stratvars,
                        clinicaldiff = clinicaldiff)

  ## OUTPUT RESULTS
  results <- list(indivsummary = res.ind, indivresults = indresults,
                  probplot = plots[[1]], forestplot = plots[[2]])
  return(results)
}

