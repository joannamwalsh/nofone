% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysisfunctions.R
\name{nof1.normal.rjags}
\alias{nof1.normal.rjags}
\title{nof1.normal.rjags}
\usage{
nof1.normal.rjags(
  nobs,
  corr.y,
  Treat.name,
  bs.trend,
  bs_df,
  beta.prior = list("dnorm", 0, 1e-06),
  hy.prior = list("dgamma", 0.001, 0.001),
  rho.prior = list("dunif", -1, 1),
  eta.prior = list("dnorm", 0, 1e-06)
)
}
\arguments{
\item{nobs}{Total number of observations in the study}

\item{corr.y}{Indicator for whether the correlation among measurements
should be modeled. The default is \code{F}.}

\item{Treat.name}{Vector of unique treatment names}

\item{bs.trend}{Indicator for whether the model should adjust for trend
using splines. The default is \code{F}.}

\item{bs_df}{The b-spline basis matrix from the bs function in the splines
package}

\item{beta.prior}{Prior for the treatment-specific intercept. It should be a
list, where the first element is the distribution and the next two are the
parameters associated with the distribution. For example,
list("dnorm", 0, 1e-6) gives a normal prior with mean 0 and standard
deviation 1e-6. If truncation is desired, the last two parameters should be
the upper and lower limits for the truncation.}

\item{hy.prior}{Prior for the heterogeneity parameter. Supports uniform,
gamma, and half normal. It should also be a list of length 3, where first
element is the distribution (one of dunif, dgamma, or dhnorm) and the next
two are the parameters associated with the distribution.}

\item{rho.prior}{Prior for the correlated error model. It should be also be
a list of length 3, where the first element is the distribution and the next
two are the parameters associated with the distribution.}

\item{eta.prior}{Prior for modelling spline terms. It should also be a list
of length 3, where first element is the distribution and the next two are
the parameters associated with the distribution.}
}
\value{
The Rjags model file code for individual analysis
}
\description{
Creates an Rjags model file code object based on the information
provided by the user that will later be used for individual analysis.
}
