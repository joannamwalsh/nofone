% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphingfunctions.R
\name{plotprobBar}
\alias{plotprobBar}
\title{plotprobBar}
\usage{
plotprobBar(
  data,
  clinicaldiff,
  title,
  title.size = 10,
  facets,
  x.colname = "ID",
  y.colname = "prob",
  bar.width = 0.7,
  bar.col = "black",
  bar.alpha = c(1, 0.5, 0.2),
  ref.col = "grey60",
  x.label = "",
  panel.y,
  label.y
)
}
\arguments{
\item{data}{A data frame containing the variables to be plotted, created by
the \code{genDfRect} function}

\item{clinicaldiff}{The lowest possible practical clinical difference between
treatment coefficients}

\item{title}{The desired title for the plot}

\item{title.size}{The desired title size for the plot}

\item{facets}{The desired facets for the plot}

\item{x.colname}{The name of the variable to be plotted on the x axis
(this is the "ID" variable)}

\item{y.colname}{The name of the variable to be plotted on the y axis
(this is the "prob" variable from the \code{genDfRect} function)}

\item{bar.width}{The desired width of the probability bars}

\item{bar.col}{The desired color of the probability bars}

\item{bar.alpha}{The alpha values for the three different bars in the
probability bar plot}

\item{ref.col}{The color of the y = 0.5 dashed reference line}

\item{x.label}{The desired label for the x axis}

\item{panel.y}{An indicator for whether the y axis facet labels should be
included in the plot}

\item{label.y}{An indicator for whether the y axis label should be included
in the plot}
}
\value{
Returns a probability bar plot with the first facet variable as a
y-axis facet, the second facet variable as an x-axis facet (if specified),
and the probabilities from the
\code{y.colname} variable plotted as stacked bar plots for each
\code{x.colname} value
}
\description{
Creates a stacked probability bar plot
}
