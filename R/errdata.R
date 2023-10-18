#' @title errdata
#'
#' @description Generates a data frame containing winsorized confidence interval
#' information that can then be passed to the \code{plot.errBar} function
#'
#' @param res A data frame containing ID and grouping variables for each
#' participant in addition to probabilities of each treatment being
#' better, worse, or the same as each other treatment under consideration. The
#' columns of this dataframe should be as follows:
#' first the ID variable, then any desired grouping variables, and then, in any
#' order,
#'      \begin{itemize}
#'          \item[] \code{comp}: the two treatments being compared in the row
#'          (in the form trt1-trt2)
#'          \item[] \code{firstbetter}: the probability that the first treatment
#'          is better than the second treatment
#'          \item[] \code{secondbetter}: the probability that the second
#'          treatment is better than the first treatment
#'          \item[] \code{neitherbetter}: 1- \code{firstbetter}-
#'          \code{secondbetter}
#'      \end{itemize}
#' @return Returns a data frame with variables representing the winsorized
#' medians and confidence intervals as well as a cutoff value representing
#' the magnitude of the y limits of the plot
#' @export

errdata <- function(res) {
  all.ci <- c(res$comp2.5, res$comp50, res$comp97.5)
  cutoff <- max(abs(unname(round(quantile(all.ci,
                                          probs = 0.025,
                                          na.rm = TRUE)))),
                abs(unname(round(quantile(all.ci,
                                          probs = 0.975,
                                          na.rm = TRUE)))))

  res$comp2.5_cut <- Winsorize(res$comp2.5, minval = -1*cutoff,
                               maxval = cutoff, probs = c(0.05, 0.95))
  res$comp50_cut <- Winsorize(res$comp50, minval = -1*cutoff,
                              maxval = cutoff, probs = c(0.05, 0.95))
  res$comp97.5_cut <- Winsorize(res$comp97.5, minval = -1*cutoff,
                                maxval = cutoff, probs = c(0.05, 0.95))

  return(list(res, cutoff))
}
