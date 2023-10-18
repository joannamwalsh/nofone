#' @title gen.DfRect
#'
#' @description Generates rectangle probability data frames from a summary
#' dataset that can then be passed to the \code{plot.probBar} function
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
#' @param n.info The number of grouping and stratifying variables
#' (including the ID variable) included in the \code{res} dataframe
#' @return Returns a data frame with variables representing the breaks and
#' regions of probability bars that can be plotted
#' with the the \code{plot.combined.probBar} function.
#' @export

gen.DfRect <- function(res, ninfo) {
  rect.prob <- res %>%
    select(1:all_of(ninfo), secondbetter, neitherbetter, firstbetter, comp)
  rect.prob <- rect.prob %>%
    mutate(brk1 = secondbetter) %>%
    mutate(brk2 = neitherbetter + secondbetter)
  rect.prob <- rect.prob %>%
    gather(region, prob, secondbetter:firstbetter)
  rect.prob <- rect.prob %>%
    mutate(ymin = ifelse(region == "secondbetter", 0, NA)) %>%
    mutate(ymin = ifelse(region == "neitherbetter", brk1, ymin)) %>%
    mutate(ymin = ifelse(region == "firstbetter", brk2, ymin))
  rect.prob <- rect.prob %>%
    mutate(ymax = ifelse(region == "secondbetter", brk1, NA)) %>%
    mutate(ymax = ifelse(region == "neitherbetter", brk2, ymax)) %>%
    mutate(ymax = ifelse(region == "firstbetter", 1, ymax))
  rect.prob <- rect.prob %>%
    mutate(region = factor(region, c("firstbetter", "neitherbetter",
                                     "secondbetter")))
  rect.prob <- rect.prob %>%
    mutate(which_comp = paste0(comp, ".", region))
  rect.prob <- rect.prob %>%
    mutate(which_comp = factor(which_comp))
  return(rect.prob)
}
