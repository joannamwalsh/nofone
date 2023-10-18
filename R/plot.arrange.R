#' @title plot.arrange
#'
#' @description Generates stratified, faceted probability bar plots and forest
#' plots and arranges them nicely
#'
#' @param res.ind A data frame generated from \code{results.table} containing ID
#' and grouping variables for each participant in addition to probabilities of
#' each treatment being better, worse, or the same as each other treatment under
#' consideration. The columns of this dataframe should be as follows:
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
#' @param groupvars The names of grouping variables that the user would like to
#' be included in the summary table (not including the stratifying variables)
#' @param stratvars The names of variables that the user would like the plots to
#'  be faceted by. A separate plot will be created for each level of the first
#' stratifying variable, and these plots will be arranged next to each other. If
#' a second stratifying variable is specified, this variable will be the facet
#' on the x axis for each of the plots.
#' @param clinicaldiff The lowest possible practical clinical difference between
#' treatment coefficients
#' @return Returns a stratified probability bar plot from the
#' \code{plot.probBar} function and a forest plot from the \code{plot.errBar}
#' function that can be viewed using \code{grid.draw()}.
#' @export

plot.arrange <- function(res.ind, groupvars, stratvars, clinicaldiff) {

  rect.prob <- gen.DfRect(res = res.ind,
                          ninfo = length(groupvars) + length(stratvars) + 1)
  ci.indSep <- errdata(res = res.ind)[[1]]

  if (length(stratvars) == 2) {
    facets = paste("comp ~", stratvars[2])
  } else {
    facets = "comp ~ ."
  }

  #create the probability bar plot and the forest plot
  p.indSep <- list()
  p.ci.indSep <- list()
  if (!is.null(stratvars)) {
    strata.title = levels(rect.prob[, stratvars[1]])
    for (i in 1:length(levels(rect.prob[, stratvars[1]]))) {
      if (i == 1) {
        panel.y = T
        label.y = T
      } else {
        panel.y = F
        label.y = F
      }

      strata <- stratvars[1]
      p.indSep[[i]] <-
        plot.probBar(data = rect.prob[rect.prob[,strata] ==
                                        levels(rect.prob[, stratvars[1]])[i],],
                     clinicaldiff = clinicaldiff,
                     title         = strata.title[i],
                     facets        = facets,
                     panel.y       = panel.y,
                     label.y       = label.y)

      p.ci.indSep[[i]] <-
        plot.errBar(data = ci.indSep[ci.indSep[,strata] ==
                                       levels(ci.indSep[, stratvars[1]])[i],],
                    title        = strata.title[i],
                    facets       = facets,
                    clinicaldiff = clinicaldiff,
                    panel.y      = panel.y,
                    label.y      = label.y,
                    cutoff       = errdata(res.ind)[[2]])
    }
  } else {
    p.indSep[[1]] <- plot.probBar(data = rect.prob,
                                  clinicaldiff = clinicaldiff,
                                  title         = "",
                                  facets        = facets,
                                  panel.y       = T,
                                  label.y       = T)


    p.ci.indSep[[1]] <- plot.errBar(data         = ci.indSep,
                                    title        = "",
                                    facets       = facets,
                                    clinicaldiff   = clinicaldiff,
                                    panel.y      = T,
                                    label.y      = T,
                                    cutoff        = errdata(res.ind)[[2]])
  }

  legend.shared <- get_legend(p.indSep[[1]])
  for (i in 1:length(p.indSep)) {
    p.indSep[[i]] <- p.indSep[[i]] + theme(legend.position = "none")
  }

  p.indSep[[(length(p.indSep) + 1)]] <- legend.shared

  g.indSep <- arrangeGrob(grobs = p.indSep, ncol = (length(p.indSep) - 1),
                          nrow = 2,
                          layout_matrix = rbind(1:(length(p.indSep) -1),
                                                rep(length(p.indSep), 3)),
                          heights = c(9.5, 0.5))

  g.ci.indSep <- arrangeGrob(grobs = p.ci.indSep, ncol = length(p.ci.indSep))
  return(list(g.indSep, g.ci.indSep))
}


