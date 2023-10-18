#' @title genDfRect
#'
#' @description Generates rectangle probability data frames from a summary
#' dataset that can then be passed to the \code{plotprobBar} function
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
#' @param ninfo The number of grouping and stratifying variables
#' (including the ID variable) included in the \code{res} dataframe
#' @return Returns a data frame with variables representing the breaks and
#' regions of probability bars that can be plotted with the
#' \code{plotprobBar} function.
#' @import tidyverse
#' @export

genDfRect <- function(res, ninfo) {
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

#' @title plotprobBar
#'
#' @description Creates a stacked probability bar plot
#'
#' @param data A data frame containing the variables to be plotted, created by
#' the \code{genDfRect} function
#' @param clinicaldiff The lowest possible practical clinical difference between
#' treatment coefficients
#' @param title The desired title for the plot
#' @param title.size The desired title size for the plot
#' @param facets The desired facets for the plot
#' @param x.colname The name of the variable to be plotted on the x axis
#' (this is the "ID" variable)
#' @param y.colname The name of the variable to be plotted on the y axis
#' (this is the "prob" variable from the \code{genDfRect} function)
#' @param bar.width The desired width of the probability bars
#' @param bar.col The desired color of the probability bars
#' @param bar.alpha The alpha values for the three different bars in the
#' probability bar plot
#' @param ref.col The color of the y = 0.5 dashed reference line
#' @param x.label The desired label for the x axis
#' @param panel.y An indicator for whether the y axis facet labels should be
#' included in the plot
#' @param label.y An indicator for whether the y axis label should be included
#' in the plot
#' @return Returns a probability bar plot with the first facet variable as a
#' y-axis facet, the second facet variable as an x-axis facet (if specified),
#' and the probabilities from the
#' \code{y.colname} variable plotted as stacked bar plots for each
#' \code{x.colname} value
#' @import tidyverse
#' @export
#'

plotprobBar <- function(data, clinicaldiff, title, title.size = 10, facets,
                         x.colname = "ID", y.colname = "prob",
                         bar.width = 0.7, bar.col = "black",
                         bar.alpha = c(1, 0.5, 0.2), ref.col = "grey60",
                         x.label = "", panel.y, label.y) {

  legend.labels <- c(paste("Probability 1st treatment better than 2nd",
                           "(by at least", clinicaldiff, "points)"),
                     paste("Probability 1st treatment same as 2nd",
                           "(less than a", clinicaldiff,
                           "-point change in either direction)"),
                     paste("Probability 1st treatment worse than 2nd",
                           "(by at least", clinicaldiff, "points)"))

  data[, x.colname] <- as.factor(data[, x.colname])
  xcolname = sym(x.colname)
  ycolname = sym(y.colname)
  p <- ggplot(data, aes(x = !!xcolname, y = !!ycolname)) +
    geom_bar(stat = "identity", position = "stack", aes(alpha = region),
             color = bar.col, fill = bar.col, width = bar.width) +
    geom_hline(yintercept = 0.5, color = ref.col, linetype = "dashed") +
    facet_grid(facets, scales = "free", space = "free") +
    ylab("Probability of estimated effect falling into different intervals") +
    xlab(x.label) +
    ggtitle(title) +
    scale_y_continuous(position = "right", sec.axis = dup_axis()) +
    scale_x_discrete(breaks = levels(data[, x.colname]),
                     labels = unique(data[, x.colname])) +
    scale_alpha_manual(name   = NULL,
                       values = bar.alpha,
                       breaks = c("secondbetter", "neitherbetter",
                                  "firstbetter"),
                       labels = legend.labels) +
    theme_bw() +
    theme(plot.title = element_text(size = title.size),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.y.right = element_blank(),
          axis.ticks.y.left  = element_blank(),
          axis.text.y.left   = element_blank(),
          legend.position    = "bottom",
          panel.grid.major   = element_blank(),
          panel.grid.minor   = element_blank(),
          strip.background.x = element_rect(fill = NA))

  if (panel.y) {
    p <- p +
      facet_grid(facets, scales = "free", space = "free", switch = "y") +
      theme(strip.background.y = element_rect(fill = NA))
  } else {
    p <- p +
      facet_grid(facets, scales = "free", space = "free") +
      theme(strip.background.y = element_blank(),
            strip.text.y       = element_blank())
  }

  if (!label.y) {
    p <- p +
      theme(axis.title.y = element_blank())
  }
  p
}

#' @title errdata
#'
#' @description Generates a data frame containing winsorized confidence interval
#' information that can then be passed to the \code{ploterrBar} function
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
#' @import DescTools
#' @import stats
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

#' @title ploterrBar
#'
#' @description Creates a forest error bar plot.
#'
#' @param data A data frame containing the confidence interval values to be
#' plotted
#' @param title The desired title for the plot
#' @param title.size The desired title size for the plot
#' @param facets The desired facets for the plot
#' @param yintercept A vector of the desired lower and upper bounds for the
#' middle shaded region of the plot
#' @param ref.col The desired color of the middle shaded region of the plot
#' @param errBar.width The desired width of the error bars
#' @param x.colname The name of the variable to be plotted on the x axis
#' (this is the "ID" variable)
#' @param y.colname A vector of the names of the variables to be plotted on the
#' y axis (these are the lower, middle, and upper confidence interval values to
#' be plotted)
#' @param x.label The desired label for the x axis
#' @param panel.y An indicator for whether the y axis facet labels should be
#' included in the plot
#' @param label.y An indicator for whether the y axis label should be included
#' in the plot
#' @param y.label The desired label for the y axis
#' @param clinicaldiff The lowest possible practical clinical difference between
#' treatment coefficients
#' @param cutoff Desired magnitude of the y limits of the plot
#' @return Returns a forest error bar plot with the first facet variable as a
#' y-axis facet, the second facet variable as an x-axis facet (if specified),
#' and the 95% confidence interval for each \code{x.colname} value.
#' @import tidyverse
#' @export


ploterrBar <- function(data,
                        title,
                        title.size = 10,
                        facets,
                        yintercept,
                        ref.col = "grey80",
                        errBar.width = 0.5,
                        x.colname ="ID",
                        y.colname = c("comp2.5_cut",
                                      "comp50_cut",
                                      "comp97.5_cut"),
                        x.label = "",
                        panel.y,
                        label.y,
                        y.label = "Average difference",
                        clinicaldiff,
                        cutoff) {

  y.limits      = c(-1*cutoff, cutoff)
  y.breaks      = c(-1*cutoff,
                    -1*floor(cutoff/2),
                    -1*clinicaldiff,
                    clinicaldiff,
                    floor(cutoff/2),
                    cutoff)
  y.breaks.labels = c(paste("\u2264-",
                            cutoff),
                      -1*floor(cutoff/2),
                      -1*clinicaldiff,
                      clinicaldiff,
                      floor(cutoff/2),
                      paste("\u2265",
                            cutoff))
  yintercept    = c(-1*clinicaldiff,
                    clinicaldiff)
  data[, x.colname] <- as.factor(data[, x.colname])
  xcolname = sym(x.colname)
  p <- ggplot(data, aes(x = !!xcolname, y = comp50_cut))

  yintercept.used <- yintercept
  p <- p +
    scale_y_continuous(limits = y.limits,
                       breaks = y.breaks,
                       labels = y.breaks.labels,
                       position = "right",
                       sec.axis = dup_axis())

  if (length(yintercept) == 1) {
    p <- p +
      geom_hline(yintercept = yintercept.used, color = ref.col,
                 linetype = "dashed")
  } else {
    p <- p +
      geom_rect(xmin = -Inf, xmax = Inf, ymin = yintercept.used[1],
                ymax = yintercept.used[2], fill = ref.col)
  }

  p <- p +
    geom_point() +
    geom_errorbar(aes(ymin = comp2.5_cut, ymax = comp97.5_cut),
                  width = errBar.width) +
    ylab(y.label) +
    xlab(x.label) +
    scale_x_discrete(breaks = levels(data[, x.colname]),
                     labels = unique(data[, x.colname])) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(size = title.size),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.y.right = element_blank(),
          axis.ticks.y.left  = element_blank(),
          axis.text.y.left   = element_blank(),
          panel.grid.major   = element_blank(),
          panel.grid.minor   = element_blank(),
          strip.background.x = element_rect(fill = NA))

  if (panel.y) {
    p <- p +
      facet_grid(facets, space = "free", switch = "y") +
      theme(strip.background.y = element_rect(fill = NA))
  } else {
    p <- p +
      facet_grid(facets, space = "free") +
      theme(strip.background.y = element_blank(),
            strip.text.y       = element_blank())
  }

  if (!label.y) {
    p <- p +
      theme(axis.title.y = element_blank())
  }

  p

}

#' @title plotarrange
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
#' \code{plotprobBar} function and a forest plot from the \code{ploterrBar}
#' function that can be viewed using \code{grid.draw()}.
#' @import tidyverse
#' @import gridExtra
#' @import cowplot
#' @export

plotarrange <- function(res.ind, groupvars, stratvars, clinicaldiff) {

  rect.prob <- genDfRect(res = res.ind,
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
        plotprobBar(data = rect.prob[rect.prob[,strata] ==
                                        levels(rect.prob[, stratvars[1]])[i],],
                     clinicaldiff = clinicaldiff,
                     title         = strata.title[i],
                     facets        = facets,
                     panel.y       = panel.y,
                     label.y       = label.y)

      p.ci.indSep[[i]] <-
        ploterrBar(data = ci.indSep[ci.indSep[,strata] ==
                                       levels(ci.indSep[, stratvars[1]])[i],],
                    title        = strata.title[i],
                    facets       = facets,
                    clinicaldiff = clinicaldiff,
                    panel.y      = panel.y,
                    label.y      = label.y,
                    cutoff       = errdata(res.ind)[[2]])
    }
  } else {
    p.indSep[[1]] <- plotprobBar(data = rect.prob,
                                  clinicaldiff = clinicaldiff,
                                  title         = "",
                                  facets        = facets,
                                  panel.y       = T,
                                  label.y       = T)


    p.ci.indSep[[1]] <- ploterrBar(data         = ci.indSep,
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
