#' @title plot.errBar
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
#' @param y.limits The limits for the y axis
#' @param y.breaks The breaks for the y axis scale
#' @param y.breaks.labels The labels for the breaks for the y axis scale
#' @return Returns a forest error bar plot with the first facet variable as a
#' y-axis facet, the second facet variable as an x-axis facet (if specified),
#' and the 95% confidence interval for each \code{x.colname} value.
#' @export


plot.errBar <- function(data,
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
