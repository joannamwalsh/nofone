#' @title plot.probBar
#'
#' @description Creates a stacked probability bar plot
#'
#' @param data A data frame containing the variables to be plotted, created by
#' the \code{gen.DfRect} function
#' @param legend.labels A vector of the desired labels for the plot legends
#' (one for each of the three probability bar colors or alphas)
#' @param title The desired title for the plot
#' @param title.size The desired title size for the plot
#' @param facets The desired facets for the plot
#' @param x.colname The name of the variable to be plotted on the x axis
#' (this is the "ID" variable)
#' @param y.colname The name of the variable to be plotted on the y axis
#' (this is the "prob" variable from the \code{gen.DfRect} function)
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
#' @export
#'

plot.probBar <- function(data, clinicaldiff, title, title.size = 10, facets,
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
