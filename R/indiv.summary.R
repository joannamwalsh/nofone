#' @title indiv.summary
#'
#' @description Creates density plots for each parameter and a summary table
#' results.table creates an overall data frame of organized results for everyone

#' Creates kernel density plots and summary tables for each parameter from an
#' individual analysis
#'
#' @param samples MCMC samples of each parameter
#' @return A list containing a list of kernel density plots for each parameter
#' and a summary table of the parameter distributions
#' @export

indiv.summary <- function(samples) {
  #create kernel density plots of each parameter
  sampledata <- as.data.frame(samples)
  plotlist <- list()
  for (col in 1:ncol(sampledata)) {
    col.name <- names(sampledata)[col]
    plotlist[[col]] <- ggplot(data = sampledata) +
      geom_density(aes(x = .data[[col.name]])) +
      geom_vline(aes(xintercept = mean(.data[[col.name]])),
                 linetype = "dashed") +
      theme_bw()
  }

  #create a parameter summary table
  extraparamtable <- data.frame(value = names(sampledata))
  for (row in 1:nrow(extraparamtable)) {
    extraparamtable$mean[row] <-
      mean(samples[, extraparamtable$value[row]])
    extraparamtable$median[row] <-
      median(samples[, extraparamtable$value[row]])
    extraparamtable$ci2.5[row] <-
      quantile(samples[, extraparamtable$value[row]], 0.025)
    extraparamtable$ci50[row] <-
      quantile(samples[, extraparamtable$value[row]], 0.5)
    extraparamtable$ci97.5[row] <-
      quantile(samples[, extraparamtable$value[row]], 0.975)
  }
  return(list(plotlist, extraparamtable))
}
